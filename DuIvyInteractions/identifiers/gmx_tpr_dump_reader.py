# -*- coding: utf-8 -*-
"""从 gmx dump 文本读取数据。

使用逐行读取方式解析，不使用正则表达式。
"""

from typing import List, Tuple, Optional
from ..core.interfaces import Reader
from ..core.data import (
    SystemData, ResidueData, AtomData, BondData, InterResidueBond
)
from ..core.constants import ATOMIC_NUMBER_TO_ELEMENT


class GmxTprDumpReader(Reader):
    """从 gmx dump 文本读取数据。"""

    @property
    def name(self) -> str:
        return "gmx_tpr_dump"

    def read(self, source: str) -> SystemData:
        """从 gmx dump 文本读取数据。

        Args:
            source: gmx dump 输出的文本文件路径

        Returns:
            SystemData 实例
        """
        # 1. 解析 molblock 段（分子数量）
        molblock_counts = self._parse_molblock(source)
        
        # 2. 解析 moltype 段（拓扑模板）
        moltypes = self._parse_dump(source)
        
        # 3. 构建 SystemData
        return self._build_system_data(moltypes, molblock_counts, source)

    def _parse_molblock(self, source: str) -> dict:
        """解析 molblock 段，获取每种分子类型的分子数量。

        Args:
            source: 文件路径

        Returns:
            {moltype_idx: num_molecules} 字典
        """
        counts = {}
        current_moltype_idx = None
        
        with open(source, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                line = line.rstrip('\n')
                
                if "molblock (" in line and "):" in line:
                    continue
                
                if "moltype" in line and "=" in line and '"' in line:
                    # 格式: '      moltype              = 0 "RBD_pro"'
                    try:
                        idx_str = line.split("=")[1].split('"')[0].strip()
                        current_moltype_idx = int(idx_str)
                    except (ValueError, IndexError):
                        current_moltype_idx = None
                    continue
                
                if "#molecules" in line and "=" in line and current_moltype_idx is not None:
                    try:
                        count = int(line.split("=")[1].strip())
                        counts[current_moltype_idx] = count
                    except (ValueError, IndexError):
                        pass
                    continue
                
                # molblock 段结束后退出（遇到非 molblock 内容）
                if counts and "molblock" not in line and "moltype" not in line and "#molecules" not in line and "posres" not in line and line.strip() and not line.startswith(" "):
                    break
        
        return counts

    def _parse_dump(self, source: str) -> list:
        """逐行解析 gmx dump 文本。

        Args:
            source: 文件路径

        Returns:
            解析后的分子类型列表
        """
        moltypes = []
        current_moltype = None
        section = ""  # "", "atoms_param", "atoms_type", "res", "ilist"
        ilist_section = ""  # Bond/Constraint/...

        with open(source, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                line = line.rstrip('\n')
                
                # 检测新 moltype
                if self._is_moltype_start(line):
                    current_moltype = self._create_moltype(line)
                    moltypes.append(current_moltype)
                    section = ""
                    continue
                
                if current_moltype is None:
                    continue
                
                # 检测 moltype 名称
                if self._is_moltype_name(line):
                    current_moltype['name'] = self._extract_moltype_name(line)
                    continue
                
                # 检测 atoms 段
                if self._is_atoms_section(line):
                    if section == "" or section == "ilist":
                        section = "atoms_param"
                    continue
                
                # 检测 type 段
                if self._is_type_section(line):
                    section = "atoms_type"
                    continue
                
                # 检测 residue 段
                if self._is_residue_section(line):
                    section = "res"
                    continue
                
                # 检测 ilist 段（Bond/Constraint等）
                if self._is_ilist_section(line):
                    ilist_section = self._extract_ilist_section(line)
                    section = "ilist"
                    continue
                
                # 根据 section 解析行内容
                if section == "atoms_param" and self._is_atom_param(line):
                    atom = self._parse_atom_param(line)
                    if atom:
                        current_moltype['atoms'].append(atom)
                
                elif section == "atoms_param" and self._is_atom_name(line):
                    self._parse_atom_name(line, current_moltype['atoms'])
                
                elif section == "atoms_type" and self._is_type_name(line):
                    self._parse_type_name(line, current_moltype['atoms'])
                
                elif section == "res" and self._is_residue_name(line):
                    residue = self._parse_residue_name(line)
                    if residue:
                        current_moltype['residues'].append(residue)
                
                elif section == "ilist" and self._is_ilist_entry(line):
                    if ilist_section == "Bond":
                        bond = self._parse_ilist_entry(line)
                        if bond:
                            current_moltype['bonds'].append((*bond, "bond"))
                    elif ilist_section == "Constraint":
                        bond = self._parse_ilist_entry(line)
                        if bond:
                            current_moltype['bonds'].append((*bond, "constrained"))
                    elif ilist_section == "Settle":
                        for bond in self._parse_settle_entry(line):
                            current_moltype['bonds'].append((*bond, "settle"))
        
        return moltypes

    def _is_moltype_start(self, line: str) -> bool:
        """检查是否是 moltype 开始行。"""
        return "moltype (" in line and "):" in line

    def _create_moltype(self, line: str) -> dict:
        """创建新的 moltype 字典。"""
        idx_str = line.split("(")[1].split(")")[0]
        return {
            'idx': int(idx_str),
            'name': "",
            'atoms': [],
            'bonds': [],      # [(a1, a2, bond_type), ...]
            'residues': []
        }

    def _is_moltype_name(self, line: str) -> bool:
        """检查是否是 moltype 名称行。"""
        return 'name="' in line and "moltype" not in line and '[' not in line

    def _extract_moltype_name(self, line: str) -> str:
        """提取 moltype 名称。"""
        # 格式: '   name="RBD_pro"'
        return line.split('"')[1]

    def _is_atoms_section(self, line: str) -> bool:
        """检查是否是 atoms 段开始。"""
        return "atoms:" in line and "atom[" not in line and "iatoms:" not in line

    def _is_type_section(self, line: str) -> bool:
        """检查是否是 type 段开始。"""
        return "type (" in line and "type[" not in line

    def _is_residue_section(self, line: str) -> bool:
        """检查是否是 residue 段开始。"""
        return "residue (" in line and "residue[" not in line

    def _is_ilist_section(self, line: str) -> bool:
        """检查是否是 ilist 段开始。"""
        ilist_sections = ["Bond:", "Constraint:", "Settle:",
                          "G96Bond:", "Angle:", "Proper Dih.:"]
        for section in ilist_sections:
            if section in line and "type=" not in line:
                return True
        return False

    def _extract_ilist_section(self, line: str) -> str:
        """提取 ilist 段名称。"""
        # 格式: "      Bond:" 或 "      Constraint:"
        return line.strip().rstrip(':')

    def _is_atom_param(self, line: str) -> bool:
        """检查是否是原子参数行。"""
        return "atom[" in line and "atomnumber=" in line

    def _parse_atom_param(self, line: str) -> Optional[dict]:
        """解析原子参数行。"""
        # 格式: '      atom[     0]={type=  0, ..., resind=    0, atomnumber=  7}'
        try:
            idx_str = line.split("[")[1].split("]")[0].strip()
            idx = int(idx_str)
            
            params = line.split("{")[1].split("}")[0]
            param_dict = {}
            for param in params.split(","):
                if "=" in param:
                    key, value = param.split("=", 1)
                    param_dict[key.strip()] = value.strip()
            
            return {
                'idx': idx,
                'type_id': int(param_dict.get('type', '0')),
                'charge': float(param_dict.get('q', '0')),
                'mass': float(param_dict.get('m', '-1')),
                'resind': int(param_dict.get('resind', '0')),
                'atomnumber': int(param_dict.get('atomnumber', '0'))
            }
        except (ValueError, IndexError):
            return None

    def _is_atom_name(self, line: str) -> bool:
        """检查是否是原子名行。"""
        return "atom[" in line and 'name="' in line

    def _parse_atom_name(self, line: str, atoms: list) -> None:
        """解析原子名行。"""
        # 格式: '      atom[     0]={name="N"}'
        try:
            idx_str = line.split("[")[1].split("]")[0].strip()
            idx = int(idx_str)
            name = line.split('"')[1]
            
            if idx < len(atoms):
                atoms[idx]['name'] = name
        except (ValueError, IndexError):
            pass

    def _is_type_name(self, line: str) -> bool:
        """检查是否是类型名行。"""
        return "type[" in line and 'name="' in line

    def _parse_type_name(self, line: str, atoms: list) -> None:
        """解析类型名行。"""
        # 格式: '      type[     0]={name="N3",nameB="N3"}'
        try:
            idx_str = line.split("[")[1].split("]")[0].strip()
            idx = int(idx_str)
            name = line.split('"')[1]
            
            if idx < len(atoms):
                atoms[idx]['type_name'] = name
        except (ValueError, IndexError):
            pass

    def _is_residue_name(self, line: str) -> bool:
        """检查是否是残基名行。"""
        return "residue[" in line and 'name="' in line

    def _parse_residue_name(self, line: str) -> Optional[Tuple[int, str, int]]:
        """解析残基名行。"""
        # 格式: '      residue[     0]={name="ASN", nr=  157, ic= ' '}'
        try:
            idx_str = line.split("[")[1].split("]")[0].strip()
            idx = int(idx_str)
            name = line.split('"')[1]
            
            # 提取 nr
            nr_str = line.split("nr=")[1].split(",")[0].strip()
            nr = int(nr_str)
            
            return (idx, name, nr)
        except (ValueError, IndexError):
            return None

    def _is_ilist_entry(self, line: str) -> bool:
        """检查是否是 ilist 条目行。"""
        return "type=" in line and "(" in line and ")" in line

    def _parse_ilist_entry(self, line: str) -> Optional[Tuple[int, int]]:
        """解析 ilist 条目行。"""
        # 格式: '         0 type=676 (BONDS)   0   4'
        try:
            parts = line.split(")")
            if len(parts) >= 2:
                atoms_str = parts[-1].strip()
                atoms = atoms_str.split()
                if len(atoms) >= 2:
                    return (int(atoms[0]), int(atoms[1]))
        except (ValueError, IndexError):
            pass
        return None

    def _parse_settle_entry(self, line: str) -> list:
        """解析 Settle 条目行，返回 3 对键。

        Settle 格式: '0 type=1525 (SETTLE)   0   1   2'
        3 个原子 a,b,c → 键 (a,b), (a,c), (b,c)
        """
        try:
            parts = line.split(")")
            if len(parts) >= 2:
                atoms_str = parts[-1].strip()
                atoms = atoms_str.split()
                if len(atoms) >= 3:
                    a, b, c = int(atoms[0]), int(atoms[1]), int(atoms[2])
                    return [(a, b), (a, c), (b, c)]
        except (ValueError, IndexError):
            pass
        return []

    def _build_system_data(self, moltypes: list, molblock_counts: dict,
                           source: str) -> SystemData:
        """构建 SystemData。"""
        residues = []
        inter_residue_bonds = []
        system_name = source.split('/')[-1].split('.')[0]
        
        # 全局偏移量
        atom_offset = 0
        residue_offset = 0
        
        for mt in moltypes:
            num_mol = molblock_counts.get(mt['idx'], 1)
            mt_res, mt_inter, atom_offset, residue_offset = self._build_moltype(
                mt, num_mol, atom_offset, residue_offset)
            residues.extend(mt_res)
            inter_residue_bonds.extend(mt_inter)
        
        return SystemData(
            system_name=system_name,
            residues=residues,
            inter_residue_bonds=inter_residue_bonds
        )

    def _build_moltype(self, mt: dict, num_mol: int,
                       atom_offset: int, residue_offset: int) -> Tuple[list, list, int, int]:
        """为一个 moltype 构建残基数据和残基间键。

        Args:
            mt: moltype 字典
            num_mol: 该类型的分子数量
            atom_offset: 全局原子索引偏移
            residue_offset: 全局残基索引偏移

        Returns:
            (residues, inter_bonds, new_atom_offset, new_residue_offset)
        """
        # 1. 模板内按 resind 分组原子
        atoms_by_resind = {}
        for atom in mt['atoms']:
            resind = atom.get('resind', 0)
            atoms_by_resind.setdefault(resind, []).append(atom)
        
        # 2. 模板内 atom_idx → (resind, local_idx)
        atom_to_local = {}
        for resind, atoms in atoms_by_resind.items():
            for local_idx, atom in enumerate(atoms):
                atom_to_local[atom['idx']] = (resind, local_idx)
        
        # 3. 模板内键分类
        template_res_global = {i: r[0] for i, r in enumerate(mt['residues'])}
        template_intra = {i: [] for i in range(len(mt['residues']))}
        template_inter = []
        
        for a1, a2, bond_type in mt['bonds']:
            r1, l1 = atom_to_local.get(a1, (None, None))
            r2, l2 = atom_to_local.get(a2, (None, None))
            if r1 is None or r2 is None:
                continue
            if r1 == r2:
                template_intra[r1].append(BondData(l1, l2, bond_type))
            else:
                template_inter.append((r1, l1, r2, l2, bond_type))
        
        # 4. 按分子数量复制
        atoms_per_mol = len(mt['atoms'])
        residues_per_mol = len(mt['residues'])
        
        all_residues = []
        all_inter_bonds = []
        
        for mol_idx in range(num_mol):
            # 原子偏移
            mol_atom_offset = atom_offset + mol_idx * atoms_per_mol
            mol_residue_offset = residue_offset + mol_idx * residues_per_mol
            
            # 复制残基
            for resind, (tpl_gidx, res_name, res_nr) in enumerate(mt['residues']):
                atoms = atoms_by_resind.get(resind, [])
                res_atoms = [
                    AtomData(
                        atom_global_idx=mol_atom_offset + a['idx'],
                        atom_idx_in_residue=local_idx,
                        atom_name=a.get('name', ''),
                        atom_type=a.get('type_name', ''),
                        atom_element=ATOMIC_NUMBER_TO_ELEMENT.get(a['atomnumber'], '?'),
                        atom_charge=a['charge'],
                        atom_mass=a['mass']
                    )
                    for local_idx, a in enumerate(atoms)
                ]
                all_residues.append(ResidueData(
                    residue_name=res_name,
                    residue_global_idx=mol_residue_offset + resind,
                    residue_idx_in_molecule=res_nr,
                    molecule_name=mt['name'],
                    atoms=res_atoms,
                    bonds=list(template_intra[resind])  # 复制键列表
                ))
            
            # 复制残基间键
            for r1, l1, r2, l2, bond_type in template_inter:
                all_inter_bonds.append(InterResidueBond(
                    mol_residue_offset + r1, l1,
                    mol_residue_offset + r2, l2, bond_type))
        
        new_atom_offset = atom_offset + num_mol * atoms_per_mol
        new_residue_offset = residue_offset + num_mol * residues_per_mol
        
        return all_residues, all_inter_bonds, new_atom_offset, new_residue_offset
