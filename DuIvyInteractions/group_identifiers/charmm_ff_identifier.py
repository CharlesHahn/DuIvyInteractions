# -*- coding: utf-8 -*-
"""CHARMM 力场基团识别器（CHARMM36 / C36m 家族）。

与 AmberFFGroupIdentifier 的差异（详见 doc/charmm_identifier_design.md §3.1）：
- A 类（力场无关）：供体、卤键、金属、金属配位、配体带电 → 继承复用
- B 类（仅换表）：受体、芳香环类型表 → 换 CHARMM 特征表
- C 类（改字典）：蛋白带电 → 换 CHARMM 残基字典（LYS/ARG/HSP + ASP/GLU/CYM）
- D 类（重写逻辑）：N 端正电（NH3+）识别 —— CHARMM N 原子带负电荷，
  父类 tertamine 单原子验证失效，需覆写为 N + 键连 H 净电荷验证
"""

from typing import List, Tuple, Dict, Set

from ..core.datas import Group, SystemData, ResidueData, AtomData
from .amber_ff_identifier import (
    AmberFFGroupIdentifier,
    METAL_IONS,
)


# ============================================================
# CHARMM 类型特征表（本地 charmm36-feb2026_cgenff-5.0.ff 实测，
# 与 CHARMM-GUI 官方 top_all36_prot.rtf 逐字一致）
# ============================================================

# H 键受体类型（蛋白 C36m + CGenFF 配体氧/氮/硫 + 卤素）
CHARMM_ACCEPTOR_TYPES = frozenset({
    # 蛋白氧
    "O", "OB", "OC", "OH1", "OS",
    # CGenFF 氧（配体扩展）
    "OG2D1", "OG2D2", "OG2D3", "OG2D4", "OG2D5",
    "OG2P1", "OG2R50", "OG301", "OG302", "OG303",
    "OG304", "OG311", "OG312", "OG3C51", "OG3C61", "OG3R60",
    # 蛋白氮
    "N", "NH1", "NH2", "NH3", "NC2", "NR1", "NR2", "NR3", "NY",
    # CGenFF 氮（配体扩展）
    "NG2D1", "NG2O1", "NG2P1", "NG2R50", "NG2R51", "NG2R52",
    "NG2R60", "NG2R61", "NG2R62", "NG2RC0", "NG2S0", "NG2S1",
    "NG2S2", "NG3N1",
    # 硫
    "S", "SM", "SS",
    # 卤素（CGenFF 命名；卤素可作 H 键受体，Lin 2017）
    "FGA1", "FGA2", "FGA3", "FGR1",
    "CLGA1", "CLGA3", "CLGR1",
    "BRGA1", "BRGA2", "BRGA3", "BRGR1",
    "IGR1",
})

# 芳香环强信号类型（蛋白 + CGenFF 配体）
CHARMM_STRONG_AROMATIC = frozenset({
    # 蛋白芳香碳
    "CA", "CAI", "CPH1", "CPH2", "CPT", "CY",
    # 蛋白芳香氮
    "NR1", "NR2", "NR3", "NY",
    # CGenFF 配体芳香碳（5/6/7 元环）
    "CG2R51", "CG2R52", "CG2R53", "CG2R57",
    "CG2R61", "CG2R62", "CG2R63", "CG2R64",
    "CG2R66", "CG2R67", "CG2R71", "CG2RC0",
    # CGenFF 配体芳香氮
    "NG2R50", "NG2R51", "NG2R52", "NG2R57",
    "NG2R60", "NG2R61", "NG2R62", "NG2R67", "NG2RC0",
    # CGenFF 杂环 O/S
    "OG2R50", "SG2R50",
})

# 芳香环兼容类型（非强信号但可参与共轭环）
CHARMM_COMPATIBLE_TYPES = frozenset({
    "C",          # 主链羰基碳（环内时）
    "CC", "CD",   # 羧基/酰胺碳（环内时）
    "CG2D1", "CG2D2", "CG2DC1", "CG2DC2", "CG2DC3",  # 共轭烯
    "CG2O1", "CG2O2", "CG2O3", "CG2O4", "CG2O5", "CG2O6",  # 环内羰基
})

# 水残基名（CHARMM 默认 TIP3/HOH，GROMACS 常规 SOL）
CHARMM_WATER_RESIDUES = frozenset({"TIP3", "HOH", "SOL", "WAT"})

# 蛋白正电残基（本地 rtp 实测：净电荷 +1）
CHARMM_POSITIVE_RESIDUES = {
    "LYS": ["NZ", "HZ1", "HZ2", "HZ3"],   # ε-铵 NH3+（同 Amber）
    "ARG": ["CZ", "NE", "NH1", "NH2", "HE", "HH11", "HH12", "HH21", "HH22"],
    # HSP 咪唑双质子化：N 为负（NR3 -0.51），必须纳入咪唑碳（CPH1/CPH2 正）
    # 才使电荷中心为正（4 原子集净 -0.14 → 7 原子集净 +0.56）
    "HSP": ["ND1", "NE2", "HD1", "HE2", "CG", "CE1", "CD2"],
}

# 蛋白负电残基（本地 rtp 实测：净电荷 -1）
CHARMM_NEGATIVE_RESIDUES = {
    "ASP": ["CG", "OD1", "OD2"],
    "GLU": ["CD", "OE1", "OE2"],
    "CYM": ["SG"],                        # 硫醇根 Cys⁻
}
# 注：中性变体 LSN/ARGN/ASPP/GLUP/HSD/HSE 不进字典。


class CharmmFFGroupIdentifier(AmberFFGroupIdentifier):
    """CHARMM 力场基团识别器（CHARMM36 / C36m）。"""

    @property
    def name(self) -> str:
        return "charmm_ff"

    # ============================================================
    # B 类：H 键受体（仅换类型表）
    # ============================================================

    def _find_acceptors(self, res: ResidueData,
                        start_id: int) -> Tuple[List[Group], int]:
        """检测 H 键受体（类型 ∈ CHARMM_ACCEPTOR_TYPES + q<0）。"""
        groups: List[Group] = []
        gid = start_id

        for atom in res.atoms:
            if atom.atom_type in CHARMM_ACCEPTOR_TYPES and atom.atom_charge < 0:
                groups.append(Group(
                    group_id=gid, group_type="H_acceptor",
                    molecule=res.molecule_name,
                    residue_name=res.residue_name,
                    residue_id=res.residue_global_idx,
                    atoms=[atom]
                ))
                gid += 1

        return groups, gid

    # ============================================================
    # B 类：芳香环（换类型表，逻辑同 amber 两条件）
    # ============================================================

    def _filter_aromatic_rings(self, rings: List[List[int]],
                               res: ResidueData) -> List[List[int]]:
        """CHARMM 芳香环两条件过滤（≥n-1 强信号 + 其余兼容）。"""
        strong = CHARMM_STRONG_AROMATIC
        compatible = CHARMM_COMPATIBLE_TYPES
        aromatic_rings = []

        for ring_atoms in rings:
            types = [res.atoms[i].atom_type for i in ring_atoms]
            strong_count = sum(t in strong for t in types)
            if strong_count < len(ring_atoms) - 1:
                continue
            non_aromatic = [t for t in types if t not in strong]
            if all(t in compatible for t in non_aromatic):
                aromatic_rings.append(ring_atoms)

        return aromatic_rings

    # ============================================================
    # C 类：蛋白带电（改字典）
    # ============================================================

    def _identify_protein_charged(self, res: ResidueData,
                                  start_id: int) -> Tuple[List[Group], int]:
        """第一层：CHARMM 残基名字典识别带电残基。"""
        groups: List[Group] = []
        gid = start_id

        for dict_, gtype in ((CHARMM_POSITIVE_RESIDUES, "charged_positive"),
                             (CHARMM_NEGATIVE_RESIDUES, "charged_negative")):
            if res.residue_name not in dict_:
                continue
            atom_names = dict_[res.residue_name]
            atoms = [a for a in res.atoms if a.atom_name in atom_names]
            if atoms:
                groups.append(self._build_charged_group(
                    atoms, gtype, res, gid, "residue_name"))
                gid += 1

        return groups, gid

    # ============================================================
    # D 类：N 端正电识别（重写官能团层验证逻辑）
    # ============================================================

    def _identify_functional_group_charged(self, res: ResidueData,
                                           bond_graph: Dict[int, Set[int]],
                                           start_id: int) -> Tuple[List[Group], int]:
        """第二层：官能团模式匹配识别带电基团（CHARMM 特有 N 端修正）。

        CHARMM N 端 NH3+ 的 N 原子带负电荷（NH3 -0.30），父类 tertamine
        单原子验证会失败；此处对满足"末端铵结构"的氮（N 只连 C/H 且
        不属于胍基），用 N + 键连 H 的净电荷验证（结构判据，不依赖类型名，
        见 _is_terminal_ammonium）。其余官能团（胍基/羧酸等）沿用父类逻辑。
        """
        groups: List[Group] = []
        gid = start_id
        atom_map = {a.atom_global_idx: a for a in res.atoms}

        for atom in res.atoms:
            neighbor_indices = bond_graph.get(atom.atom_global_idx, set())
            neighbors = [atom_map[idx] for idx in neighbor_indices if idx in atom_map]

            if self._is_quartamine(atom, neighbors):
                gid = self._append_charge_group(
                    groups, [atom], "charged_positive", res, gid, "quartamine")
            elif self._is_tertamine(atom, neighbors):
                charge_atoms = [atom]
                if self._is_terminal_ammonium(atom, neighbors,
                                              bond_graph, atom_map):
                    charge_atoms = [atom] + [
                        n for n in neighbors if n.atom_element == "H"]
                gid = self._append_charge_group(
                    groups, charge_atoms, "charged_positive", res, gid, "tertamine")
            elif self._is_guanidine(atom, neighbors, bond_graph, atom_map):
                n_atoms = [n for n in neighbors if n.atom_element == 'N']
                gid = self._append_charge_group(
                    groups, [atom] + n_atoms, "charged_positive", res, gid, "guanidine")
            elif self._is_sulfonium(atom, neighbors):
                gid = self._append_charge_group(
                    groups, [atom], "charged_positive", res, gid, "sulfonium")
            elif self._is_phosphate(atom, neighbors):
                o_atoms = [n for n in neighbors if n.atom_element == 'O']
                gid = self._append_charge_group(
                    groups, [atom] + o_atoms, "charged_negative", res, gid, "phosphate")
            elif self._is_sulfonicacid(atom, neighbors):
                o_atoms = [n for n in neighbors if n.atom_element == 'O']
                gid = self._append_charge_group(
                    groups, [atom] + o_atoms, "charged_negative", res, gid, "sulfonicacid")
            elif self._is_sulfate(atom, neighbors):
                o_atoms = [n for n in neighbors if n.atom_element == 'O']
                gid = self._append_charge_group(
                    groups, [atom] + o_atoms, "charged_negative", res, gid, "sulfate")
            elif self._is_carboxylate(atom, neighbors):
                o_atoms = [n for n in neighbors if n.atom_element == 'O']
                gid = self._append_charge_group(
                    groups, [atom] + o_atoms, "charged_negative", res, gid, "carboxylate")

        return groups, gid

    def _append_charge_group(self, groups: List[Group],
                             atoms: List[AtomData], group_type: str,
                             res: ResidueData, gid: int,
                             func_group: str) -> int:
        """验证并追加一个官能团带电基团，返回新 gid。"""
        g, gid = self._verify_and_build_group(
            atoms, group_type, res, gid, func_group)
        if g:
            groups.append(g)
        return gid

    # ============================================================
    # D 类：N 端正电识别（结构判据 + 重写验证逻辑）
    # ============================================================

    @staticmethod
    def _is_terminal_ammonium(atom: AtomData, neighbors: List[AtomData],
                              bond_graph: Dict[int, Set[int]],
                              atom_map: Dict[int, AtomData]) -> bool:
        """末端铵结构判据：N 邻居∈{C,H} 且重原子邻居不连 ≥2 个 N。"""
        h_neighbors = [n for n in neighbors if n.atom_element == "H"]
        heavy = [n for n in neighbors if n.atom_element != "H"]
        if not h_neighbors or not all(n.atom_element == "C" for n in heavy):
            return False
        for c in heavy:
            cn = [atom_map[i] for i in bond_graph.get(c.atom_global_idx, set())
                  if i in atom_map]
            if sum(1 for n in cn if n.atom_element == "N") >= 2:
                return False
        return True

    # ============================================================
    # B 类：水（换残基名集合）
    # ============================================================

    def _find_water(self, res: ResidueData,
                    start_id: int) -> Tuple[List[Group], int]:
        """检测水分子（残基名 ∈ CHARMM_WATER_RESIDUES）。"""
        groups: List[Group] = []
        gid = start_id

        if res.residue_name in CHARMM_WATER_RESIDUES:
            groups.append(Group(
                group_id=gid, group_type="water",
                molecule=res.molecule_name,
                residue_name=res.residue_name,
                residue_id=res.residue_global_idx,
                atoms=res.atoms
            ))
            gid += 1

        return groups, gid