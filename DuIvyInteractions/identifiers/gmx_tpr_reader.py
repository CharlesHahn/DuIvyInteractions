# -*- coding: utf-8 -*-
"""从 tpr 二进制文件直接读取数据（使用 MDAnalysis）。"""

from typing import List, Tuple, Dict
from ..core.interfaces import Reader
from ..core.data import (
    SystemData, ResidueData, AtomData, BondData, InterResidueBond
)


class GmxTprReader(Reader):
    """从 tpr 二进制文件直接读取数据。"""

    @property
    def name(self) -> str:
        return "gmx_tpr"

    def read(self, source: str) -> SystemData:
        """从 tpr 文件读取数据。

        Args:
            source: tpr 文件路径

        Returns:
            SystemData 实例
        """
        import MDAnalysis as mda
        u = mda.Universe(source)
        return self._build_system_data(u, source)

    def _build_system_data(self, u, source: str) -> SystemData:
        """从 MDAnalysis Universe 构建 SystemData。"""
        system_name = source.split('/')[-1].split('.')[0]

        # 1. 建立 atom.index → (residue_local_idx, global_res_idx) 映射
        atom_to_res_local, atom_to_res_global = self._build_atom_maps(u)

        # 2. 构建残基数据
        residues = self._build_residues(u, atom_to_res_local)

        # 3. 构建残基间键
        inter_bonds = self._build_inter_bonds(u, atom_to_res_local,
                                                atom_to_res_global)

        return SystemData(
            system_name=system_name,
            residues=residues,
            inter_residue_bonds=inter_bonds
        )

    def _build_atom_maps(self, u) -> Tuple[Dict[int, int], Dict[int, int]]:
        """建立原子索引到残基局部索引和全局残基索引的映射。

        Returns:
            (atom_to_res_local, atom_to_res_global)
            - atom_to_res_local: atom_index → atom_idx_in_residue
            - atom_to_res_global: atom_index → residue_global_idx
        """
        atom_to_res_local = {}
        atom_to_res_global = {}

        global_res_idx = 0
        for res in u.residues:
            for local_idx, atom in enumerate(res.atoms):
                atom_to_res_local[atom.index] = local_idx
                atom_to_res_global[atom.index] = global_res_idx
            global_res_idx += 1

        return atom_to_res_local, atom_to_res_global

    def _build_residues(self, u, atom_to_res_local: Dict[int, int]
                        ) -> List[ResidueData]:
        """构建所有残基数据。"""
        residues = []
        global_res_idx = 0

        for seg in u.segments:
            for res in seg.residues:
                atoms = self._build_atoms(res, atom_to_res_local)
                bonds = self._build_intra_bonds(res, atom_to_res_local)

                residues.append(ResidueData(
                    residue_name=res.resname,
                    residue_global_idx=global_res_idx,
                    residue_idx_in_molecule=res.resid,
                    molecule_name=seg.segid,
                    atoms=atoms,
                    bonds=bonds
                ))
                global_res_idx += 1

        return residues

    def _build_atoms(self, res, atom_to_res_local: Dict[int, int]
                     ) -> List[AtomData]:
        """构建残基内的原子数据。"""
        atoms = []
        for atom in res.atoms:
            atoms.append(AtomData(
                atom_global_idx=atom.index,
                atom_idx_in_residue=atom_to_res_local[atom.index],
                atom_name=atom.name,
                atom_type=atom.type,
                atom_element=atom.element,
                atom_charge=atom.charge,
                atom_mass=atom.mass
            ))
        return atoms

    def _build_intra_bonds(self, res, atom_to_res_local: Dict[int, int]
                           ) -> List[BondData]:
        """构建残基内的键数据。"""
        bonds = []
        res_atoms = set(atom.index for atom in res.atoms)

        for bond in res.atoms.bonds:
            idx0, idx1 = bond.indices
            if idx0 in res_atoms and idx1 in res_atoms:
                bonds.append(BondData(
                    atom1_idx_in_residue=atom_to_res_local[idx0],
                    atom2_idx_in_residue=atom_to_res_local[idx1],
                    bond_type="bond"
                ))

        return bonds

    def _build_inter_bonds(self, u, atom_to_res_local: Dict[int, int],
                           atom_to_res_global: Dict[int, int]
                           ) -> List[InterResidueBond]:
        """构建残基间键数据。"""
        inter_bonds = []

        for bond in u.bonds:
            idx0, idx1 = bond.indices
            res0 = atom_to_res_global[idx0]
            res1 = atom_to_res_global[idx1]

            if res0 != res1:
                inter_bonds.append(InterResidueBond(
                    residue1_global_idx=res0,
                    atom_idx_in_residue1=atom_to_res_local[idx0],
                    residue2_global_idx=res1,
                    atom_idx_in_residue2=atom_to_res_local[idx1],
                    bond_type="bond"
                ))

        return inter_bonds
