# -*- coding: utf-8 -*-
"""Amber 力场基团识别器。

从 SystemData 的原子类型、键连接、电荷等信息识别化学基团。
"""

from typing import List, Dict, Set, Tuple
from collections import defaultdict

from ..core.interfaces import GroupIdentifier
from ..core.data import Group, SystemData, ResidueData, AtomData


# 强芳香类型（直接由类型名确定的芳香原子）
STRONG_AROMATIC = frozenset({
    "ca", "cg", "ch", "cm", "cn", "cp", "cq", "c1",
    "n1", "n2", "na", "nb", "nh", "ni", "nj",
    "CA", "CB", "CC", "CK", "CM", "C5", "C6", "C7", "C*", "CW", "CR", "CN", "CV", "CQ",
    "NA", "NB", "NC", "N*"
})

# H 键供体的 D 原子原子序数（N, O, S, F）
DONOR_ATOMIC_NUMBERS = frozenset({7, 8, 16, 9})

# H 键受体的原子序数（N, O, F, S）
ACCEPTOR_ATOMIC_NUMBERS = frozenset({7, 8, 9, 16})

# 卤素原子序数
HALOGEN_ATOMIC_NUMBERS = frozenset({9, 17, 35, 53})

# 金属原子序数（常见生物金属）
METAL_ATOMIC_NUMBERS = frozenset({3, 11, 12, 19, 20, 25, 26, 29, 30})

# 水分子残基名
WATER_RESIDUES = frozenset({"SOL", "HOH", "WAT"})
ACCEPTOR_TYPES = frozenset({
    "o", "o2", "oh", "os", "oe", "o1", "ow",  # GAFF 氧
    "n", "n2", "n3",                            # GAFF 氮（酰胺/胺）
    "nb", "ni", "nj", "nc", "ne", "nf", "nk",  # GAFF 芳香/胺氮（排除 na, nh）
    "s", "ss", "sh", "sx", "s2",              # GAFF 硫
    "f", "cl", "br", "i",                      # 卤素
    "O", "OH", "O2", "OS", "OW",              # Amber 蛋白氧 + 水氧
    "N", "N2", "N3",                           # Amber 蛋白氮（酰胺/胺）
    "NA", "NB", "N*", "NC",                    # Amber 芳香氮
    "S", "SH",                                 # Amber 硫
})


# 金属离子（来自 PLIP config.py）
METAL_IONS = frozenset({
    "Ca", "Co", "Mg", "Mn", "Fe", "Cu", "Zn",
    "Li", "Na", "K", "Rb", "Sr", "Cs", "Ba",
    "Cr", "Ni", "Ru", "Rh", "Pd", "Ag", "Cd",
    "La", "W", "Os", "Ir", "Pt", "Au", "Hg",
    "Ce", "Pr", "Sm", "Eu", "Gd", "Tb", "Yb", "Lu",
    "Al", "Ga", "In", "Sb", "Tl", "Pb"
})


class AmberFFGroupIdentifier(GroupIdentifier):
    """Amber 力场基团识别器。"""

    @property
    def name(self) -> str:
        return "amber_ff"

    def identify(self, system_data: SystemData) -> List[Group]:
        """从 SystemData 识别基团。"""
        bond_graph = self._build_bond_graph(system_data)
        groups: List[Group] = []
        group_id = 0

        for res in system_data.residues:
            res_groups, group_id = self._identify_residue(
                res, bond_graph, group_id)
            groups.extend(res_groups)

        # 跨残基键的 H_donor 检测
        inter_donors, group_id = self._find_inter_residue_donors(
            system_data, group_id)
        groups.extend(inter_donors)

        return groups

    def _build_bond_graph(self, sd: SystemData) -> Dict[int, Set[int]]:
        """构建全局原子连接图。"""
        graph: Dict[int, Set[int]] = defaultdict(set)
        res_lookup = {r.residue_global_idx: r for r in sd.residues}

        for res in sd.residues:
            for b in res.bonds:
                g1 = res.atoms[b.atom1_idx_in_residue].atom_global_idx
                g2 = res.atoms[b.atom2_idx_in_residue].atom_global_idx
                graph[g1].add(g2)
                graph[g2].add(g1)

        for ib in sd.inter_residue_bonds:
            r1 = res_lookup[ib.residue1_global_idx]
            r2 = res_lookup[ib.residue2_global_idx]
            g1 = r1.atoms[ib.atom_idx_in_residue1].atom_global_idx
            g2 = r2.atoms[ib.atom_idx_in_residue2].atom_global_idx
            graph[g1].add(g2)
            graph[g2].add(g1)

        return graph

    def _identify_residue(self, res: ResidueData,
                          bond_graph: Dict[int, Set[int]],
                          start_id: int) -> Tuple[List[Group], int]:
        """识别一个残基内的所有基团。"""
        groups: List[Group] = []
        gid = start_id

        # 芳香环
        rings, gid = self._find_aromatic_rings(res, bond_graph, gid)
        groups.extend(rings)

        # H 键供体/受体
        donors, gid = self._find_donors(res, bond_graph, gid)
        groups.extend(donors)

        acceptors, gid = self._find_acceptors(res, gid)
        groups.extend(acceptors)

        # 带电基团
        charged, gid = self._find_charged(res, gid)
        groups.extend(charged)

        # 卤素
        halogens, gid = self._find_halogens(res, gid)
        groups.extend(halogens)

        # 金属
        metals, gid = self._find_metals(res, gid)
        groups.extend(metals)

        # 水
        water, gid = self._find_water(res, gid)
        groups.extend(water)

        return groups, gid

    def _find_aromatic_rings(self, res: ResidueData,
                             bond_graph: Dict[int, Set[int]],
                             start_id: int) -> Tuple[List[Group], int]:
        """检测芳香环。"""
        rings = self._detect_rings(res, bond_graph)
        groups: List[Group] = []
        gid = start_id

        for ring_atoms in rings:
            types = [res.atoms[i].atom_type for i in ring_atoms]
            strong = sum(1 for t in types if t in STRONG_AROMATIC)
            aromatic = strong >= len(ring_atoms) - 1

            if not aromatic:
                continue

            global_indices = [res.atoms[i].atom_global_idx
                              for i in ring_atoms]
            groups.append(Group(
                group_id=gid, group_type="aromatic_ring",
                molecule=res.molecule_name,
                residue_name=res.residue_name,
                residue_id=res.residue_global_idx,
                atom_indices=global_indices,
                atom_types=[res.atoms[i].atom_type for i in ring_atoms],
                elements=[res.atoms[i].atom_element for i in ring_atoms],
                charges=[res.atoms[i].atom_charge for i in ring_atoms]
            ))
            gid += 1

        return groups, gid

    def _detect_rings(self, res: ResidueData,
                      bond_graph: Dict[int, Set[int]]) -> List[List[int]]:
        """检测残基内的环（边删除 + BFS）。"""
        # 构建残基内局部连接图
        local_graph: Dict[int, Set[int]] = defaultdict(set)
        global_to_local = {}
        for i, atom in enumerate(res.atoms):
            global_to_local[atom.atom_global_idx] = i

        for b in res.bonds:
            i1 = b.atom1_idx_in_residue
            i2 = b.atom2_idx_in_residue
            local_graph[i1].add(i2)
            local_graph[i2].add(i1)

        rings: List[List[int]] = []
        seen: Set[Tuple[int, ...]] = set()
        edges = [(i, j) for i in local_graph for j in local_graph[i] if i < j]

        for u, v in edges:
            found = self._bfs_ring(local_graph, u, v, max_size=8)
            for ring in found:
                key = tuple(sorted(ring))
                if key not in seen:
                    seen.add(key)
                    rings.append(ring)

        return self._remove_redundant_rings(rings)

    def _bfs_ring(self, graph: Dict[int, Set[int]],
                  u: int, v: int, max_size: int) -> List[List[int]]:
        """从边 (u,v) 出发找环。"""
        from collections import deque
        rings: List[List[int]] = []
        q = deque([(v, [v])])
        visited = {v}

        while q:
            node, path = q.popleft()
            if len(path) > max_size:
                continue
            for nb in graph[node]:
                if (node == u and nb == v) or (node == v and nb == u):
                    continue
                if nb == u and len(path) >= 2:
                    rings.append(path + [u])
                    continue
                if nb not in visited:
                    visited.add(nb)
                    q.append((nb, path + [nb]))

        return rings

    def _remove_redundant_rings(self, rings: List[List[int]]
                                ) -> List[List[int]]:
        """去除稠合冗余环。"""
        def ring_edges(r):
            n = len(r)
            return {(min(r[i], r[(i+1) % n]), max(r[i], r[(i+1) % n]))
                    for i in range(n)}

        rings_sorted = sorted(rings, key=len)
        final: List[List[int]] = []
        coverage: Set[Tuple[int, int]] = set()

        for r in rings_sorted:
            e = ring_edges(r)
            if not e.issubset(coverage):
                final.append(r)
                coverage |= e

        return final

    def _find_donors(self, res: ResidueData,
                     bond_graph: Dict[int, Set[int]],
                     start_id: int) -> Tuple[List[Group], int]:
        """检测 H 键供体（D-H 键且 q(H)>0）。"""
        groups: List[Group] = []
        gid = start_id

        for b in res.bonds:
            a1 = res.atoms[b.atom1_idx_in_residue]
            a2 = res.atoms[b.atom2_idx_in_residue]
            d_atom, h_atom = self._classify_dh_pair(a1, a2)
            if d_atom is None:
                continue
            groups.append(Group(
                group_id=gid, group_type="H_donor",
                molecule=res.molecule_name,
                residue_name=res.residue_name,
                residue_id=res.residue_global_idx,
                atom_indices=[d_atom.atom_global_idx],
                atom_types=[d_atom.atom_type],
                elements=[d_atom.atom_element],
                charges=[d_atom.atom_charge],
                metadata={"h_atom": h_atom.atom_global_idx}
            ))
            gid += 1

        return groups, gid

    def _find_inter_residue_donors(self, sd: SystemData,
                                   start_id: int) -> Tuple[List[Group], int]:
        """检测跨残基键的 H 键供体。"""
        groups: List[Group] = []
        gid = start_id
        res_lookup = {r.residue_global_idx: r for r in sd.residues}

        for ib in sd.inter_residue_bonds:
            r1 = res_lookup[ib.residue1_global_idx]
            r2 = res_lookup[ib.residue2_global_idx]
            a1 = r1.atoms[ib.atom_idx_in_residue1]
            a2 = r2.atoms[ib.atom_idx_in_residue2]
            d_atom, h_atom = self._classify_dh_pair(a1, a2)
            if d_atom is None:
                continue
            d_res = r1 if d_atom is a1 else r2
            groups.append(Group(
                group_id=gid, group_type="H_donor",
                molecule=d_res.molecule_name,
                residue_name=d_res.residue_name,
                residue_id=d_res.residue_global_idx,
                atom_indices=[d_atom.atom_global_idx],
                atom_types=[d_atom.atom_type],
                elements=[d_atom.atom_element],
                charges=[d_atom.atom_charge],
                metadata={"h_atom": h_atom.atom_global_idx}
            ))
            gid += 1

        return groups, gid

    def _classify_dh_pair(self, a1: AtomData, a2: AtomData
                          ) -> Tuple[AtomData, AtomData]:
        """判断两个原子是否构成 D-H 对。"""
        is_h1 = a1.atom_element == 'H'
        is_h2 = a2.atom_element == 'H'

        if is_h1 and not is_h2:
            h, d = a1, a2
        elif is_h2 and not is_h1:
            h, d = a2, a1
        else:
            return None, None

        # D 必须是 N/O/S/F，H 必须带正电
        if d.atom_element in ('N', 'O', 'S', 'F') and h.atom_charge > 0:
            return d, h

        return None, None

    def _find_acceptors(self, res: ResidueData,
                        start_id: int) -> Tuple[List[Group], int]:
        """检测 H 键受体。"""
        groups: List[Group] = []
        gid = start_id

        for atom in res.atoms:
            if atom.atom_type in ACCEPTOR_TYPES and atom.atom_charge < 0:
                groups.append(Group(
                    group_id=gid, group_type="H_acceptor",
                    molecule=res.molecule_name,
                    residue_name=res.residue_name,
                    residue_id=res.residue_global_idx,
                    atom_indices=[atom.atom_global_idx],
                    atom_types=[atom.atom_type],
                    elements=[atom.atom_element],
                    charges=[atom.atom_charge]
                ))
                gid += 1

        return groups, gid

    def _find_charged(self, res: ResidueData,
                      start_id: int) -> Tuple[List[Group], int]:
        """检测带电基团。"""
        groups: List[Group] = []
        gid = start_id

        for atom in res.atoms:
            if abs(atom.atom_charge) <= 0.3:
                continue
            if atom.atom_element == 'H':
                continue
            sign = "charged_positive" if atom.atom_charge > 0 else "charged_negative"
            groups.append(Group(
                group_id=gid, group_type=sign,
                molecule=res.molecule_name,
                residue_name=res.residue_name,
                residue_id=res.residue_global_idx,
                atom_indices=[atom.atom_global_idx],
                atom_types=[atom.atom_type],
                elements=[atom.atom_element],
                charges=[atom.atom_charge]
            ))
            gid += 1

        return groups, gid

    def _find_halogens(self, res: ResidueData,
                       start_id: int) -> Tuple[List[Group], int]:
        """检测卤素。"""
        groups: List[Group] = []
        gid = start_id

        for atom in res.atoms:
            if atom.atom_element in ('F', 'Cl', 'Br', 'I'):
                groups.append(Group(
                    group_id=gid, group_type="halogen",
                    molecule=res.molecule_name,
                    residue_name=res.residue_name,
                    residue_id=res.residue_global_idx,
                    atom_indices=[atom.atom_global_idx],
                    atom_types=[atom.atom_type],
                    elements=[atom.atom_element],
                    charges=[atom.atom_charge]
                ))
                gid += 1

        return groups, gid

    def _find_metals(self, res: ResidueData,
                     start_id: int) -> Tuple[List[Group], int]:
        """检测金属离子。"""
        groups: List[Group] = []
        gid = start_id

        for atom in res.atoms:
            if atom.atom_element in METAL_IONS:
                groups.append(Group(
                    group_id=gid, group_type="metal",
                    molecule=res.molecule_name,
                    residue_name=res.residue_name,
                    residue_id=res.residue_global_idx,
                    atom_indices=[atom.atom_global_idx],
                    atom_types=[atom.atom_type],
                    elements=[atom.atom_element],
                    charges=[atom.atom_charge]
                ))
                gid += 1

        return groups, gid

    def _find_water(self, res: ResidueData,
                    start_id: int) -> Tuple[List[Group], int]:
        """检测水分子。"""
        groups: List[Group] = []
        gid = start_id

        if res.residue_name in WATER_RESIDUES:
            groups.append(Group(
                group_id=gid, group_type="water",
                molecule=res.molecule_name,
                residue_name=res.residue_name,
                residue_id=res.residue_global_idx,
                atom_indices=[a.atom_global_idx for a in res.atoms],
                atom_types=[a.atom_type for a in res.atoms],
                elements=[a.atom_element for a in res.atoms],
                charges=[a.atom_charge for a in res.atoms]
            ))
            gid += 1

        return groups, gid
