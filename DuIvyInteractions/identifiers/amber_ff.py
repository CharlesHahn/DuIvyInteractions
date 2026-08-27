# -*- coding: utf-8 -*-
"""Amber 力场基团识别器。

从 SystemData 的原子类型、键连接、电荷等信息识别化学基团。
"""

from typing import List, Dict, Set, Tuple
from collections import defaultdict
import numpy as np

from ..core.interfaces import GroupIdentifier
from ..core.data import Group, SystemData, ResidueData, AtomData


# 芳香类型（直接由类型名确定的芳香原子）
STRONG_AROMATIC = frozenset({
    # GAFF 芳香碳
    "ca", "cg", "ch", "cm", "cn", "cp", "cq", "c1",
    # GAFF 芳香氮
    "na", "nb", "nh", "ni", "nj", "n1", "n2",
    # GAFF 芳香磷
    "pb",
    # Amber 蛋白芳香碳
    "CA", "CB", "CC", "CK", "CM", "C5", "C6", "C7", "C*", "CW", "CR", "CN", "CV", "CQ",
    # Amber 蛋白芳香氮
    "NA", "NB", "NC", "N*",
})

# 兼容原子类型（非芳香，但在 n-1 个芳香原子"强制"下可参与共轭）
COMPATIBLE_TYPES = frozenset({
    "C", "N",           # Amber14sb 歧义类型
    "os", "ss",         # GAFF 呋喃 O / 噻吩 S
    "cc", "cd",         # GAFF 非纯芳香共轭环碳
    "pc", "pd",         # GAFF 共轭环内 sp2 磷
})

# 平面性阈值（来源：PLIP AROMATIC_PLANARITY）
AROMATIC_PLANARITY = 5.0  # 度

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

# 电荷验证阈值
CHARGE_THRESHOLD = 0.1

# 蛋白正电残基（第一层：残基名字典）
POSITIVE_RESIDUES = {
    "ARG": ["CZ", "NE", "NH1", "NH2", "HE", "HH11", "HH12", "HH21", "HH22"],
    "LYS": ["NZ", "HZ1", "HZ2", "HZ3"],
    "HIP": ["ND1", "NE2", "HD1", "HE2"],
    "ORN": ["NE", "HE1", "HE2", "HE3"],
    "DAB": ["ND", "HD1", "HD2", "HD3"],
    "M3L": ["NZ", "CM1", "CM2", "CM3",
            "HM11", "HM12", "HM13", "HM21", "HM22", "HM23",
            "HM31", "HM32", "HM33"],
    "MLY": ["NZ", "CH1", "CH2",
            "HH11", "HH12", "HH13", "HH21", "HH22", "HH23"],
}

# 蛋白负电残基（第一层：残基名字典）
NEGATIVE_RESIDUES = {
    "ASP": ["CG", "OD1", "OD2"],
    "GLU": ["CD", "OE1", "OE2"],
    "CYM": ["SG"],
    "KCX": ["NZ", "CX", "OQ1", "OQ2", "HZ"],
    "PCA": ["CA", "C", "O", "N", "CD", "CG"],
    "SEP": ["OG", "P", "O1P", "O2P", "O3P"],
    "TPO": ["OG1", "P", "O1P", "O2P", "O3P"],
    "PTR": ["OH", "P", "O1P", "O2P", "O3P"],
}


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
        charged, gid = self._find_charged(res, bond_graph, gid)
        groups.extend(charged)

        # 卤键供体
        hal_donors, gid = self._find_halogen_donors(res, bond_graph, gid)
        groups.extend(hal_donors)

        # 卤键受体
        hal_acceptors, gid = self._find_halogen_acceptors(res, bond_graph, gid)
        groups.extend(hal_acceptors)

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
        """检测芳香环（三条件：n-1芳香 + 兼容原子 + 平面性）。"""
        rings = self._detect_rings(res, bond_graph)
        aromatic_rings = self._filter_aromatic_rings(rings, res)
        aromatic_rings = self._deduplicate_aromatic_rings(aromatic_rings)

        groups: List[Group] = []
        gid = start_id

        for ring_atoms in aromatic_rings:
            ring_atom_data = [res.atoms[i] for i in ring_atoms]
            groups.append(Group(
                group_id=gid, group_type="aromatic_ring",
                molecule=res.molecule_name,
                residue_name=res.residue_name,
                residue_id=res.residue_global_idx,
                atoms=ring_atom_data
            ))
            gid += 1

        return groups, gid

    def _filter_aromatic_rings(self, rings: List[List[int]],
                               res: ResidueData) -> List[List[int]]:
        """过滤出满足芳香性三条件的环。"""
        # 构建残基内局部连接图
        local_graph: Dict[int, Set[int]] = defaultdict(set)
        for b in res.bonds:
            i1, i2 = b.atom1_idx_in_residue, b.atom2_idx_in_residue
            local_graph[i1].add(i2)
            local_graph[i2].add(i1)

        aromatic_rings = []
        for ring_atoms in rings:
            n = len(ring_atoms)
            types = [res.atoms[i].atom_type for i in ring_atoms]

            # 条件1：至少 n-1 个原子在 STRONG_AROMIC 中
            strong_count = sum(1 for t in types if t in STRONG_AROMIC)
            if strong_count < n - 1:
                continue

            # 条件2：不在 STRONG_AROMIC 中的原子必须在 COMPATIBLE 中
            non_aromatic = [t for t in types if t not in STRONG_AROMIC]
            if not all(t in COMPATIBLE_TYPES for t in non_aromatic):
                continue

            # 条件3：平面性（暂不检查，需要坐标数据）
            # TODO: 当有坐标数据时，调用 _is_ring_planar 检查平面性

            aromatic_rings.append(ring_atoms)

        return aromatic_rings

    def _detect_rings(self, res: ResidueData,
                      bond_graph: Dict[int, Set[int]]) -> List[List[int]]:
        """检测残基内的所有环（BFS，无大小限制）。"""
        # 构建残基内局部连接图
        local_graph: Dict[int, Set[int]] = defaultdict(set)
        for b in res.bonds:
            i1, i2 = b.atom1_idx_in_residue, b.atom2_idx_in_residue
            local_graph[i1].add(i2)
            local_graph[i2].add(i1)

        rings: List[List[int]] = []
        seen: Set[Tuple[int, ...]] = set()
        edges = [(i, j) for i in local_graph for j in local_graph[i] if i < j]

        for u, v in edges:
            found = self._bfs_ring(local_graph, u, v)
            for ring in found:
                key = tuple(sorted(ring))
                if key not in seen:
                    seen.add(key)
                    rings.append(ring)

        return rings

    def _bfs_ring(self, graph: Dict[int, Set[int]],
                  u: int, v: int) -> List[List[int]]:
        """从边 (u,v) 出发找环（无大小限制）。"""
        from collections import deque
        rings: List[List[int]] = []
        q = deque([(v, [v])])
        visited = {v}

        while q:
            node, path = q.popleft()
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

    @staticmethod
    def _is_ring_planar(ring_atoms: List[int],
                        local_graph: Dict[int, Set[int]],
                        coordinates: np.ndarray) -> bool:
        """检查环是否平面（法向量夹角阈值 5°）。

        Args:
            ring_atoms: 环内原子的残基内索引列表
            local_graph: 残基内局部连接图
            coordinates: 原子坐标数组 (N, 3)

        Returns:
            True 如果环是平面的
        """
        normals = []
        ring_set = set(ring_atoms)

        for idx in ring_atoms:
            neighbors_in_ring = [nb for nb in local_graph[idx]
                                 if nb in ring_set]
            if len(neighbors_in_ring) < 2:
                continue
            n1, n2 = neighbors_in_ring[0], neighbors_in_ring[1]
            a_pos = coordinates[idx]
            v1 = coordinates[n1] - a_pos
            v2 = coordinates[n2] - a_pos
            normal = np.cross(v1, v2)
            norm = np.linalg.norm(normal)
            if norm > 1e-10:
                normals.append(normal / norm)

        if len(normals) < 3:
            return False

        for i in range(len(normals)):
            for j in range(i + 1, len(normals)):
                cos_angle = np.clip(np.dot(normals[i], normals[j]), -1.0, 1.0)
                angle = np.degrees(np.arccos(cos_angle))
                if AROMATIC_PLANARITY < angle < 180.0 - AROMATIC_PLANARITY:
                    return False
        return True

    def _deduplicate_aromatic_rings(self, rings: List[List[int]]
                                    ) -> List[List[int]]:
        """去重：按原子数排序，排除被小环覆盖的大环。"""
        sorted_rings = sorted(rings, key=len)
        accepted: List[List[int]] = []
        covered_atoms: Set[int] = set()

        for ring in sorted_rings:
            ring_set = set(ring)
            if not ring_set.issubset(covered_atoms):
                accepted.append(ring)
                covered_atoms.update(ring_set)

        return accepted

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
                atoms=[d_atom],
                metadata={"h_atom": h_atom}
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
                atoms=[d_atom],
                metadata={"h_atom": h_atom}
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
                    atoms=[atom]
                ))
                gid += 1

        return groups, gid

    def _find_charged(self, res: ResidueData,
                      bond_graph: Dict[int, Set[int]],
                      start_id: int) -> Tuple[List[Group], int]:
        """检测带电基团（三层递进）。

        第一层：残基名字典（蛋白残基 R 基团）
        第二层：官能团模式匹配（N/C 端 + 非蛋白残基）
        第三层：部分电荷交叉验证
        """
        groups: List[Group] = []
        gid = start_id

        # 第一层：残基名字典
        layer1_groups, gid = self._identify_protein_charged(res, gid)
        groups.extend(layer1_groups)

        # 第二层：官能团模式匹配
        layer2_groups, gid = self._identify_functional_group_charged(
            res, bond_graph, gid)
        groups.extend(layer2_groups)

        # 去重：同一残基+同类型，优先保留第一层结果
        groups = self._deduplicate_charged(groups)

        return groups, gid

    def _deduplicate_charged(self, groups: List[Group]) -> List[Group]:
        """去重：基于原子集合去重，优先保留第一层结果。"""
        seen = {}
        result = []

        for g in groups:
            # 用 frozenset(atom_indices) 作为去重 key
            key = frozenset(g.atom_indices)
            if key not in seen:
                seen[key] = g
                result.append(g)
            else:
                # 已存在，优先保留 residue_name 来源
                existing = seen[key]
                if existing.metadata.get('source') != 'residue_name' and g.metadata.get('source') == 'residue_name':
                    result.remove(existing)
                    seen[key] = g
                    result.append(g)

        return result

    def _identify_protein_charged(self, res: ResidueData,
                                  start_id: int) -> Tuple[List[Group], int]:
        """第一层：残基名字典识别带电残基。"""
        groups: List[Group] = []
        gid = start_id

        # 正电残基
        if res.residue_name in POSITIVE_RESIDUES:
            atom_names = POSITIVE_RESIDUES[res.residue_name]
            atoms = [a for a in res.atoms if a.atom_name in atom_names]
            if atoms:
                groups.append(self._build_charged_group(
                    atoms, "charged_positive", res, gid, "residue_name"))
                gid += 1

        # 负电残基
        if res.residue_name in NEGATIVE_RESIDUES:
            atom_names = NEGATIVE_RESIDUES[res.residue_name]
            atoms = [a for a in res.atoms if a.atom_name in atom_names]
            if atoms:
                groups.append(self._build_charged_group(
                    atoms, "charged_negative", res, gid, "residue_name"))
                gid += 1

        return groups, gid

    def _identify_functional_group_charged(self, res: ResidueData,
                                           bond_graph: Dict[int, Set[int]],
                                           start_id: int) -> Tuple[List[Group], int]:
        """第二层：官能团模式匹配识别带电基团。"""
        groups: List[Group] = []
        gid = start_id

        # 构建 atom_global_idx → atom 映射
        atom_map = {a.atom_global_idx: a for a in res.atoms}

        for atom in res.atoms:
            # 获取邻居
            neighbor_indices = bond_graph.get(atom.atom_global_idx, set())
            neighbors = [atom_map[idx] for idx in neighbor_indices if idx in atom_map]

            # 正电官能团
            if self._is_quartamine(atom, neighbors):
                g, gid = self._verify_and_build_group(
                    [atom], "charged_positive", res, gid, "quartamine")
                if g:
                    groups.append(g)
            elif self._is_tertamine(atom, neighbors):
                g, gid = self._verify_and_build_group(
                    [atom], "charged_positive", res, gid, "tertamine")
                if g:
                    groups.append(g)
            elif self._is_guanidine(atom, neighbors, bond_graph, atom_map):
                n_atoms = [n for n in neighbors if n.atom_element == 'N']
                g, gid = self._verify_and_build_group(
                    [atom] + n_atoms, "charged_positive", res, gid, "guanidine")
                if g:
                    groups.append(g)
            elif self._is_sulfonium(atom, neighbors):
                g, gid = self._verify_and_build_group(
                    [atom], "charged_positive", res, gid, "sulfonium")
                if g:
                    groups.append(g)

            # 负电官能团
            elif self._is_phosphate(atom, neighbors):
                o_atoms = [n for n in neighbors if n.atom_element == 'O']
                g, gid = self._verify_and_build_group(
                    [atom] + o_atoms, "charged_negative", res, gid, "phosphate")
                if g:
                    groups.append(g)
            elif self._is_sulfonicacid(atom, neighbors):
                o_atoms = [n for n in neighbors if n.atom_element == 'O']
                g, gid = self._verify_and_build_group(
                    [atom] + o_atoms, "charged_negative", res, gid, "sulfonicacid")
                if g:
                    groups.append(g)
            elif self._is_sulfate(atom, neighbors):
                o_atoms = [n for n in neighbors if n.atom_element == 'O']
                g, gid = self._verify_and_build_group(
                    [atom] + o_atoms, "charged_negative", res, gid, "sulfate")
                if g:
                    groups.append(g)
            elif self._is_carboxylate(atom, neighbors):
                o_atoms = [n for n in neighbors if n.atom_element == 'O']
                g, gid = self._verify_and_build_group(
                    [atom] + o_atoms, "charged_negative", res, gid, "carboxylate")
                if g:
                    groups.append(g)

        return groups, gid

    def _verify_and_build_group(self, atoms: List[AtomData],
                                group_type: str, res: ResidueData,
                                gid: int, func_group: str
                                ) -> Tuple[Group, int]:
        """第三层：电荷验证 + 构建 Group。"""
        # 计算净电荷
        net_charge = sum(a.atom_charge for a in atoms)

        # 验证电荷方向
        if group_type == "charged_positive" and net_charge <= CHARGE_THRESHOLD:
            return None, gid
        if group_type == "charged_negative" and net_charge >= -CHARGE_THRESHOLD:
            return None, gid

        # 构建 Group
        group = Group(
            group_id=gid,
            group_type=group_type,
            molecule=res.molecule_name,
            residue_name=res.residue_name,
            residue_id=res.residue_global_idx,
            atoms=atoms,
            metadata={
                "source": "functional_group",
                "func_group": func_group
            }
        )
        return group, gid + 1

    def _build_charged_group(self, atoms: List[AtomData],
                             group_type: str, res: ResidueData,
                             gid: int, source: str) -> Group:
        """构建带电基团 Group。"""
        return Group(
            group_id=gid,
            group_type=group_type,
            molecule=res.molecule_name,
            residue_name=res.residue_name,
            residue_id=res.residue_global_idx,
            atoms=atoms,
            metadata={
                "source": source,
                "func_group": None
            }
        )

    # 第二层官能团模式匹配函数（严格参照 PLIP is_functional_group）

    def _is_quartamine(self, atom: AtomData, neighbors: List[AtomData]) -> bool:
        """季铵：N 有 4 个邻居，且无 H 邻居。"""
        return (atom.atom_element == 'N'
                and len(neighbors) == 4
                and all(n.atom_element != 'H' for n in neighbors))

    def _is_tertamine(self, atom: AtomData, neighbors: List[AtomData]) -> bool:
        """叔胺：N 有 ≥3 个邻居（包括 H）。"""
        return (atom.atom_element == 'N'
                and len(neighbors) >= 3)

    def _is_guanidine(self, atom: AtomData, neighbors: List[AtomData],
                      bond_graph: Dict[int, Set[int]],
                      atom_map: Dict[int, AtomData]) -> bool:
        """胍基：C 有 3 个 N 邻居，且至少一个 N 只连该 C（非 H 邻居）。"""
        if atom.atom_element != 'C' or len(neighbors) != 3:
            return False
        if not all(n.atom_element == 'N' for n in neighbors):
            return False
        for n in neighbors:
            n_all = bond_graph.get(n.atom_global_idx, set())
            n_heavy = sum(1 for idx in n_all
                          if idx in atom_map and atom_map[idx].atom_element != 'H')
            if n_heavy == 1:
                return True
        return False

    def _is_sulfonium(self, atom: AtomData, neighbors: List[AtomData]) -> bool:
        """锍：S 有 3 个邻居，且无 H 邻居。"""
        return (atom.atom_element == 'S'
                and len(neighbors) == 3
                and all(n.atom_element != 'H' for n in neighbors))

    def _is_phosphate(self, atom: AtomData, neighbors: List[AtomData]) -> bool:
        """磷酸盐：P 的邻居全是 O。"""
        return (atom.atom_element == 'P'
                and all(n.atom_element == 'O' for n in neighbors))

    def _is_sulfonicacid(self, atom: AtomData, neighbors: List[AtomData]) -> bool:
        """磺酸：S 有 3 个 O 邻居。"""
        o_count = sum(1 for n in neighbors if n.atom_element == 'O')
        return (atom.atom_element == 'S' and o_count == 3)

    def _is_sulfate(self, atom: AtomData, neighbors: List[AtomData]) -> bool:
        """硫酸盐：S 有 4 个 O 邻居。"""
        o_count = sum(1 for n in neighbors if n.atom_element == 'O')
        return (atom.atom_element == 'S' and o_count == 4)

    def _is_carboxylate(self, atom: AtomData, neighbors: List[AtomData]) -> bool:
        """羧酸盐：C 有 2 个 O + 恰好 1 个 C。"""
        if atom.atom_element != 'C':
            return False
        o_count = sum(1 for n in neighbors if n.atom_element == 'O')
        c_count = sum(1 for n in neighbors if n.atom_element == 'C')
        return o_count == 2 and c_count == 1

    def _find_halogen_donors(self, res: ResidueData,
                       bond_graph: Dict[int, Set[int]],
                       start_id: int) -> Tuple[List[Group], int]:
        """检测卤键供体（卤素连接到碳）。"""
        groups: List[Group] = []
        gid = start_id
        local_to_global = {a.atom_idx_in_residue: a.atom_global_idx for a in res.atoms}

        for atom in res.atoms:
            if atom.atom_element not in ('F', 'Cl', 'Br', 'I'):
                continue
            if not self._is_bonded_to_element(atom, res, bond_graph, 'C'):
                continue
            groups.append(Group(
                group_id=gid, group_type="halogen_donor",
                molecule=res.molecule_name,
                residue_name=res.residue_name,
                residue_id=res.residue_global_idx,
                atoms=[atom]
            ))
            gid += 1

        return groups, gid

    def _find_halogen_acceptors(self, res: ResidueData,
                                bond_graph: Dict[int, Set[int]],
                                start_id: int) -> Tuple[List[Group], int]:
        """检测卤键受体（C/P/S 连接到 O/P/N/S）。"""
        groups: List[Group] = []
        gid = start_id
        ACCEPTOR_ELEMENTS = ('C', 'P', 'S')
        NEIGHBOR_ELEMENTS = ('O', 'P', 'N', 'S')

        for atom in res.atoms:
            if atom.atom_element not in ACCEPTOR_ELEMENTS:
                continue
            if not self._is_bonded_to_any_element(atom, res, bond_graph, NEIGHBOR_ELEMENTS):
                continue
            groups.append(Group(
                group_id=gid, group_type="halogen_acceptor",
                molecule=res.molecule_name,
                residue_name=res.residue_name,
                residue_id=res.residue_global_idx,
                atoms=[atom]
            ))
            gid += 1

        return groups, gid

    def _is_bonded_to_element(self, atom: AtomData, res: ResidueData,
                              bond_graph: Dict[int, Set[int]],
                              target_elem: str) -> bool:
        """检查原子是否连接到指定元素。"""
        g_idx = atom.atom_global_idx
        neighbors = bond_graph.get(g_idx, set())
        for n_idx in neighbors:
            for a in res.atoms:
                if a.atom_global_idx == n_idx and a.atom_element == target_elem:
                    return True
        return False

    def _is_bonded_to_any_element(self, atom: AtomData, res: ResidueData,
                                  bond_graph: Dict[int, Set[int]],
                                  target_elems: Tuple[str, ...]) -> bool:
        """检查原子是否连接到任一指定元素。"""
        g_idx = atom.atom_global_idx
        neighbors = bond_graph.get(g_idx, set())
        for n_idx in neighbors:
            for a in res.atoms:
                if a.atom_global_idx == n_idx and a.atom_element in target_elems:
                    return True
        return False

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
                    atoms=[atom]
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
                atoms=res.atoms
            ))
            gid += 1

        return groups, gid
