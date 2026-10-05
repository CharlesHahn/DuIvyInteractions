# -*- coding: utf-8 -*-
"""GROMOS 力场基团识别器（53A6 / 54A7 家族）。

与 AmberFFGroupIdentifier 的差异（详见 doc/gromos_identifier_design.md §3.1）：
- A 类（力场无关）：供体、卤键、金属、水、金属配位、配体带电 → 继承复用
- B 类（仅换表）：H 键受体 → 换 GROMOS_ACCEPTOR_TYPES
- C 类（改字典）：蛋白带电 → 换 GROMOS 残基字典（LYS 无 HZ3、只含 5 残基）
- D 类（重写逻辑）：芳香环（分级判定）、疏水（反列举）→ 覆写
"""

from typing import List, Tuple, Dict, Set

from ..core.interfaces import GroupIdentifier
from ..core.datas import Group, SystemData, ResidueData, AtomData
from .amber_ff_identifier import (
    AmberFFGroupIdentifier,
    METAL_IONS,
    WATER_RESIDUES,
)


# ============================================================
# GROMOS 类型特征表（本地 gromos54a7.ff 实测，53A6/54A7 通用）
# ============================================================

# H 键受体类型（atomtypes.atp：O/OM/OA/OE/OW、N/NT/NL/NR/NZ/NE、S、F/CL/BR）
GROMOS_ACCEPTOR_TYPES = frozenset({
    "O", "OM", "OA", "OE", "OW",       # 氧：羰基/羧基/羟基/醚/水
    "N", "NT", "NL", "NR", "NZ", "NE", # 氮：肽胺/终端/芳香/胍基
    "S",                                # 硫
    "F", "CL", "BR",                   # 卤素
})

# 蛋白正电残基（GROMOS rtp + hdb 实测：净电荷 +1 的残基）
# 质子化变体约定：GROMOS 将 LYS 拆为 中性 LYS(NH2, 电荷 0) / 质子化 LYSH(NH3+, +1)，
# 后者才带正电；对应的 Amber 残基为 LYS（Amber 默认恒为 NH3+）。
GROMOS_POSITIVE_RESIDUES = {
    "ARG": ["CZ", "NE", "NH1", "NH2", "HE", "HH11", "HH12", "HH21", "HH22"],
    "LYSH": ["NZ", "HZ1", "HZ2", "HZ3"],   # 质子化 LYS（ε-铵 NH3+，净电荷 +1）
    "HISH": ["ND1", "NE2", "HD1", "HE2"],  # 质子化 His（咪唑双质子化，净电荷 +1）
}

# 蛋白负电残基（GROMOS rtp 实测：净电荷 -1）
GROMOS_NEGATIVE_RESIDUES = {
    "ASP": ["CG", "OD1", "OD2"],
    "GLU": ["CD", "OE1", "OE2"],
}
# 注：LYS（中性 ε-胺 NH2）/ ARGN / ASPH / GLUH / HISA/B/1/2 等中性变体不进字典；
# CYS（-0.5 硫醇根）强度不足、PLIP 通行不含，排除。

# 芳香环分级判定用集合
GROMOS_RING_TYPES = frozenset({"C", "CR1", "NR"})        # GROMOS 环碳/环氮全集
GROMOS_AROMATIC_STRONG = frozenset({"NR"})               # NR：芳香氮，强信号
GROMOS_AROMATIC_RESIDUES = frozenset({                    # 弱信号：全 C 环的残基名白名单
    "PHE", "TYR", "TRP", "HISA", "HISB", "HISH",
    "HIS1", "HIS2",
})

# 疏水反列举
GROMOS_HYDROPHOBIC_TYPES = frozenset({
    "C", "CH0", "CH1", "CH2", "CH3", "CH4", "CH2r",
    # 不含 CH3p（极性，胆碱 N+）、CR1（芳香/烯 sp2）
})
GROMOS_HYDROPHOBIC_EXCLUDED = frozenset({"O", "N", "S"})


class GromosFFGroupIdentifier(AmberFFGroupIdentifier):
    """GROMOS 力场基团识别器。"""

    @property
    def name(self) -> str:
        return "gromos_ff"

    # ============================================================
    # B 类：H 键受体（仅换类型表，逻辑同 amber）
    # ============================================================

    def _find_acceptors(self, res: ResidueData,
                        start_id: int,
                        exclude_atoms: Set[int] = None) -> Tuple[List[Group], int]:
        """检测 H 键受体（类型 ∈ GROMOS_ACCEPTOR_TYPES + q<0，排除无孤对原子）。"""
        exclude_atoms = exclude_atoms or set()
        groups: List[Group] = []
        gid = start_id

        for atom in res.atoms:
            if atom.atom_type in GROMOS_ACCEPTOR_TYPES and atom.atom_charge < 0:
                if atom.atom_global_idx in exclude_atoms:
                    continue
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
    # C 类：蛋白带电（改字典，逻辑同 amber）
    # ============================================================

    def _deduplicate_charged(self, groups: List[Group]) -> List[Group]:
        """带电基团去重：父类（原子集完全相同）+ GROMOS 子集去重。

        GROMOS 的质子化残基正电 N 原子本身带正电荷（LYSH NZ +0.129、
        HISH ND1 +0.38 / NE2 +0.31），tertamine 官能团层会对每个 N 单原子
        额外产出一个正电基团（如 {NZ}、{ND1}、{NE2}），与字典层的完整基团
        （{NZ,HZ1,HZ2,HZ3}、{ND1,NE2,HD1,HE2}）构成「子集」关系。

        这里删除被同类型更完整基团严格包含的基团，只保留最完整者；
        Amber 正电 N 为负电荷、无此问题，故仅 GROMOS 需要此覆写。
        """
        # 1. 父类：原子集完全相同的去重（优先保留 residue_name 来源）
        result = super()._deduplicate_charged(groups)

        # 2. 子集去重：同类型下，若 A 的原子集严格包含于 B，删 A 保 B
        kept: List[Group] = []
        for i, g in enumerate(result):
            gi = set(g.atom_indices)
            is_duplicate = False
            for j, other in enumerate(result):
                if i == j or g.group_type != other.group_type:
                    continue
                if gi < set(other.atom_indices):   # 严格真子集
                    is_duplicate = True
                    break
            if not is_duplicate:
                kept.append(g)
        return kept

    def _identify_protein_charged(self, res: ResidueData,
                                  start_id: int) -> Tuple[List[Group], int]:
        """第一层：GROMOS 残基名字典识别带电残基。"""
        groups: List[Group] = []
        gid = start_id

        if res.residue_name in GROMOS_POSITIVE_RESIDUES:
            atom_names = GROMOS_POSITIVE_RESIDUES[res.residue_name]
            atoms = [a for a in res.atoms if a.atom_name in atom_names]
            if atoms:
                groups.append(self._build_charged_group(
                    atoms, "charged_positive", res, gid, "residue_name"))
                gid += 1

        if res.residue_name in GROMOS_NEGATIVE_RESIDUES:
            atom_names = GROMOS_NEGATIVE_RESIDUES[res.residue_name]
            atoms = [a for a in res.atoms if a.atom_name in atom_names]
            if atoms:
                groups.append(self._build_charged_group(
                    atoms, "charged_negative", res, gid, "residue_name"))
                gid += 1

        return groups, gid

    # ============================================================
    # D 类：芳香环分级判定（重写逻辑）
    # ============================================================

    def _filter_aromatic_rings(self, rings: List[List[int]],
                               res: ResidueData) -> List[List[int]]:
        """GROMOS 芳香环分级判定。

        强信号（NR）直接判芳香；全 C / 含 CR1 环回退残基名白名单。
        CR1 也用于链状烯烃（MEBMT），但环检测先行保证此处是闭合环，
        环内 CR1 + 环成分过滤即可接受。
        """
        aromatic_rings = []

        for ring_atoms in rings:
            types = [res.atoms[i].atom_type for i in ring_atoms]

            # 条件 1：环原子类型必须是 GROMOS 环碳/环氮全集
            if not all(t in GROMOS_RING_TYPES for t in types):
                continue

            # 条件 2：强信号——环内含 NR（芳香氮）
            if any(t in GROMOS_AROMATIC_STRONG for t in types):
                aromatic_rings.append(ring_atoms)
                continue

            # 条件 3：弱信号——全 C 环（或含 CR1 无 NR），仅标准蛋白芳香残基
            if res.residue_name in GROMOS_AROMATIC_RESIDUES:
                aromatic_rings.append(ring_atoms)

        return aromatic_rings

    # ============================================================
    # D 类：疏水反列举（重写逻辑）
    # ============================================================

    def _find_hydrophobic(self, res: ResidueData,
                          bond_graph: Dict[int, Set[int]],
                          start_id: int,
                          atom_map: Dict[int, AtomData]) -> Tuple[List[Group], int]:
        """检测疏水原子（GROMOS UA 碳 + 邻居不含 O/N/S）。

        与 amber 的"C 且邻居∈{C,H}"不同：GROMOS 无碳上 H 邻居，
        改为反列举排除极性取代（邻居 ∉ {O,N,S}）。
        """
        groups: List[Group] = []
        gid = start_id

        for atom in res.atoms:
            # 条件 1：类型必须是 GROMOS 非极性碳（UA 碳）
            if atom.atom_type not in GROMOS_HYDROPHOBIC_TYPES:
                continue

            # 条件 2：所有邻居不含极性元素 O/N/S
            neighbor_indices = bond_graph.get(atom.atom_global_idx, set())
            neighbors = [atom_map[idx] for idx in neighbor_indices
                         if idx in atom_map]
            if any(n.atom_element in GROMOS_HYDROPHOBIC_EXCLUDED
                   for n in neighbors):
                continue

            groups.append(Group(
                group_id=gid, group_type="hydrophobic",
                molecule=res.molecule_name,
                residue_name=res.residue_name,
                residue_id=res.residue_global_idx,
                atoms=[atom]
            ))
            gid += 1

        return groups, gid