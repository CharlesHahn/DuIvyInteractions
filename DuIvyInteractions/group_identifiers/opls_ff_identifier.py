# -*- coding: utf-8 -*-
"""OPLS 力场基团识别器（OPLS-AA/L，GROMACS oplsaa.ff）。

与 AmberFFGroupIdentifier 的差异（详见 doc/opls_identifier_design.md §3.1）：
- A 类（力场无关）：供体、卤键、金属、金属配位、配体带电、疏水 → 继承复用
- B 类（仅换表）：受体、芳香环类型表 → 换 OPLS 特征表（opls_XXX 编号）
- C 类（改字典）：蛋白带电 → 换 OPLS 残基字典（ARG/LYSH/HISH + ASP/GLU）
- D 类（重写逻辑）：N 端/侧链铵（LYSH）识别 —— OPLS 铵氮带负电荷，
  复用 CHARMM 的 _is_terminal_ammonium 结构判据（N 邻居∈{C,H} + 排除胍基）
"""

from typing import List, Tuple, Dict, Set

from ..core.datas import Group, ResidueData, AtomData
from .amber_ff_identifier import AmberFFGroupIdentifier


# ============================================================
# OPLS 类型特征表（蛋白 rtp 实测 + ffnonbonded 符号映射，
# 版本 = OPLS-AA/L 2001，符号语义见 Towhee/Robertson）
# ============================================================

# H 键受体类型（OPLS 蛋白实际使用的 O/N/S + 卤素，q<0 时作受体）
OPLS_ACCEPTOR_TYPES = frozenset({
    # 氧（O/O2/OH/O_3）
    "opls_236", "opls_272", "opls_154", "opls_167", "opls_268", "opls_269",
    # 氮（保留受体类型，依据 OPLS-AA rtp 实测 + 文献）
    # 剔除 opls_238(主链 N)/opls_237(侧链酰胺 N)/opls_287(LYSH 铵)/
    #      opls_300,303(Arg 胍基)/opls_503(带 H 吡咯: HISD ND1, TRP NE1)/
    #      opls_512(双质子化 His N)/opls_749,750,751(ARGN 中性胍)
    #     —— 均为孤被共振占用/无孤对，非受体（Eildal 2013; 教科书）
    # 保留 opls_239(Pro N，无 H，受体，Deepak 2016)、opls_511(无 H 吡啶型 His N，受体)、
    #      opls_900(LYS 中性胺 NZ，受体，Luisi 1998; Baik 2003)
    "opls_239", "opls_511", "opls_900",
    # 硫（S/SH）
    "opls_202", "opls_200",
    # 卤素（离子 F-/Cl-/Br-/I- + 有机卤素；卤素可作 H 键受体，Lin 2017；
    # 与 amber ACCEPTOR_TYPES 含 f/cl/br/i 同理念）
    "opls_400", "opls_401", "opls_402", "opls_403",   # 离子卤素
    "opls_123", "opls_151", "opls_164", "opls_226", "opls_264", "opls_709",
    "opls_719", "opls_721", "opls_722", "opls_726", "opls_728",
    "opls_730", "opls_732", "opls_786", "opls_956", "opls_965",  # 有机卤素
})

# 芳香环强信号类型（CA/C*/CB/CN/CW/CV/CR/CX + NA/NB）
OPLS_STRONG_AROMATIC = frozenset({
    # CA 系（苯环/芳香碳）
    "opls_145", "opls_166", "opls_302", "opls_752",
    # 5 元环：C*/CB/CN/CW/CV/CR/CX
    "opls_500", "opls_501", "opls_502", "opls_506", "opls_507",
    "opls_508", "opls_509", "opls_510", "opls_514",
    # 芳香氮：NA/NB
    "opls_503", "opls_511", "opls_512",
})

# 芳香环兼容类型（非强信号但可参与共轭环：环内羰基碳）
OPLS_COMPATIBLE_TYPES = frozenset({
    "opls_235",    # C 主链羰基碳（环内时）
    "opls_267",    # C 变体
    "opls_271",    # C_3 羧酸碳（环内时）
})

# 水残基名（OPLS: HOH/SPC、HO4/TIP4P、HO5/TIP5P + 兼容 SOL/WAT）
OPLS_WATER_RESIDUES = frozenset({"HOH", "HO4", "HO5", "SOL", "WAT"})

# 蛋白正电残基（rtp 净电荷 +1；OPLS LYS 中性、LYSH 质子化）
OPLS_POSITIVE_RESIDUES = {
    "ARG": ["CZ", "NE", "NH1", "NH2", "HE", "HH11", "HH12", "HH21", "HH22"],
    "LYSH": ["NZ", "HZ1", "HZ2", "HZ3"],   # 质子化 ε-铵（NH3+）
    "HISH": ["ND1", "NE2", "HD1", "HE2", "CG", "CE1", "CD2"],  # 质子化 His（含咪唑碳）
}

# 蛋白负电残基（rtp 净电荷 -1）
OPLS_NEGATIVE_RESIDUES = {
    "ASP": ["CG", "OD1", "OD2"],
    "GLU": ["CD", "OE1", "OE2"],
}
# 注：LYS（中性）/ARGN/ASPH/GLUH/HISD/HISE 不进字典。


class OplsFFGroupIdentifier(AmberFFGroupIdentifier):
    """OPLS 力场基团识别器（OPLS-AA/L）。"""

    # 水残基名（OPLS: HOH/SPC、HO4/TIP4P、HO5/TIP5P + SOL/WAT）
    WATER_RESIDUES = OPLS_WATER_RESIDUES

    @property
    def name(self) -> str:
        return "opls_ff"

    # ============================================================
    # B 类：H 键受体（仅换类型表）
    # ============================================================

    def _find_acceptors(self, res: ResidueData,
                        start_id: int,
                        exclude_atoms: Set[int] = None) -> Tuple[List[Group], int]:
        """检测 H 键受体（类型 ∈ OPLS_ACCEPTOR_TYPES + q<0，排除无孤对原子）。"""
        exclude_atoms = exclude_atoms or set()
        groups: List[Group] = []
        gid = start_id

        for atom in res.atoms:
            if atom.atom_type in OPLS_ACCEPTOR_TYPES and atom.atom_charge < 0:
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
    # B 类：芳香环（仅换类型表，两条件逻辑同 amber）
    # ============================================================

    def _filter_aromatic_rings(self, rings: List[List[int]],
                               res: ResidueData) -> List[List[int]]:
        """OPLS 芳香环两条件过滤（≥n-1 强信号 + 其余兼容）。"""
        strong = OPLS_STRONG_AROMATIC
        compatible = OPLS_COMPATIBLE_TYPES
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
        """第一层：OPLS 残基名字典识别带电残基。"""
        groups: List[Group] = []
        gid = start_id

        for dict_, gtype in ((OPLS_POSITIVE_RESIDUES, "charged_positive"),
                             (OPLS_NEGATIVE_RESIDUES, "charged_negative")):
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
    # D 类：N 端/侧链铵识别（自实现结构判据）
    # ============================================================

    @staticmethod
    def _is_terminal_ammonium(atom: AtomData, neighbors: List[AtomData],
                              bond_graph: Dict[int, Set[int]],
                              atom_map: Dict[int, AtomData]) -> bool:
        """末端铵结构判据：N 邻居∈{C,H} 且重原子邻居不连 ≥2 个 N。

        用于 OPLS N 端 NH3+/LYSH 侧链铵识别：铵氮带负电荷，
        须 N+键连 H 一起验证，但仅限"末端铵"结构（胍基/酰胺除外）。
        力场无关的纯结构规则（与 CHARMM 实现逻辑相同，独立维护）。
        """
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

    def _identify_functional_group_charged(self, res: ResidueData,
                                           bond_graph: Dict[int, Set[int]],
                                           start_id: int) -> Tuple[List[Group], int]:
        """第二层：官能团模式匹配识别带电基团（OPLS 铵氮负电荷修正）。

        对满足"末端铵结构"的氮（N 邻居∈{C,H} 且不属于胍基），
        用 N + 键连 H 的净电荷验证（LYSH NZ/N 端）。
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

    # （_find_water / _find_metal_binding 继承父类，经 self.WATER_RESIDUES 取用）