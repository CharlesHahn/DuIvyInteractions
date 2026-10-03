# -*- coding: utf-8 -*-
"""GROMOS 识别器单元测试（合成 SystemData，不依赖真实 tpr）。

原子/键/电荷数据取自本地 GROMACS 2019 `gromos54a7.ff/aminoacids.rtp`（与 53A6 一致），
验证 GromosFFGroupIdentifier 的基团识别逻辑：
- 芳香环分级判定（NR 强信号 / 全 C 弱信号白名单 / PRO 假环排除）
- 带电残基字典（LYSH/HISH 质子化 +1；LYS 中性不判）
- LYSH/HISH 子集去重（不产生双识别正电基团）
- 供体/受体/疏水/水
"""

import pytest

from DuIvyInteractions.core.datas import (
    SystemData, ResidueData, AtomData, BondData,
)
from DuIvyInteractions.group_identifiers.gromos_ff_identifier import (
    GromosFFGroupIdentifier,
)

# ============================================================
# 残基模板（GROMOS 54A7 aminoacids.rtp 实测）
# 每项: (atom_name, atom_type, element, charge)
# ============================================================

_PHE = [
    ("N", "N", "N", -0.310), ("H", "H", "H", 0.310),
    ("CA", "CH1", "C", 0.000), ("CB", "CH2", "C", 0.000),
    ("CG", "C", "C", 0.000), ("CD1", "C", "C", -0.140), ("HD1", "HC", "H", 0.140),
    ("CD2", "C", "C", -0.140), ("HD2", "HC", "H", 0.140),
    ("CE1", "C", "C", -0.140), ("HE1", "HC", "H", 0.140),
    ("CE2", "C", "C", -0.140), ("HE2", "HC", "H", 0.140),
    ("CZ", "C", "C", -0.140), ("HZ", "HC", "H", 0.140),
    ("C", "C", "C", 0.450), ("O", "O", "O", -0.450),
]
_PHE_BONDS = [
    ("N", "H"), ("N", "CA"), ("CA", "CB"), ("CA", "C"), ("CB", "CG"),
    ("CG", "CD1"), ("CG", "CD2"), ("CD1", "HD1"), ("CD1", "CE1"),
    ("CD2", "HD2"), ("CD2", "CE2"), ("CE1", "HE1"), ("CE1", "CZ"),
    ("CE2", "HE2"), ("CE2", "CZ"), ("CZ", "HZ"), ("C", "O"),
]

_TRP = [
    ("N", "N", "N", -0.310), ("H", "H", "H", 0.310),
    ("CA", "CH1", "C", 0.000), ("CB", "CH2", "C", 0.000),
    ("CG", "C", "C", -0.210), ("CD1", "C", "C", -0.140), ("HD1", "HC", "H", 0.140),
    ("CD2", "C", "C", 0.000), ("NE1", "NR", "N", -0.100), ("HE1", "H", "H", 0.310),
    ("CE2", "C", "C", 0.000), ("CE3", "C", "C", -0.140), ("HE3", "HC", "H", 0.140),
    ("CZ2", "C", "C", -0.140), ("HZ2", "HC", "H", 0.140),
    ("CZ3", "C", "C", -0.140), ("HZ3", "HC", "H", 0.140),
    ("CH2", "C", "C", -0.140), ("HH2", "HC", "H", 0.140),
    ("C", "C", "C", 0.450), ("O", "O", "O", -0.450),
]
_TRP_BONDS = [
    ("N", "H"), ("N", "CA"), ("CA", "CB"), ("CA", "C"), ("CB", "CG"),
    ("CG", "CD1"), ("CG", "CD2"), ("CD1", "HD1"), ("CD1", "NE1"),
    ("CD2", "CE2"), ("CD2", "CE3"), ("NE1", "HE1"), ("NE1", "CE2"),
    ("CE2", "CZ2"), ("CE3", "HE3"), ("CE3", "CZ3"), ("CZ2", "HZ2"),
    ("CZ2", "CH2"), ("CZ3", "HZ3"), ("CZ3", "CH2"), ("CH2", "HH2"),
    ("C", "O"),
]

_PRO = [
    ("N", "N", "N", 0.000), ("CA", "CH1", "C", 0.000), ("CB", "CH2r", "C", 0.000),
    ("CG", "CH2r", "C", 0.000), ("CD", "CH2r", "C", 0.000),
    ("C", "C", "C", 0.450), ("O", "O", "O", -0.450),
]
_PRO_BONDS = [
    ("N", "CA"), ("N", "CD"), ("CA", "CB"), ("CA", "C"),
    ("CB", "CG"), ("CG", "CD"), ("C", "O"),
]

_HISH = [
    ("N", "N", "N", -0.310), ("H", "H", "H", 0.310),
    ("CA", "CH1", "C", 0.000), ("CB", "CH2", "C", 0.000),
    ("CG", "C", "C", -0.050), ("ND1", "NR", "N", 0.380), ("HD1", "H", "H", 0.300),
    ("CD2", "C", "C", -0.100), ("HD2", "HC", "H", 0.100),
    ("CE1", "C", "C", -0.340), ("HE1", "HC", "H", 0.100),
    ("NE2", "NR", "N", 0.310), ("HE2", "H", "H", 0.300),
    ("C", "C", "C", 0.450), ("O", "O", "O", -0.450),
]
_HISH_BONDS = [
    ("N", "H"), ("N", "CA"), ("CA", "CB"), ("CA", "C"), ("CB", "CG"),
    ("CG", "ND1"), ("CG", "CD2"), ("ND1", "HD1"), ("ND1", "CE1"),
    ("CD2", "HD2"), ("CD2", "NE2"), ("CE1", "HE1"), ("CE1", "NE2"),
    ("NE2", "HE2"), ("C", "O"),
]

_HISA = [
    ("N", "N", "N", -0.310), ("H", "H", "H", 0.310),
    ("CA", "CH1", "C", 0.000), ("CB", "CH2", "C", 0.000),
    ("CG", "C", "C", 0.000), ("ND1", "NR", "N", -0.050), ("HD1", "H", "H", 0.310),
    ("CD2", "C", "C", 0.000), ("HD2", "HC", "H", 0.140),
    ("CE1", "C", "C", 0.000), ("HE1", "HC", "H", 0.140),
    ("NE2", "NR", "N", -0.540),
    ("C", "C", "C", 0.450), ("O", "O", "O", -0.450),
]
_HISA_BONDS = [
    ("N", "H"), ("N", "CA"), ("CA", "CB"), ("CA", "C"), ("CB", "CG"),
    ("CG", "ND1"), ("CG", "CD2"), ("ND1", "HD1"), ("ND1", "CE1"),
    ("CD2", "HD2"), ("CD2", "NE2"), ("CE1", "HE1"), ("CE1", "NE2"),
    ("C", "O"),
]

_LYSH = [
    ("N", "N", "N", -0.310), ("H", "H", "H", 0.310),
    ("CA", "CH1", "C", 0.000), ("CB", "CH2", "C", 0.000),
    ("CG", "CH2", "C", 0.000), ("CD", "CH2", "C", 0.000),
    ("CE", "CH2", "C", 0.127), ("NZ", "NL", "N", 0.129),
    ("HZ1", "H", "H", 0.248), ("HZ2", "H", "H", 0.248), ("HZ3", "H", "H", 0.248),
    ("C", "C", "C", 0.450), ("O", "O", "O", -0.450),
]
_LYSH_BONDS = [
    ("N", "H"), ("N", "CA"), ("CA", "CB"), ("CA", "C"), ("CB", "CG"),
    ("CG", "CD"), ("CD", "CE"), ("CE", "NZ"),
    ("NZ", "HZ1"), ("NZ", "HZ2"), ("NZ", "HZ3"), ("C", "O"),
]

_LYS = [
    ("N", "N", "N", -0.310), ("H", "H", "H", 0.310),
    ("CA", "CH1", "C", 0.000), ("CB", "CH2", "C", 0.000),
    ("CG", "CH2", "C", 0.000), ("CD", "CH2", "C", 0.000),
    ("CE", "CH2", "C", -0.240), ("NZ", "NT", "N", -0.640),
    ("HZ1", "H", "H", 0.440), ("HZ2", "H", "H", 0.440),
    ("C", "C", "C", 0.450), ("O", "O", "O", -0.450),
]
_LYS_BONDS = [
    ("N", "H"), ("N", "CA"), ("CA", "CB"), ("CA", "C"), ("CB", "CG"),
    ("CG", "CD"), ("CD", "CE"), ("CE", "NZ"),
    ("NZ", "HZ1"), ("NZ", "HZ2"), ("C", "O"),
]

_ARG = [
    ("N", "N", "N", -0.310), ("H", "H", "H", 0.310),
    ("CA", "CH1", "C", 0.000), ("CB", "CH2", "C", 0.000),
    ("CG", "CH2", "C", 0.000), ("CD", "CH2", "C", 0.090),
    ("NE", "NE", "N", -0.110), ("HE", "H", "H", 0.240),
    ("CZ", "C", "C", 0.340), ("NH1", "NZ", "N", -0.260),
    ("HH11", "H", "H", 0.240), ("HH12", "H", "H", 0.240),
    ("NH2", "NZ", "N", -0.260), ("HH21", "H", "H", 0.240), ("HH22", "H", "H", 0.240),
    ("C", "C", "C", 0.450), ("O", "O", "O", -0.450),
]
_ARG_BONDS = [
    ("N", "H"), ("N", "CA"), ("CA", "CB"), ("CA", "C"), ("CB", "CG"),
    ("CG", "CD"), ("CD", "NE"), ("NE", "HE"), ("NE", "CZ"),
    ("CZ", "NH1"), ("CZ", "NH2"), ("NH1", "HH11"), ("NH1", "HH12"),
    ("NH2", "HH21"), ("NH2", "HH22"), ("C", "O"),
]

_ASP = [
    ("N", "N", "N", -0.310), ("H", "H", "H", 0.310),
    ("CA", "CH1", "C", 0.000), ("CB", "CH2", "C", 0.000),
    ("CG", "C", "C", 0.270), ("OD1", "OM", "O", -0.635), ("OD2", "OM", "O", -0.635),
    ("C", "C", "C", 0.450), ("O", "O", "O", -0.450),
]
_ASP_BONDS = [
    ("N", "H"), ("N", "CA"), ("CA", "CB"), ("CA", "C"), ("CB", "CG"),
    ("CG", "OD1"), ("CG", "OD2"), ("C", "O"),
]

_GLU = [
    ("N", "N", "N", -0.310), ("H", "H", "H", 0.310),
    ("CA", "CH1", "C", 0.000), ("CB", "CH2", "C", 0.000),
    ("CG", "CH2", "C", 0.000), ("CD", "C", "C", 0.270),
    ("OE1", "OM", "O", -0.635), ("OE2", "OM", "O", -0.635),
    ("C", "C", "C", 0.450), ("O", "O", "O", -0.450),
]
_GLU_BONDS = [
    ("N", "H"), ("N", "CA"), ("CA", "CB"), ("CA", "C"), ("CB", "CG"),
    ("CG", "CD"), ("CD", "OE1"), ("CD", "OE2"), ("C", "O"),
]

_SOL = [
    ("OW", "OW", "O", -0.82), ("HW1", "H", "H", 0.41), ("HW2", "H", "H", 0.41),
]
_SOL_BONDS = [("OW", "HW1"), ("OW", "HW2")]

_TEMPLATES = {
    "PHE": (_PHE, _PHE_BONDS),
    "TRP": (_TRP, _TRP_BONDS),
    "PRO": (_PRO, _PRO_BONDS),
    "HISH": (_HISH, _HISH_BONDS),
    "HISA": (_HISA, _HISA_BONDS),
    "LYSH": (_LYSH, _LYSH_BONDS),
    "LYS": (_LYS, _LYS_BONDS),
    "ARG": (_ARG, _ARG_BONDS),
    "ASP": (_ASP, _ASP_BONDS),
    "GLU": (_GLU, _GLU_BONDS),
    "SOL": (_SOL, _SOL_BONDS),
}


def _make_residue(resname: str, resid: int, atom_offset: int) -> ResidueData:
    """按模板构造一个残基（原子全局索引连续分配）。"""
    atoms_spec, bonds_spec = _TEMPLATES[resname]
    atoms = [
        AtomData(
            atom_global_idx=atom_offset + i,
            atom_idx_in_residue=i,
            atom_name=name,
            atom_type=atype,
            atom_element=elem,
            atom_charge=charge,
            atom_mass=0.0,
        )
        for i, (name, atype, elem, charge) in enumerate(atoms_spec)
    ]
    name2local = {name: i for i, (name, _, _, _) in enumerate(atoms_spec)}
    bonds = [
        BondData(name2local[a], name2local[b], "bond")
        for a, b in bonds_spec
    ]
    return ResidueData(
        residue_name=resname,
        residue_global_idx=resid,
        residue_idx_in_molecule=resid + 1,
        molecule_name="testmol",
        atoms=atoms,
        bonds=bonds,
    )


def _make_system(residues: list) -> SystemData:
    """构造合成体系（残基列表 → SystemData），全局原子索引连续。"""
    sd_residues = []
    atom_offset = 0
    for resid, resname in enumerate(residues):
        res = _make_residue(resname, resid, atom_offset)
        sd_residues.append(res)
        atom_offset += len(res.atoms)
    return SystemData(
        system_name="synthetic_gromos",
        residues=sd_residues,
        inter_residue_bonds=[],
    )


@pytest.fixture(scope="module")
def identifier():
    return GromosFFGroupIdentifier()


@pytest.fixture(scope="module")
def groups(identifier):
    """识别合成体系（含多种代表残基）的全部基团。"""
    sd = _make_system(
        ["PHE", "TRP", "PRO", "HISH", "HISA", "LYSH", "LYS", "ARG", "ASP", "GLU", "SOL"]
    )
    return identifier.identify(sd)


def _groups_of(resname, groups):
    return [g for g in groups if g.residue_name == resname]


# ============================================================
# 芳香环识别
# ============================================================

class TestAromaticRings:

    def test_phe_single_ring(self, groups):
        rings = [g for g in _groups_of("PHE", groups)
                 if g.group_type == "aromatic_ring"]
        assert len(rings) == 1
        assert len(rings[0].atoms) == 6   # 全 C 环（弱信号白名单 PHE）

    def test_trp_two_rings(self, groups):
        rings = [g for g in _groups_of("TRP", groups)
                 if g.group_type == "aromatic_ring"]
        assert len(rings) == 2           # 六元全 C 环 + 五元含 NR 环
        assert sorted(len(r.atoms) for r in rings) == [5, 6]

    def test_trp_ring_contains_nr(self, groups):
        """TRP 五元环含 NE1（NR 强信号）→ 不依赖残基名白名单。"""
        rings = [g for g in _groups_of("TRP", groups)
                 if g.group_type == "aromatic_ring"]
        five = [r for r in rings if len(r.atoms) == 5][0]
        types = {a.atom_type for a in five.atoms}
        assert "NR" in types

    def test_hish_imidazole_ring(self, groups):
        rings = [g for g in _groups_of("HISH", groups)
                 if g.group_type == "aromatic_ring"]
        assert len(rings) == 1
        assert len(rings[0].atoms) == 5  # 咪唑 5 元环（含 NR 强信号）

    def test_hisa_imidazole_ring(self, groups):
        rings = [g for g in _groups_of("HISA", groups)
                 if g.group_type == "aromatic_ring"]
        assert len(rings) == 1
        assert len(rings[0].atoms) == 5

    def test_pro_ring_not_aromatic(self, groups):
        """PRO 五元环（N/CH2r）不得误判为芳香环。"""
        rings = [g for g in _groups_of("PRO", groups)
                 if g.group_type == "aromatic_ring"]
        assert rings == []


# ============================================================
# 带电基团
# ============================================================

class TestChargedGroups:

    def test_lysh_positive_single(self, groups):
        """LYSH 质子化 ε-铵 → 恰好 1 个正电基团（子集去重生效，无双识别）。"""
        pos = [g for g in _groups_of("LYSH", groups)
               if g.group_type == "charged_positive"]
        assert len(pos) == 1
        # 完整原子集（NZ + 3H），不是 tertamine 单原子 {NZ}
        assert pos[0].num_atoms == 4

    def test_lysh_net_charge_positive(self, groups):
        pos = [g for g in _groups_of("LYSH", groups)
               if g.group_type == "charged_positive"]
        assert pos[0].net_charge > 0.1

    def test_lys_neutral_not_in_positive_dict(self, groups):
        """中性 LYS（NH2）不得产出正电基团。"""
        pos = [g for g in _groups_of("LYS", groups)
               if g.group_type == "charged_positive"]
        assert pos == []

    def test_hish_positive_single(self, groups):
        """HISH → 恰好 1 个正电基团（ND1/NE2 双 tertamine 被去重）。"""
        pos = [g for g in _groups_of("HISH", groups)
               if g.group_type == "charged_positive"]
        assert len(pos) == 1

    def test_arg_positive_single(self, groups):
        pos = [g for g in _groups_of("ARG", groups)
               if g.group_type == "charged_positive"]
        assert len(pos) == 1
        assert pos[0].net_charge > 0.1

    def test_asp_glu_negative(self, groups):
        for resname in ("ASP", "GLU"):
            neg = [g for g in _groups_of(resname, groups)
                   if g.group_type == "charged_negative"]
            assert len(neg) == 1, resname
            assert neg[0].net_charge < -0.1, resname


# ============================================================
# 供体 / 受体 / 水 / 疏水
# ============================================================

class TestDonorAcceptor:

    def test_donors_present(self, groups):
        """主链 N-H 供体应被识别（q(H)>0）。"""
        phe_donors = [g for g in _groups_of("PHE", groups)
                      if g.group_type == "H_donor"]
        assert len(phe_donors) == 1
        assert phe_donors[0].atoms[0].atom_element == "N"
        assert phe_donors[0].atoms[1].atom_element == "H"

    def test_tyr_like_o_h_donor(self, groups):
        """TRP 无 O-H；TYR 未建模——用 TRP NE1-H 验证芳香 N-H 供体。"""
        trp_donors = [g for g in _groups_of("TRP", groups)
                      if g.group_type == "H_donor"]
        # 主链 N-H + 吲哚 NE1-H
        n_donors = [d for d in trp_donors
                    if d.atoms[0].atom_name in ("N", "NE1")]
        assert len(n_donors) == 2

    def test_acceptors_negative_charged(self, groups):
        """受体应带负电荷（q<0）；ASP 含 OD1/OD2 羧基氧受体。"""
        asp_acc = [g for g in _groups_of("ASP", groups)
                   if g.group_type == "H_acceptor"]
        # OD1、OD2 必须为受体；主链 O、主链 N（q<0）也会是受体
        acc_names = {a.atoms[0].atom_name for a in asp_acc}
        assert "OD1" in acc_names and "OD2" in acc_names
        assert all(a.atoms[0].atom_charge < 0 for a in asp_acc)

    def test_water_group(self, groups):
        waters = [g for g in _groups_of("SOL", groups)
                  if g.group_type == "water"]
        assert len(waters) == 1
        assert waters[0].num_atoms == 3   # OW, HW1, HW2


class TestHydrophobic:

    def test_aliphatic_carbon_hydrophobic(self, groups):
        """LYS 侧链 CH2（无极性邻居）→ 疏水。"""
        lys_hydro = [g for g in _groups_of("LYS", groups)
                     if g.group_type == "hydrophobic"]
        assert len(lys_hydro) >= 1        # CB/CG/CD/CE

    def test_polar_adjacent_carbon_not_hydrophobic(self, groups):
        """ASP CG 连 O → 不得判为疏水。"""
        asp_hydro = [g for g in _groups_of("ASP", groups)
                     if g.group_type == "hydrophobic"]
        asp_hydro_names = {g.atoms[0].atom_name for g in asp_hydro}
        assert "CG" not in asp_hydro_names

    def test_lysh_nz_not_hydrophobic(self, groups):
        """LYSH NZ 连 N → 不得判为疏水。"""
        lysh_hydro = [g for g in _groups_of("LYSH", groups)
                      if g.group_type == "hydrophobic"]
        lysh_hydro_names = {g.atoms[0].atom_name for g in lysh_hydro}
        assert "NZ" not in lysh_hydro_names


# ============================================================
# 注册表 / 接口
# ============================================================

class TestRegistry:

    def test_name(self, identifier):
        assert identifier.name == "gromos_ff"

    def test_registry_contains_gromos(self):
        from DuIvyInteractions.group_identifiers import IDENTIFIER_CLASSES
        assert "gromos" in IDENTIFIER_CLASSES
        assert IDENTIFIER_CLASSES["gromos"] is GromosFFGroupIdentifier

    def test_pipeline_make_identifier(self):
        from DuIvyInteractions.pipeline import Pipeline
        assert isinstance(Pipeline("gromos")._make_identifier(),
                          GromosFFGroupIdentifier)

    def test_filter_groups_excludes_water_for_saltbridge(self):
        """盐桥检测器过滤：水分子应被排除。"""
        from DuIvyInteractions.pipeline import Pipeline
        det = Pipeline("gromos", "two_pass")._make_detector("salt_bridge")
        sd = _make_system(["ARG", "ASP", "SOL"])
        gs = GromosFFGroupIdentifier().identify(sd)
        filtered = Pipeline._filter_groups(gs, det)
        types = {g.group_type for g in filtered}
        assert types == {"charged_positive", "charged_negative"}
        assert all(g.residue_name != "SOL" for g in filtered)