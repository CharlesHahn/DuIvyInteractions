# -*- coding: utf-8 -*-
"""OPLS 识别器单元测试（合成 SystemData，不依赖真实 tpr）。

原子/键/电荷数据取自本地 `oplsaa.ff`（OPLS-AA/L 2001，aminoacids.rtp +
n.tdb/c.tdb），验证 OplsFFGroupIdentifier：
- N 端 NH3+ / LYSH 侧链铵（结构判据 N+H 验证）
- 芳香环（opls_145 CA 系 + C*/NA/NB）
- 带电残基字典（ARG/LYSH/HISH + ASP/GLU）
- 供体/受体/水/疏水
"""

import pytest

from DuIvyInteractions.core.datas import (
    SystemData, ResidueData, AtomData, BondData,
)
from DuIvyInteractions.group_identifiers.opls_ff_identifier import (
    OplsFFGroupIdentifier,
)

# ============================================================
# 残基模板（OPLS-AA/L oplsaa.ff 实测）
# ============================================================

_GLY_NTER = [
    ("N", "opls_287", "N", -0.30), ("H1", "opls_290", "H", 0.33),
    ("H2", "opls_290", "H", 0.33), ("H3", "opls_290", "H", 0.33),
    ("CA", "opls_292B", "C", 0.19), ("HA2", "opls_140", "H", 0.06),
    ("C", "opls_235", "C", 0.50), ("O", "opls_236", "O", -0.50),
]
_GLY_NTER_BONDS = [
    ("N", "H1"), ("N", "H2"), ("N", "H3"), ("N", "CA"),
    ("CA", "C"), ("C", "O"), ("CA", "HA2"),
]

_PHE = [
    ("N", "opls_238", "N", -0.50), ("H", "opls_241", "H", 0.30),
    ("CA", "opls_224B", "C", 0.14), ("HA", "opls_140", "H", 0.06),
    ("CB", "opls_149", "C", -0.005), ("HB1", "opls_140", "H", 0.06), ("HB2", "opls_140", "H", 0.06),
    ("CG", "opls_145", "C", -0.115), ("CD1", "opls_145", "C", -0.115), ("HD1", "opls_146", "H", 0.115),
    ("CD2", "opls_145", "C", -0.115), ("HD2", "opls_146", "H", 0.115),
    ("CE1", "opls_145", "C", -0.115), ("HE1", "opls_146", "H", 0.115),
    ("CE2", "opls_145", "C", -0.115), ("HE2", "opls_146", "H", 0.115),
    ("CZ", "opls_145", "C", -0.115), ("HZ", "opls_146", "H", 0.115),
    ("C", "opls_235", "C", 0.50), ("O", "opls_236", "O", -0.50),
]
_PHE_BONDS = [
    ("N", "H"), ("N", "CA"), ("CA", "C"), ("C", "O"),
    ("CA", "CB"), ("CB", "CG"), ("CA", "HA"), ("CB", "HB1"), ("CB", "HB2"),
    ("CG", "CD1"), ("CG", "CD2"), ("CD1", "HD1"), ("CD1", "CE1"),
    ("CD2", "HD2"), ("CD2", "CE2"), ("CE1", "HE1"), ("CE1", "CZ"),
    ("CE2", "HE2"), ("CE2", "CZ"), ("CZ", "HZ"),
]

_TRP = [
    ("N", "opls_238", "N", -0.50), ("H", "opls_241", "H", 0.30),
    ("CA", "opls_224B", "C", 0.14), ("HA", "opls_140", "H", 0.06),
    ("CB", "opls_149", "C", -0.005), ("HB1", "opls_140", "H", 0.06), ("HB2", "opls_140", "H", 0.06),
    ("CG", "opls_500", "C", 0.075), ("CD1", "opls_514", "C", -0.115), ("HD1", "opls_146", "H", 0.115),
    ("NE1", "opls_503", "N", -0.57), ("HE1", "opls_504", "H", 0.41),
    ("CE2", "opls_502", "C", 0.13), ("CD2", "opls_501", "C", -0.055),
    ("CE3", "opls_145", "C", -0.115), ("HE3", "opls_146", "H", 0.115),
    ("CZ3", "opls_145", "C", -0.115), ("HZ3", "opls_146", "H", 0.115),
    ("CZ2", "opls_145", "C", -0.115), ("HZ2", "opls_146", "H", 0.115),
    ("CH2", "opls_145", "C", -0.115), ("HH2", "opls_146", "H", 0.115),
    ("C", "opls_235", "C", 0.50), ("O", "opls_236", "O", -0.50),
]
_TRP_BONDS = [
    ("N", "H"), ("N", "CA"), ("CA", "C"), ("C", "O"),
    ("CA", "CB"), ("CB", "CG"), ("CA", "HA"), ("CB", "HB1"), ("CB", "HB2"),
    ("CG", "CD1"), ("CG", "CD2"), ("CD1", "HD1"), ("CD1", "NE1"),
    ("NE1", "HE1"), ("NE1", "CE2"), ("CE2", "CD2"), ("CE2", "CZ2"),
    ("CD2", "CE3"), ("CE3", "HE3"), ("CE3", "CZ3"), ("CZ3", "HZ3"),
    ("CZ3", "CH2"), ("CZ2", "HZ2"), ("CZ2", "CH2"), ("CH2", "HH2"),
]

_HISH = [
    ("N", "opls_238", "N", -0.50), ("H", "opls_241", "H", 0.30),
    ("CA", "opls_224B", "C", 0.14), ("HA", "opls_140", "H", 0.06),
    ("CB", "opls_505", "C", -0.005), ("HB1", "opls_140", "H", 0.06), ("HB2", "opls_140", "H", 0.06),
    ("CG", "opls_510", "C", 0.215), ("ND1", "opls_512", "N", -0.54), ("HD1", "opls_513", "H", 0.46),
    ("CD2", "opls_510", "C", 0.215), ("HD2", "opls_146", "H", 0.115),
    ("CE1", "opls_509", "C", 0.385), ("HE1", "opls_146", "H", 0.115),
    ("NE2", "opls_512", "N", -0.54), ("HE2", "opls_513", "H", 0.46),
    ("C", "opls_235", "C", 0.50), ("O", "opls_236", "O", -0.50),
]
_HISH_BONDS = [
    ("N", "H"), ("N", "CA"), ("CA", "C"), ("C", "O"),
    ("CA", "CB"), ("CB", "CG"), ("CA", "HA"), ("CB", "HB1"), ("CB", "HB2"),
    ("CG", "ND1"), ("CG", "CD2"), ("ND1", "HD1"), ("ND1", "CE1"),
    ("CD2", "HD2"), ("CD2", "NE2"), ("CE1", "HE1"), ("CE1", "NE2"),
    ("NE2", "HE2"),
]

_LYSH = [
    ("N", "opls_238", "N", -0.50), ("H", "opls_241", "H", 0.30),
    ("CA", "opls_224B", "C", 0.14), ("HA", "opls_140", "H", 0.06),
    ("CB", "opls_136", "C", -0.06), ("HB1", "opls_140", "H", 0.06), ("HB2", "opls_140", "H", 0.06),
    ("CG", "opls_136", "C", -0.06), ("HG1", "opls_140", "H", 0.06), ("HG2", "opls_140", "H", 0.06),
    ("CD", "opls_136", "C", -0.06), ("HD1", "opls_140", "H", 0.06), ("HD2", "opls_140", "H", 0.06),
    ("CE", "opls_292", "C", 0.19), ("HE1", "opls_140", "H", 0.06), ("HE2", "opls_140", "H", 0.06),
    ("NZ", "opls_287", "N", -0.30), ("HZ1", "opls_290", "H", 0.33),
    ("HZ2", "opls_290", "H", 0.33), ("HZ3", "opls_290", "H", 0.33),
    ("C", "opls_235", "C", 0.50), ("O", "opls_236", "O", -0.50),
]
_LYSH_BONDS = [
    ("N", "H"), ("N", "CA"), ("CA", "C"), ("C", "O"),
    ("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "CE"), ("CE", "NZ"),
    ("NZ", "HZ1"), ("NZ", "HZ2"), ("NZ", "HZ3"),
    ("CA", "HA"), ("CB", "HB1"), ("CB", "HB2"), ("CG", "HG1"), ("CG", "HG2"),
    ("CD", "HD1"), ("CD", "HD2"), ("CE", "HE1"), ("CE", "HE2"),
]

_ARG = [
    ("N", "opls_238", "N", -0.50), ("H", "opls_241", "H", 0.30),
    ("CA", "opls_224B", "C", 0.14), ("HA", "opls_140", "H", 0.06),
    ("CB", "opls_308", "C", -0.05), ("HB1", "opls_140", "H", 0.06), ("HB2", "opls_140", "H", 0.06),
    ("CG", "opls_308", "C", -0.05), ("HG1", "opls_140", "H", 0.06), ("HG2", "opls_140", "H", 0.06),
    ("CD", "opls_307", "C", 0.19), ("HD1", "opls_140", "H", 0.06), ("HD2", "opls_140", "H", 0.06),
    ("NE", "opls_303", "N", -0.70), ("HE", "opls_304", "H", 0.44),
    ("CZ", "opls_302", "C", 0.64), ("NH1", "opls_300", "N", -0.80),
    ("HH11", "opls_301", "H", 0.46), ("HH12", "opls_301", "H", 0.46),
    ("NH2", "opls_300", "N", -0.80), ("HH21", "opls_301", "H", 0.46),
    ("HH22", "opls_301", "H", 0.46),
    ("C", "opls_235", "C", 0.50), ("O", "opls_236", "O", -0.50),
]
_ARG_BONDS = [
    ("N", "H"), ("N", "CA"), ("CA", "C"), ("C", "O"),
    ("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "NE"),
    ("NE", "HE"), ("NE", "CZ"), ("CZ", "NH1"), ("CZ", "NH2"),
    ("NH1", "HH11"), ("NH1", "HH12"), ("NH2", "HH21"), ("NH2", "HH22"),
    ("CA", "HA"), ("CB", "HB1"), ("CB", "HB2"), ("CG", "HG1"), ("CG", "HG2"),
    ("CD", "HD1"), ("CD", "HD2"),
]

_ASP = [
    ("N", "opls_238", "N", -0.50), ("H", "opls_241", "H", 0.30),
    ("CA", "opls_224B", "C", 0.14), ("HA", "opls_140", "H", 0.06),
    ("CB", "opls_135", "C", -0.03), ("HB1", "opls_140", "H", 0.06), ("HB2", "opls_140", "H", 0.06),
    ("CG", "opls_271", "C", 0.70), ("OD1", "opls_272", "O", -0.80), ("OD2", "opls_272", "O", -0.80),
    ("C", "opls_235", "C", 0.50), ("O", "opls_236", "O", -0.50),
]
_ASP_BONDS = [
    ("N", "H"), ("N", "CA"), ("CA", "C"), ("C", "O"),
    ("CA", "CB"), ("CB", "CG"), ("CG", "OD1"), ("CG", "OD2"),
    ("CA", "HA"), ("CB", "HB1"), ("CB", "HB2"),
]

_GLU = [
    ("N", "opls_238", "N", -0.50), ("H", "opls_241", "H", 0.30),
    ("CA", "opls_224B", "C", 0.14), ("HA", "opls_140", "H", 0.06),
    ("CB", "opls_274", "C", -0.22), ("HB1", "opls_140", "H", 0.06), ("HB2", "opls_140", "H", 0.06),
    ("CG", "opls_274", "C", -0.22), ("HG1", "opls_140", "H", 0.06), ("HG2", "opls_140", "H", 0.06),
    ("CD", "opls_271", "C", 0.70), ("OE1", "opls_272", "O", -0.80), ("OE2", "opls_272", "O", -0.80),
    ("C", "opls_235", "C", 0.50), ("O", "opls_236", "O", -0.50),
]
_GLU_BONDS = [
    ("N", "H"), ("N", "CA"), ("CA", "C"), ("C", "O"),
    ("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "OE1"), ("CD", "OE2"),
    ("CA", "HA"), ("CB", "HB1"), ("CB", "HB2"), ("CG", "HG1"), ("CG", "HG2"),
]

_SOL = [
    ("OW", "opls_116", "O", -0.82), ("HW1", "opls_117", "H", 0.41), ("HW2", "opls_117", "H", 0.41),
]
_SOL_BONDS = [("OW", "HW1"), ("OW", "HW2")]

_TEMPLATES = {
    "GLY_NTER": (_GLY_NTER, _GLY_NTER_BONDS),
    "PHE": (_PHE, _PHE_BONDS),
    "TRP": (_TRP, _TRP_BONDS),
    "HISH": (_HISH, _HISH_BONDS),
    "LYSH": (_LYSH, _LYSH_BONDS),
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
        system_name="synthetic_opls",
        residues=sd_residues,
        inter_residue_bonds=[],
    )


@pytest.fixture(scope="module")
def identifier():
    return OplsFFGroupIdentifier()


@pytest.fixture(scope="module")
def groups(identifier):
    """识别合成体系（含多种代表残基）的全部基团。"""
    sd = _make_system(
        ["GLY_NTER", "PHE", "TRP", "HISH", "LYSH", "ARG", "ASP", "GLU", "SOL"]
    )
    return identifier.identify(sd)


def _groups_of(resname, groups):
    return [g for g in groups if g.residue_name == resname]


# ============================================================
# N 端 / 带电残基
# ============================================================

class TestChargedGroups:

    def test_nterm_positive(self, groups):
        """N 端 NH3+ → 1 个正电基团（结构判据 N+H 验证）。"""
        pos = [g for g in _groups_of("GLY_NTER", groups)
               if g.group_type == "charged_positive"]
        assert len(pos) == 1
        assert pos[0].num_atoms == 4
        assert pos[0].net_charge > 0.1

    def test_lysh_positive_single(self, groups):
        """LYSH（NH3+）→ 1 个正电基团（无双识别）。"""
        pos = [g for g in _groups_of("LYSH", groups)
               if g.group_type == "charged_positive"]
        assert len(pos) == 1
        assert pos[0].num_atoms == 4

    def test_arg_positive_single(self, groups):
        pos = [g for g in _groups_of("ARG", groups)
               if g.group_type == "charged_positive"]
        assert len(pos) == 1
        assert pos[0].net_charge > 0.1

    def test_hish_positive_single(self, groups):
        """HISH（质子化 His）→ 1 个正电基团（含咪唑碳）。"""
        pos = [g for g in _groups_of("HISH", groups)
               if g.group_type == "charged_positive"]
        assert len(pos) == 1
        assert pos[0].num_atoms == 7          # ND1/NE2/HD1/HE2 + CG/CE1/CD2
        assert pos[0].net_charge > 0.1

    def test_asp_glu_negative(self, groups):
        for resname in ("ASP", "GLU"):
            neg = [g for g in _groups_of(resname, groups)
                   if g.group_type == "charged_negative"]
            assert len(neg) == 1, resname
            assert neg[0].net_charge < -0.1, resname


# ============================================================
# 芳香环
# ============================================================

class TestAromaticRings:

    def test_phe_single_ring(self, groups):
        rings = [g for g in _groups_of("PHE", groups)
                 if g.group_type == "aromatic_ring"]
        assert len(rings) == 1
        assert len(rings[0].atoms) == 6

    def test_trp_two_rings(self, groups):
        rings = [g for g in _groups_of("TRP", groups)
                 if g.group_type == "aromatic_ring"]
        assert len(rings) == 2
        assert sorted(len(r.atoms) for r in rings) == [5, 6]

    def test_trp_ring_types(self, groups):
        """TRP 环应含 C*/CN/NA 等 OPLS 芳香类型。"""
        rings = [g for g in _groups_of("TRP", groups)
                 if g.group_type == "aromatic_ring"]
        all_types = {a.atom_type for r in rings for a in r.atoms}
        assert {"opls_500", "opls_502", "opls_503"} <= all_types

    def test_hish_imidazole_ring(self, groups):
        rings = [g for g in _groups_of("HISH", groups)
                 if g.group_type == "aromatic_ring"]
        assert len(rings) == 1
        assert len(rings[0].atoms) == 5


# ============================================================
# 供体 / 受体 / 水
# ============================================================

class TestDonorAcceptor:

    def test_donors_present(self, groups):
        phe_donors = [g for g in _groups_of("PHE", groups)
                      if g.group_type == "H_donor"]
        assert len(phe_donors) == 1
        assert phe_donors[0].atoms[0].atom_element == "N"
        assert phe_donors[0].atoms[1].atom_element == "H"

    def test_trp_nh_donors(self, groups):
        """TRP 主链 N-H + 吲哚 NE1-H 两个供体。"""
        trp_donors = [g for g in _groups_of("TRP", groups)
                      if g.group_type == "H_donor"]
        assert len(trp_donors) == 2

    def test_asp_acceptors(self, groups):
        """ASP 羧基氧（opls_272，q<0）为受体。"""
        asp_acc = [g for g in _groups_of("ASP", groups)
                   if g.group_type == "H_acceptor"]
        acc_names = {a.atoms[0].atom_name for a in asp_acc}
        assert "OD1" in acc_names and "OD2" in acc_names
        assert all(a.atoms[0].atom_charge < 0 for a in asp_acc)

    def test_water_group(self, groups):
        waters = [g for g in _groups_of("SOL", groups)
                  if g.group_type == "water"]
        assert len(waters) == 1
        assert waters[0].num_atoms == 3


class TestHydrophobic:

    def test_aliphatic_carbon_hydrophobic(self, groups):
        """LYSH 脂肪碳（CT，邻居 C/H）→ 疏水。"""
        lysh_hydro = [g for g in _groups_of("LYSH", groups)
                      if g.group_type == "hydrophobic"]
        assert len(lysh_hydro) >= 1

    def test_polar_adjacent_carbon_not_hydrophobic(self, groups):
        """ASP CG（opls_271，连 O）→ 不得判为疏水。"""
        asp_hydro = [g for g in _groups_of("ASP", groups)
                     if g.group_type == "hydrophobic"]
        asp_hydro_names = {g.atoms[0].atom_name for g in asp_hydro}
        assert "CG" not in asp_hydro_names


# ============================================================
# 注册表 / 接口
# ============================================================

class TestRegistry:

    def test_name(self, identifier):
        assert identifier.name == "opls_ff"

    def test_registry_contains_opls(self):
        from DuIvyInteractions.group_identifiers import IDENTIFIER_CLASSES
        assert "opls" in IDENTIFIER_CLASSES
        assert IDENTIFIER_CLASSES["opls"] is OplsFFGroupIdentifier

    def test_pipeline_make_identifier(self):
        from DuIvyInteractions.pipeline import Pipeline
        assert isinstance(Pipeline("opls")._make_identifier(),
                          OplsFFGroupIdentifier)

    def test_unknown_ff_error_lists_opls(self):
        from DuIvyInteractions.pipeline import Pipeline
        try:
            Pipeline("nope")._make_identifier()
            raise AssertionError("should raise")
        except ValueError as e:
            assert "opls" in str(e)