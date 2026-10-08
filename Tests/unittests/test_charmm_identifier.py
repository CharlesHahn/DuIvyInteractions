# -*- coding: utf-8 -*-
"""CHARMM 识别器单元测试（合成 SystemData，不依赖真实 tpr）。

原子/键/电荷数据取自本地 `charmm36-feb2026_cgenff-5.0.ff`（aminoacids.rtp、
n.tdb/c.tdb patch），验证 CharmmFFGroupIdentifier 的基团识别逻辑：
- N 端 NH3+ 正电识别（CHARMM 特有：N 负电荷但 N+H 净电荷通过）
- 芳香环识别（CA/CPH*/CPT/CY/NR* 类型表）
- 带电残基字典（LYS/ARG/HSP + ASP/GLU/CYM）
- 供体/受体/水/疏水
"""

import pytest

from DuIvyInteractions.core.datas import (
    SystemData, ResidueData, AtomData, BondData,
)
from DuIvyInteractions.group_identifiers.charmm_ff_identifier import (
    CharmmFFGroupIdentifier,
)

# ============================================================
# 残基模板（CHARMM36 5.0 aminoacids.rtp / n.tdb / c.tdb 实测）
# 每项: (atom_name, atom_type, element, charge)
# ============================================================

# N 端 GLY（NH3+ patch：N NH3 -0.30, 3H +0.33）
_GLY_NTER = [
    ("N", "NH3", "N", -0.30), ("HT1", "HC", "H", 0.33),
    ("HT2", "HC", "H", 0.33), ("HT3", "HC", "H", 0.33),
    ("CA", "CT2", "C", 0.13), ("HA2", "HB2", "H", 0.09),
    ("C", "C", "C", 0.51), ("O", "O", "O", -0.51),
]
_GLY_NTER_BONDS = [
    ("N", "HT1"), ("N", "HT2"), ("N", "HT3"), ("N", "CA"),
    ("CA", "C"), ("C", "O"), ("CA", "HA2"),
]

_PHE = [
    ("N", "NH1", "N", -0.47), ("HN", "H", "H", 0.31),
    ("CA", "CT1", "C", 0.07), ("HA", "HB1", "H", 0.09),
    ("CB", "CT2", "C", -0.18), ("HB1", "HA2", "H", 0.09), ("HB2", "HA2", "H", 0.09),
    ("CG", "CA", "C", 0.00), ("CD1", "CA", "C", -0.115), ("HD1", "HP", "H", 0.115),
    ("CD2", "CA", "C", -0.115), ("HD2", "HP", "H", 0.115),
    ("CE1", "CA", "C", -0.115), ("HE1", "HP", "H", 0.115),
    ("CE2", "CA", "C", -0.115), ("HE2", "HP", "H", 0.115),
    ("CZ", "CA", "C", -0.115), ("HZ", "HP", "H", 0.115),
    ("C", "C", "C", 0.51), ("O", "O", "O", -0.51),
]
_PHE_BONDS = [
    ("N", "HN"), ("N", "CA"), ("CA", "C"), ("C", "O"),
    ("CA", "CB"), ("CB", "CG"), ("CA", "HA"),
    ("CB", "HB1"), ("CB", "HB2"),
    ("CG", "CD1"), ("CG", "CD2"), ("CD1", "HD1"), ("CD1", "CE1"),
    ("CD2", "HD2"), ("CD2", "CE2"), ("CE1", "HE1"), ("CE1", "CZ"),
    ("CE2", "HE2"), ("CE2", "CZ"), ("CZ", "HZ"),
]

_TRP = [
    ("N", "NH1", "N", -0.47), ("HN", "H", "H", 0.31),
    ("CA", "CT1", "C", 0.07), ("HA", "HB1", "H", 0.09),
    ("CB", "CT2", "C", -0.18), ("HB1", "HA2", "H", 0.09), ("HB2", "HA2", "H", 0.09),
    ("CG", "CY", "C", -0.03), ("CD1", "CA", "C", -0.15), ("HD1", "HP", "H", 0.22),
    ("NE1", "NY", "N", -0.51), ("HE1", "H", "H", 0.37),
    ("CE2", "CPT", "C", 0.24), ("CD2", "CPT", "C", 0.11),
    ("CE3", "CAI", "C", -0.25), ("HE3", "HP", "H", 0.17),
    ("CZ3", "CA", "C", -0.20), ("HZ3", "HP", "H", 0.14),
    ("CZ2", "CAI", "C", -0.27), ("HZ2", "HP", "H", 0.16),
    ("CH2", "CA", "C", -0.14), ("HH2", "HP", "H", 0.14),
    ("C", "C", "C", 0.51), ("O", "O", "O", -0.51),
]
_TRP_BONDS = [
    ("N", "HN"), ("N", "CA"), ("CA", "C"), ("C", "O"),
    ("CA", "CB"), ("CB", "CG"), ("CA", "HA"),
    ("CB", "HB1"), ("CB", "HB2"),
    ("CG", "CD1"), ("CG", "CD2"), ("CD1", "HD1"), ("CD1", "NE1"),
    ("NE1", "HE1"), ("NE1", "CE2"), ("CE2", "CD2"), ("CE2", "CZ2"),
    ("CD2", "CE3"), ("CE3", "HE3"), ("CE3", "CZ3"),
    ("CZ3", "HZ3"), ("CZ3", "CH2"), ("CZ2", "HZ2"), ("CZ2", "CH2"),
    ("CH2", "HH2"),
]

_HSP = [
    ("N", "NH1", "N", -0.47), ("HN", "H", "H", 0.31),
    ("CA", "CT1", "C", 0.07), ("HA", "HB1", "H", 0.09),
    ("CB", "CT2A", "C", -0.05), ("HB1", "HA2", "H", 0.09), ("HB2", "HA2", "H", 0.09),
    ("CD2", "CPH1", "C", 0.19), ("HD2", "HR1", "H", 0.13),
    ("CG", "CPH1", "C", 0.19), ("NE2", "NR3", "N", -0.51), ("HE2", "H", "H", 0.44),
    ("ND1", "NR3", "N", -0.51), ("HD1", "H", "H", 0.44),
    ("CE1", "CPH2", "C", 0.32), ("HE1", "HR2", "H", 0.18),
    ("C", "C", "C", 0.51), ("O", "O", "O", -0.51),
]
_HSP_BONDS = [
    ("N", "HN"), ("N", "CA"), ("CA", "C"), ("C", "O"),
    ("CA", "CB"), ("CB", "CG"), ("CA", "HA"),
    ("CB", "HB1"), ("CB", "HB2"),
    ("CG", "ND1"), ("CG", "CD2"), ("ND1", "HD1"), ("ND1", "CE1"),
    ("CD2", "HD2"), ("CD2", "NE2"), ("CE1", "HE1"), ("CE1", "NE2"),
    ("NE2", "HE2"),
]

_LYS = [
    ("N", "NH1", "N", -0.47), ("HN", "H", "H", 0.31),
    ("CA", "CT1", "C", 0.07), ("HA", "HB1", "H", 0.09),
    ("CB", "CT2", "C", -0.18), ("HB1", "HA2", "H", 0.09), ("HB2", "HA2", "H", 0.09),
    ("CG", "CT2", "C", -0.18), ("HG1", "HA2", "H", 0.09), ("HG2", "HA2", "H", 0.09),
    ("CD", "CT2", "C", -0.18), ("HD1", "HA2", "H", 0.09), ("HD2", "HA2", "H", 0.09),
    ("CE", "CT2", "C", -0.02), ("HE1", "HA2", "H", 0.09), ("HE2", "HA2", "H", 0.09),
    ("NZ", "NH3", "N", -0.30), ("HZ1", "H", "H", 0.33),
    ("HZ2", "H", "H", 0.33), ("HZ3", "H", "H", 0.33),
    ("C", "C", "C", 0.51), ("O", "O", "O", -0.51),
]
_LYS_BONDS = [
    ("N", "HN"), ("N", "CA"), ("CA", "C"), ("C", "O"),
    ("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "CE"), ("CE", "NZ"),
    ("NZ", "HZ1"), ("NZ", "HZ2"), ("NZ", "HZ3"),
    ("CA", "HA"), ("CB", "HB1"), ("CB", "HB2"),
    ("CG", "HG1"), ("CG", "HG2"), ("CD", "HD1"), ("CD", "HD2"),
    ("CE", "HE1"), ("CE", "HE2"),
]

_ARG = [
    ("N", "NH1", "N", -0.47), ("HN", "H", "H", 0.31),
    ("CA", "CT1", "C", 0.07), ("HA", "HB1", "H", 0.09),
    ("CB", "CT2", "C", -0.18), ("HB1", "HA2", "H", 0.09), ("HB2", "HA2", "H", 0.09),
    ("CG", "CT2", "C", -0.18), ("HG1", "HA2", "H", 0.09), ("HG2", "HA2", "H", 0.09),
    ("CD", "CT2", "C", 0.20), ("HD1", "HA2", "H", 0.09), ("HD2", "HA2", "H", 0.09),
    ("NE", "NC2", "N", -0.70), ("HE", "HC", "H", 0.44),
    ("CZ", "C", "C", 0.64), ("NH1", "NC2", "N", -0.80),
    ("HH11", "HC", "H", 0.46), ("HH12", "HC", "H", 0.46),
    ("NH2", "NC2", "N", -0.80), ("HH21", "HC", "H", 0.46),
    ("HH22", "HC", "H", 0.46),
    ("C", "C", "C", 0.51), ("O", "O", "O", -0.51),
]
_ARG_BONDS = [
    ("N", "HN"), ("N", "CA"), ("CA", "C"), ("C", "O"),
    ("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "NE"),
    ("NE", "HE"), ("NE", "CZ"), ("CZ", "NH1"), ("CZ", "NH2"),
    ("NH1", "HH11"), ("NH1", "HH12"), ("NH2", "HH21"), ("NH2", "HH22"),
    ("CA", "HA"), ("CB", "HB1"), ("CB", "HB2"),
    ("CG", "HG1"), ("CG", "HG2"), ("CD", "HD1"), ("CD", "HD2"),
]

_ASP = [
    ("N", "NH1", "N", -0.47), ("HN", "H", "H", 0.31),
    ("CA", "CT1", "C", 0.07), ("HA", "HB1", "H", 0.09),
    ("CB", "CT2", "C", -0.28), ("HB1", "HA2", "H", 0.09), ("HB2", "HA2", "H", 0.09),
    ("CG", "CC", "C", 0.62), ("OD1", "OC", "O", -0.76), ("OD2", "OC", "O", -0.76),
    ("C", "C", "C", 0.51), ("O", "O", "O", -0.51),
]
_ASP_BONDS = [
    ("N", "HN"), ("N", "CA"), ("CA", "C"), ("C", "O"),
    ("CA", "CB"), ("CB", "CG"), ("CG", "OD1"), ("CG", "OD2"),
    ("CA", "HA"), ("CB", "HB1"), ("CB", "HB2"),
]

_GLU = [
    ("N", "NH1", "N", -0.47), ("HN", "H", "H", 0.31),
    ("CA", "CT1", "C", 0.07), ("HA", "HB1", "H", 0.09),
    ("CB", "CT2A", "C", -0.18), ("HB1", "HA2", "H", 0.09), ("HB2", "HA2", "H", 0.09),
    ("CG", "CT2", "C", -0.28), ("HG1", "HA2", "H", 0.09), ("HG2", "HA2", "H", 0.09),
    ("CD", "CC", "C", 0.62), ("OE1", "OC", "O", -0.76), ("OE2", "OC", "O", -0.76),
    ("C", "C", "C", 0.51), ("O", "O", "O", -0.51),
]
_GLU_BONDS = [
    ("N", "HN"), ("N", "CA"), ("CA", "C"), ("C", "O"),
    ("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "OE1"), ("CD", "OE2"),
    ("CA", "HA"), ("CB", "HB1"), ("CB", "HB2"),
    ("CG", "HG1"), ("CG", "HG2"),
]

_CYM = [
    ("N", "NH1", "N", -0.47), ("HN", "H", "H", 0.31),
    ("CA", "CT1", "C", 0.07), ("HA", "HB1", "H", 0.09),
    ("CB", "CT2", "C", 0.05), ("HB1", "HA2", "H", 0.09), ("HB2", "HA2", "H", 0.09),
    ("SG", "SS", "S", -0.80),
    ("C", "C", "C", 0.51), ("O", "O", "O", -0.51),
]
_CYM_BONDS = [
    ("N", "HN"), ("N", "CA"), ("CA", "C"), ("C", "O"),
    ("CA", "CB"), ("CB", "SG"), ("CA", "HA"), ("CB", "HB1"), ("CB", "HB2"),
]

_SOL = [
    ("OH2", "OT", "O", -0.834), ("H1", "HT", "H", 0.417), ("H2", "HT", "H", 0.417),
]
_SOL_BONDS = [("OH2", "H1"), ("OH2", "H2")]

_TEMPLATES = {
    "GLY_NTER": (_GLY_NTER, _GLY_NTER_BONDS),
    "PHE": (_PHE, _PHE_BONDS),
    "TRP": (_TRP, _TRP_BONDS),
    "HSP": (_HSP, _HSP_BONDS),
    "LYS": (_LYS, _LYS_BONDS),
    "ARG": (_ARG, _ARG_BONDS),
    "ASP": (_ASP, _ASP_BONDS),
    "GLU": (_GLU, _GLU_BONDS),
    "CYM": (_CYM, _CYM_BONDS),
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
        system_name="synthetic_charmm",
        residues=sd_residues,
        inter_residue_bonds=[],
    )


@pytest.fixture(scope="module")
def identifier():
    return CharmmFFGroupIdentifier()


@pytest.fixture(scope="module")
def groups(identifier):
    """识别合成体系（含多种代表残基）的全部基团。"""
    sd = _make_system(
        ["GLY_NTER", "PHE", "TRP", "HSP", "LYS", "ARG", "ASP", "GLU", "CYM", "SOL"]
    )
    return identifier.identify(sd)


def _groups_of(resname, groups):
    return [g for g in groups if g.residue_name == resname]


# ============================================================
# N 端 / 带电残基
# ============================================================

class TestChargedGroups:

    def test_nterm_positive(self, groups):
        """N 端 NH3+ → 恰好 1 个正电基团（N+H 净电荷验证生效）。"""
        pos = [g for g in _groups_of("GLY_NTER", groups)
               if g.group_type == "charged_positive"]
        assert len(pos) == 1
        assert pos[0].num_atoms == 4          # N + 3H
        assert pos[0].net_charge > 0.1

    def test_lys_positive_single(self, groups):
        """LYS（NH3+）→ 1 个正电基团（无 tertamine 双识别）。"""
        pos = [g for g in _groups_of("LYS", groups)
               if g.group_type == "charged_positive"]
        assert len(pos) == 1
        assert pos[0].num_atoms == 4          # NZ + 3H

    def test_arg_positive_single(self, groups):
        pos = [g for g in _groups_of("ARG", groups)
               if g.group_type == "charged_positive"]
        assert len(pos) == 1
        assert pos[0].net_charge > 0.1

    def test_hsp_positive_single(self, groups):
        """HSP（质子化 His）→ 1 个正电基团（含咪唑碳，净电荷为正）。"""
        pos = [g for g in _groups_of("HSP", groups)
               if g.group_type == "charged_positive"]
        assert len(pos) == 1
        assert pos[0].num_atoms == 7         # ND1/NE2/HD1/HE2 + CG/CE1/CD2
        assert pos[0].net_charge > 0.1

    def test_asp_glu_cym_negative(self, groups):
        for resname in ("ASP", "GLU", "CYM"):
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
        assert len(rings[0].atoms) == 6       # 全 CA 环

    def test_trp_two_rings(self, groups):
        rings = [g for g in _groups_of("TRP", groups)
                 if g.group_type == "aromatic_ring"]
        assert len(rings) == 2
        assert sorted(len(r.atoms) for r in rings) == [5, 6]

    def test_trp_ring_types(self, groups):
        """TRP 环应含 CY/CPT/NY 等 CHARMM 芳香类型。"""
        rings = [g for g in _groups_of("TRP", groups)
                 if g.group_type == "aromatic_ring"]
        all_types = {a.atom_type for r in rings for a in r.atoms}
        assert {"CY", "CPT", "NY"} <= all_types

    def test_hsp_imidazole_ring(self, groups):
        rings = [g for g in _groups_of("HSP", groups)
                 if g.group_type == "aromatic_ring"]
        assert len(rings) == 1
        assert len(rings[0].atoms) == 5       # CPH1/CPH2 咪唑


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

    def test_trp_nh_donors(self, groups):
        """TRP 主链 N-H + 吲哚 NE1-H 两个供体。"""
        trp_donors = [g for g in _groups_of("TRP", groups)
                      if g.group_type == "H_donor"]
        assert len(trp_donors) == 2

    def test_asp_acceptors(self, groups):
        """ASP 羧基氧（OC，q<0）为受体。"""
        asp_acc = [g for g in _groups_of("ASP", groups)
                   if g.group_type == "H_acceptor"]
        acc_names = {a.atoms[0].atom_name for a in asp_acc}
        assert "OD1" in acc_names and "OD2" in acc_names
        assert all(a.atoms[0].atom_charge < 0 for a in asp_acc)

    def test_backbone_N_not_acceptor(self, groups):
        """主链肽键 N（NH1=peptide nitrogen，带 H）不应是受体（Eildal 2013, JACS）。

        回归保护：NH1/NH2/NH3/NC2/NY/NR1/NR3 已从 CHARMM_ACCEPTOR_TYPES 剔除
        （孤对被共振占用/无孤对）；仅保留 N(Pro)/NR2 作为蛋白氮受体。
        """
        acc_names = {a.atoms[0].atom_name for a in
                     [g for g in _groups_of("PHE", groups)
                      if g.group_type == "H_acceptor"]}
        assert "N" not in acc_names, "主链肽键 N 不应是受体"

    def test_water_group(self, groups):
        waters = [g for g in _groups_of("SOL", groups)
                  if g.group_type == "water"]
        assert len(waters) == 1
        assert waters[0].num_atoms == 3       # OH2, H1, H2


class TestHydrophobic:

    def test_aliphatic_carbon_hydrophobic(self, groups):
        """LYS 脂肪碳（CT2，邻居 C/H）→ 疏水。"""
        lys_hydro = [g for g in _groups_of("LYS", groups)
                     if g.group_type == "hydrophobic"]
        assert len(lys_hydro) >= 1

    def test_polar_adjacent_carbon_not_hydrophobic(self, groups):
        """ASP CG（CC，连 O）→ 不得判为疏水。"""
        asp_hydro = [g for g in _groups_of("ASP", groups)
                     if g.group_type == "hydrophobic"]
        asp_hydro_names = {g.atoms[0].atom_name for g in asp_hydro}
        assert "CG" not in asp_hydro_names


# ============================================================
# 注册表 / 接口
# ============================================================

class TestRegistry:

    def test_name(self, identifier):
        assert identifier.name == "charmm_ff"

    def test_registry_contains_charmm(self):
        from DuIvyInteractions.group_identifiers import IDENTIFIER_CLASSES
        assert "charmm" in IDENTIFIER_CLASSES
        assert IDENTIFIER_CLASSES["charmm"] is CharmmFFGroupIdentifier

    def test_pipeline_make_identifier(self):
        from DuIvyInteractions.pipeline import Pipeline
        assert isinstance(Pipeline("charmm")._make_identifier(),
                          CharmmFFGroupIdentifier)

    def test_unknown_ff_error_lists_charmm(self):
        from DuIvyInteractions.pipeline import Pipeline
        try:
            Pipeline("nope")._make_identifier()
            raise AssertionError("should raise")
        except ValueError as e:
            assert "charmm" in str(e)