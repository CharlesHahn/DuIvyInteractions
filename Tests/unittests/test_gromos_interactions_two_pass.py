# -*- coding: utf-8 -*-
"""GROMOS 53A6 真实体系相互作用测试（策略三：TwoPass）。

覆盖：氢键 / 疏水（蛋白-配体间）；盐桥 / π-π（蛋白内）。
使用 `Tests/test_MD_case_gromos/` 真实 tpr+xtc。

⚠️ 注意：本测试的全部数值断言（pair 数：氢键 379、疏水 161、
盐桥 18、π-π 8）来自程序自动探测，【尚未经人工核验】
（未与 gmx 工具或人工分子分析交叉确认）。数值可能与真实化学
不符，请谨慎引用，待人工核验后更新。
"""

import pytest
from pathlib import Path

import MDAnalysis as mda

from DuIvyInteractions.system_readers import GmxTprReader
from DuIvyInteractions.group_identifiers import GromosFFGroupIdentifier
from DuIvyInteractions.interaction_detectors import (
    HydrogenBondDetectorTwoPass,
    HydrophobicDetectorTwoPass,
    SaltBridgeDetectorTwoPass,
    PiStackingDetectorTwoPass,
)


TPR_FILE = Path(__file__).parent.parent / "test_MD_case_gromos" / "gromos53a6_md.tpr"
XTC_FILE = Path(__file__).parent.parent / "test_MD_case_gromos" / "gromos53a6_md10ns.xtc"
LIGAND_RESIDUES = {"1ZIN", "2ZIN", "3ZIN", "4ZIN", "5ZIN", "6ZIN"}


def _protein_ligand_filter(gt):
    """蛋白-配体间配对（ZIN 与非 ZIN 残基）。"""
    return ("ZIN" in gt[0].residue_name) != ("ZIN" in gt[1].residue_name)


@pytest.fixture(scope="module")
def groups():
    sd = GmxTprReader().read(str(TPR_FILE))
    return GromosFFGroupIdentifier().identify(sd)


@pytest.fixture(scope="module")
def trajectory():
    u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    return u.trajectory


@pytest.fixture(scope="module")
def hydrogen_bonds(groups, trajectory):
    d = HydrogenBondDetectorTwoPass()
    return d.detect(groups, trajectory=trajectory,
                    tuple_filter=_protein_ligand_filter)


@pytest.fixture(scope="module")
def hydrophobics(groups, trajectory):
    d = HydrophobicDetectorTwoPass()
    return d.detect(groups, trajectory=trajectory,
                    tuple_filter=_protein_ligand_filter)


@pytest.fixture(scope="module")
def salt_bridges(groups, trajectory):
    d = SaltBridgeDetectorTwoPass()
    return d.detect(groups, trajectory=trajectory)


@pytest.fixture(scope="module")
def pi_stackings(groups, trajectory):
    d = PiStackingDetectorTwoPass()
    return d.detect(groups, trajectory=trajectory)


# ============================================================
# 氢键（蛋白-配体）
# ============================================================

class TestHydrogenBondTwoPass:

    def test_has_results(self, hydrogen_bonds):
        assert len(hydrogen_bonds) > 0

    def test_n_pairs(self, hydrogen_bonds):
        assert hydrogen_bonds[0].n_pairs == 379

    def test_all_inter_protein_ligand(self, hydrogen_bonds):
        it = hydrogen_bonds[0]
        for g1, g2 in it.groups:
            assert ("ZIN" in g1.residue_name) != ("ZIN" in g2.residue_name)

    def test_no_nan_metrics(self, hydrogen_bonds):
        it = hydrogen_bonds[0]
        assert it.metrics["distance"].ndim == 2
        assert it.metrics["angle"].ndim == 2


# ============================================================
# 疏水（蛋白-配体）
# ============================================================

class TestHydrophobicTwoPass:

    def test_has_results(self, hydrophobics):
        assert len(hydrophobics) > 0

    def test_n_pairs(self, hydrophobics):
        assert hydrophobics[0].n_pairs == 161


# ============================================================
# 盐桥（蛋白内）
# ============================================================

class TestSaltBridgeTwoPass:

    def test_has_results(self, salt_bridges):
        assert len(salt_bridges) > 0

    def test_n_pairs(self, salt_bridges):
        assert salt_bridges[0].n_pairs == 18

    def test_all_opposite_charge(self, salt_bridges):
        it = salt_bridges[0]
        for g1, g2 in it.groups:
            assert g1.group_type != g2.group_type


# ============================================================
# π-π（蛋白内 PHE）
# ============================================================

class TestPiStackingTwoPass:

    def test_has_results(self, pi_stackings):
        assert len(pi_stackings) > 0

    def test_n_pairs(self, pi_stackings):
        assert pi_stackings[0].n_pairs == 8

    def test_all_phe(self, pi_stackings):
        it = pi_stackings[0]
        for g1, g2 in it.groups:
            assert g1.residue_name == "PHE" and g2.residue_name == "PHE"