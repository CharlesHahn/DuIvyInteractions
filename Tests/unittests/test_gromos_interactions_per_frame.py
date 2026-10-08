# -*- coding: utf-8 -*-
"""GROMOS 53A6 真实体系相互作用测试（策略二：PerFrame）。

覆盖：氢键 / 疏水（蛋白-配体间）；盐桥 / π-π（蛋白内）。
使用 `Tests/test_MD_case_gromos/` 真实 tpr+xtc。

⚠️ 注意：本测试的全部数值断言（pair 数：氢键 290、疏水 16、
盐桥 16、π-π 8）来自程序自动探测，【尚未经人工核验】
（未与 gmx 工具或人工分子分析交叉确认）。数值可能与真实化学
不符，请谨慎引用，待人工核验后更新。
"""

import pytest
from pathlib import Path

import MDAnalysis as mda

from DuIvyInteractions.system_readers import GmxTprReader
from DuIvyInteractions.group_identifiers import GromosFFGroupIdentifier
from DuIvyInteractions.interaction_detectors import (
    HydrogenBondDetectorPerFrame,
    HydrophobicDetectorPerFrame,
    SaltBridgeDetectorPerFrame,
    PiStackingDetectorPerFrame,
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
    d = HydrogenBondDetectorPerFrame()
    return d.detect(groups, trajectory=trajectory,
                    tuple_filter=_protein_ligand_filter)


@pytest.fixture(scope="module")
def hydrophobics(groups, trajectory):
    d = HydrophobicDetectorPerFrame()
    return d.detect(groups, trajectory=trajectory,
                    tuple_filter=_protein_ligand_filter)


@pytest.fixture(scope="module")
def salt_bridges(groups, trajectory):
    d = SaltBridgeDetectorPerFrame()
    return d.detect(groups, trajectory=trajectory)


@pytest.fixture(scope="module")
def pi_stackings(groups, trajectory):
    d = PiStackingDetectorPerFrame()
    return d.detect(groups, trajectory=trajectory)


# ============================================================
# 氢键（蛋白-配体）
# ============================================================

class TestHydrogenBondPerFrame:

    def test_has_results(self, hydrogen_bonds):
        assert len(hydrogen_bonds) > 0

    def test_n_pairs(self, hydrogen_bonds):
        assert hydrogen_bonds[0].n_pairs == 213  # 修复后：剔除带 H 非受体 N，H 键 pair 减少

    def test_all_inter_protein_ligand(self, hydrogen_bonds):
        it = hydrogen_bonds[0]
        for g1, g2 in it.groups:
            assert ("ZIN" in g1.residue_name) != ("ZIN" in g2.residue_name)

    def test_metric_ndim(self, hydrogen_bonds):
        assert hydrogen_bonds[0].metrics["distance"].ndim == 2


# ============================================================
# 疏水（蛋白-配体）
# ============================================================

class TestHydrophobicPerFrame:

    def test_has_results(self, hydrophobics):
        assert len(hydrophobics) > 0

    def test_n_pairs(self, hydrophobics):
        assert hydrophobics[0].n_pairs == 16


# ============================================================
# 盐桥（蛋白内）
# ============================================================

class TestSaltBridgePerFrame:

    def test_has_results(self, salt_bridges):
        assert len(salt_bridges) > 0

    def test_n_pairs(self, salt_bridges):
        assert salt_bridges[0].n_pairs == 16


# ============================================================
# π-π（蛋白内 PHE）
# ============================================================

class TestPiStackingPerFrame:

    def test_has_results(self, pi_stackings):
        assert len(pi_stackings) > 0

    def test_n_pairs(self, pi_stackings):
        assert pi_stackings[0].n_pairs == 8

    def test_all_phe(self, pi_stackings):
        it = pi_stackings[0]
        for g1, g2 in it.groups:
            assert g1.residue_name == "PHE" and g2.residue_name == "PHE"