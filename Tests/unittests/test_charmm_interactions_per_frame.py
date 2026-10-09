# -*- coding: utf-8 -*-
"""CHARMM36 真实体系相互作用测试（策略二：PerFrame）。

SMO-BST（Smoothened + β-谷甾醇），`Tests/test_MD_case_charmm36/` 真实 tpr+xtc
（GROMACS 2018.1，Omar et al. 2020, *Data in Brief*；Mendeley v94vzbwzf3）。

覆盖：氢键 / 疏水（蛋白-配体间）；盐桥 / π-π / π-阳离子（蛋白内）；
卤键 / 金属配位（体系无卤素/金属，断言 0 对）。

⚠️ 断言基线（pair 数：氢键 3、疏水 25、盐桥 58、π-π 47）来自本策略实现输出，
【尚未经人工核验】（未与 gmx 工具或人工分子分析交叉确认），请谨慎引用。
"""

import pytest
from pathlib import Path

import MDAnalysis as mda

from DuIvyInteractions.system_readers import GmxTprReader
from DuIvyInteractions.group_identifiers.charmm_ff_identifier import CharmmFFGroupIdentifier
from DuIvyInteractions.interaction_detectors import (
    HydrogenBondDetectorPerFrame,
    HydrophobicDetectorPerFrame,
    SaltBridgeDetectorPerFrame,
    PiStackingDetectorPerFrame,
    PiCationDetectorPerFrame,
    HalogenBondDetectorPerFrame,
    MetalCoordinationDetectorPerFrame,
)


TPR_FILE = Path(__file__).parent.parent / "test_MD_case_charmm36" / "SMO-BST_md_complex.tpr"
XTC_FILE = Path(__file__).parent.parent / "test_MD_case_charmm36" / "SMO-BST_md_complex.xtc"


def _protein_ligand_filter(gt):
    """蛋白-配体间配对（BST 配体 vs 蛋白残基）。"""
    return (gt[0].residue_name == "BST") != (gt[1].residue_name == "BST")


@pytest.fixture(scope="module")
def groups():
    sd = GmxTprReader().read(str(TPR_FILE))
    return CharmmFFGroupIdentifier().identify(sd)


@pytest.fixture(scope="module")
def trajectory():
    u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    return u.trajectory


@pytest.fixture(scope="module")
def hydrogen_bonds(groups, trajectory):
    return HydrogenBondDetectorPerFrame().detect(
        groups, trajectory=trajectory, tuple_filter=_protein_ligand_filter)


@pytest.fixture(scope="module")
def hydrophobics(groups, trajectory):
    return HydrophobicDetectorPerFrame().detect(
        groups, trajectory=trajectory, tuple_filter=_protein_ligand_filter)


@pytest.fixture(scope="module")
def salt_bridges(groups, trajectory):
    return SaltBridgeDetectorPerFrame().detect(groups, trajectory=trajectory)


@pytest.fixture(scope="module")
def pi_stackings(groups, trajectory):
    return PiStackingDetectorPerFrame().detect(groups, trajectory=trajectory)


@pytest.fixture(scope="module")
def pi_cations(groups, trajectory):
    return PiCationDetectorPerFrame().detect(groups, trajectory=trajectory)


@pytest.fixture(scope="module")
def halogen_bonds(groups, trajectory):
    return HalogenBondDetectorPerFrame().detect(groups, trajectory=trajectory)


@pytest.fixture(scope="module")
def metal_coordinations(groups, trajectory):
    return MetalCoordinationDetectorPerFrame().detect(groups, trajectory=trajectory)


# ============================================================
# 氢键（蛋白-配体）
# ============================================================

class TestHydrogenBondPerFrame:

    def test_has_results(self, hydrogen_bonds):
        assert len(hydrogen_bonds) > 0

    def test_n_pairs(self, hydrogen_bonds):
        assert hydrogen_bonds[0].n_pairs == 3

    def test_all_protein_ligand(self, hydrogen_bonds):
        it = hydrogen_bonds[0]
        for g1, g2 in it.groups:
            assert (g1.residue_name == "BST") != (g2.residue_name == "BST")

    def test_metric_ndim(self, hydrogen_bonds):
        assert hydrogen_bonds[0].metrics["distance"].ndim == 2


# ============================================================
# 疏水（蛋白-配体）
# ============================================================

class TestHydrophobicPerFrame:

    def test_has_results(self, hydrophobics):
        assert len(hydrophobics) > 0

    def test_n_pairs(self, hydrophobics):
        assert hydrophobics[0].n_pairs == 25

    def test_all_protein_ligand(self, hydrophobics):
        it = hydrophobics[0]
        for g1, g2 in it.groups:
            assert (g1.residue_name == "BST") != (g2.residue_name == "BST")


# ============================================================
# 盐桥（蛋白内）
# ============================================================

class TestSaltBridgePerFrame:

    def test_has_results(self, salt_bridges):
        assert len(salt_bridges) > 0

    def test_n_pairs(self, salt_bridges):
        assert salt_bridges[0].n_pairs == 58


# ============================================================
# π-π（蛋白内）
# ============================================================

class TestPiStackingPerFrame:

    def test_has_results(self, pi_stackings):
        assert len(pi_stackings) > 0

    def test_n_pairs(self, pi_stackings):
        assert pi_stackings[0].n_pairs == 47

    def test_all_aromatic_residues(self, pi_stackings):
        it = pi_stackings[0]
        for g1, g2 in it.groups:
            assert g1.group_type == "aromatic_ring"
            assert g2.group_type == "aromatic_ring"

# ============================================================
# π-阳离子（蛋白内）
# ============================================================

class TestPiCationPerFrame:

    def test_has_results(self, pi_cations):
        assert len(pi_cations) > 0

    def test_n_pairs(self, pi_cations):
        assert pi_cations[0].n_pairs == 24


# ============================================================
# 卤键 / 金属配位（体系无对应基团 → 0 对）
# ============================================================

class TestHalogenBondPerFrame:

    def test_no_results(self, halogen_bonds):
        """SMO/BST 无卤素（无 σ-hole 供体） → 无卤键。"""
        assert len(halogen_bonds) == 0


class TestMetalCoordinationPerFrame:

    def test_no_results(self, metal_coordinations):
        """SMO/BST 无金属中心 → 无金属配位。"""
        assert len(metal_coordinations) == 0
