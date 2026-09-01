# -*- coding: utf-8 -*-
"""氢键 TwoPass 检测器集成测试。

验证策略三（两轮遍历 + 稀疏存储）的氢键检测结果。
"""

import pytest
from pathlib import Path

import numpy as np
import MDAnalysis as mda

from DuIvyInteractions.input_readers import GmxTprReader
from DuIvyInteractions.group_identifiers import AmberFFGroupIdentifier
from DuIvyInteractions.interaction_detectors import (
    HydrogenBondDetectorTwoPass, HydrogenBondDetectorPerFrame)


TPR_FILE = Path(__file__).parent.parent / "test_MD_case" / "md.tpr"
XTC_FILE = Path(__file__).parent.parent / "test_MD_case" / "md1ns.xtc"

MOL_RBD = "seg_0_RBD_pro"
MOL_KRAS = "seg_2_KRAS_pro"
WATER_RESIDUES = {"SOL", "HOH", "WAT"}
PROTEIN_MOLECULES = {MOL_RBD, MOL_KRAS}


@pytest.fixture(scope="module")
def system_data():
    reader = GmxTprReader()
    return reader.read(str(TPR_FILE))


@pytest.fixture(scope="module")
def groups(system_data):
    return AmberFFGroupIdentifier().identify(system_data)


@pytest.fixture(scope="module")
def protein_groups(groups):
    return [g for g in groups
            if g.residue_name not in WATER_RESIDUES
            and g.molecule in PROTEIN_MOLECULES]


@pytest.fixture(scope="module")
def detector():
    return HydrogenBondDetectorTwoPass()


@pytest.fixture(scope="module")
def sparse_result(detector, protein_groups):
    u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    return detector.detect_pass1_only(
        protein_groups, u.trajectory,
        tuple_filter=lambda gt: gt[0].molecule != gt[1].molecule), detector


@pytest.fixture(scope="module")
def hydrogen_bonds(detector, protein_groups):
    u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    return detector.detect(
        protein_groups, trajectory=u.trajectory,
        tuple_filter=lambda gt: gt[0].molecule != gt[1].molecule)


# ============================================================
# Pass1 only
# ============================================================

class TestPass1Only:

    def test_returns_dict(self, sparse_result):
        sparse, _ = sparse_result
        assert isinstance(sparse, dict)

    def test_pair_count(self, sparse_result):
        sparse, _ = sparse_result
        assert len(sparse) == 119

    def test_each_pair_has_metrics(self, sparse_result):
        sparse, _ = sparse_result
        for idx, data in sparse.items():
            assert "distance" in data["metrics"]
            assert "angle" in data["metrics"]
            assert len(data["metrics"]["distance"]) == len(data["frames"])


# ============================================================
# detect：Pass1 + Pass2
# ============================================================

class TestDetect:

    def test_has_results(self, hydrogen_bonds):
        assert len(hydrogen_bonds) > 0

    def test_interaction_type(self, hydrogen_bonds):
        assert hydrogen_bonds[0].interaction_type == "hydrogen_bond"

    def test_n_pairs(self, hydrogen_bonds):
        assert hydrogen_bonds[0].n_pairs == 119

    def test_all_inter_protein(self, hydrogen_bonds):
        for g1, g2 in hydrogen_bonds[0].groups:
            assert g1.molecule != g2.molecule

    def test_top_pair(self, hydrogen_bonds):
        it = hydrogen_bonds[0]
        occ = it.occupancy()
        top = np.argmax(occ)
        g1, g2 = it.groups[top]
        assert g1.residue_name == "ARG" and g1.residue_id == 73
        assert g2.residue_name == "ASP" and g2.residue_id == 175
        assert occ[top] > 0.8

    def test_no_nan(self, hydrogen_bonds):
        assert np.sum(np.isnan(hydrogen_bonds[0].metrics["distance"])) == 0
        assert np.sum(np.isnan(hydrogen_bonds[0].metrics["angle"])) == 0

    def test_active_distance_range(self, hydrogen_bonds):
        it = hydrogen_bonds[0]
        for i in range(it.n_pairs):
            active = it.metrics["distance"][i][it.existence[i]]
            if len(active) > 0:
                assert np.all(active <= 4.1)

    def test_active_angle_range(self, hydrogen_bonds):
        it = hydrogen_bonds[0]
        for i in range(it.n_pairs):
            active = it.metrics["angle"][i][it.existence[i]]
            if len(active) > 0:
                assert np.all(active >= 100.0)


# ============================================================
# 与 PerFrame 交叉验证
# ============================================================

class TestVsPerFrame:

    def test_pair_count_matches(self, hydrogen_bonds, protein_groups):
        u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
        per_frame = HydrogenBondDetectorPerFrame().detect(
            protein_groups, trajectory=u.trajectory,
            tuple_filter=lambda gt: gt[0].molecule != gt[1].molecule)
        assert hydrogen_bonds[0].n_pairs == per_frame[0].n_pairs

    def test_pair_set_matches(self, hydrogen_bonds, protein_groups):
        u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
        per_frame = HydrogenBondDetectorPerFrame().detect(
            protein_groups, trajectory=u.trajectory,
            tuple_filter=lambda gt: gt[0].molecule != gt[1].molecule)
        two = {(g1.group_id, g2.group_id) for g1, g2 in hydrogen_bonds[0].groups}
        per = {(g1.group_id, g2.group_id) for g1, g2 in per_frame[0].groups}
        assert two == per
