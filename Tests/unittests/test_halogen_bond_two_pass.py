# -*- coding: utf-8 -*-
"""卤键 TwoPass 检测器集成测试。

从真实 tpr + xtc 输入，验证检测结果。
与 PerFrame 版本对比，验证结果一致性。
"""

import pytest
from pathlib import Path

import numpy as np

from DuIvyInteractions.input_readers import GmxTprReader
from DuIvyInteractions.group_identifiers import AmberFFGroupIdentifier
from DuIvyInteractions.interaction_detectors import (
    HalogenBondDetectorTwoPass, HalogenBondDetectorPerFrame)
import MDAnalysis as mda


TPR_FILE = Path(__file__).parent.parent / "test_MD_case" / "md.tpr"
XTC_FILE = Path(__file__).parent.parent / "test_MD_case" / "md1ns.xtc"

MOL_D927 = "seg_1_D927"


@pytest.fixture(scope="module")
def system_data():
    reader = GmxTprReader()
    return reader.read(str(TPR_FILE))


@pytest.fixture(scope="module")
def groups(system_data):
    identifier = AmberFFGroupIdentifier()
    return identifier.identify(system_data)


@pytest.fixture(scope="module")
def relevant_groups(groups):
    return [g for g in groups
            if g.group_type in ("halogen_donor", "halogen_acceptor")]


@pytest.fixture(scope="module")
def halogen_bonds(relevant_groups):
    u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    detector = HalogenBondDetectorTwoPass()
    return detector.detect(relevant_groups, trajectory=u.trajectory)


@pytest.fixture(scope="module")
def halogen_bonds_per_frame(relevant_groups):
    u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    detector = HalogenBondDetectorPerFrame()
    return detector.detect(relevant_groups, trajectory=u.trajectory)


# ============================================================
# 基本功能测试
# ============================================================

class TestHalogenBondTwoPassBasic:

    def test_has_results(self, halogen_bonds):
        assert len(halogen_bonds) > 0

    def test_interaction_type(self, halogen_bonds):
        assert halogen_bonds[0].interaction_type == "halogen_bond"

    def test_n_pairs(self, halogen_bonds):
        assert halogen_bonds[0].n_pairs > 0

    def test_metrics_keys(self, halogen_bonds):
        for interaction in halogen_bonds:
            assert "distance" in interaction.metrics
            assert "don_angle" in interaction.metrics
            assert "acc_angle" in interaction.metrics


# ============================================================
# 输出格式测试
# ============================================================

class TestHalogenBondTwoPassFormat:

    def test_groups_are_donor_acceptor(self, halogen_bonds):
        it = halogen_bonds[0]
        for g1, g2 in it.groups:
            assert g1.group_type == "halogen_donor"
            assert g2.group_type == "halogen_acceptor"

    def test_existence_is_2d(self, halogen_bonds):
        for interaction in halogen_bonds:
            assert interaction.existence.ndim == 2

    def test_metrics_shape(self, halogen_bonds):
        for interaction in halogen_bonds:
            n_frames = interaction.n_frames
            for key in ["distance", "don_angle", "acc_angle"]:
                assert interaction.metrics[key].shape == (
                    interaction.n_pairs, n_frames)


# ============================================================
# 阈值测试
# ============================================================

class TestHalogenBondTwoPassThreshold:

    def test_active_distance_range(self, halogen_bonds):
        it = halogen_bonds[0]
        for i in range(it.n_pairs):
            active = it.metrics["distance"][i][it.existence[i]]
            if len(active) > 0:
                assert np.all(active <= 4.0)

    def test_active_don_angle_range(self, halogen_bonds):
        it = halogen_bonds[0]
        for i in range(it.n_pairs):
            active = it.metrics["don_angle"][i][it.existence[i]]
            if len(active) > 0:
                assert np.all(active >= 135.0)
                assert np.all(active <= 195.0)

    def test_active_acc_angle_range(self, halogen_bonds):
        it = halogen_bonds[0]
        for i in range(it.n_pairs):
            active = it.metrics["acc_angle"][i][it.existence[i]]
            if len(active) > 0:
                assert np.all(active >= 90.0)
                assert np.all(active <= 150.0)


# ============================================================
# 与 PerFrame 对比测试
# ============================================================

class TestHalogenBondTwoPassVsPerFrame:

    def test_pair_count_matches(self, halogen_bonds, halogen_bonds_per_frame):
        assert halogen_bonds[0].n_pairs == halogen_bonds_per_frame[0].n_pairs

    def test_pair_set_matches(self, halogen_bonds, halogen_bonds_per_frame):
        two_pairs = {(g1.group_id, g2.group_id)
                     for g1, g2 in halogen_bonds[0].groups}
        per_pairs = {(g1.group_id, g2.group_id)
                     for g1, g2 in halogen_bonds_per_frame[0].groups}
        assert two_pairs == per_pairs

    def test_occupancy_matches(self, halogen_bonds, halogen_bonds_per_frame):
        two_occ = halogen_bonds[0].occupancy()
        per_occ = halogen_bonds_per_frame[0].occupancy()
        np.testing.assert_array_almost_equal(
            np.sort(two_occ), np.sort(per_occ), decimal=3)


# ============================================================
# 检测器元信息测试
# ============================================================

class TestHalogenBondTwoPassDetectorMeta:

    def test_name(self):
        assert HalogenBondDetectorTwoPass().name == "halogen_bond"

    def test_required_group_types(self):
        det = HalogenBondDetectorTwoPass()
        assert "halogen_donor" in det.required_group_types
        assert "halogen_acceptor" in det.required_group_types

    def test_metric_names(self):
        det = HalogenBondDetectorTwoPass()
        assert "distance" in det.metric_names
        assert "don_angle" in det.metric_names
        assert "acc_angle" in det.metric_names
