# -*- coding: utf-8 -*-
"""水桥 TwoPass 检测器集成测试。

从真实 tpr + xtc 输入，验证检测结果。
"""

import pytest
from pathlib import Path

import numpy as np

from DuIvyInteractions.system_readers import GmxTprReader
from DuIvyInteractions.group_identifiers import AmberFFGroupIdentifier
from DuIvyInteractions.interaction_detectors import WaterBridgeDetectorTwoPass
import MDAnalysis as mda


TPR_FILE = Path(__file__).parent.parent / "test_MD_case_amber" / "md.tpr"
XTC_FILE = Path(__file__).parent.parent / "test_MD_case_amber" / "md1ns.xtc"


@pytest.fixture(scope="module")
def system_data():
    reader = GmxTprReader()
    return reader.read(str(TPR_FILE))


@pytest.fixture(scope="module")
def groups(system_data):
    identifier = AmberFFGroupIdentifier()
    return identifier.identify(system_data)


@pytest.fixture(scope="module")
def water_bridges(groups):
    u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    detector = WaterBridgeDetectorTwoPass()
    return detector.detect(groups, trajectory=u.trajectory)


# ============================================================
# 基本功能测试
# ============================================================

class TestWaterBridgeTwoPassBasic:

    def test_returns_list(self, water_bridges):
        assert isinstance(water_bridges, list)

    def test_interaction_type(self, water_bridges):
        if water_bridges:
            assert water_bridges[0].interaction_type == "water_bridge"

    def test_metrics_keys(self, water_bridges):
        for interaction in water_bridges:
            for key in ["dist_h_ow", "angle_dh_ow", "dist_d_ow",
                        "dist_hw_a", "angle_ow_hw_a", "dist_ow_a", "angle_h_ow_a"]:
                assert key in interaction.metrics


# ============================================================
# 输出格式测试
# ============================================================

class TestWaterBridgeTwoPassFormat:

    def test_groups_are_triples(self, water_bridges):
        """每个 groups 元素应为 (donor, water, acceptor) 三元组。"""
        for interaction in water_bridges:
            for triple in interaction.groups:
                assert len(triple) == 3
                donor, water, acceptor = triple
                assert donor.group_type == "H_donor"
                assert water.group_type == "water"
                assert acceptor.group_type == "H_acceptor"

    def test_existence_is_2d(self, water_bridges):
        for interaction in water_bridges:
            assert interaction.existence.ndim == 2

    def test_metrics_shape(self, water_bridges):
        for interaction in water_bridges:
            n_frames = interaction.n_frames
            for key in ["dist_h_ow", "angle_dh_ow", "dist_d_ow",
                        "dist_hw_a", "angle_ow_hw_a", "dist_ow_a", "angle_h_ow_a"]:
                assert interaction.metrics[key].shape == (
                    interaction.n_pairs, n_frames)


# ============================================================
# 结果数据测试
# ============================================================

class TestWaterBridgeTwoPassResults:

    def test_has_results(self, water_bridges):
        assert len(water_bridges) > 0

    def test_n_pairs_candidates(self, water_bridges):
        """候选三元组总数（TwoPass 保留全部候选，含从未存在者）。

        新判据（显式 H + Baker-Hubbard + 水在中间 + 三层配额）：3130 个候选
        （原 PLIP 宽松判据为 14598）。
        """
        assert water_bridges[0].n_pairs == 3130

    def test_active_pairs(self, water_bridges):
        """有存在帧的对数（排除从未存在的候选）。"""
        it = water_bridges[0]
        active = np.any(it.existence.astype(bool), axis=1).sum()
        assert active == 2913

    def test_occupancy_range(self, water_bridges):
        occ = water_bridges[0].occupancy()
        assert np.all(occ >= 0) and np.all(occ <= 1)

    def test_top_pair(self, water_bridges):
        it = water_bridges[0]
        occ = it.occupancy()
        top = np.argmax(occ)
        g1, g2, g3 = it.groups[top]
        assert occ[top] > 0.8
        assert g2.group_type == "water"

    def test_top_pair_identity(self, water_bridges):
        """top pair（LEU 供体 → 水 → PHE 受体），新判据下帧占有率 0.861。"""
        it = water_bridges[0]
        occ = it.occupancy()
        top = np.argmax(occ)
        g1, g2, g3 = it.groups[top]
        assert round(float(occ[top]), 3) == 0.861
        assert g1.residue_name == "LEU" and g1.group_type == "H_donor"
        assert g2.residue_name == "SOL" and g2.group_type == "water"
        assert g3.residue_name == "PHE" and g3.group_type == "H_acceptor"


# ============================================================
# 阈值测试
# ============================================================

class TestWaterBridgeTwoPassThreshold:

    def test_active_dist_h_ow_range(self, water_bridges):
        """活跃帧：供体 H···水 O 应在 (0.5, 2.5] Å（Baker-Hubbard）。"""
        for interaction in water_bridges:
            for i in range(interaction.n_pairs):
                active = interaction.metrics["dist_h_ow"][i][interaction.existence[i]]
                if len(active) > 0:
                    assert np.all(active > 0.5)
                    assert np.all(active <= 2.5)

    def test_active_dist_hw_a_range(self, water_bridges):
        """活跃帧：水 H···受体应在 (0.5, 2.5] Å。"""
        for interaction in water_bridges:
            for i in range(interaction.n_pairs):
                active = interaction.metrics["dist_hw_a"][i][interaction.existence[i]]
                if len(active) > 0:
                    assert np.all(active > 0.5)
                    assert np.all(active <= 2.5)

    def test_active_angle_dh_ow_range(self, water_bridges):
        """活跃帧：D-H···Ow 角度应 ≥ 120°。"""
        for interaction in water_bridges:
            for i in range(interaction.n_pairs):
                active = interaction.metrics["angle_dh_ow"][i][interaction.existence[i]]
                if len(active) > 0:
                    assert np.all(active >= 120.0)

    def test_active_angle_ow_hw_a_range(self, water_bridges):
        """活跃帧：Ow-Hw···A 角度应 ≥ 120°。"""
        for interaction in water_bridges:
            for i in range(interaction.n_pairs):
                active = interaction.metrics["angle_ow_hw_a"][i][interaction.existence[i]]
                if len(active) > 0:
                    assert np.all(active >= 120.0)

    def test_active_omega_range(self, water_bridges):
        """活跃帧：H(供体)-Ow···A 角应在 71°~140°（"水在中间"约束）。"""
        for interaction in water_bridges:
            for i in range(interaction.n_pairs):
                active = interaction.metrics["angle_h_ow_a"][i][interaction.existence[i]]
                if len(active) > 0:
                    assert np.all(active >= 71.0)
                    assert np.all(active <= 140.0)


# ============================================================
# 检测器元信息测试
# ============================================================

class TestWaterBridgeTwoPassDetectorMeta:

    def test_name(self):
        detector = WaterBridgeDetectorTwoPass()
        assert detector.name == "water_bridge"

    def test_required_group_types(self):
        detector = WaterBridgeDetectorTwoPass()
        assert "H_donor" in detector.required_group_types
        assert "water" in detector.required_group_types
        assert "H_acceptor" in detector.required_group_types

    def test_metric_names(self):
        detector = WaterBridgeDetectorTwoPass()
        names = detector.metric_names
        for key in ["dist_h_ow", "angle_dh_ow", "dist_d_ow",
                        "dist_hw_a", "angle_ow_hw_a", "dist_ow_a", "angle_h_ow_a"]:
            assert key in names
