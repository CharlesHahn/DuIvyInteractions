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


TPR_FILE = Path(__file__).parent.parent / "test_MD_case" / "md.tpr"
XTC_FILE = Path(__file__).parent.parent / "test_MD_case" / "md1ns.xtc"


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
            assert "dist_dw" in interaction.metrics
            assert "dist_wa" in interaction.metrics
            assert "theta" in interaction.metrics
            assert "omega" in interaction.metrics


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
            for key in ["dist_dw", "dist_wa", "theta", "omega"]:
                assert interaction.metrics[key].shape == (
                    interaction.n_pairs, n_frames)


# ============================================================
# 阈值测试
# ============================================================

class TestWaterBridgeTwoPassThreshold:

    def test_active_dist_dw_range(self, water_bridges):
        """活跃帧的 dist_dw 应 < 4.1 Å。"""
        for interaction in water_bridges:
            for i in range(interaction.n_pairs):
                active = interaction.metrics["dist_dw"][i][interaction.existence[i]]
                if len(active) > 0:
                    assert np.all(active < 4.1)

    def test_active_dist_wa_range(self, water_bridges):
        """活跃帧的 dist_wa 应 < 4.1 Å。"""
        for interaction in water_bridges:
            for i in range(interaction.n_pairs):
                active = interaction.metrics["dist_wa"][i][interaction.existence[i]]
                if len(active) > 0:
                    assert np.all(active < 4.1)

    def test_active_theta_range(self, water_bridges):
        """活跃帧的 theta 应 ≥ 100°。"""
        for interaction in water_bridges:
            for i in range(interaction.n_pairs):
                active = interaction.metrics["theta"][i][interaction.existence[i]]
                if len(active) > 0:
                    assert np.all(active >= 100.0)

    def test_active_omega_range(self, water_bridges):
        """活跃帧的 omega 应在 71°~140°。"""
        for interaction in water_bridges:
            for i in range(interaction.n_pairs):
                active = interaction.metrics["omega"][i][interaction.existence[i]]
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
        assert "dist_dw" in names
        assert "dist_wa" in names
        assert "theta" in names
        assert "omega" in names
