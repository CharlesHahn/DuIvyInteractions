# -*- coding: utf-8 -*-
"""水桥检测器集成测试：从真实 tpr + xtc 输入，验证完整的检测结果。

使用 D927 体系的真实数据作为测试基准。
"""

import pytest
from pathlib import Path

import numpy as np

from DuIvyInteractions.input_readers import GmxTprReader
from DuIvyInteractions.group_identifiers import AmberFFGroupIdentifier
from DuIvyInteractions.interaction_detectors import WaterBridgeDetector
import MDAnalysis as mda


# 测试数据路径
TPR_FILE = Path(__file__).parent.parent / "test_MD_case" / "md.tpr"
XTC_FILE = Path(__file__).parent.parent / "test_MD_case" / "md1ns.xtc"


@pytest.fixture(scope="module")
def system_data():
    """读取 tpr 文件，返回 SystemData。"""
    reader = GmxTprReader()
    return reader.read(str(TPR_FILE))


@pytest.fixture(scope="module")
def groups(system_data):
    """运行基团识别，返回 Group 列表。"""
    identifier = AmberFFGroupIdentifier()
    return identifier.identify(system_data)


@pytest.fixture(scope="module")
def water_bridges(groups):
    """运行水桥检测，返回 Interaction 列表。"""
    u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    detector = WaterBridgeDetector()
    return detector.detect(groups, trajectory=u.trajectory)


# ============================================================
# 基本功能测试
# ============================================================

class TestWaterBridgeBasic:
    """验证水桥检测的基本功能。"""

    def test_returns_list(self, water_bridges):
        assert isinstance(water_bridges, list)

    def test_interaction_type(self, water_bridges):
        for interaction in water_bridges:
            assert interaction.interaction_type == "water_bridge"

    def test_metrics_keys(self, water_bridges):
        for interaction in water_bridges:
            assert "dist_dw" in interaction.metrics
            assert "dist_wa" in interaction.metrics
            assert "theta" in interaction.metrics
            assert "omega" in interaction.metrics


# ============================================================
# 输出格式测试
# ============================================================

class TestWaterBridgeFormat:
    """验证输出格式符合设计。"""

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
        """existence 应为二维数组 (n_pairs, n_frames)。"""
        for interaction in water_bridges:
            assert interaction.existence.ndim == 2

    def test_metrics_shape(self, water_bridges):
        """metrics 各数组应为 (n_pairs, n_frames)。"""
        for interaction in water_bridges:
            n_frames = interaction.n_frames
            for key in ["dist_dw", "dist_wa", "theta", "omega"]:
                assert interaction.metrics[key].shape == (interaction.n_pairs, n_frames)


# ============================================================
# 阈值测试
# ============================================================

class TestWaterBridgeThreshold:
    """验证阈值判定正确。"""

    def test_active_dist_dw_range(self, water_bridges):
        """活跃帧的 dist_dw 应在 2.5~4.1 Å。"""
        for interaction in water_bridges:
            for i in range(interaction.n_pairs):
                active = interaction.metrics["dist_dw"][i][interaction.existence[i]]
                if len(active) > 0:
                    assert np.all(active > 2.5)
                    assert np.all(active < 4.1)

    def test_active_dist_wa_range(self, water_bridges):
        """活跃帧的 dist_wa 应在 2.5~4.1 Å。"""
        for interaction in water_bridges:
            for i in range(interaction.n_pairs):
                active = interaction.metrics["dist_wa"][i][interaction.existence[i]]
                if len(active) > 0:
                    assert np.all(active > 2.5)
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

class TestWaterBridgeDetectorMeta:
    """验证检测器的元信息。"""

    def test_name(self):
        detector = WaterBridgeDetector()
        assert detector.name == "water_bridge"

    def test_required_group_types(self):
        detector = WaterBridgeDetector()
        assert "H_donor" in detector.required_group_types
        assert "water" in detector.required_group_types
        assert "H_acceptor" in detector.required_group_types

    def test_metric_names(self):
        detector = WaterBridgeDetector()
        names = detector.metric_names
        assert "dist_dw" in names
        assert "dist_wa" in names
        assert "theta" in names
        assert "omega" in names
