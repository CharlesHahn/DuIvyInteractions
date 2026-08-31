# -*- coding: utf-8 -*-
"""盐桥 TwoPass 检测器集成测试：从真实 tpr + xtc 输入，验证检测结果。

与 PerFrame 版本对比，验证结果一致性。
"""

import pytest
import time
from pathlib import Path

import numpy as np

from DuIvyInteractions.input_readers import GmxTprReader
from DuIvyInteractions.group_identifiers import AmberFFGroupIdentifier
from DuIvyInteractions.interaction_detectors import (
    SaltBridgeDetectorTwoPass, SaltBridgeDetectorPerFrame)
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
def filtered_groups(groups):
    """只保留正电和负电基团。"""
    return [g for g in groups
            if g.group_type in ("charged_positive", "charged_negative")]


@pytest.fixture(scope="module")
def salt_bridges(filtered_groups):
    """运行 TwoPass 盐桥检测，返回 Interaction 列表。"""
    u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    detector = SaltBridgeDetectorTwoPass()
    return detector.detect(filtered_groups, trajectory=u.trajectory)


@pytest.fixture(scope="module")
def salt_bridges_per_frame(filtered_groups):
    """运行 PerFrame 盐桥检测，用于对比。"""
    u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    detector = SaltBridgeDetectorPerFrame()
    return detector.detect(filtered_groups, trajectory=u.trajectory)


# ============================================================
# 基本功能测试
# ============================================================

class TestSaltBridgeTwoPassBasic:
    """验证 TwoPass 盐桥检测的基本功能。"""

    def test_has_results(self, salt_bridges):
        assert len(salt_bridges) > 0

    def test_interaction_type(self, salt_bridges):
        assert salt_bridges[0].interaction_type == "salt_bridge"

    def test_n_pairs(self, salt_bridges):
        """应检测到 47 对盐桥（与 PerFrame 一致）。"""
        assert salt_bridges[0].n_pairs == 47

    def test_all_opposite_charge(self, salt_bridges):
        it = salt_bridges[0]
        for g1, g2 in it.groups:
            assert g1.group_type == "charged_positive"
            assert g2.group_type == "charged_negative"


# ============================================================
# 输出格式测试
# ============================================================

class TestSaltBridgeTwoPassFormat:
    """验证输出格式符合设计。"""

    def test_existence_is_2d(self, salt_bridges):
        for interaction in salt_bridges:
            assert interaction.existence.ndim == 2

    def test_metrics_shape(self, salt_bridges):
        for interaction in salt_bridges:
            n_frames = interaction.n_frames
            assert interaction.metrics["distance"].shape == (interaction.n_pairs, n_frames)

    def test_occupancy_range(self, salt_bridges):
        occupancy = salt_bridges[0].occupancy()
        assert np.all(occupancy >= 0)
        assert np.all(occupancy <= 1)


# ============================================================
# 阈值测试
# ============================================================

class TestSaltBridgeTwoPassThreshold:
    """验证阈值判定正确。"""

    def test_active_distance_range(self, salt_bridges):
        """活跃帧的电荷中心距离应 ≤ 5.5 Å。"""
        it = salt_bridges[0]
        for i in range(it.n_pairs):
            active_dist = it.metrics["distance"][i][it.existence[i]]
            if len(active_dist) > 0:
                assert np.all(active_dist <= 5.5)


# ============================================================
# 与 PerFrame 对比测试
# ============================================================

class TestSaltBridgeTwoPassVsPerFrame:
    """验证 TwoPass 与 PerFrame 结果一致性。"""

    def test_pair_count_matches(self, salt_bridges, salt_bridges_per_frame):
        """两策略应检测到相同数量的盐桥对。"""
        assert salt_bridges[0].n_pairs == salt_bridges_per_frame[0].n_pairs

    def test_pair_set_matches(self, salt_bridges, salt_bridges_per_frame):
        """两策略应检测到完全相同的盐桥对。"""
        two_pairs = {(g1.group_id, g2.group_id)
                     for g1, g2 in salt_bridges[0].groups}
        per_pairs = {(g1.group_id, g2.group_id)
                     for g1, g2 in salt_bridges_per_frame[0].groups}
        assert two_pairs == per_pairs

    def test_occupancy_matches(self, salt_bridges, salt_bridges_per_frame):
        """两策略的占位率应一致。"""
        two_occ = salt_bridges[0].occupancy()
        per_occ = salt_bridges_per_frame[0].occupancy()
        # 排序后比较（pair 顺序可能不同）
        np.testing.assert_array_almost_equal(
            np.sort(two_occ), np.sort(per_occ), decimal=3)


# ============================================================
# 检测器元信息测试
# ============================================================

class TestSaltBridgeTwoPassDetectorMeta:
    """验证检测器的元信息。"""

    def test_name(self):
        detector = SaltBridgeDetectorTwoPass()
        assert detector.name == "salt_bridge"

    def test_required_group_types(self):
        detector = SaltBridgeDetectorTwoPass()
        assert "charged_positive" in detector.required_group_types
        assert "charged_negative" in detector.required_group_types

    def test_metric_names(self):
        detector = SaltBridgeDetectorTwoPass()
        assert "distance" in detector.metric_names
