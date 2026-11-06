# -*- coding: utf-8 -*-
"""水桥 PerFrame 检测器集成测试。

验证 WaterBridgeDetectorPerFrame 的结果正确性。
"""

import pytest
from pathlib import Path

import numpy as np
import MDAnalysis as mda

from DuIvyInteractions.system_readers import GmxTprReader
from DuIvyInteractions.group_identifiers import AmberFFGroupIdentifier
from DuIvyInteractions.interaction_detectors import WaterBridgeDetectorPerFrame


TPR_FILE = Path(__file__).parent.parent / "test_MD_case_amber" / "md.tpr"
XTC_FILE = Path(__file__).parent.parent / "test_MD_case_amber" / "md1ns.xtc"


@pytest.fixture(scope="module")
def system_data():
    reader = GmxTprReader()
    return reader.read(str(TPR_FILE))


@pytest.fixture(scope="module")
def groups(system_data):
    return AmberFFGroupIdentifier().identify(system_data)


@pytest.fixture(scope="module")
def water_bridges(groups):
    u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    detector = WaterBridgeDetectorPerFrame()
    return detector.detect(groups, trajectory=u.trajectory)


# ============================================================
# 基本功能
# ============================================================

class TestWaterBridgePerFrameBasic:

    def test_returns_list(self, water_bridges):
        assert isinstance(water_bridges, list)

    def test_interaction_type(self, water_bridges):
        for it in water_bridges:
            assert it.interaction_type == "water_bridge"

    def test_metrics_keys(self, water_bridges):
        for it in water_bridges:
            for key in ["dist_h_ow", "angle_dh_ow", "dist_d_ow",
                        "dist_hw_a", "angle_ow_hw_a", "dist_ow_a", "angle_h_ow_a"]:
                assert key in it.metrics


# ============================================================
# 输出格式
# ============================================================

class TestWaterBridgePerFrameFormat:

    def test_groups_are_triples(self, water_bridges):
        for it in water_bridges:
            for triple in it.groups:
                assert len(triple) == 3
                donor, water, acceptor = triple
                assert donor.group_type == "H_donor"
                assert water.group_type == "water"
                assert acceptor.group_type == "H_acceptor"

    def test_existence_is_2d(self, water_bridges):
        for it in water_bridges:
            assert it.existence.ndim == 2

    def test_metrics_shape(self, water_bridges):
        for it in water_bridges:
            for key in ["dist_h_ow", "angle_dh_ow", "dist_d_ow",
                        "dist_hw_a", "angle_ow_hw_a", "dist_ow_a", "angle_h_ow_a"]:
                assert it.metrics[key].shape == (it.n_pairs, it.n_frames)


# ============================================================
# 阈值
# ============================================================

class TestWaterBridgePerFrameThreshold:

    def test_active_dist_h_ow_range(self, water_bridges):
        """活跃帧：供体 H···水 O 应在 (0.5, 2.5] Å（Baker-Hubbard）。"""
        for it in water_bridges:
            for i in range(it.n_pairs):
                active = it.metrics["dist_h_ow"][i][it.existence[i]]
                if len(active) > 0:
                    assert np.all(active > 0.5) and np.all(active <= 2.5)

    def test_active_dist_hw_a_range(self, water_bridges):
        """活跃帧：水 H···受体应在 (0.5, 2.5] Å。"""
        for it in water_bridges:
            for i in range(it.n_pairs):
                active = it.metrics["dist_hw_a"][i][it.existence[i]]
                if len(active) > 0:
                    assert np.all(active > 0.5) and np.all(active <= 2.5)

    def test_active_angle_dh_ow_range(self, water_bridges):
        """活跃帧：D-H···Ow 角度应 ≥ 120°。"""
        for it in water_bridges:
            for i in range(it.n_pairs):
                active = it.metrics["angle_dh_ow"][i][it.existence[i]]
                if len(active) > 0:
                    assert np.all(active >= 120.0)

    def test_active_angle_ow_hw_a_range(self, water_bridges):
        """活跃帧：Ow-Hw···A 角度应 ≥ 120°。"""
        for it in water_bridges:
            for i in range(it.n_pairs):
                active = it.metrics["angle_ow_hw_a"][i][it.existence[i]]
                if len(active) > 0:
                    assert np.all(active >= 120.0)

    def test_active_omega_range(self, water_bridges):
        """活跃帧：H(供体)-Ow···A 角应在 71°~140°（"水在中间"约束）。"""
        for it in water_bridges:
            for i in range(it.n_pairs):
                active = it.metrics["angle_h_ow_a"][i][it.existence[i]]
                if len(active) > 0:
                    assert np.all(active >= 71.0)
                    assert np.all(active <= 140.0)


# ============================================================
# 结果数据
# ============================================================

class TestWaterBridgePerFrameResults:

    def test_has_results(self, water_bridges):
        assert len(water_bridges) > 0

    def test_n_pairs(self, water_bridges):
        # 新判据（显式 H + Baker-Hubbard + 水在中间 + 三层配额）：86 对（原 PLIP 1698）。
        # 注：PerFrame 候选只基于第一帧（设计原则），故全轨迹对数少于 TwoPass 的累积候选。
        assert water_bridges[0].n_pairs == 86

    def test_occupancy_range(self, water_bridges):
        occ = water_bridges[0].occupancy()
        assert np.all(occ >= 0) and np.all(occ <= 1)

    def test_top_pair(self, water_bridges):
        it = water_bridges[0]
        occ = it.occupancy()
        top = np.argmax(occ)
        g1, g2, g3 = it.groups[top]
        assert occ[top] > 0.7
        assert g2.group_type == "water"
