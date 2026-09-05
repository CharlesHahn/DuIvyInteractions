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


TPR_FILE = Path(__file__).parent.parent / "test_MD_case" / "md.tpr"
XTC_FILE = Path(__file__).parent.parent / "test_MD_case" / "md1ns.xtc"


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
            for key in ["dist_dw", "dist_wa", "theta", "omega"]:
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
            for key in ["dist_dw", "dist_wa", "theta", "omega"]:
                assert it.metrics[key].shape == (it.n_pairs, it.n_frames)


# ============================================================
# 阈值
# ============================================================

class TestWaterBridgePerFrameThreshold:

    def test_active_dist_dw_range(self, water_bridges):
        for it in water_bridges:
            for i in range(it.n_pairs):
                active = it.metrics["dist_dw"][i][it.existence[i]]
                if len(active) > 0:
                    assert np.all(active > 2.5) and np.all(active < 4.1)

    def test_active_dist_wa_range(self, water_bridges):
        for it in water_bridges:
            for i in range(it.n_pairs):
                active = it.metrics["dist_wa"][i][it.existence[i]]
                if len(active) > 0:
                    assert np.all(active > 2.5) and np.all(active < 4.1)

    def test_active_theta_range(self, water_bridges):
        for it in water_bridges:
            for i in range(it.n_pairs):
                active = it.metrics["theta"][i][it.existence[i]]
                if len(active) > 0:
                    assert np.all(active >= 100.0)

    def test_active_omega_range(self, water_bridges):
        for it in water_bridges:
            for i in range(it.n_pairs):
                active = it.metrics["omega"][i][it.existence[i]]
                if len(active) > 0:
                    assert np.all(active >= 71.0) and np.all(active <= 140.0)


# ============================================================
# 结果数据
# ============================================================

class TestWaterBridgePerFrameResults:

    def test_has_results(self, water_bridges):
        assert len(water_bridges) > 0

    def test_n_pairs(self, water_bridges):
        assert water_bridges[0].n_pairs == 2541

    def test_occupancy_range(self, water_bridges):
        occ = water_bridges[0].occupancy()
        assert np.all(occ >= 0) and np.all(occ <= 1)

    def test_top_pair(self, water_bridges):
        it = water_bridges[0]
        occ = it.occupancy()
        top = np.argmax(occ)
        g1, g2, g3 = it.groups[top]
        assert occ[top] > 0.9
        assert g2.group_type == "water"
