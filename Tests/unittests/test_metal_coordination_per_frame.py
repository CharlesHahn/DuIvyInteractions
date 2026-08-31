# -*- coding: utf-8 -*-
"""金属配位 PerFrame 检测器集成测试。

验证 MetalCoordinationDetectorPerFrame 的结果正确性。
"""

import pytest
from pathlib import Path

import numpy as np
import MDAnalysis as mda

from DuIvyInteractions.input_readers import GmxTprReader
from DuIvyInteractions.group_identifiers import AmberFFGroupIdentifier
from DuIvyInteractions.interaction_detectors import MetalCoordinationDetectorPerFrame


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
def metal_coordination(groups):
    u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    detector = MetalCoordinationDetectorPerFrame()
    return detector.detect(groups, trajectory=u.trajectory)


# ============================================================
# 基本功能
# ============================================================

class TestMetalCoordinationPerFrameBasic:

    def test_has_results(self, metal_coordination):
        assert len(metal_coordination) > 0

    def test_interaction_type(self, metal_coordination):
        assert metal_coordination[0].interaction_type == "metal_coordination"

    def test_metrics_keys(self, metal_coordination):
        assert "distance" in metal_coordination[0].metrics


# ============================================================
# 结果数据
# ============================================================

class TestMetalCoordinationPerFrameResults:

    def test_n_pairs(self, metal_coordination):
        assert metal_coordination[0].n_pairs == 6

    def test_occupancy_range(self, metal_coordination):
        occ = metal_coordination[0].occupancy()
        assert np.all(occ >= 0) and np.all(occ <= 1)

    def test_top_pair_occupancy(self, metal_coordination):
        occ = metal_coordination[0].occupancy()
        assert np.max(occ) == 1.0

    def test_active_distance_range(self, metal_coordination):
        it = metal_coordination[0]
        for i in range(it.n_pairs):
            active = it.metrics["distance"][i][it.existence[i]]
            if len(active) > 0:
                assert np.all(active < 3.0)

    def test_metrics_shape(self, metal_coordination):
        it = metal_coordination[0]
        assert it.metrics["distance"].shape == (it.n_pairs, it.n_frames)
