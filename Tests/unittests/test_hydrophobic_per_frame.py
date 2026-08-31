# -*- coding: utf-8 -*-
"""疏水 PerFrame 检测器集成测试。

验证 HydrophobicDetectorPerFrame 的结果正确性。
"""

import pytest
from pathlib import Path

import numpy as np
import MDAnalysis as mda

from DuIvyInteractions.input_readers import GmxTprReader
from DuIvyInteractions.group_identifiers import AmberFFGroupIdentifier
from DuIvyInteractions.interaction_detectors import HydrophobicDetectorPerFrame


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
def hydrophobic(groups):
    u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    detector = HydrophobicDetectorPerFrame()
    return detector.detect(groups, trajectory=u.trajectory)


# ============================================================
# 基本功能
# ============================================================

class TestHydrophobicPerFrameBasic:

    def test_has_results(self, hydrophobic):
        assert len(hydrophobic) > 0

    def test_interaction_type(self, hydrophobic):
        assert hydrophobic[0].interaction_type == "hydrophobic"

    def test_metrics_keys(self, hydrophobic):
        assert "distance" in hydrophobic[0].metrics


# ============================================================
# 结果数据
# ============================================================

class TestHydrophobicPerFrameResults:

    def test_n_pairs(self, hydrophobic):
        assert hydrophobic[0].n_pairs == 705

    def test_occupancy_range(self, hydrophobic):
        occ = hydrophobic[0].occupancy()
        assert np.all(occ >= 0) and np.all(occ <= 1)

    def test_active_distance_range(self, hydrophobic):
        it = hydrophobic[0]
        for i in range(it.n_pairs):
            active = it.metrics["distance"][i][it.existence[i]]
            if len(active) > 0:
                assert np.all(active > 0.5) and np.all(active < 4.0)

    def test_metrics_shape(self, hydrophobic):
        it = hydrophobic[0]
        assert it.metrics["distance"].shape == (it.n_pairs, it.n_frames)

    def test_pairs_are_different(self, hydrophobic):
        """每个对的两个基团不应相同。"""
        it = hydrophobic[0]
        for g1, g2 in it.groups:
            assert g1.group_id != g2.group_id
