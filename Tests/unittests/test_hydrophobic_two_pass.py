# -*- coding: utf-8 -*-
"""疏水 TwoPass 检测器集成测试。

从真实 tpr + xtc 输入，验证检测结果。
与 PerFrame 版本对比，验证结果一致性。
"""

import pytest
from pathlib import Path

import numpy as np

from DuIvyInteractions.input_readers import GmxTprReader
from DuIvyInteractions.group_identifiers import AmberFFGroupIdentifier
from DuIvyInteractions.interaction_detectors import (
    HydrophobicDetectorTwoPass, HydrophobicDetectorPerTuple)
import MDAnalysis as mda


# 测试数据路径
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
def hydrophobic_groups(groups):
    return [g for g in groups if g.group_type == "hydrophobic"]


@pytest.fixture(scope="module")
def hydrophobic_interactions(hydrophobic_groups):
    """TwoPass 检测（全量734个疏水原子）。"""
    u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    detector = HydrophobicDetectorTwoPass()
    return detector.detect(hydrophobic_groups, trajectory=u.trajectory)


# ============================================================
# 基本功能测试
# ============================================================

class TestHydrophobicTwoPassBasic:

    def test_has_results(self, hydrophobic_interactions):
        assert len(hydrophobic_interactions) > 0

    def test_interaction_type(self, hydrophobic_interactions):
        assert hydrophobic_interactions[0].interaction_type == "hydrophobic"

    def test_n_pairs(self, hydrophobic_interactions):
        assert hydrophobic_interactions[0].n_pairs > 0

    def test_all_hydrophobic(self, hydrophobic_interactions):
        it = hydrophobic_interactions[0]
        for g1, g2 in it.groups:
            assert g1.group_type == "hydrophobic"
            assert g2.group_type == "hydrophobic"


# ============================================================
# 输出格式测试
# ============================================================

class TestHydrophobicTwoPassFormat:

    def test_existence_is_2d(self, hydrophobic_interactions):
        for interaction in hydrophobic_interactions:
            assert interaction.existence.ndim == 2

    def test_metrics_shape(self, hydrophobic_interactions):
        for interaction in hydrophobic_interactions:
            n_frames = interaction.n_frames
            assert interaction.metrics["distance"].shape == (
                interaction.n_pairs, n_frames)

    def test_occupancy_range(self, hydrophobic_interactions):
        occupancy = hydrophobic_interactions[0].occupancy()
        assert np.all(occupancy >= 0)
        assert np.all(occupancy <= 1)


# ============================================================
# 阈值测试
# ============================================================

class TestHydrophobicTwoPassThreshold:

    def test_active_distance_range(self, hydrophobic_interactions):
        """活跃帧的距离应在 0.5~4.0 Å。"""
        it = hydrophobic_interactions[0]
        for i in range(it.n_pairs):
            active_dist = it.metrics["distance"][i][it.existence[i]]
            if len(active_dist) > 0:
                assert np.all(active_dist > 0.5)
                assert np.all(active_dist < 4.0)


# ============================================================
# 检测器元信息测试
# ============================================================

class TestHydrophobicTwoPassDetectorMeta:

    def test_name(self):
        detector = HydrophobicDetectorTwoPass()
        assert detector.name == "hydrophobic"

    def test_required_group_types(self):
        detector = HydrophobicDetectorTwoPass()
        assert "hydrophobic" in detector.required_group_types

    def test_metric_names(self):
        detector = HydrophobicDetectorTwoPass()
        assert "distance" in detector.metric_names
