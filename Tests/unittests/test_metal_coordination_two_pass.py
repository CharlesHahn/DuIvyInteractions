# -*- coding: utf-8 -*-
"""金属配位 TwoPass 检测器集成测试。

从真实 tpr + xtc 输入，验证检测结果。
与 PerFrame 版本对比，验证结果一致性。
"""

import pytest
from pathlib import Path

import numpy as np

from DuIvyInteractions.system_readers import GmxTprReader
from DuIvyInteractions.group_identifiers import AmberFFGroupIdentifier
from DuIvyInteractions.interaction_detectors import (
    MetalCoordinationDetectorTwoPass, MetalCoordinationDetectorPerFrame)
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
def relevant_groups(groups):
    return [g for g in groups
            if g.group_type in ("metal", "metal_binding")]


@pytest.fixture(scope="module")
def metal_coordination(relevant_groups):
    u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    detector = MetalCoordinationDetectorTwoPass()
    return detector.detect(relevant_groups, trajectory=u.trajectory)


@pytest.fixture(scope="module")
def metal_coordination_per_frame(relevant_groups):
    u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    detector = MetalCoordinationDetectorPerFrame()
    return detector.detect(relevant_groups, trajectory=u.trajectory)


# ============================================================
# 基本功能测试
# ============================================================

class TestMetalCoordinationTwoPassBasic:

    def test_has_results(self, metal_coordination):
        assert len(metal_coordination) > 0

    def test_interaction_type(self, metal_coordination):
        assert metal_coordination[0].interaction_type == "metal_coordination"

    def test_n_pairs(self, metal_coordination):
        assert metal_coordination[0].n_pairs > 0

    def test_all_metal_binding(self, metal_coordination):
        it = metal_coordination[0]
        for g1, g2 in it.groups:
            assert g1.group_type == "metal"
            assert g2.group_type == "metal_binding"


# ============================================================
# 输出格式测试
# ============================================================

class TestMetalCoordinationTwoPassFormat:

    def test_existence_is_2d(self, metal_coordination):
        for interaction in metal_coordination:
            assert interaction.existence.ndim == 2

    def test_metrics_shape(self, metal_coordination):
        for interaction in metal_coordination:
            n_frames = interaction.n_frames
            assert interaction.metrics["distance"].shape == (
                interaction.n_pairs, n_frames)


# ============================================================
# 阈值测试
# ============================================================

class TestMetalCoordinationTwoPassThreshold:

    def test_active_distance_range(self, metal_coordination):
        """活跃帧的距离应 < 3.0 Å。"""
        it = metal_coordination[0]
        for i in range(it.n_pairs):
            active_dist = it.metrics["distance"][i][it.existence[i]]
            if len(active_dist) > 0:
                assert np.all(active_dist < 3.0)


# ============================================================
# 与 PerFrame 对比测试
# ============================================================

class TestMetalCoordinationTwoPassVsPerFrame:

    def test_two_pass_finds_at_least_as_many(self, metal_coordination, metal_coordination_per_frame):
        """TwoPass 应找到至少与 PerFrame 一样多的 pair（零遗漏）。"""
        assert metal_coordination[0].n_pairs >= metal_coordination_per_frame[0].n_pairs

    def test_per_frame_pairs_are_subset(self, metal_coordination, metal_coordination_per_frame):
        """PerFrame 的 pair 应是 TwoPass 的子集。"""
        two_pairs = {(g1.group_id, g2.group_id)
                     for g1, g2 in metal_coordination[0].groups}
        per_pairs = {(g1.group_id, g2.group_id)
                     for g1, g2 in metal_coordination_per_frame[0].groups}
        assert per_pairs.issubset(two_pairs)


# ============================================================
# 检测器元信息测试
# ============================================================

class TestMetalCoordinationTwoPassDetectorMeta:

    def test_name(self):
        assert MetalCoordinationDetectorTwoPass().name == "metal_coordination"

    def test_required_group_types(self):
        det = MetalCoordinationDetectorTwoPass()
        assert "metal" in det.required_group_types
        assert "metal_binding" in det.required_group_types

    def test_metric_names(self):
        assert "distance" in MetalCoordinationDetectorTwoPass().metric_names
