# -*- coding: utf-8 -*-
"""金属配位检测器集成测试：从真实 tpr + xtc 输入，验证完整的检测结果。

使用 D927 体系的真实数据作为测试基准。
"""

import pytest
from pathlib import Path

import numpy as np

from DuIvyInteractions.input_readers import GmxTprReader
from DuIvyInteractions.group_identifiers import AmberFFGroupIdentifier
from DuIvyInteractions.interaction_detectors import MetalCoordinationDetector
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
def metal_coordination(groups):
    """运行金属配位检测，返回 Interaction 列表。"""
    u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    detector = MetalCoordinationDetector()
    return detector.detect(groups, trajectory=u.trajectory)


# ============================================================
# 基本功能测试
# ============================================================

class TestMetalCoordinationBasic:
    """验证金属配位检测的基本功能。"""

    def test_returns_list(self, metal_coordination):
        assert isinstance(metal_coordination, list)

    def test_has_results(self, metal_coordination):
        """应检测到金属配位。"""
        assert len(metal_coordination) > 0

    def test_interaction_type(self, metal_coordination):
        for interaction in metal_coordination:
            assert interaction.interaction_type == "metal_coordination"

    def test_metrics_keys(self, metal_coordination):
        for interaction in metal_coordination:
            assert "distance" in interaction.metrics


# ============================================================
# 输出格式测试
# ============================================================

class TestMetalCoordinationFormat:
    """验证输出格式符合设计。"""

    def test_groups_are_metal_atom_pairs(self, metal_coordination):
        """每个 groups 元素应为 (metal, coordination_atom) 元组。"""
        for interaction in metal_coordination:
            for pair in interaction.groups:
                assert len(pair) == 2
                metal, atom = pair
                assert metal.group_type == "metal"
                assert atom.group_type == "metal_binding"

    def test_existence_is_2d(self, metal_coordination):
        """existence 应为二维数组 (n_pairs, n_frames)。"""
        for interaction in metal_coordination:
            assert interaction.existence.ndim == 2

    def test_metrics_shape(self, metal_coordination):
        """metrics['distance'] 应为 (n_pairs, n_frames)。"""
        for interaction in metal_coordination:
            dist = interaction.metrics["distance"]
            assert dist.shape == (interaction.n_pairs, interaction.n_frames)


# ============================================================
# 阈值测试
# ============================================================

class TestMetalCoordinationThreshold:
    """验证阈值判定正确。"""

    def test_active_distance_range(self, metal_coordination):
        """活跃帧的距离应 < 3.0 Å。"""
        for interaction in metal_coordination:
            for i in range(interaction.n_pairs):
                active_dist = interaction.metrics["distance"][i][interaction.existence[i]]
                if len(active_dist) > 0:
                    assert np.all(active_dist < 3.0)


# ============================================================
# 检测器元信息测试
# ============================================================

class TestMetalCoordinationDetectorMeta:
    """验证检测器的元信息。"""

    def test_name(self):
        detector = MetalCoordinationDetector()
        assert detector.name == "metal_coordination"

    def test_required_group_types(self):
        detector = MetalCoordinationDetector()
        assert "metal" in detector.required_group_types
        assert "metal_binding" in detector.required_group_types

    def test_metric_names(self):
        detector = MetalCoordinationDetector()
        assert "distance" in detector.metric_names
