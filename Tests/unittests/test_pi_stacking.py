# -*- coding: utf-8 -*-
"""π-π 堆积检测器集成测试：从真实 tpr + xtc 输入，验证完整的 π-π 堆积检测结果。

使用 D927 体系的真实数据作为测试基准。
"""

import pytest
from pathlib import Path

import numpy as np

from DuIvyInteractions.input_readers import GmxTprReader
from DuIvyInteractions.group_identifiers import AmberFFGroupIdentifier
from DuIvyInteractions.interaction_detectors import PiStackingDetector
import MDAnalysis as mda


# 测试数据路径
TPR_FILE = Path(__file__).parent.parent / "test_MD_case" / "md.tpr"
XTC_FILE = Path(__file__).parent.parent / "test_MD_case" / "md1ns.xtc"

# GmxTprReader 的分子名前缀
MOL_RBD = "seg_0_RBD_pro"
MOL_D927 = "seg_1_D927"
MOL_KRAS = "seg_2_KRAS_pro"


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
def aromatic_groups(groups):
    """只保留芳香环基团。"""
    return [g for g in groups if g.group_type == "aromatic_ring"]


@pytest.fixture(scope="module")
def pi_stackings(aromatic_groups):
    """运行 π-π 堆积检测（不检查平面性），返回 Interaction 列表。"""
    u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    detector = PiStackingDetector(check_planarity=False)
    return detector.detect(aromatic_groups, trajectory=u.trajectory)


@pytest.fixture(scope="module")
def pi_stackings_with_planarity(aromatic_groups):
    """运行 π-π 堆积检测（检查平面性），返回 Interaction 列表。"""
    u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    detector = PiStackingDetector(check_planarity=True)
    return detector.detect(aromatic_groups, trajectory=u.trajectory)


# ============================================================
# 基本功能测试
# ============================================================

class TestPiStackingBasic:
    """验证 π-π 堆积检测的基本功能。"""

    def test_returns_list(self, pi_stackings):
        assert isinstance(pi_stackings, list)

    def test_interaction_type(self, pi_stackings):
        for interaction in pi_stackings:
            assert interaction.interaction_type == "pi_stacking"

    def test_metrics_keys(self, pi_stackings):
        for interaction in pi_stackings:
            assert "distance" in interaction.metrics
            assert "angle" in interaction.metrics
            assert "offset" in interaction.metrics
            assert "pistacking_type" in interaction.metrics

    def test_no_planarity_by_default(self, pi_stackings):
        for interaction in pi_stackings:
            assert "planarity_ring1" not in interaction.metrics
            assert "planarity_ring2" not in interaction.metrics


# ============================================================
# 平面性参数测试
# ============================================================

class TestPiStackingPlanarity:
    """验证平面性参数。"""

    def test_planarity_metrics_present(self, pi_stackings_with_planarity):
        for interaction in pi_stackings_with_planarity:
            assert "planarity_ring1" in interaction.metrics
            assert "planarity_ring2" in interaction.metrics

    def test_planarity_shape(self, pi_stackings_with_planarity):
        for interaction in pi_stackings_with_planarity:
            n_pairs = interaction.n_pairs
            n_frames = interaction.n_frames
            assert interaction.metrics["planarity_ring1"].shape == (n_pairs, n_frames)
            assert interaction.metrics["planarity_ring2"].shape == (n_pairs, n_frames)


# ============================================================
# 指标值范围测试
# ============================================================

class TestPiStackingMetrics:
    """验证指标值的合理范围。"""

    def test_distance_positive(self, pi_stackings):
        for interaction in pi_stackings:
            assert np.all(interaction.metrics["distance"] >= 0)

    def test_angle_range(self, pi_stackings):
        for interaction in pi_stackings:
            angle = interaction.metrics["angle"]
            assert np.all(angle >= 0)
            assert np.all(angle <= 90)  # min(θ, 180-θ) 范围是 0~90°

    def test_offset_positive(self, pi_stackings):
        for interaction in pi_stackings:
            assert np.all(interaction.metrics["offset"] >= 0)

    def test_pistacking_type_values(self, pi_stackings):
        valid_types = {'N', 'P', 'T'}
        for interaction in pi_stackings:
            types = set(interaction.metrics["pistacking_type"].flatten())
            assert types.issubset(valid_types)

    def test_existence_matches_type(self, pi_stackings):
        """existence=True 的帧，pistacking_type 不应为 'N'。"""
        for interaction in pi_stackings:
            existence = interaction.existence
            types = interaction.metrics["pistacking_type"]
            for i in range(interaction.n_pairs):
                active_frames = existence[i]
                active_types = types[i][active_frames]
                assert np.all(active_types != 'N')


# ============================================================
# P 型 / T 型分类测试
# ============================================================

class TestPiStackingClassification:
    """验证 P 型和 T 型分类的正确性。"""

    def test_p_type_angle(self, pi_stackings):
        """P 型堆积的 angle 应 ≤ 30°。"""
        for interaction in pi_stackings:
            types = interaction.metrics["pistacking_type"]
            angle = interaction.metrics["angle"]
            p_mask = types == 'P'
            if np.any(p_mask):
                assert np.all(angle[p_mask] <= 30.0)

    def test_t_type_angle(self, pi_stackings):
        """T 型堆积的 angle 应 ≥ 60°。"""
        for interaction in pi_stackings:
            types = interaction.metrics["pistacking_type"]
            angle = interaction.metrics["angle"]
            t_mask = types == 'T'
            if np.any(t_mask):
                assert np.all(angle[t_mask] >= 60.0)

    def test_active_offset(self, pi_stackings):
        """活跃帧的 offset 应 < 2.0 Å。"""
        for interaction in pi_stackings:
            existence = interaction.existence
            offset = interaction.metrics["offset"]
            for i in range(interaction.n_pairs):
                active_frames = existence[i]
                if np.any(active_frames):
                    assert np.all(offset[i][active_frames] < 2.0)


# ============================================================
# 候选对测试
# ============================================================

class TestPiStackingPairs:
    """验证候选对的生成。"""

    def test_pairs_are_different_rings(self, pi_stackings):
        """每个堆积对的两个环不应相同。"""
        for interaction in pi_stackings:
            for ring1, ring2 in interaction.groups:
                assert ring1.group_id != ring2.group_id

    def test_pairs_are_aromatic(self, pi_stackings):
        """每个堆积对的两个环都应是芳香环。"""
        for interaction in pi_stackings:
            for ring1, ring2 in interaction.groups:
                assert ring1.group_type == "aromatic_ring"
                assert ring2.group_type == "aromatic_ring"


# ============================================================
# 检测器元信息测试
# ============================================================

class TestPiStackingDetectorMeta:
    """验证检测器的元信息。"""

    def test_name(self):
        detector = PiStackingDetector()
        assert detector.name == "pi_stacking"

    def test_required_group_types(self):
        detector = PiStackingDetector()
        assert detector.required_group_types == ["aromatic_ring"]

    def test_metric_names_without_planarity(self):
        detector = PiStackingDetector(check_planarity=False)
        names = detector.metric_names
        assert "distance" in names
        assert "angle" in names
        assert "offset" in names
        assert "pistacking_type" in names
        assert "planarity_ring1" not in names
        assert "planarity_ring2" not in names

    def test_metric_names_with_planarity(self):
        detector = PiStackingDetector(check_planarity=True)
        names = detector.metric_names
        assert "planarity_ring1" in names
        assert "planarity_ring2" in names
