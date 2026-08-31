# -*- coding: utf-8 -*-
"""π-阳离子相互作用检测器集成测试：从真实 tpr + xtc 输入，验证完整的检测结果。

使用 D927 体系的真实数据作为测试基准。
"""

import pytest
from pathlib import Path

import numpy as np

from DuIvyInteractions.input_readers import GmxTprReader
from DuIvyInteractions.group_identifiers import AmberFFGroupIdentifier
from DuIvyInteractions.interaction_detectors import PiCationDetectorPerTuple
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
def relevant_groups(groups):
    """只保留芳香环和正电荷基团。"""
    return [g for g in groups
            if g.group_type in ("aromatic_ring", "charged_positive")]


@pytest.fixture(scope="module")
def pi_cations(relevant_groups):
    """运行 π-阳离子检测，返回 Interaction 列表。"""
    u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    detector = PiCationDetectorPerTuple()
    return detector.detect(relevant_groups, trajectory=u.trajectory)


# ============================================================
# 基本功能测试
# ============================================================

class TestPiCationBasic:
    """验证 π-阳离子检测的基本功能。"""

    def test_returns_list(self, pi_cations):
        assert isinstance(pi_cations, list)

    def test_interaction_type(self, pi_cations):
        for interaction in pi_cations:
            assert interaction.interaction_type == "pi_cation"

    def test_metrics_keys(self, pi_cations):
        for interaction in pi_cations:
            assert "distance" in interaction.metrics
            assert "offset" in interaction.metrics


# ============================================================
# 指标值范围测试
# ============================================================

class TestPiCationMetrics:
    """验证指标值的合理范围。"""

    def test_distance_positive(self, pi_cations):
        for interaction in pi_cations:
            assert np.all(interaction.metrics["distance"] >= 0)

    def test_offset_positive(self, pi_cations):
        for interaction in pi_cations:
            assert np.all(interaction.metrics["offset"] >= 0)

    def test_metrics_shape(self, pi_cations):
        """metrics 数组形状应为 (n_pairs, n_frames)。"""
        for interaction in pi_cations:
            n_frames = interaction.n_frames
            assert interaction.metrics["distance"].shape == (interaction.n_pairs, n_frames)
            assert interaction.metrics["offset"].shape == (interaction.n_pairs, n_frames)


# ============================================================
# 活跃帧值范围测试
# ============================================================

class TestPiCationActiveFrames:
    """验证活跃帧的指标值范围。"""

    def test_active_distance_range(self, pi_cations):
        """活跃帧的距离应 < 6.0 Å。"""
        for interaction in pi_cations:
            for i in range(interaction.n_pairs):
                active_dist = interaction.metrics["distance"][i][interaction.existence[i]]
                if len(active_dist) > 0:
                    assert np.all(active_dist < 6.0)

    def test_active_offset_range(self, pi_cations):
        """活跃帧的 offset 应 < 2.0 Å。"""
        for interaction in pi_cations:
            for i in range(interaction.n_pairs):
                active_offset = interaction.metrics["offset"][i][interaction.existence[i]]
                if len(active_offset) > 0:
                    assert np.all(active_offset < 2.0)


# ============================================================
# 候选对测试
# ============================================================

class TestPiCationPairs:
    """验证候选对的生成。"""

    def test_pairs_have_correct_types(self, pi_cations):
        """每个对的两个基团类型应正确。"""
        for interaction in pi_cations:
            for ring, pos in interaction.groups:
                assert ring.group_type == "aromatic_ring"
                assert pos.group_type == "charged_positive"


# ============================================================
# 检测器元信息测试
# ============================================================

class TestPiCationDetectorPerTupleMeta:
    """验证检测器的元信息。"""

    def test_name(self):
        detector = PiCationDetectorPerTuple()
        assert detector.name == "pi_cation"

    def test_required_group_types(self):
        detector = PiCationDetectorPerTuple()
        assert "aromatic_ring" in detector.required_group_types
        assert "charged_positive" in detector.required_group_types

    def test_metric_names(self):
        detector = PiCationDetectorPerTuple()
        names = detector.metric_names
        assert "distance" in names
        assert "offset" in names


# ============================================================
# 结果数据测试
# ============================================================

class TestPiCationResults:
    """验证基于真实数据的检测结果。"""

    def test_has_results(self, pi_cations):
        """应检测到 π-阳离子相互作用。"""
        assert len(pi_cations) > 0

    def test_n_pairs(self, pi_cations):
        """应检测到 8 对 π-阳离子相互作用。"""
        assert pi_cations[0].n_pairs == 8

    def test_aromatic_ring_count(self, groups):
        """体系应有 38 个芳香环。"""
        aromatic = [g for g in groups if g.group_type == "aromatic_ring"]
        assert len(aromatic) == 38

    def test_positive_charge_count(self, groups):
        """体系应有 42 个正电荷基团。"""
        pos = [g for g in groups if g.group_type == "charged_positive"]
        assert len(pos) == 42

    def test_positive_charge_by_residue(self, groups):
        """正电荷应主要来自 ARG 和 LYS。"""
        from collections import Counter
        pos = [g for g in groups if g.group_type == "charged_positive"]
        res_counts = Counter(g.residue_name for g in pos)
        assert res_counts["ARG"] == 15
        assert res_counts["LYS"] == 25

    def test_top_pair(self, pi_cations):
        """最高占位率的 π-阳离子应为 RBD:TRP38 ↔ RBD:LYS47。"""
        it = pi_cations[0]
        occ = it.occupancy()
        top_idx = np.argmax(occ)
        r, p = it.groups[top_idx]
        assert r.residue_name == "TRP"
        assert r.residue_id == 38
        assert p.residue_name == "LYS"
        assert p.residue_id == 47
        assert occ[top_idx] > 0.6

    def test_d927_participation(self, pi_cations):
        """D927 配体应参与 π-阳离子相互作用。"""
        it = pi_cations[0]
        d927_pairs = 0
        for r, p in it.groups:
            if r.molecule == MOL_D927 or p.molecule == MOL_D927:
                d927_pairs += 1
        assert d927_pairs >= 1

    def test_lys_dominant(self, pi_cations):
        """LYS 应是 π-阳离子中主要的正电荷来源。"""
        it = pi_cations[0]
        lys_count = sum(1 for _, p in it.groups if p.residue_name == "LYS")
        arg_count = sum(1 for _, p in it.groups if p.residue_name == "ARG")
        assert lys_count > arg_count

    def test_occupancy_range(self, pi_cations):
        """占位率应在 [0, 1] 范围内。"""
        occ = pi_cations[0].occupancy()
        assert np.all(occ >= 0)
        assert np.all(occ <= 1)

    def test_metrics_shape(self, pi_cations):
        """metrics 数组形状应为 (n_pairs, n_frames)。"""
        it = pi_cations[0]
        n_frames = it.n_frames
        assert it.metrics["distance"].shape == (it.n_pairs, n_frames)
        assert it.metrics["offset"].shape == (it.n_pairs, n_frames)
