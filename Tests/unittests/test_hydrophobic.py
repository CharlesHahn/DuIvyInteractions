# -*- coding: utf-8 -*-
"""疏水相互作用检测器集成测试：从真实 tpr + xtc 输入，验证完整的检测结果。

使用 D927 体系的真实数据作为测试基准。
"""

import pytest
from pathlib import Path
from collections import Counter

import numpy as np

from DuIvyInteractions.input_readers import GmxTprReader
from DuIvyInteractions.group_identifiers import AmberFFGroupIdentifier
from DuIvyInteractions.interaction_detectors import HydrophobicDetectorPerTuple
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
def hydrophobic_groups(groups):
    """只保留疏水基团。"""
    return [g for g in groups if g.group_type == "hydrophobic"]


@pytest.fixture(scope="module")
def hydrophobic_interactions(hydrophobic_groups):
    """运行疏水检测（并行），返回 Interaction 列表。"""
    detector = HydrophobicDetectorPerTuple()

    # 只取前 100 个疏水原子，减少候选对数量
    subset = hydrophobic_groups[:100]

    return detector.detect(
        subset,
        n_workers=16,
        topology_path=str(TPR_FILE),
        trajectory_path=str(XTC_FILE)
    )


# ============================================================
# 基本功能测试
# ============================================================

class TestHydrophobicBasic:
    """验证疏水检测的基本功能。"""

    def test_returns_list(self, hydrophobic_interactions):
        assert isinstance(hydrophobic_interactions, list)

    def test_interaction_type(self, hydrophobic_interactions):
        for interaction in hydrophobic_interactions:
            assert interaction.interaction_type == "hydrophobic"

    def test_metrics_keys(self, hydrophobic_interactions):
        for interaction in hydrophobic_interactions:
            assert "distance" in interaction.metrics


# ============================================================
# 指标值范围测试
# ============================================================

class TestHydrophobicMetrics:
    """验证指标值的合理范围。"""

    def test_distance_positive(self, hydrophobic_interactions):
        for interaction in hydrophobic_interactions:
            assert np.all(interaction.metrics["distance"] >= 0)

    def test_metrics_shape(self, hydrophobic_interactions):
        """metrics 数组形状应为 (n_pairs, n_frames)。"""
        for interaction in hydrophobic_interactions:
            n_frames = interaction.n_frames
            assert interaction.metrics["distance"].shape == (interaction.n_pairs, n_frames)


# ============================================================
# 活跃帧值范围测试
# ============================================================

class TestHydrophobicActiveFrames:
    """验证活跃帧的指标值范围。"""

    def test_active_distance_range(self, hydrophobic_interactions):
        """活跃帧的距离应 < 4.0 Å。"""
        for interaction in hydrophobic_interactions:
            for i in range(interaction.n_pairs):
                active_dist = interaction.metrics["distance"][i][interaction.existence[i]]
                if len(active_dist) > 0:
                    assert np.all(active_dist < 4.0)


# ============================================================
# 候选对测试
# ============================================================

class TestHydrophobicPairs:
    """验证候选对的生成。"""

    def test_pairs_are_different_atoms(self, hydrophobic_interactions):
        """每个疏水对的两个基团不应相同。"""
        for interaction in hydrophobic_interactions:
            for g1, g2 in interaction.groups:
                assert g1.group_id != g2.group_id

    def test_pairs_are_hydrophobic(self, hydrophobic_interactions):
        """每个对的两个基团都应是疏水基团。"""
        for interaction in hydrophobic_interactions:
            for g1, g2 in interaction.groups:
                assert g1.group_type == "hydrophobic"
                assert g2.group_type == "hydrophobic"


# ============================================================
# 检测器元信息测试
# ============================================================

class TestHydrophobicDetectorPerTupleMeta:
    """验证检测器的元信息。"""

    def test_name(self):
        detector = HydrophobicDetectorPerTuple()
        assert detector.name == "hydrophobic"

    def test_required_group_types(self):
        detector = HydrophobicDetectorPerTuple()
        assert detector.required_group_types == ["hydrophobic"]

    def test_metric_names(self):
        detector = HydrophobicDetectorPerTuple()
        assert "distance" in detector.metric_names


# ============================================================
# 结果数据测试
# ============================================================

class TestHydrophobicResults:
    """验证基于真实数据的检测结果。"""

    def test_has_results(self, hydrophobic_interactions):
        """应检测到疏水相互作用。"""
        assert len(hydrophobic_interactions) > 0

    def test_n_pairs_reasonable(self, hydrophobic_interactions):
        """100 个原子子集，去重后对数应 < 4950（C(100,2)）。"""
        n_pairs = hydrophobic_interactions[0].n_pairs
        assert n_pairs > 0
        assert n_pairs < 4950

    def test_hydrophobic_count(self, hydrophobic_groups):
        """体系应有 734 个疏水原子。"""
        assert len(hydrophobic_groups) == 734

    def test_deduplication(self, hydrophobic_interactions):
        """去重后，同一原子不应与同残基多个原子同时保留。"""
        inter = hydrophobic_interactions[0]
        # 收集每个 (g1_id, r2_id) 对
        seen = {}
        for i, (g1, g2) in enumerate(inter.groups):
            key1 = (g1.group_id, g2.residue_id)
            key2 = (g2.group_id, g1.residue_id)
            # 每个 (atom, residue) 对应最多一个接触
            assert key1 not in seen, f"重复: g1={g1}, g2={g2}"
            assert key2 not in seen, f"重复: g2={g2}, g1={g1}"
            seen[key1] = i
            seen[key2] = i

    def test_occupancy_range(self, hydrophobic_interactions):
        """占位率应在 [0, 1] 范围内。"""
        occ = hydrophobic_interactions[0].occupancy()
        assert np.all(occ >= 0)
        assert np.all(occ <= 1)

    def test_same_molecule_pairs_exist(self, hydrophobic_interactions):
        """应存在同分子的疏水相互作用。"""
        inter = hydrophobic_interactions[0]
        same_mol = sum(1 for g1, g2 in inter.groups if g1.molecule == g2.molecule)
        assert same_mol > 0
