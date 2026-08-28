# -*- coding: utf-8 -*-
"""盐桥检测器集成测试：从真实 tpr + xtc 输入，验证完整的盐桥检测结果。

使用 D927 体系的真实数据作为测试基准。
"""

import pytest
from pathlib import Path

import numpy as np

from DuIvyInteractions.input_readers import GmxTprReader
from DuIvyInteractions.group_identifiers import AmberFFGroupIdentifier
from DuIvyInteractions.interaction_detectors import SaltBridgeDetector
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
def salt_bridges(groups):
    """运行盐桥检测，返回 Interaction 列表。"""
    filtered = [g for g in groups
                if g.group_type in ("charged_positive", "charged_negative")]
    detector = SaltBridgeDetector()
    return detector.detect(
        filtered,
        n_workers=32,
        topology_path=str(TPR_FILE),
        trajectory_path=str(XTC_FILE)
    )


# ============================================================
# 盐桥检测结果测试
# ============================================================

class TestSaltBridgeDetection:
    """验证盐桥检测结果。"""

    def test_has_results(self, salt_bridges):
        """应检测到盐桥。"""
        assert len(salt_bridges) > 0

    def test_interaction_type(self, salt_bridges):
        """相互作用类型应为 salt_bridge。"""
        assert salt_bridges[0].interaction_type == "salt_bridge"

    def test_n_pairs(self, salt_bridges):
        """应检测到 47 对盐桥。"""
        assert salt_bridges[0].n_pairs == 47

    def test_all_opposite_charge(self, salt_bridges):
        """所有盐桥应为正电-负电对。"""
        it = salt_bridges[0]
        for g1, g2 in it.groups:
            assert g1.group_type == "charged_positive"
            assert g2.group_type == "charged_negative"

    def test_occupancy_range(self, salt_bridges):
        """占位率应在 [0, 1] 范围内。"""
        occupancy = salt_bridges[0].occupancy()
        assert np.all(occupancy >= 0)
        assert np.all(occupancy <= 1)

    def test_high_occupancy_pairs(self, salt_bridges):
        """应有占位率为 1.0 的盐桥（如 ARG303↔ASP189）。"""
        it = salt_bridges[0]
        occupancy = it.occupancy()
        assert np.sum(occupancy == 1.0) > 0

    def test_top_pair(self, salt_bridges):
        """最高占位率的盐桥应为 ARG/ASP 或 ARG/GLU 对。"""
        it = salt_bridges[0]
        occupancy = it.occupancy()
        top_idx = np.argmax(occupancy)
        g1, g2 = it.groups[top_idx]
        assert g1.group_type == "charged_positive"
        assert g2.group_type == "charged_negative"
        assert occupancy[top_idx] == 1.0

    def test_metrics_shape(self, salt_bridges):
        """metrics 数组形状应为 (n_pairs, n_frames)。"""
        it = salt_bridges[0]
        n_frames = it.n_frames
        assert it.metrics["distance"].shape == (it.n_pairs, n_frames)

    def test_distance_range(self, salt_bridges):
        """活跃帧的电荷中心距离应 ≤ 5.5 Å。"""
        it = salt_bridges[0]
        for i in range(it.n_pairs):
            active_dist = it.metrics["distance"][i][it.existence[i]]
            if len(active_dist) > 0:
                assert np.all(active_dist <= 5.5)

    def test_residue_coverage(self, salt_bridges):
        """盐桥应覆盖多个残基类型（ARG, LYS, ASP, GLU）。"""
        it = salt_bridges[0]
        residue_names = set()
        for g1, g2 in it.groups:
            residue_names.add(g1.residue_name)
            residue_names.add(g2.residue_name)
        assert "ARG" in residue_names
        assert "LYS" in residue_names
        assert "ASP" in residue_names
        assert "GLU" in residue_names
