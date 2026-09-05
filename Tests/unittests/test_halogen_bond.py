# -*- coding: utf-8 -*-
"""卤键检测器集成测试：从真实 tpr + xtc 输入，验证完整的卤键检测结果。

使用 D927 体系的真实数据作为测试基准。
"""

import pytest
from pathlib import Path

import numpy as np

from DuIvyInteractions.system_readers import GmxTprReader
from DuIvyInteractions.group_identifiers import AmberFFGroupIdentifier
from DuIvyInteractions.interaction_detectors import HalogenBondDetectorPerTuple
import MDAnalysis as mda


# 测试数据路径
TPR_FILE = Path(__file__).parent.parent / "test_MD_case" / "md.tpr"
XTC_FILE = Path(__file__).parent.parent / "test_MD_case" / "md1ns.xtc"

# GmxTprReader 的分子名前缀
MOL_D927 = "seg_1_D927"


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
def halogen_bonds(groups):
    """运行卤键检测，返回 Interaction 列表。"""
    filtered = [g for g in groups
                if g.group_type in ("halogen_donor", "halogen_acceptor")]
    u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    detector = HalogenBondDetectorPerTuple()
    return detector.detect(filtered, trajectory=u.trajectory)


# ============================================================
# 卤键检测结果测试
# ============================================================

class TestHalogenBondDetection:
    """验证卤键检测结果。"""

    def test_has_results(self, halogen_bonds):
        """应检测到卤键。"""
        assert len(halogen_bonds) > 0

    def test_interaction_type(self, halogen_bonds):
        """相互作用类型应为 halogen_bond。"""
        assert halogen_bonds[0].interaction_type == "halogen_bond"

    def test_n_pairs(self, halogen_bonds):
        """应检测到 1 对卤键（D927 内部）。"""
        assert halogen_bonds[0].n_pairs == 1

    def test_metrics_keys(self, halogen_bonds):
        """metrics 应包含 distance, don_angle, acc_angle。"""
        it = halogen_bonds[0]
        assert set(it.metrics.keys()) == {"distance", "don_angle", "acc_angle"}

    def test_metrics_shape(self, halogen_bonds):
        """metrics 数组形状应为 (n_pairs, n_frames)。"""
        it = halogen_bonds[0]
        n_frames = it.n_frames
        assert it.metrics["distance"].shape == (it.n_pairs, n_frames)
        assert it.metrics["don_angle"].shape == (it.n_pairs, n_frames)
        assert it.metrics["acc_angle"].shape == (it.n_pairs, n_frames)

    def test_occupancy_range(self, halogen_bonds):
        """占位率应在 [0, 1] 范围内。"""
        occupancy = halogen_bonds[0].occupancy()
        assert np.all(occupancy >= 0)
        assert np.all(occupancy <= 1)

    def test_active_pair_is_d927(self, halogen_bonds):
        """活跃的卤键应来自 D927 分子。"""
        it = halogen_bonds[0]
        occupancy = it.occupancy()
        active_idx = np.where(occupancy > 0)[0]
        assert len(active_idx) == 1
        g1, g2 = it.groups[active_idx[0]]
        assert g1.molecule == MOL_D927

    def test_distance_range(self, halogen_bonds):
        """活跃帧的 X-A 距离应 ≤ 4.0 Å。"""
        it = halogen_bonds[0]
        for i in range(it.n_pairs):
            active_dist = it.metrics["distance"][i][it.existence[i]]
            if len(active_dist) > 0:
                assert np.all(active_dist <= 4.0)

    def test_don_angle_range(self, halogen_bonds):
        """活跃帧的 C-X···A 角度应在 [135°, 195°]。"""
        it = halogen_bonds[0]
        for i in range(it.n_pairs):
            active_angle = it.metrics["don_angle"][i][it.existence[i]]
            if len(active_angle) > 0:
                assert np.all(active_angle >= 135.0)
                assert np.all(active_angle <= 195.0)

    def test_acc_angle_range(self, halogen_bonds):
        """活跃帧的 X···A-R 角度应在 [90°, 150°]。"""
        it = halogen_bonds[0]
        for i in range(it.n_pairs):
            active_angle = it.metrics["acc_angle"][i][it.existence[i]]
            if len(active_angle) > 0:
                assert np.all(active_angle >= 90.0)
                assert np.all(active_angle <= 150.0)
