# -*- coding: utf-8 -*-
"""PBC（周期边界）处理正确性测试：MIC 距离/角度、无盒退化。

设计文档：doc/pbc_design.md
被测对象：core/pbc.get_box + 检测器距离计算的 PBC 化
（疏水 two_pass/per_frame、金属 two_pass/per_frame）。
"""

import numpy as np
import pytest

from DuIvyInteractions.core.pbc import get_box
from MDAnalysis.lib.distances import calc_bonds, self_capped_distance, capped_distance


# ============================================================
# get_box 防护
# ============================================================

class TestGetBox:
    """无盒/退化盒返回 None（PBC 退化欧氏）。"""

    class _FakeTs:
        def __init__(self, dims):
            self.dimensions = dims

    def test_normal_box(self):
        ts = self._FakeTs(np.array([10.0, 10.0, 10.0, 90, 90, 90]))
        box = get_box(ts)
        assert box is not None and box[0] == 10.0

    def test_no_dimensions(self):
        ts = self._FakeTs(None)
        assert get_box(ts) is None

    def test_zero_box(self):
        ts = self._FakeTs(np.array([0.0, 0.0, 0.0, 90, 90, 90]))
        assert get_box(ts) is None


# ============================================================
# MIC 距离正确性（跨盒近邻）
# ============================================================

class TestPbcDistance:
    """最小镜像约定：跨盒原子对距离取最近镜像。"""

    def test_cross_box_pair(self):
        """x=0.1 与 x=L-0.1（L=10）→ MIC 距离 ≈ 0.2Å（欧氏 ≈ 9.8Å）。"""
        box = np.array([10.0, 10.0, 10.0, 90, 90, 90])
        c1 = np.array([[0.1, 5.0, 5.0]])
        c2 = np.array([[9.9, 5.0, 5.0]])
        mic = calc_bonds(c1, c2, box=box)
        euclid = np.linalg.norm(c1 - c2, axis=1)
        assert mic[0] == pytest.approx(0.2, abs=1e-6)
        assert euclid[0] == pytest.approx(9.8, abs=1e-6)

    def test_cross_box_hydrophobic_detected(self):
        """跨盒疏水对：MIC 距离落在 (0.5, 4.0) → 应被检测；欧氏 > 4.0 → 旧逻辑漏检。

        构造：x=0.25 与 x=8.75（L=10）→ 欧氏 8.5Å（>4.0 漏检），
              MIC = 10 - 8.5 = 1.5Å（>0.5 且 <4.0 → 应检测）。
        """
        box = np.array([10.0, 10.0, 10.0, 90, 90, 90])
        coords = np.array([[0.25, 5.0, 5.0], [8.75, 5.0, 5.0], [5.0, 5.0, 5.0]])
        pairs, dists = self_capped_distance(
            coords, 4.0, min_cutoff=0.5, box=box, method='pkdtree')
        # 跨盒对 (0,1) 应出现且距离 1.5
        pair_map = dict(zip([tuple(p) for p in pairs], dists))
        assert (0, 1) in pair_map
        assert pair_map[(0, 1)] == pytest.approx(1.5, abs=1e-5)

    def test_no_box_degrades_to_euclid(self):
        """无盒（box=None）：与欧氏距离一致。"""
        coords = np.array([[0.1, 5.0, 5.0], [9.9, 5.0, 5.0]])
        pairs, dists = self_capped_distance(
            coords, 10.0, box=None, method='pkdtree')
        assert dists[0] == pytest.approx(9.8, abs=1e-6)

    def test_min_dist_lower_bound(self):
        """距离 < 0.5 的对被 min_cutoff 排除（下界不破）。"""
        coords = np.array([[0.0, 0.0, 0.0], [0.2, 0.0, 0.0]])
        pairs, dists = self_capped_distance(
            coords, 4.0, min_cutoff=0.5, box=None, method='pkdtree')
        assert len(pairs) == 0

    def test_metal_cross_box(self):
        """金属与配位原子跨盒（MIC < 3.0）→ 存在性判定正确。"""
        box = np.array([10.0, 10.0, 10.0, 90, 90, 90])
        m = np.array([[0.2, 5.0, 5.0]])   # 金属
        b = np.array([[9.9, 5.0, 5.0]])   # 配位原子（跨盒）
        mic = calc_bonds(m, b, box=box)
        assert mic[0] == pytest.approx(0.3, abs=1e-6)
        assert mic[0] < 3.0  # METAL_DIST_MAX

    def test_hbond_cross_box_distance(self):
        """氢键跨盒：D 与 A 跨盒（MIC ≤4.1）→ 应判存在；欧氏 >4.1 → 旧逻辑漏检。"""
        box = np.array([10.0, 10.0, 10.0, 90, 90, 90])
        d = np.array([[0.5, 5.0, 5.0]])   # 供体（跨盒侧）
        a = np.array([[9.5, 5.0, 5.0]])   # 受体（盒另一侧）
        mic = calc_bonds(d, a, box=box)
        euclid = np.linalg.norm(d - a, axis=1)
        assert mic[0] == pytest.approx(1.0, abs=1e-6)   # MIC = 10-9 = 1.0 Å
        assert euclid[0] == pytest.approx(9.0, abs=1e-6)
        assert mic[0] <= 4.1  # HBOND_DIST_MAX（闭区间）

    def test_hbond_cross_box_angle(self):
        """氢键跨盒角度：H 处 D-H···A 角经 MIC 正确（跨盒向量回绕）。"""
        from MDAnalysis.lib.distances import calc_angles
        box = np.array([10.0, 10.0, 10.0, 90, 90, 90])
        d = np.array([[0.9, 5.0, 5.0]])   # D
        h = np.array([[0.0, 5.0, 5.0]])   # H（顶点）
        a = np.array([[9.9, 5.0, 5.0]])   # A（跨盒，h→a 欧氏 9.9 > 半盒 5）
        mic = np.degrees(calc_angles(d, h, a, box=box))[0]
        euclid = np.degrees(calc_angles(d, h, a, box=None))[0]
        # h→d=(0.9,0,0)；h→a MIC 后 = (-0.1,0,0)（a-10）→ 夹角 180°
        assert mic == pytest.approx(180.0, abs=1e-4)
        assert euclid == pytest.approx(0.0, abs=1e-4)  # 无盒手写给出错误 0°


# ============================================================
# self_capped_distance 去重与返回结构（对抗性审查锁定）
# ============================================================

class TestSelfCappedStructure:
    """去重 i<j、无自对；return_distances=False 只返回 pairs（不解包）。"""

    def test_dedup_no_self(self):
        coords = np.array([[0, 0, 0], [2.0, 0, 0], [3.0, 0, 0]])
        pairs, dists = self_capped_distance(
            coords, 4.0, method='pkdtree', return_distances=True)
        assert pairs.tolist() == [[0, 1], [0, 2], [1, 2]]
        assert dists.tolist() == [2.0, 3.0, 1.0]

    def test_no_distances_returns_single_array(self):
        coords = np.array([[0, 0, 0], [2.0, 0, 0], [3.0, 0, 0]])
        pairs = self_capped_distance(
            coords, 4.0, method='pkdtree', return_distances=False)
        assert isinstance(pairs, np.ndarray)
        assert pairs.shape == (3, 2)

    def test_empty_pairs_shape(self):
        """self_capped_distance 空对返回 pairs.shape=(0,2)、dists.shape=(0,)（实测）。"""
        coords = np.array([[0, 0, 0], [20.0, 0, 0]])
        pairs, dists = self_capped_distance(
            coords, 3.0, method='pkdtree', return_distances=True)
        assert pairs.shape == (0, 2)
        assert dists.shape == (0,)