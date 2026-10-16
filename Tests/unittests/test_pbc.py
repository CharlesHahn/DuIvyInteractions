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

    def test_halogen_cross_box_distance(self):
        """卤键跨盒：X 与 A 跨盒（MIC ≤4.0）→ 应判存在；欧氏 >4.0 → 旧逻辑漏检。"""
        box = np.array([10.0, 10.0, 10.0, 90, 90, 90])
        x = np.array([[0.5, 5.0, 5.0]])   # 卤素供体（跨盒侧）
        a = np.array([[9.5, 5.0, 5.0]])   # 受体（盒另一侧）
        mic = calc_bonds(x, a, box=box)
        euclid = np.linalg.norm(x - a, axis=1)
        assert mic[0] == pytest.approx(1.0, abs=1e-6)   # MIC = 10-9 = 1.0 Å
        assert euclid[0] == pytest.approx(9.0, abs=1e-6)
        assert mic[0] <= 4.0  # HALOGEN_DIST_MAX（闭区间）

    def test_halogen_cross_box_don_angle(self):
        """卤键跨盒 C-X···A 角：顶点 X 经 MIC 正确。"""
        from MDAnalysis.lib.distances import calc_angles
        box = np.array([10.0, 10.0, 10.0, 90, 90, 90])
        c = np.array([[0.8, 5.0, 5.0]])   # C
        x = np.array([[0.0, 5.0, 5.0]])   # X（顶点）
        a = np.array([[9.9, 5.0, 5.0]])   # A（跨盒，x→a 欧氏 9.9 > 半盒 5）
        mic = np.degrees(calc_angles(c, x, a, box=box))[0]
        euclid = np.degrees(calc_angles(c, x, a, box=None))[0]
        # x→c=(0.8,0,0)；x→a MIC 后 = (-0.1,0,0)（a-10）→ 夹角 180°
        assert mic == pytest.approx(180.0, abs=1e-4)
        assert euclid == pytest.approx(0.0, abs=1e-4)  # 无盒欧氏给出错误 0°

    def test_halogen_acc_angle_expand_matches_hand(self):
        """卤键 acc_angle 展开-还原（calc_angles 批量 MIC）与逐 R 手算一致。"""
        from MDAnalysis.lib.distances import calc_angles
        rng = np.random.default_rng(3)
        box = np.array([100.0, 100.0, 100.0, 90, 90, 90])
        n_pairs, max_atoms = 4, 5
        x = rng.random((n_pairs, 3)) * 50
        a = rng.random((n_pairs, 3)) * 50
        r = rng.random((n_pairs, max_atoms, 3)) * 50
        r_only = r[:, 1:, :]              # 跳过第 0 列（A 自身）
        mr = r_only.shape[1]

        # 手写（原实现逐 R）
        ax = x - a
        ar = r_only - a[:, None, :]
        cos = np.sum(ax[:, None, :] * ar, axis=2) / (
            np.linalg.norm(ax, axis=1, keepdims=True) * np.linalg.norm(ar, axis=2))
        hand_all = np.degrees(np.arccos(np.clip(cos, -1, 1)))
        hand_best = np.argmin(np.abs(hand_all - 120.0), axis=1)
        hand_acc = hand_all[np.arange(n_pairs), hand_best]

        # 展开-还原（新实现）
        x_rep = np.repeat(x, mr, axis=0)
        a_rep = np.repeat(a, mr, axis=0)
        all_acc = np.degrees(calc_angles(x_rep, a_rep, r_only.reshape(-1, 3), box=box))
        new_all = all_acc.reshape(n_pairs, mr)
        new_best = np.argmin(np.abs(new_all - 120.0), axis=1)
        new_acc = new_all[np.arange(n_pairs), new_best]

        assert np.array_equal(hand_best, new_best)      # argmin 选择一致
        assert np.allclose(hand_acc, new_acc, atol=1e-4)  # 最终 acc_angle 一致

    def test_saltbridge_cross_box_charge_center_wrap(self):
        """盐桥跨盒残基：电荷中心以参考原子 MIC 回绕（保基团连续）→ 中心正确。"""
        from MDAnalysis.lib.distances import minimize_vectors
        box = np.array([10.0, 10.0, 10.0, 90, 90, 90])
        # 残基原子跨盒: 参考 1.0、另两原子 9.9/9.5（L=10）→ 应回绕到参考附近
        coords = np.array([[[1.0, 5.0, 5.0], [9.9, 5.0, 5.0], [9.5, 5.0, 5.0]]])
        n, ma = coords.shape[:2]
        ref = coords[:, 0:1, :]
        rel = coords - ref
        rel_mic = minimize_vectors(rel.reshape(-1, 3), box).reshape(n, ma, 3)
        wrapped = ref + rel_mic
        # 期望: 9.9 → 参考 1.0 的最近像（9.9-10=-0.1），9.5 → -0.5
        assert wrapped[0, 1, 0] == pytest.approx(-0.1, abs=1e-5)
        assert wrapped[0, 2, 0] == pytest.approx(-0.5, abs=1e-5)
        # 不回绕的欧氏中心（错误）: (1.0+9.9+9.5)/3=6.8 vs 回绕中心 (1.0-0.1-0.5)/3≈0.133
        eu_center = coords[0].mean(axis=0)[0]
        wrap_center = wrapped[0].mean(axis=0)[0]
        assert eu_center == pytest.approx(6.8, abs=1e-5)
        assert wrap_center == pytest.approx(0.1333, abs=1e-3)

    def test_saltbridge_cross_box_distance_mic(self):
        """盐桥跨盒：正/负电荷中心跨盒（MIC ≤5.5）→ 判存在；欧氏 >5.5 旧漏检。"""
        box = np.array([10.0, 10.0, 10.0, 90, 90, 90])
        pc = np.array([[0.2, 5.0, 5.0]])   # 正电中心
        nc = np.array([[9.8, 5.0, 5.0]])   # 负电中心（跨盒）
        mic = calc_bonds(pc, nc, box=box)
        euclid = np.linalg.norm(pc - nc, axis=1)
        assert mic[0] == pytest.approx(0.4, abs=1e-6)
        assert euclid[0] == pytest.approx(9.6, abs=1e-6)
        assert mic[0] <= 5.5  # SALTBRIDGE_DIST_MAX（闭区间）

    def test_saltbridge_normal_residue_wrap_no_effect(self):
        """盐桥正常残基（原子间距 <半盒）：回绕中心 = 欧氏中心（MIC 不触发）。"""
        from MDAnalysis.lib.distances import minimize_vectors
        box = np.array([10.0, 10.0, 10.0, 90, 90, 90])
        coords = np.array([[[1.0, 5.0, 5.0], [2.0, 5.0, 5.0], [1.5, 5.0, 5.0]]])
        n, ma = coords.shape[:2]
        ref = coords[:, 0:1, :]
        rel = coords - ref
        rel_mic = minimize_vectors(rel.reshape(-1, 3), box).reshape(n, ma, 3)
        wrapped = ref + rel_mic
        assert np.allclose(wrapped, coords)  # 回绕不改变（不跨盒）
        wrap_center = wrapped[0].mean(axis=0)
        eu_center = coords[0].mean(axis=0)
        assert np.allclose(wrap_center, eu_center)

    def test_saltbridge_reference_atom_independence(self):
        """盐桥参考原子选择无关性：不同参考原子回绕 → 中心间 MIC 距离一致。"""
        from MDAnalysis.lib.distances import minimize_vectors
        box = np.array([10.0, 10.0, 10.0, 90, 90, 90])

        def center_with_ref(coords, ref_idx):
            n, ma = coords.shape[:2]
            ref = coords[:, ref_idx:ref_idx+1, :]
            rel = coords - ref
            rel_mic = minimize_vectors(rel.reshape(-1, 3), box).reshape(n, ma, 3)
            return (ref + rel_mic)[0].mean(axis=0)

        cP = np.array([[0.2, 5.0, 5.0], [9.9, 5.0, 5.0], [9.5, 5.0, 5.0]])  # 跨盒残基
        cN = np.array([[5.2, 5.0, 5.0], [5.8, 5.0, 5.0]])                    # 正常残基
        pc0 = center_with_ref(cP[None], 0)
        pc1 = center_with_ref(cP[None], 1)
        nc = center_with_ref(cN[None], 0)
        d0 = calc_bonds(pc0[None], nc[None], box=box)[0]
        d1 = calc_bonds(pc1[None], nc[None], box=box)[0]
        assert d0 == pytest.approx(d1, abs=1e-6)  # 判定不依赖参考原子选择

    def test_pi_cation_cross_box_ring_wrap(self):
        """π-阳离子跨盒环：环原子跨盒（部分 wrap 到 L 侧）→ 回绕后环心/法向量正确。"""
        from MDAnalysis.lib.distances import minimize_vectors
        box = np.array([10.0, 10.0, 10.0, 90, 90, 90])
        # 真实苯环（xy 平面 z=1, 环心 x=0.2 跨盒），部分原子 x<0 wrap 到 L 侧
        R = 1.4
        angles = np.arange(6) * np.pi / 3
        ring_true = np.array([[0.2 + R*np.cos(a), 5 + R*np.sin(a), 1.0]
                              for a in angles])
        ring_wrap = ring_true.copy()
        ring_wrap[:, 0][ring_wrap[:, 0] < 0] += 10   # wrap 到 [0,L)
        # 回绕（以原子 0 为参考）
        ref = ring_wrap[0:1]
        rel = ring_wrap - ref
        rel_mic = minimize_vectors(rel.reshape(-1, 3), box).reshape(-1, 3)
        ring_back = ref + rel_mic
        # 环心
        c_back = ring_back.mean(axis=0)
        c_true = ring_true.mean(axis=0)
        assert np.linalg.norm(c_back - c_true) < 1e-5   # 环心正确（差整盒或 0）
        # 法向量
        def normal_of(r):
            v1 = r - np.roll(r, 1, axis=0)
            v2 = r - np.roll(r, -1, axis=0)
            n = np.cross(v1, v2).sum(axis=0)
            return n / np.linalg.norm(n)
        n_back = normal_of(ring_back)
        n_true = normal_of(ring_true)
        cos_ang = np.clip(np.dot(n_back, n_true), -1, 1)
        assert np.degrees(np.arccos(cos_ang)) < 1e-3   # 法向量一致

    def test_pi_cation_cross_box_distance(self):
        """π-阳离子跨盒：电荷中心与环心跨盒（MIC 0.5<d<6.0）→ 判存在；欧氏 >6.0 漏检。"""
        box = np.array([10.0, 10.0, 10.0, 90, 90, 90])
        ring_c = np.array([[0.5, 5.0, 5.0]])    # 环心
        cat_c = np.array([[9.5, 5.0, 5.0]])     # 电荷中心（跨盒，MIC=1.0 ∈ (0.5,6.0)）
        mic = calc_bonds(ring_c, cat_c, box=box)
        euclid = np.linalg.norm(ring_c - cat_c, axis=1)
        assert mic[0] == pytest.approx(1.0, abs=1e-6)
        assert euclid[0] == pytest.approx(9.0, abs=1e-6)
        assert 0.5 < mic[0] < 6.0   # PICATION_DIST_MAX（开区间）

    def test_pi_cation_projection_offset_mic(self):
        """π-阳离子投影偏移：跨盒时欧氏 t 偏移错（>2.0），MIC t 偏移对（<2.0）。"""
        from MDAnalysis.lib.distances import minimize_vectors
        box = np.array([10.0, 10.0, 10.0, 90, 90, 90])
        # 环在 xy 平面 z=1（法向量 z），环心 x=0.2；电荷中心在环平面内、距环心 0.3Å 但跨盒
        ring_c = np.array([[0.2, 5.0, 1.0]])
        cat_c = np.array([[9.9, 5.0, 1.0]])    # MIC 距环心 0.3Å（x 方向）
        n_hat = np.array([[0.0, 0.0, 1.0]])    # 环法向量

        def proj_old(n, p, tgt):
            d1 = np.linalg.norm(tgt - (p + n), axis=1)
            d2 = np.linalg.norm(tgt - (p - n), axis=1)
            sign = np.where(d1 < d2, 1.0, -1.0)[:, None]
            on = n * sign
            t = tgt - p
            pd = np.sum(t * on, axis=1)
            pp = tgt - pd[:, None] * on
            return np.linalg.norm(pp - p, axis=1)

        def proj_mic(n, p, tgt, box):
            d1 = np.linalg.norm(tgt - (p + n), axis=1)
            d2 = np.linalg.norm(tgt - (p - n), axis=1)
            sign = np.where(d1 < d2, 1.0, -1.0)[:, None]
            on = n * sign
            t = minimize_vectors((tgt - p).reshape(-1, 3), box).reshape(-1, 3)
            pd = np.sum(t * on, axis=1)
            pp = p + (t - pd[:, None] * on)
            return np.linalg.norm(pp - p, axis=1)

        off_old = proj_old(n_hat, ring_c, cat_c)[0]
        off_mic = proj_mic(n_hat, ring_c, cat_c, box)[0]
        # 物理正确: 电荷中心在环平面内距环心 0.3Å → offset=0.3（<2.0 判存在）
        assert off_mic == pytest.approx(0.3, abs=1e-3)
        assert off_old > 2.0   # 欧氏 t 错误（电荷中心在 x=9.9 欧氏距 9.7 > 2.0）
        assert off_mic < 2.0   # MIC t 正确 → 判定翻转（修复）


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