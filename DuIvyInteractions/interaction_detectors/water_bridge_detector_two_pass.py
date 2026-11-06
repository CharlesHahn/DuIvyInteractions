# -*- coding: utf-8 -*-
"""水桥检测器（两轮遍历 + 稀疏存储）。

判据：水分子位于供体和受体之间，距离和角度满足 PLIP 定义。
参考：PLIP, Jiang et al. 2005。

Pass1：KDTree 发现三元组 + 内联计算角度（PBC 感知）
Pass2：从 InteractionSparse 重建索引 + 向量化计算（PBC 感知）

PBC：候选发现 capped_distance(box)、距离 calc_bonds(box)、角度 calc_angles(box)，
见 doc/pbc_design.md §6B。
"""

import numpy as np
from typing import List, Tuple, Dict
from collections import defaultdict

from ..core.interfaces import InteractionDetectorTwoPass
from ..core.datas import Group, Interaction, InteractionSparse
from ..core.pbc import get_box


# 氢键几何判据（Baker & Hubbard 1984；设计见 doc/water_bridge_design.md）
WATER_BRIDGE_HB_DIST_MAX = 2.5      # Å，H···A 距离上限
WATER_BRIDGE_HB_ANGLE_MIN = 120.0   # °，D-H···A 角度下限
WATER_BRIDGE_MIN_DIST = 0.5         # Å，防原子重叠下界
WATER_BRIDGE_MAX_DONORS = 2         # 每水每帧最多供体数（2 对孤对）
WATER_BRIDGE_MAX_BRIDGES = 2        # 每水每帧最多桥数（2 个 H）
WATER_BRIDGE_PREFILTER = 3.0        # Å，候选预筛半径（H···O / H_w···A）
WATER_BRIDGE_OMEGA_MIN = 71.0       # °，A-Ow···H(供体) 最小角（"水在中间"约束）
WATER_BRIDGE_OMEGA_MAX = 140.0      # °，A-Ow···H(供体) 最大角


class WaterBridgeDetectorTwoPass(InteractionDetectorTwoPass):
    """水桥检测器（两轮遍历 + 稀疏存储）。"""

    @property
    def name(self) -> str:
        return "water_bridge"

    @property
    def required_group_types(self) -> List[str]:
        return ["H_donor", "water", "H_acceptor"]

    @property
    def metric_names(self) -> List[str]:
        return ["dist_h_ow", "angle_dh_ow", "dist_d_ow",
                "dist_hw_a", "angle_ow_hw_a", "dist_ow_a", "angle_h_ow_a"]

    # ==================== Pass1：KDTree 发现 ====================

    def run_pass1(self, groups: List[Group], trajectory,
                  tuple_filter=None) -> InteractionSparse:
        """执行 Pass1：KDTree 发现三元组 + 内联计算角度。"""

        # 分组：供体/受体排除水分子自身（按 water 组原子集合，跨力场通用，不依赖残基名）
        waters = [g for g in groups if g.group_type == "water"]
        water_atoms = set()
        for w in waters:
            water_atoms.update(w.atom_indices)
        donors = [g for g in groups
                  if g.group_type == "H_donor" and not (set(g.atom_indices) & water_atoms)]
        acceptors = [g for g in groups
                     if g.group_type == "H_acceptor" and not (set(g.atom_indices) & water_atoms)]

        if not donors or not waters or not acceptors:
            return InteractionSparse(interaction_type=self.name, data={})

        # 退化检查：判据依赖水的显式 H（见设计文档 §6）
        if any(sum(1 for a in w.atoms if a.atom_element == "H") < 2 for w in waters):
            print("[WARN] water_bridge: water has no explicit H, skipped")
            return InteractionSparse(interaction_type=self.name, data={})

        # 构建索引
        self._donor_d_idx = np.array([g.atoms[0].atom_global_idx for g in donors])
        self._donor_h_idx = np.array([g.atoms[1].atom_global_idx for g in donors])
        self._water_ow_idx = np.array([g.atoms[0].atom_global_idx for g in waters])
        self._acceptor_a_idx = np.array([g.atoms[0].atom_global_idx for g in acceptors])
        # 水的两个 H（形状 (n_water, 2)）
        self._water_hw_idx = np.array(
            [[w.atoms[1].atom_global_idx, w.atoms[2].atom_global_idx] for w in waters])

        sparse_data: Dict[Tuple[int, ...], dict] = {}

        for f, ts in enumerate(trajectory):
            positions = ts.positions
            triples = self._kdtree_discover(positions, get_box(ts))

            for idx in range(len(triples[0])):
                di = triples[0][idx]
                wi = triples[1][idx]
                ai = triples[2][idx]
                donor = donors[di]
                water = waters[wi]
                acceptor = acceptors[ai]

                if tuple_filter is not None and not tuple_filter((donor, water, acceptor)):
                    continue

                group_ids = (donor.group_id, water.group_id, acceptor.group_id)

                if group_ids not in sparse_data:
                    sparse_data[group_ids] = {
                        "groups": (donor, water, acceptor),
                        "frames": [],
                        "metrics": {name: [] for name in self.metric_names}
                    }

                sparse_data[group_ids]["frames"].append(f)
                # triples[3:] 依次对应 metric_names（见 _kdtree_discover）
                for k, name in enumerate(self.metric_names, start=3):
                    sparse_data[group_ids]["metrics"][name].append(
                        float(triples[k][idx]))

        return InteractionSparse(interaction_type=self.name, data=sparse_data)

    # ==================== Pass2：从 InteractionSparse 重建 ====================

    def _get_box(self, ts):
        """PBC：返回当前帧盒尺寸（None=无盒退化欧氏）。"""
        return get_box(ts)

    # ==================== compute_pair_metrics ====================

    def compute_pair_metrics(self, group_tuples, all_positions, box=None):
        """向量化计算 6 个指标 + hw_id（使用缓存索引）。

        受体侧对水的两个 H 各算一次，取角度更优者记录 hw_id。
        box: 盒尺寸（MIC）；None = 欧氏距离（无 PBC）。
        """
        if not hasattr(self, '_cached_d_idx'):
            raise RuntimeError(
                "compute_pair_metrics requires cached indices. "
                "Call _build_indices_from_sparse first via run_pass2.")

        from MDAnalysis.lib.distances import calc_bonds, calc_angles

        d = all_positions[self._cached_d_idx]
        h = all_positions[self._cached_h_idx]
        ow = all_positions[self._cached_ow_idx]
        a = all_positions[self._cached_a_idx]

        # 供体侧
        dist_h_ow = calc_bonds(h, ow, box=box)
        dist_d_ow = calc_bonds(d, ow, box=box)
        angle_dh_ow = np.degrees(calc_angles(d, h, ow, box=box))

        # 受体侧：HW1 / HW2 取角度更优者
        hw1 = all_positions[self._cached_hw1_idx]
        hw2 = all_positions[self._cached_hw2_idx]
        dist_hw1 = calc_bonds(hw1, a, box=box)
        dist_hw2 = calc_bonds(hw2, a, box=box)
        ang1 = np.degrees(calc_angles(ow, hw1, a, box=box))
        ang2 = np.degrees(calc_angles(ow, hw2, a, box=box))
        use2 = ang2 > ang1

        return {
            "dist_h_ow": dist_h_ow,
            "angle_dh_ow": angle_dh_ow,
            "dist_d_ow": dist_d_ow,
            "dist_hw_a": np.where(use2, dist_hw2, dist_hw1),
            "angle_ow_hw_a": np.where(use2, ang2, ang1),
            "dist_ow_a": calc_bonds(ow, a, box=box),
            "hw_id": use2.astype(float),
            "angle_h_ow_a": np.degrees(calc_angles(h, ow, a, box=box)),
        }

    def apply_threshold(self, metrics):
        """H···A ≤ 2.5 Å 且 D-H···A ≥ 120°（两侧各一，见设计文档 §2.2）。"""
        dw_ok = (metrics["dist_h_ow"] > WATER_BRIDGE_MIN_DIST) & \
                (metrics["dist_h_ow"] <= WATER_BRIDGE_HB_DIST_MAX)
        wa_ok = (metrics["dist_hw_a"] > WATER_BRIDGE_MIN_DIST) & \
                (metrics["dist_hw_a"] <= WATER_BRIDGE_HB_DIST_MAX)
        theta_ok = metrics["angle_dh_ow"] >= WATER_BRIDGE_HB_ANGLE_MIN
        omega_ok = metrics["angle_ow_hw_a"] >= WATER_BRIDGE_HB_ANGLE_MIN
        mid_ok = (metrics["angle_h_ow_a"] >= WATER_BRIDGE_OMEGA_MIN) & \
                 (metrics["angle_h_ow_a"] <= WATER_BRIDGE_OMEGA_MAX)
        return dw_ok & wa_ok & theta_ok & omega_ok & mid_ok

    def _post_threshold(self, ok, frame_metrics, group_tuples):
        """逐帧三层配额（TwoPass 独立实现，与 PerFrame 逻辑等价）。

        第 2 层按供体去重后，第 3 层以「受体侧 (H_w,A)」×「供体侧 D」组合择优取 2，
        从而保留"同一供体 + 水的两个 H 指向两个受体"的 2 个合法桥。
        """
        water_uid = np.array([gt[1].group_id for gt in group_tuples])
        donor_uid = np.array([gt[0].group_id for gt in group_tuples])
        idx = np.where(ok)[0]
        if len(idx) == 0:
            return ok
        ang = frame_metrics["angle_ow_hw_a"]
        ang_d = frame_metrics["angle_dh_ow"]
        dh = frame_metrics["dist_h_ow"]
        hw = frame_metrics["hw_id"]
        keep = np.zeros_like(ok)
        for w in np.unique(water_uid[idx]):
            cand = idx[water_uid[idx] == w]
            by_h = {}
            for i in cand:
                h = hw[i]
                if h not in by_h or ang[i] > ang[by_h[h]]:
                    by_h[h] = i
            by_d = {}
            for i in cand:
                d = donor_uid[i]
                if d not in by_d or dh[i] < dh[by_d[d]]:
                    by_d[d] = i
            donors_kept = sorted(by_d.values(), key=lambda i: dh[i])[:WATER_BRIDGE_MAX_DONORS]
            combos = [(i, j) for i in by_h.values() for j in donors_kept]
            combos.sort(key=lambda p: -(ang[p[0]] + ang_d[p[1]]))
            for i, _ in combos[:WATER_BRIDGE_MAX_BRIDGES]:
                keep[i] = True
        return keep

    # ==================== 内部辅助方法 ====================

    def _kdtree_discover(self, positions, box=None):
        """PBC 感知发现三元组 + 向量化计算指标（基于显式 H）。

        box: 盒尺寸（MIC）；None = 欧氏距离（无 PBC）。

        Returns:
            (di, wi, ai, dist_h_ow, angle_dh_ow, dist_d_ow,
             dist_hw_a, angle_ow_hw_a, dist_ow_a, hw_id)
            每个都是 numpy 数组，长度 = 满足阈值的三元组数量。
        """
        from MDAnalysis.lib.distances import (
            capped_distance, calc_bonds, calc_angles)

        w_coords = positions[self._water_ow_idx]
        h_coords = positions[self._donor_h_idx]
        a_coords = positions[self._acceptor_a_idx]
        hw_coords = positions[self._water_hw_idx.ravel()]   # 2n 个水 H

        # 供体 H ↔ 水 O；水 H ↔ 受体 A
        ho_pairs = capped_distance(
            h_coords, w_coords, WATER_BRIDGE_PREFILTER,
            box=box, method='pkdtree', return_distances=False)
        hwa_pairs = capped_distance(
            hw_coords, a_coords, WATER_BRIDGE_PREFILTER,
            box=box, method='pkdtree', return_distances=False)

        n_out = 10
        if len(ho_pairs) == 0 or len(hwa_pairs) == 0:
            return tuple(np.array([]) for _ in range(n_out))

        # 按水分组：供体侧（H→O_w 的水索引 = pairs[:,1]）、受体侧（H_w→A 的水索引）
        d_by_w = defaultdict(set)
        for idx in range(len(ho_pairs)):
            d_by_w[ho_pairs[idx, 1]].add(ho_pairs[idx, 0])
        a_by_w = defaultdict(set)
        for idx in range(len(hwa_pairs)):
            a_by_w[hwa_pairs[idx, 0] // 2].add(hwa_pairs[idx, 1])

        common_w = set(d_by_w) & set(a_by_w)
        if not common_w:
            return tuple(np.array([]) for _ in range(n_out))

        all_di, all_wi, all_ai = [], [], []
        for wi in common_w:
            dl = np.array(sorted(d_by_w[wi]))
            al = np.array(sorted(a_by_w[wi]))
            n_d, n_a = len(dl), len(al)
            all_di.append(np.repeat(dl, n_a))
            all_wi.append(np.full(n_d * n_a, wi))
            all_ai.append(np.tile(al, n_d))

        di_arr = np.concatenate(all_di)
        wi_arr = np.concatenate(all_wi)
        ai_arr = np.concatenate(all_ai)

        # 坐标
        d = positions[self._donor_d_idx[di_arr]]
        h = h_coords[di_arr]
        ow = w_coords[wi_arr]
        a = a_coords[ai_arr]

        # 供体侧
        dist_h_ow = calc_bonds(h, ow, box=box)
        dist_d_ow = calc_bonds(d, ow, box=box)
        angle_dh_ow = np.degrees(calc_angles(d, h, ow, box=box))

        # 受体侧：HW1 / HW2 取角度更优者
        hw1 = positions[self._water_hw_idx[wi_arr, 0]]
        hw2 = positions[self._water_hw_idx[wi_arr, 1]]
        dist_hw1 = calc_bonds(hw1, a, box=box)
        dist_hw2 = calc_bonds(hw2, a, box=box)
        ang1 = np.degrees(calc_angles(ow, hw1, a, box=box))
        ang2 = np.degrees(calc_angles(ow, hw2, a, box=box))
        use2 = ang2 > ang1
        dist_hw_a = np.where(use2, dist_hw2, dist_hw1)
        angle_ow_hw_a = np.where(use2, ang2, ang1)
        dist_ow_a = calc_bonds(ow, a, box=box)
        hw_id = use2.astype(float)
        angle_h_ow_a = np.degrees(calc_angles(h, ow, a, box=box))

        # 阈值筛选
        mask = (dist_h_ow > WATER_BRIDGE_MIN_DIST) & \
               (dist_h_ow <= WATER_BRIDGE_HB_DIST_MAX) & \
               (dist_hw_a > WATER_BRIDGE_MIN_DIST) & \
               (dist_hw_a <= WATER_BRIDGE_HB_DIST_MAX) & \
               (angle_dh_ow >= WATER_BRIDGE_HB_ANGLE_MIN) & \
               (angle_ow_hw_a >= WATER_BRIDGE_HB_ANGLE_MIN) & \
               (angle_h_ow_a >= WATER_BRIDGE_OMEGA_MIN) & \
               (angle_h_ow_a <= WATER_BRIDGE_OMEGA_MAX)

        # 返回顺序：3 个索引 + 7 个 metric（与 metric_names 顺序一致）+ hw_id（末尾，非 metric）
        return (di_arr[mask], wi_arr[mask], ai_arr[mask],
                dist_h_ow[mask], angle_dh_ow[mask], dist_d_ow[mask],
                dist_hw_a[mask], angle_ow_hw_a[mask], dist_ow_a[mask],
                angle_h_ow_a[mask], hw_id[mask])

    def _build_indices_from_sparse(self, sparse: InteractionSparse):
        """从 InteractionSparse 重建缓存索引（Pass2 用）。"""
        group_tuples = [entry["groups"] for entry in sparse.data.values()]

        self._cached_d_idx = np.array(
            [gt[0].atoms[0].atom_global_idx for gt in group_tuples])
        self._cached_h_idx = np.array(
            [gt[0].atoms[1].atom_global_idx for gt in group_tuples])
        self._cached_ow_idx = np.array(
            [gt[1].atoms[0].atom_global_idx for gt in group_tuples])
        self._cached_a_idx = np.array(
            [gt[2].atoms[0].atom_global_idx for gt in group_tuples])
        # 水的两个 H（供/受体侧均需要）
        self._cached_hw1_idx = np.array(
            [gt[1].atoms[1].atom_global_idx for gt in group_tuples])
        self._cached_hw2_idx = np.array(
            [gt[1].atoms[2].atom_global_idx for gt in group_tuples])
