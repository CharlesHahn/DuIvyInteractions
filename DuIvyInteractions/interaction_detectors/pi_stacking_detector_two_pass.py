# -*- coding: utf-8 -*-
"""π-π 堆积检测器（两轮遍历 + 稀疏存储）。

判据：环心距离 ≤ 5.5 Å，法向量夹角满足 P 型或 T 型，投影偏移 ≤ 2.0 Å。
参考：PLIP, McGaughey 1998。
"""

import numpy as np
from typing import List, Tuple, Dict

from ..core.interfaces import InteractionDetectorTwoPass
from ..core.datas import Group, Interaction, InteractionSparse


# PLIP 阈值
PISTACK_DIST_MAX = 5.5       # Å
PISTACK_ANG_DEV = 30.0       # °
PISTACK_OFFSET_MAX = 2.0     # Å
PISTACK_PLANARITY = 5.0      # °
PISTACK_MIN_DIST = 0.5       # Å


class PiStackingDetectorTwoPass(InteractionDetectorTwoPass):
    """π-π 堆积检测器（两轮遍历 + 稀疏存储）。"""

    def __init__(self, check_planarity: bool = False):
        self.check_planarity = check_planarity

    @property
    def name(self) -> str:
        return "pi_stacking"

    @property
    def required_group_types(self) -> List[str]:
        return ["aromatic_ring"]

    @property
    def metric_names(self) -> List[str]:
        names = ["distance", "angle", "offset", "pistacking_type"]
        if self.check_planarity:
            names.extend(["planarity_ring1", "planarity_ring2"])
        return names

    # ==================== initialize_candidates ====================

    def initialize_candidates(self, groups, trajectory, tuple_filter=None):
        """筛选环 → 环形 padding → 组合 C(n,2) → tuple_filter → 设置缓存索引。"""
        rings = [g for g in groups if g.group_type == "aromatic_ring"]
        if len(rings) < 2:
            return []

        # 环形 padding
        ring_idx, ring_prev, ring_next, ring_valid, ring_n = \
            self._build_circular_padding(rings)

        # 组合 C(n,2)
        n = len(rings)
        r1_list, r2_list = [], []
        for i in range(n):
            for j in range(i + 1, n):
                r1_list.append(i)
                r2_list.append(j)
        pair_r1 = np.array(r1_list)
        pair_r2 = np.array(r2_list)

        # tuple_filter
        if tuple_filter is not None:
            mask = np.array([
                tuple_filter((rings[r1], rings[r2]))
                for r1, r2 in zip(pair_r1, pair_r2)
            ])
            pair_r1 = pair_r1[mask]
            pair_r2 = pair_r2[mask]

        if len(pair_r1) == 0:
            return []

        # 设置缓存索引
        self._cached_ring_idx = ring_idx
        self._cached_ring_prev = ring_prev
        self._cached_ring_next = ring_next
        self._cached_ring_valid = ring_valid
        self._cached_ring_n = ring_n
        self._cached_pair_r1 = pair_r1
        self._cached_pair_r2 = pair_r2

        return [(rings[pair_r1[i]], rings[pair_r2[i]])
                for i in range(len(pair_r1))]

    # ==================== compute_pair_metrics ====================

    def compute_pair_metrics(self, group_tuples, all_positions):
        """向量化计算距离、角度、投影偏移、分类（使用缓存索引）。"""
        if not hasattr(self, '_cached_ring_idx'):
            raise RuntimeError(
                "compute_pair_metrics requires cached indices. "
                "Call initialize_candidates or _build_indices_from_sparse first.")

        r1 = self._cached_pair_r1
        r2 = self._cached_pair_r2

        # 环心
        coords = all_positions[self._cached_ring_idx]
        masked = coords * self._cached_ring_valid[:, :, None]
        centers = np.sum(masked, axis=1) / \
                  np.sum(self._cached_ring_valid, axis=1)[:, None]
        c1 = centers[r1]
        c2 = centers[r2]

        # 环心距
        distance = np.linalg.norm(c1 - c2, axis=1)

        # 法向量
        prev_coords = all_positions[self._cached_ring_prev]
        next_coords = all_positions[self._cached_ring_next]
        atom_normals = np.cross(coords - prev_coords, coords - next_coords)
        n_masked = atom_normals * self._cached_ring_valid[:, :, None]
        normal = np.sum(n_masked, axis=1) / \
                 np.sum(self._cached_ring_valid, axis=1)[:, None]
        norm = np.linalg.norm(normal, axis=1, keepdims=True)
        normals = normal / norm

        n1 = normals[r1]
        n2 = normals[r2]

        # 法向量夹角
        cos_angle = np.sum(n1 * n2, axis=1)
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        raw_angle = np.degrees(np.arccos(cos_angle))
        angle = np.minimum(raw_angle, 180.0 - raw_angle)

        # 投影偏移（双向取较小值）
        offset1 = self._projection_distance(n2, c2, c1)
        offset2 = self._projection_distance(n1, c1, c2)
        offset = np.minimum(offset1, offset2)

        # 分类
        pistacking_type = self._classify(angle, offset)

        result = {
            "distance": distance,
            "angle": angle,
            "offset": offset,
            "pistacking_type": pistacking_type,
        }

        # 平面性
        if self.check_planarity:
            planarity = self._ring_planarity(atom_normals)
            result["planarity_ring1"] = planarity[r1]
            result["planarity_ring2"] = planarity[r2]

        return result

    def apply_threshold(self, metrics):
        """距离 + 分类 + 可选平面性。"""
        base = (metrics["distance"] > PISTACK_MIN_DIST) & \
               (metrics["distance"] <= PISTACK_DIST_MAX) & \
               (metrics["pistacking_type"] != 'N')
        if self.check_planarity:
            planar = (metrics["planarity_ring1"] <= PISTACK_PLANARITY) & \
                     (metrics["planarity_ring2"] <= PISTACK_PLANARITY)
            return base & planar
        return base

    # ==================== Pass2（覆盖基类） ====================

    def run_pass2(self, sparse: InteractionSparse,
                  trajectory) -> List[Interaction]:
        """执行 Pass2：重建缓存索引 + 处理 pistacking_type 的 char dtype。"""
        if not sparse.data:
            return []

        self._build_indices_from_sparse(sparse)

        group_tuples = [entry["groups"] for entry in sparse.data.values()]
        n_groups = len(group_tuples)
        n_frames = trajectory.n_frames

        # 预分配（pistacking_type 用 'U1' + 'N'，其他用 float + NaN）
        existence = np.zeros((n_groups, n_frames), dtype=bool)
        metrics = {}
        for name in self.metric_names:
            if name == "pistacking_type":
                metrics[name] = np.full((n_groups, n_frames), 'N', dtype='U1')
            else:
                metrics[name] = np.full((n_groups, n_frames), np.nan)

        # 逐帧计算
        for f, ts in enumerate(trajectory):
            frame_metrics = self.compute_pair_metrics(group_tuples, ts.positions)
            existence[:, f] = self.apply_threshold(frame_metrics)
            for name in self.metric_names:
                metrics[name][:, f] = frame_metrics[name]

        # 构建 results
        results = [
            (entry["groups"], existence[i],
             {k: v[i] for k, v in metrics.items()})
            for i, entry in enumerate(sparse.data.values())
        ]

        results = self._post_process(results)
        return self._build_interaction(results)

    # ==================== 内部辅助方法 ====================

    def _build_indices_from_sparse(self, sparse: InteractionSparse):
        """从 InteractionSparse 重建缓存索引（Pass2 用）。"""
        group_tuples = [entry["groups"] for entry in sparse.data.values()]
        rings = [gt[0] for gt in group_tuples] + [gt[1] for gt in group_tuples]
        # 去重保持顺序
        seen = set()
        unique_rings = []
        for r in rings:
            if r.group_id not in seen:
                seen.add(r.group_id)
                unique_rings.append(r)

        ring_idx, ring_prev, ring_next, ring_valid, ring_n = \
            self._build_circular_padding(unique_rings)
        ring_id_to_local = {r.group_id: i for i, r in enumerate(unique_rings)}

        pair_r1 = np.array([ring_id_to_local[gt[0].group_id] for gt in group_tuples])
        pair_r2 = np.array([ring_id_to_local[gt[1].group_id] for gt in group_tuples])

        self._cached_ring_idx = ring_idx
        self._cached_ring_prev = ring_prev
        self._cached_ring_next = ring_next
        self._cached_ring_valid = ring_valid
        self._cached_ring_n = ring_n
        self._cached_pair_r1 = pair_r1
        self._cached_pair_r2 = pair_r2

    @staticmethod
    def _build_circular_padding(rings: List[Group]):
        """构建环形 padding + 显式邻居索引。"""
        n = len(rings)
        max_atoms = max(len(r.atoms) for r in rings)
        indices = np.zeros((n, max_atoms), dtype=int)
        prev_idx = np.zeros((n, max_atoms), dtype=int)
        next_idx = np.zeros((n, max_atoms), dtype=int)
        valid = np.zeros((n, max_atoms), dtype=bool)
        n_atoms = np.zeros(n, dtype=int)

        for i, ring in enumerate(rings):
            na = len(ring.atoms)
            n_atoms[i] = na
            atom_indices = ring.atom_indices
            for j in range(max_atoms):
                indices[i, j] = atom_indices[j % na]
            for j in range(max_atoms):
                prev_idx[i, j] = atom_indices[(j - 1) % na]
                next_idx[i, j] = atom_indices[(j + 1) % na]
            valid[i, :na] = True

        return indices, prev_idx, next_idx, valid, n_atoms

    @staticmethod
    def _projection_distance(normal, plane_point, target_point):
        """投影距离。"""
        d1 = np.linalg.norm(target_point - (plane_point + normal), axis=1)
        d2 = np.linalg.norm(target_point - (plane_point - normal), axis=1)
        sign = np.where(d1 < d2, 1.0, -1.0)[:, None]
        oriented_normal = normal * sign
        t = target_point - plane_point
        proj_dist = np.sum(t * oriented_normal, axis=1)
        proj_point = target_point - proj_dist[:, None] * oriented_normal
        return np.linalg.norm(proj_point - plane_point, axis=1)

    @staticmethod
    def _classify(angle, offset):
        """分类 P 型 / T 型 / N 型。"""
        pistacking_type = np.full(angle.shape[0], 'N', dtype='U1')
        p_mask = (angle <= PISTACK_ANG_DEV) & (offset < PISTACK_OFFSET_MAX)
        t_mask = (angle >= 90 - PISTACK_ANG_DEV) & (offset < PISTACK_OFFSET_MAX)
        pistacking_type[p_mask] = 'P'
        pistacking_type[t_mask] = 'T'
        return pistacking_type

    def _ring_planarity(self, normals_all):
        """计算每环的平面性（法向量两两最大夹角）。"""
        n_rings = normals_all.shape[0]
        max_atoms = normals_all.shape[1]
        max_angle = np.zeros(n_rings)
        valid = self._cached_ring_valid

        for i in range(max_atoms):
            for j in range(i + 1, max_atoms):
                both_valid = valid[:, i] & valid[:, j]
                if not np.any(both_valid):
                    continue
                cos = np.sum(normals_all[:, i, :] * normals_all[:, j, :], axis=1)
                norm_i = np.linalg.norm(normals_all[:, i, :], axis=1)
                norm_j = np.linalg.norm(normals_all[:, j, :], axis=1)
                cos = np.clip(cos / (norm_i * norm_j), -1.0, 1.0)
                ang = np.degrees(np.arccos(cos))
                max_angle = np.where(both_valid, np.maximum(max_angle, ang), max_angle)

        return max_angle
