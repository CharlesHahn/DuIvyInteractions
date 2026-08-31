# -*- coding: utf-8 -*-
"""π-π 堆积检测器（逐帧策略，环形 padding 向量化）。

判据：环心距离 ≤ 5.5 Å，法向量夹角满足 P 型或 T 型，投影偏移 ≤ 2.0 Å。
参考：PLIP, McGaughey 1998。
"""

import numpy as np
from typing import List, Tuple, Dict

from ..core.interfaces import InteractionDetectorPerFrame
from ..core.datas import Group, Interaction


# PLIP 阈值
PISTACK_DIST_MAX = 5.5       # Å，环心最大距离
PISTACK_ANG_DEV = 30.0       # °，角度偏差上限
PISTACK_OFFSET_MAX = 2.0     # Å，最大投影偏移
PISTACK_PLANARITY = 5.0      # °，平面性检验阈值
PISTACK_MIN_DIST = 0.5       # Å，最小距离（排除自身）
PREFILTER_CUTOFF = PISTACK_DIST_MAX * 3  # 16.5 Å


class PiStackingDetectorPerFrame(InteractionDetectorPerFrame):
    """π-π 堆积检测器（逐帧策略，环形 padding 向量化）。"""

    def __init__(self, check_planarity: bool = False):
        self.check_planarity = check_planarity
        self._ring_indices = None
        self._ring_prev = None
        self._ring_next = None
        self._ring_valid = None
        self._ring_n = None
        self._pair_r1 = None
        self._pair_r2 = None

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

    # ==================== 抽象方法（基类要求，此处不使用） ====================

    def get_candidate_tuples(self, groups, coordinates=None):
        raise NotImplementedError("detect() is overridden, this method is not used")

    def compute_metrics_for_frame(self, tuples, all_positions, frame):
        raise NotImplementedError("detect() is overridden, this method is not used")

    def apply_threshold(self, metrics):
        base = (metrics["distance"] > PISTACK_MIN_DIST) & \
               (metrics["distance"] <= PISTACK_DIST_MAX) & \
               (metrics["pistacking_type"] != 'N')
        if self.check_planarity:
            planar = (metrics["planarity_ring1"] <= PISTACK_PLANARITY) & \
                     (metrics["planarity_ring2"] <= PISTACK_PLANARITY)
            return base & planar
        return base

    # ==================== 重写 detect ====================

    def detect(self, groups, trajectory=None,
               n_workers: int = 1,
               topology_path: str = None,
               trajectory_path: str = None,
               tuple_filter=None) -> List[Interaction]:
        """π-π 堆积检测主流程。"""
        if trajectory is None:
            raise ValueError("trajectory is required")
        if n_workers > 1:
            raise NotImplementedError("PerFrame 检测器暂不支持并行")

        # 1. 筛选芳香环
        rings = [g for g in groups if g.group_type == "aromatic_ring"]
        if len(rings) < 2:
            return []

        # 2. 构建环形 padding
        self._ring_indices, self._ring_prev, self._ring_next, self._ring_valid, self._ring_n = \
            self._build_circular_padding(rings)

        # 3. 候选环对（组合）
        n_rings = len(rings)
        r1_list, r2_list = [], []
        for i in range(n_rings):
            for j in range(i + 1, n_rings):
                r1_list.append(i)
                r2_list.append(j)
        pair_r1 = np.array(r1_list)
        pair_r2 = np.array(r2_list)

        # 4. tuple_filter
        if tuple_filter is not None:
            mask = np.array([
                tuple_filter((rings[r1], rings[r2]))
                for r1, r2 in zip(pair_r1, pair_r2)
            ])
            pair_r1 = pair_r1[mask]
            pair_r2 = pair_r2[mask]

        if len(pair_r1) == 0:
            return []

        # 5. 第一帧：距离预过滤
        first_pos = trajectory[0].positions
        centers = self._ring_centers(first_pos)
        dist = np.linalg.norm(centers[pair_r1] - centers[pair_r2], axis=1)
        mask = (dist > PISTACK_MIN_DIST) & (dist < PREFILTER_CUTOFF)

        # 平面性预过滤（可选）
        if self.check_planarity:
            normals_all = self._ring_normals(first_pos)
            planarity = self._ring_planarity(normals_all)
            mask = mask & (planarity[pair_r1] <= PISTACK_PLANARITY) & \
                          (planarity[pair_r2] <= PISTACK_PLANARITY)

        pair_r1 = pair_r1[mask]
        pair_r2 = pair_r2[mask]

        if len(pair_r1) == 0:
            return []

        self._pair_r1 = pair_r1
        self._pair_r2 = pair_r2

        n_pairs = len(pair_r1)
        n_frames = trajectory.n_frames

        # 6. 预分配结果数组
        existence = np.zeros((n_pairs, n_frames), dtype=bool)
        distance = np.zeros((n_pairs, n_frames))
        angle = np.zeros((n_pairs, n_frames))
        offset = np.zeros((n_pairs, n_frames))
        pistacking_type = np.full((n_pairs, n_frames), 'N', dtype='U1')
        if self.check_planarity:
            planarity1 = np.zeros((n_pairs, n_frames))
            planarity2 = np.zeros((n_pairs, n_frames))

        # 7. 逐帧计算
        for f, ts in enumerate(trajectory):
            positions = ts.positions
            metrics = self._compute_metrics(positions)
            distance[:, f] = metrics["distance"]
            angle[:, f] = metrics["angle"]
            offset[:, f] = metrics["offset"]
            pistacking_type[:, f] = metrics["pistacking_type"]
            existence[:, f] = self.apply_threshold(metrics)
            if self.check_planarity:
                planarity1[:, f] = metrics["planarity_ring1"]
                planarity2[:, f] = metrics["planarity_ring2"]

        # 8. 过滤从未存在的 pair
        has_any = np.any(existence, axis=1)
        if not np.any(has_any):
            return []

        result_metrics = {
            "distance": distance[has_any],
            "angle": angle[has_any],
            "offset": offset[has_any],
            "pistacking_type": pistacking_type[has_any],
        }
        if self.check_planarity:
            result_metrics["planarity_ring1"] = planarity1[has_any]
            result_metrics["planarity_ring2"] = planarity2[has_any]

        tuples = [(rings[pair_r1[i]], rings[pair_r2[i]])
                  for i in range(n_pairs) if has_any[i]]

        return [Interaction(
            interaction_type=self.name,
            groups=tuples,
            existence=existence[has_any],
            metrics=result_metrics
        )]

    # ==================== 环形 padding 构建 ====================

    @staticmethod
    def _build_circular_padding(rings: List[Group]):
        """构建环形 padding 索引矩阵 + 显式邻居索引。

        对每个环，循环填充原子索引到 max_atoms，
        同时构建 prev/next 邻居的全局原子索引。

        Returns:
            indices: (n_rings, max_atoms) int
            prev_idx: (n_rings, max_atoms) int — 前邻居的全局原子索引
            next_idx: (n_rings, max_atoms) int — 后邻居的全局原子索引
            valid: (n_rings, max_atoms) bool
            n_atoms: (n_rings,) int
        """
        n_rings = len(rings)
        max_atoms = max(len(ring.atoms) for ring in rings)
        indices = np.zeros((n_rings, max_atoms), dtype=int)
        prev_idx = np.zeros((n_rings, max_atoms), dtype=int)
        next_idx = np.zeros((n_rings, max_atoms), dtype=int)
        valid = np.zeros((n_rings, max_atoms), dtype=bool)
        n_atoms = np.zeros(n_rings, dtype=int)

        for i, ring in enumerate(rings):
            n = len(ring.atoms)
            n_atoms[i] = n
            atom_indices = ring.atom_indices
            # 循环填充原子索引
            for j in range(max_atoms):
                indices[i, j] = atom_indices[j % n]
            # 显式邻居索引（全局原子索引）
            for j in range(max_atoms):
                prev_j = (j - 1) % n
                next_j = (j + 1) % n
                prev_idx[i, j] = atom_indices[prev_j]
                next_idx[i, j] = atom_indices[next_j]
            valid[i, :n] = True

        return indices, prev_idx, next_idx, valid, n_atoms

    # ==================== 向量化计算 ====================

    def _ring_centers(self, positions):
        """计算全部环的几何中心。

        Returns:
            centers: (n_rings, 3)
        """
        coords = positions[self._ring_indices]  # (n_rings, max_atoms, 3)
        masked = coords * self._ring_valid[:, :, None]
        centers = np.sum(masked, axis=1) / np.sum(self._ring_valid, axis=1)[:, None]
        return centers

    def _ring_normals(self, positions):
        """计算全部环的法向量。

        使用显式邻居索引（非 np.roll），对任意环大小正确。

        Returns:
            normals: (n_rings, 3) 单位向量
        """
        coords = positions[self._ring_indices]   # (n_rings, max_atoms, 3)
        prev_coords = positions[self._ring_prev]  # (n_rings, max_atoms, 3)
        next_coords = positions[self._ring_next]  # (n_rings, max_atoms, 3)
        atom_normals = np.cross(coords - prev_coords, coords - next_coords)
        masked = atom_normals * self._ring_valid[:, :, None]
        normal = np.sum(masked, axis=1) / np.sum(self._ring_valid, axis=1)[:, None]
        norm = np.linalg.norm(normal, axis=1, keepdims=True)
        return normal / norm

    def _ring_planarity(self, normals_all):
        """计算每环的平面性（法向量两两最大夹角）。

        Args:
            normals_all: (n_rings, max_atoms, 3) 每个位置的法向量

        Returns:
            planarity: (n_rings,) 最大夹角（°）
        """
        n_rings = normals_all.shape[0]
        max_atoms = normals_all.shape[1]
        max_angle = np.zeros(n_rings)
        valid = self._ring_valid  # (n_rings, max_atoms)

        for i in range(max_atoms):
            for j in range(i + 1, max_atoms):
                both_valid = valid[:, i] & valid[:, j]
                if not np.any(both_valid):
                    continue
                cos = np.sum(normals_all[:, i, :] * normals_all[:, j, :], axis=1)
                norm_i = np.linalg.norm(normals_all[:, i, :], axis=1)
                norm_j = np.linalg.norm(normals_all[:, j, :], axis=1)
                cos = cos / (norm_i * norm_j)
                cos = np.clip(cos, -1.0, 1.0)
                ang = np.degrees(np.arccos(cos))
                max_angle = np.where(both_valid, np.maximum(max_angle, ang), max_angle)

        return max_angle

    def _compute_metrics(self, positions):
        """向量化计算全部候选环对的指标。"""
        r1 = self._pair_r1
        r2 = self._pair_r2

        # 环心
        centers = self._ring_centers(positions)  # (n_rings, 3)
        c1 = centers[r1]  # (n_pairs, 3)
        c2 = centers[r2]

        # 环心距
        distance = np.linalg.norm(c1 - c2, axis=1)

        # 法向量
        normals = self._ring_normals(positions)  # (n_rings, 3)
        n1 = normals[r1]
        n2 = normals[r2]

        # 法向量夹角（取最小角，处理方向不确定性）
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
            coords = positions[self._ring_indices]
            prev_coords = positions[self._ring_prev]
            next_coords = positions[self._ring_next]
            atom_normals = np.cross(coords - prev_coords, coords - next_coords)
            planarity_all = self._ring_planarity(atom_normals)
            result["planarity_ring1"] = planarity_all[r1]
            result["planarity_ring2"] = planarity_all[r2]

        return result

    # ==================== 辅助方法 ====================

    @staticmethod
    def _projection_distance(normal, plane_point, target_point):
        """将 target_point 投影到 (normal, plane_point) 平面，返回投影点到 plane_point 的距离。

        Args:
            normal: (n, 3) 平面法向量
            plane_point: (n, 3) 平面上一点
            target_point: (n, 3) 待投影的点
        Returns:
            dist: (n,)
        """
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
