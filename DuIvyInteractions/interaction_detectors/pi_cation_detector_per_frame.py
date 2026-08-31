# -*- coding: utf-8 -*-
"""π-阳离子相互作用检测器（逐帧策略，向量化）。

判据：电荷中心到环心距离 < 6.0 Å，投影偏移 < 2.0 Å。
参考：PLIP, Gallivan and Dougherty 1999。
"""

import numpy as np
from typing import List, Tuple, Dict

from ..core.interfaces import InteractionDetectorPerFrame
from ..core.datas import Group, Interaction


# PLIP 阈值
PICATION_DIST_MAX = 6.0    # Å，电荷中心到环心最大距离
PICATION_OFFSET_MAX = 2.0  # Å，最大投影偏移
PICATION_MIN_DIST = 0.5    # Å，最小距离
PREFILTER_CUTOFF = PICATION_DIST_MAX * 3  # 18.0 Å


class PiCationDetectorPerFrame(InteractionDetectorPerFrame):
    """π-阳离子相互作用检测器（逐帧策略，向量化）。"""

    def __init__(self):
        self._ring_indices = None
        self._ring_prev = None
        self._ring_next = None
        self._ring_valid = None
        self._cation_indices = None
        self._cation_charges = None
        self._cation_valid = None
        self._pair_ring = None
        self._pair_cation = None

    @property
    def name(self) -> str:
        return "pi_cation"

    @property
    def required_group_types(self) -> List[str]:
        return ["aromatic_ring", "charged_positive"]

    @property
    def metric_names(self) -> List[str]:
        return ["distance", "offset"]

    # ==================== 抽象方法（基类要求，此处不使用） ====================

    def get_candidate_tuples(self, groups, coordinates=None):
        raise NotImplementedError("detect() is overridden, this method is not used")

    def compute_metrics_for_frame(self, tuples, all_positions, frame):
        raise NotImplementedError("detect() is overridden, this method is not used")

    def apply_threshold(self, metrics):
        return (metrics["distance"] > PICATION_MIN_DIST) & \
               (metrics["distance"] < PICATION_DIST_MAX) & \
               (metrics["offset"] < PICATION_OFFSET_MAX)

    # ==================== 重写 detect ====================

    def detect(self, groups, trajectory=None,
               n_workers: int = 1,
               topology_path: str = None,
               trajectory_path: str = None,
               tuple_filter=None) -> List[Interaction]:
        """π-阳离子检测主流程。"""
        if trajectory is None:
            raise ValueError("trajectory is required")
        if n_workers > 1:
            raise NotImplementedError("PerFrame 检测器暂不支持并行")

        # 1. 分组
        rings = [g for g in groups if g.group_type == "aromatic_ring"]
        cations = [g for g in groups if g.group_type == "charged_positive"]

        if not rings or not cations:
            return []

        # 2. 构建数据结构
        # 环侧：环形 padding + 显式邻居索引（同 PiStacking）
        self._ring_indices, self._ring_prev, self._ring_next, self._ring_valid, _ = \
            self._build_circular_padding(rings)

        # 阳离子侧：padding 矩阵（同 SaltBridge）
        self._cation_indices, self._cation_charges, self._cation_valid = \
            self._build_padding(cations)

        # 3. 笛卡尔积
        n_rings = len(rings)
        n_cations = len(cations)
        r_grid, c_grid = np.meshgrid(np.arange(n_rings), np.arange(n_cations), indexing='ij')
        pair_ring = r_grid.ravel()
        pair_cation = c_grid.ravel()

        # 4. tuple_filter
        if tuple_filter is not None:
            mask = np.array([
                tuple_filter((rings[ri], cations[ci]))
                for ri, ci in zip(pair_ring, pair_cation)
            ])
            pair_ring = pair_ring[mask]
            pair_cation = pair_cation[mask]

        if len(pair_ring) == 0:
            return []

        # 5. 第一帧：距离预过滤
        first_pos = trajectory[0].positions
        ring_centers = self._ring_centers(first_pos)
        cation_centers = self._cation_centers(first_pos)
        dist = np.linalg.norm(ring_centers[pair_ring] - cation_centers[pair_cation], axis=1)
        mask = (dist > PICATION_MIN_DIST) & (dist < PREFILTER_CUTOFF)

        pair_ring = pair_ring[mask]
        pair_cation = pair_cation[mask]

        if len(pair_ring) == 0:
            return []

        self._pair_ring = pair_ring
        self._pair_cation = pair_cation

        n_pairs = len(pair_ring)
        n_frames = trajectory.n_frames

        # 6. 预分配结果数组
        existence = np.zeros((n_pairs, n_frames), dtype=bool)
        distance = np.zeros((n_pairs, n_frames))
        offset = np.zeros((n_pairs, n_frames))

        # 7. 逐帧计算
        for f, ts in enumerate(trajectory):
            positions = ts.positions
            metrics = self._compute_metrics(positions)
            distance[:, f] = metrics["distance"]
            offset[:, f] = metrics["offset"]
            existence[:, f] = self.apply_threshold(metrics)

        # 8. 过滤从未存在的 pair
        has_any = np.any(existence, axis=1)
        if not np.any(has_any):
            return []

        tuples = [(rings[pair_ring[i]], cations[pair_cation[i]])
                  for i in range(n_pairs) if has_any[i]]

        return [Interaction(
            interaction_type=self.name,
            groups=tuples,
            existence=existence[has_any],
            metrics={
                "distance": distance[has_any],
                "offset": offset[has_any],
            }
        )]

    # ==================== 数据结构构建 ====================

    @staticmethod
    def _build_circular_padding(rings: List[Group]):
        """构建环形 padding + 显式邻居索引。同 PiStacking。"""
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
            for j in range(max_atoms):
                indices[i, j] = atom_indices[j % n]
            for j in range(max_atoms):
                prev_idx[i, j] = atom_indices[(j - 1) % n]
                next_idx[i, j] = atom_indices[(j + 1) % n]
            valid[i, :n] = True

        return indices, prev_idx, next_idx, valid, n_atoms

    @staticmethod
    def _build_padding(groups: List[Group]):
        """构建 padding 索引矩阵。同 SaltBridge。"""
        n = len(groups)
        max_atoms = max(len(g.atoms) for g in groups)
        indices = np.zeros((n, max_atoms), dtype=int)
        charges = np.zeros((n, max_atoms))
        valid = np.zeros((n, max_atoms), dtype=bool)
        for i, g in enumerate(groups):
            na = len(g.atoms)
            indices[i, :na] = g.atom_indices
            charges[i, :na] = [a.atom_charge for a in g.atoms]
            valid[i, :na] = True
        return indices, charges, valid

    # ==================== 向量化计算 ====================

    def _ring_centers(self, positions):
        """计算全部环的几何中心。"""
        coords = positions[self._ring_indices]
        masked = coords * self._ring_valid[:, :, None]
        return np.sum(masked, axis=1) / np.sum(self._ring_valid, axis=1)[:, None]

    def _ring_normals(self, positions):
        """计算全部环的法向量（显式邻居索引）。"""
        coords = positions[self._ring_indices]
        prev_coords = positions[self._ring_prev]
        next_coords = positions[self._ring_next]
        atom_normals = np.cross(coords - prev_coords, coords - next_coords)
        masked = atom_normals * self._ring_valid[:, :, None]
        normal = np.sum(masked, axis=1) / np.sum(self._ring_valid, axis=1)[:, None]
        norm = np.linalg.norm(normal, axis=1, keepdims=True)
        return normal / norm

    def _cation_centers(self, positions):
        """计算全部阳离子的电荷中心。"""
        coords = positions[self._cation_indices]
        weighted = coords * self._cation_charges[:, :, None] * self._cation_valid[:, :, None]
        q_sum = np.sum(self._cation_charges * self._cation_valid, axis=1)
        return np.sum(weighted, axis=1) / q_sum[:, None]

    def _compute_metrics(self, positions):
        """向量化计算全部候选对的指标。"""
        ri = self._pair_ring
        ci = self._pair_cation

        r_centers = self._ring_centers(positions)
        r_normals = self._ring_normals(positions)
        c_centers = self._cation_centers(positions)

        c1 = r_centers[ri]
        c2 = c_centers[ci]
        n1 = r_normals[ri]

        distance = np.linalg.norm(c1 - c2, axis=1)
        offset = self._projection_distance(n1, c1, c2)

        return {"distance": distance, "offset": offset}

    # ==================== 辅助方法 ====================

    @staticmethod
    def _projection_distance(normal, plane_point, target_point):
        """投影距离。同 PiStacking。"""
        d1 = np.linalg.norm(target_point - (plane_point + normal), axis=1)
        d2 = np.linalg.norm(target_point - (plane_point - normal), axis=1)
        sign = np.where(d1 < d2, 1.0, -1.0)[:, None]
        oriented_normal = normal * sign
        t = target_point - plane_point
        proj_dist = np.sum(t * oriented_normal, axis=1)
        proj_point = target_point - proj_dist[:, None] * oriented_normal
        return np.linalg.norm(proj_point - plane_point, axis=1)
