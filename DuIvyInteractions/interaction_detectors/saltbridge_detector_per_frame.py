# -*- coding: utf-8 -*-
"""盐桥检测器（逐帧策略，padding 向量化）。

判据：两个相反电荷中心距离 ≤ SALTBRIDGE_DIST_MAX。
参考：PLIP, Barlow and Thornton 1983（+1.5 Å 扩展）。
"""

import numpy as np
from typing import List, Tuple, Dict

from ..core.interfaces import InteractionDetectorPerFrame
from ..core.datas import Group, Interaction


# PLIP 阈值：5.5 Å
SALTBRIDGE_DIST_MAX = 5.5  # Å
PREFILTER_CUTOFF = SALTBRIDGE_DIST_MAX * 3  # 16.5 Å


class SaltBridgeDetectorPerFrame(InteractionDetectorPerFrame):
    """盐桥检测器（逐帧策略，padding 向量化）。"""

    def __init__(self):
        self._pos_indices = None
        self._neg_indices = None
        self._pos_idx = None
        self._neg_idx = None
        self._pos_q = None
        self._neg_q = None
        self._pos_valid = None
        self._neg_valid = None

    @property
    def name(self) -> str:
        return "salt_bridge"

    @property
    def required_group_types(self) -> List[str]:
        return ["charged_positive", "charged_negative"]

    @property
    def metric_names(self) -> List[str]:
        return ["distance"]

    # ==================== 抽象方法（基类要求，此处不使用） ====================

    def get_candidate_tuples(self, groups, coordinates=None):
        raise NotImplementedError("detect() is overridden, this method is not used")

    def compute_metrics_for_frame(self, tuples, all_positions, frame):
        raise NotImplementedError("detect() is overridden, this method is not used")

    def apply_threshold(self, metrics):
        return metrics["distance"] <= SALTBRIDGE_DIST_MAX

    # ==================== 重写 detect ====================

    def detect(self, groups, trajectory=None,
               n_workers: int = 1,
               topology_path: str = None,
               trajectory_path: str = None,
               tuple_filter=None) -> List[Interaction]:
        """盐桥检测主流程。

        流程：
        1. 分组：正电 / 负电
        2. 分别构建 padding 矩阵
        3. 笛卡尔积生成全部 tuple → tuple_filter 过滤
        4. 第一帧：padding 向量化算全部电荷中心 + 距离矩阵 → 距离过滤
        5. 后续帧：只对有效 tuple 算距离
        6. 构建 Interaction
        """
        if trajectory is None:
            raise ValueError("trajectory is required")
        if n_workers > 1:
            raise NotImplementedError("PerFrame 检测器暂不支持并行")

        # 1. 分组
        pos_groups = [g for g in groups if g.group_type == "charged_positive"]
        neg_groups = [g for g in groups if g.group_type == "charged_negative"]

        if not pos_groups or not neg_groups:
            return []

        # 2. 构建 padding
        pos_idx, pos_q, pos_valid = self._build_padding(pos_groups)
        neg_idx, neg_q, neg_valid = self._build_padding(neg_groups)

        # 3. 笛卡尔积 + tuple_filter
        n_pos = len(pos_groups)
        n_neg = len(neg_groups)
        pos_grid, neg_grid = np.meshgrid(np.arange(n_pos), np.arange(n_neg), indexing='ij')
        pos_flat = pos_grid.ravel()
        neg_flat = neg_grid.ravel()

        if tuple_filter is not None:
            pair_mask = np.array([
                tuple_filter((pos_groups[pi], neg_groups[ni]))
                for pi, ni in zip(pos_flat, neg_flat)
            ])
            pos_flat = pos_flat[pair_mask]
            neg_flat = neg_flat[pair_mask]

        if len(pos_flat) == 0:
            return []

        # 4. 第一帧：向量化过滤
        first_pos = trajectory[0].positions
        pos_centers = self._charge_centers(first_pos, pos_idx, pos_q, pos_valid)
        neg_centers = self._charge_centers(first_pos, neg_idx, neg_q, neg_valid)
        all_dist = np.linalg.norm(
            pos_centers[pos_flat] - neg_centers[neg_flat], axis=1)

        mask = all_dist < PREFILTER_CUTOFF
        if not np.any(mask):
            return []

        self._pos_indices = pos_flat[mask]
        self._neg_indices = neg_flat[mask]

        # 存 padding（逐帧复用）
        self._pos_idx = pos_idx
        self._neg_idx = neg_idx
        self._pos_q = pos_q
        self._neg_q = neg_q
        self._pos_valid = pos_valid
        self._neg_valid = neg_valid

        n_tuples = int(np.sum(mask))
        n_frames = trajectory.n_frames
        existence = np.zeros((n_tuples, n_frames), dtype=bool)
        distance = np.zeros((n_tuples, n_frames))

        # 5. 逐帧计算
        for f, ts in enumerate(trajectory):
            positions = ts.positions
            pos_c = self._charge_centers(positions, pos_idx, pos_q, pos_valid)
            neg_c = self._charge_centers(positions, neg_idx, neg_q, neg_valid)
            dist = np.linalg.norm(
                pos_c[self._pos_indices] - neg_c[self._neg_indices], axis=1)
            distance[:, f] = dist
            existence[:, f] = dist <= SALTBRIDGE_DIST_MAX

        # 6. 过滤从未存在的 tuple
        has_any = np.any(existence, axis=1)
        if not np.any(has_any):
            return []

        tuples = [(pos_groups[self._pos_indices[i]], neg_groups[self._neg_indices[i]])
                  for i in range(n_tuples) if has_any[i]]
        existence = existence[has_any]
        distance = distance[has_any]

        return [Interaction(
            interaction_type=self.name,
            groups=tuples,
            existence=existence,
            metrics={"distance": distance}
        )]

    # ==================== 内部辅助方法 ====================

    @staticmethod
    def _build_padding(groups: List[Group]):
        """构建 padding 索引矩阵。

        Returns:
            indices: (n_groups, max_atoms) int — 全局原子索引
            charges: (n_groups, max_atoms) float — 原子电荷
            valid: (n_groups, max_atoms) bool — 有效位
        """
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

    @staticmethod
    def _charge_centers(positions, group_idx, group_q, group_valid):
        """padding 向量化计算全部 group 的电荷中心。

        Args:
            positions: (n_atoms_total, 3)
            group_idx: (n_groups, max_atoms) int
            group_q: (n_groups, max_atoms) float
            group_valid: (n_groups, max_atoms) bool

        Returns:
            centers: (n_groups, 3)
        """
        coords = positions[group_idx]                              # (n_groups, max_atoms, 3)
        weighted = coords * group_q[:, :, None] * group_valid[:, :, None]
        q_sum = np.sum(group_q * group_valid, axis=1)              # (n_groups,)
        centers = np.sum(weighted, axis=1) / q_sum[:, None]        # (n_groups, 3)
        return centers
