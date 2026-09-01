# -*- coding: utf-8 -*-
"""盐桥检测器（两轮遍历 + 稀疏存储）。

判据：两个相反电荷中心距离 ≤ SALTBRIDGE_DIST_MAX。
参考：PLIP, Barlow and Thornton 1983（+1.5 Å 扩展）。
"""

import numpy as np
from typing import List, Tuple, Dict

from ..core.interfaces import InteractionDetectorTwoPass
from ..core.datas import Group, Interaction


# PLIP 阈值：5.5 Å
SALTBRIDGE_DIST_MAX = 5.5  # Å


class SaltBridgeDetectorTwoPass(InteractionDetectorTwoPass):
    """盐桥检测器（两轮遍历 + 稀疏存储）。"""

    @property
    def name(self) -> str:
        return "salt_bridge"

    @property
    def required_group_types(self) -> List[str]:
        return ["charged_positive", "charged_negative"]

    @property
    def metric_names(self) -> List[str]:
        return ["distance"]

    # ==================== setup ====================

    def setup(self, groups, trajectory, tuple_filter=None):
        """分组 → padding → 笛卡尔积 → tuple_filter。"""
        pos_groups = [g for g in groups if g.group_type == "charged_positive"]
        neg_groups = [g for g in groups if g.group_type == "charged_negative"]

        if not pos_groups or not neg_groups:
            return []

        # 构建 padding
        self._pos_idx, self._pos_q, self._pos_valid = self._build_padding(pos_groups)
        self._neg_idx, self._neg_q, self._neg_valid = self._build_padding(neg_groups)

        # 笛卡尔积
        n_pos = len(pos_groups)
        n_neg = len(neg_groups)
        pos_grid, neg_grid = np.meshgrid(
            np.arange(n_pos), np.arange(n_neg), indexing='ij')
        pos_flat = pos_grid.ravel()
        neg_flat = neg_grid.ravel()

        # tuple_filter
        if tuple_filter is not None:
            mask = np.array([
                tuple_filter((pos_groups[pi], neg_groups[ni]))
                for pi, ni in zip(pos_flat, neg_flat)
            ])
            pos_flat = pos_flat[mask]
            neg_flat = neg_flat[mask]

        if len(pos_flat) == 0:
            return []

        # 存入 self（detect_frame 用）
        self._pos_flat = pos_flat
        self._neg_flat = neg_flat

        # 返回 items 列表（基类用于 sparse["groups"]）
        items = [(pos_groups[pos_flat[i]], neg_groups[neg_flat[i]])
                 for i in range(len(pos_flat))]
        return items

    # ==================== detect_frame ====================

    def detect_frame(self, all_positions, frame):
        """padding 向量化算全部 pair 距离，返回 active 的 pair。"""
        pos_c = self._charge_centers(
            all_positions, self._pos_idx, self._pos_q, self._pos_valid)
        neg_c = self._charge_centers(
            all_positions, self._neg_idx, self._neg_q, self._neg_valid)
        dist = np.linalg.norm(
            pos_c[self._pos_flat] - neg_c[self._neg_flat], axis=1)

        mask = dist <= SALTBRIDGE_DIST_MAX
        active_indices = np.where(mask)[0].tolist()
        active_distances = dist[mask].tolist()

        return active_indices, {"distance": active_distances}

    # ==================== 内部辅助方法 ====================

    @staticmethod
    def _build_padding(groups: List[Group]):
        """构建 padding 索引矩阵。"""
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
        """padding 向量化计算全部 group 的电荷中心。"""
        coords = positions[group_idx]
        weighted = coords * group_q[:, :, None] * group_valid[:, :, None]
        q_sum = np.sum(group_q * group_valid, axis=1)
        centers = np.sum(weighted, axis=1) / q_sum[:, None]
        return centers
