# -*- coding: utf-8 -*-
"""盐桥检测器（两轮遍历 + 稀疏存储）。

判据：两个相反电荷中心距离 ≤ SALTBRIDGE_DIST_MAX。
参考：PLIP, Barlow and Thornton 1983（+1.5 Å 扩展）。
"""

import numpy as np
from typing import List, Tuple, Dict

from ..core.interfaces import InteractionDetectorTwoPass
from ..core.datas import Group, Interaction, InteractionSparse


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

    # ==================== 子类实现 ====================

    def initialize_candidates(self, groups, trajectory, tuple_filter=None):
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

        # 存入 self
        self._pos_flat = pos_flat
        self._neg_flat = neg_flat

        # 返回 group_tuples 列表
        group_tuples = [(pos_groups[pos_flat[i]], neg_groups[neg_flat[i]])
                        for i in range(len(pos_flat))]
        return group_tuples

    def compute_pair_metrics(self, group_tuples, all_positions):
        """向量化计算电荷中心距离。"""
        pos_c = self._charge_centers(
            all_positions, self._pos_idx, self._pos_q, self._pos_valid)
        neg_c = self._charge_centers(
            all_positions, self._neg_idx, self._neg_q, self._neg_valid)
        return {"distance": np.linalg.norm(
            pos_c[self._pos_flat] - neg_c[self._neg_flat], axis=1)}

    def apply_threshold(self, metrics):
        """距离 ≤ SALTBRIDGE_DIST_MAX。"""
        return metrics["distance"] <= SALTBRIDGE_DIST_MAX

    def run_pass2(self, sparse: InteractionSparse,
                  trajectory) -> List[Interaction]:
        """执行 Pass2：重新构建 padding，然后向量化计算。"""
        if not sparse.data:
            return []

        # 从 InteractionSparse 重新构建 padding
        self._build_padding_from_sparse(sparse)

        # 调用父类的 run_pass2
        return super().run_pass2(sparse, trajectory)

    # ==================== 内部辅助方法 ====================

    def _build_padding_from_sparse(self, sparse: InteractionSparse):
        """从 InteractionSparse 重新构建 padding。"""
        group_tuples = [entry["groups"] for entry in sparse.data.values()]
        pos_groups = [gt[0] for gt in group_tuples]
        neg_groups = [gt[1] for gt in group_tuples]

        # 重新构建 padding
        self._pos_idx, self._pos_q, self._pos_valid = self._build_padding(pos_groups)
        self._neg_idx, self._neg_q, self._neg_valid = self._build_padding(neg_groups)

        # 设置 _pos_flat 和 _neg_flat（顺序索引，长度 = discovered pairs）
        self._pos_flat = np.arange(len(pos_groups))
        self._neg_flat = np.arange(len(neg_groups))

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
