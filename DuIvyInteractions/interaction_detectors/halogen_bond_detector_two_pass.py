# -*- coding: utf-8 -*-
"""卤键检测器（两轮遍历 + 稀疏存储）。

判据：X-A ≤ 4.0 Å，C-X···A 在 165°±30°，X···A-R 在 120°±30°。
参考：PLIP, Auffinger et al.
"""

import numpy as np
from typing import List, Tuple, Dict

from ..core.interfaces import InteractionDetectorTwoPass
from ..core.datas import Group, Interaction, InteractionSparse


# PLIP 阈值
HALOGEN_DIST_MAX = 4.0        # Å
HALOGEN_DON_ANGLE = 165.0     # °，C-X···A 最优角度
HALOGEN_ACC_ANGLE = 120.0     # °，X···A-R 最优角度
HALOGEN_ANGLE_DEV = 30.0      # °，角度偏差上限


class HalogenBondDetectorTwoPass(InteractionDetectorTwoPass):
    """卤键检测器（两轮遍历 + 稀疏存储）。"""

    @property
    def name(self) -> str:
        return "halogen_bond"

    @property
    def required_group_types(self) -> List[str]:
        return ["halogen_donor", "halogen_acceptor"]

    @property
    def metric_names(self) -> List[str]:
        return ["distance", "don_angle", "acc_angle"]

    # ==================== initialize_candidates ====================

    def initialize_candidates(self, groups, trajectory, tuple_filter=None):
        """分组 → padding → 笛卡尔积 → tuple_filter → 设置缓存索引。"""
        donors = [g for g in groups if g.group_type == "halogen_donor"]
        acceptors = [g for g in groups if g.group_type == "halogen_acceptor"]

        if not donors or not acceptors:
            return []

        # 构建索引
        donor_c_idx = np.array([g.atoms[0].atom_global_idx for g in donors])
        donor_x_idx = np.array([g.atoms[1].atom_global_idx for g in donors])
        acceptor_a_idx = np.array([g.atoms[0].atom_global_idx for g in acceptors])
        acceptor_r_padding = self._build_acceptor_padding(acceptors)

        # 笛卡尔积
        n_d = len(donors)
        n_a = len(acceptors)
        d_grid, a_grid = np.meshgrid(np.arange(n_d), np.arange(n_a), indexing='ij')
        pair_donor = d_grid.ravel()
        pair_acceptor = a_grid.ravel()

        # tuple_filter
        if tuple_filter is not None:
            mask = np.array([
                tuple_filter((donors[di], acceptors[ai]))
                for di, ai in zip(pair_donor, pair_acceptor)
            ])
            pair_donor = pair_donor[mask]
            pair_acceptor = pair_acceptor[mask]

        if len(pair_donor) == 0:
            return []

        # 设置缓存索引（compute_pair_metrics 使用）
        self._cached_c_idx = donor_c_idx[pair_donor]
        self._cached_x_idx = donor_x_idx[pair_donor]
        self._cached_a_idx = acceptor_a_idx[pair_acceptor]
        self._cached_r_padding = acceptor_r_padding[pair_acceptor]

        # 返回 group_tuples 列表
        return [(donors[pair_donor[i]], acceptors[pair_acceptor[i]])
                for i in range(len(pair_donor))]

    # ==================== compute_pair_metrics ====================

    def compute_pair_metrics(self, group_tuples, all_positions):
        """向量化计算距离和角度（使用缓存索引）。"""
        if not hasattr(self, '_cached_c_idx'):
            raise RuntimeError(
                "compute_pair_metrics requires cached indices. "
                "Call initialize_candidates or _build_indices_from_sparse first.")

        c = all_positions[self._cached_c_idx]
        x = all_positions[self._cached_x_idx]
        a = all_positions[self._cached_a_idx]
        r = all_positions[self._cached_r_padding]  # (n_pairs, max_r, 3)

        # X-A 距离
        distance = np.linalg.norm(x - a, axis=1)

        # C-X···A 角度
        don_angle = self._vec_angle(c - x, a - x)

        # X···A-R 角度（跳过 padding 第 0 列，那是 A 本身）
        ax = x - a
        r_only = r[:, 1:, :]
        ar = r_only - a[:, None, :]
        cos = np.sum(ax[:, None, :] * ar, axis=2)
        norm_ax = np.linalg.norm(ax, axis=1, keepdims=True)
        norm_ar = np.linalg.norm(ar, axis=2)
        cos = np.clip(cos / (norm_ax * norm_ar), -1.0, 1.0)
        all_acc_angles = np.degrees(np.arccos(cos))
        best_r = np.argmin(np.abs(all_acc_angles - HALOGEN_ACC_ANGLE), axis=1)
        acc_angle = all_acc_angles[np.arange(len(best_r)), best_r]

        return {"distance": distance, "don_angle": don_angle, "acc_angle": acc_angle}

    def apply_threshold(self, metrics):
        """距离 ≤ 4.0Å，don_angle ∈ [135°,195°]，acc_angle ∈ [90°,150°]。"""
        dist_ok = metrics["distance"] <= HALOGEN_DIST_MAX
        don_ok = (metrics["don_angle"] >= HALOGEN_DON_ANGLE - HALOGEN_ANGLE_DEV) & \
                 (metrics["don_angle"] <= HALOGEN_DON_ANGLE + HALOGEN_ANGLE_DEV)
        acc_ok = (metrics["acc_angle"] >= HALOGEN_ACC_ANGLE - HALOGEN_ANGLE_DEV) & \
                 (metrics["acc_angle"] <= HALOGEN_ACC_ANGLE + HALOGEN_ANGLE_DEV)
        return dist_ok & don_ok & acc_ok

    # ==================== Pass2 ====================

    def run_pass2(self, sparse: InteractionSparse,
                  trajectory) -> List[Interaction]:
        """执行 Pass2：从 InteractionSparse 重建缓存索引。"""
        if not sparse.data:
            return []

        self._build_indices_from_sparse(sparse)
        return super().run_pass2(sparse, trajectory)

    # ==================== 内部辅助方法 ====================

    def _build_indices_from_sparse(self, sparse: InteractionSparse):
        """从 InteractionSparse 重建缓存索引（Pass2 用）。"""
        group_tuples = [entry["groups"] for entry in sparse.data.values()]
        donors = [gt[0] for gt in group_tuples]
        acceptors = [gt[1] for gt in group_tuples]

        self._cached_c_idx = np.array([g.atoms[0].atom_global_idx for g in donors])
        self._cached_x_idx = np.array([g.atoms[1].atom_global_idx for g in donors])
        self._cached_a_idx = np.array([g.atoms[0].atom_global_idx for g in acceptors])
        self._cached_r_padding = self._build_acceptor_padding(acceptors)

    @staticmethod
    def _build_acceptor_padding(acceptors: List[Group]):
        """构建 acceptor R 部分的环形 padding 索引矩阵。"""
        n = len(acceptors)
        max_atoms = max(len(g.atoms) for g in acceptors)
        indices = np.zeros((n, max_atoms), dtype=int)

        for i, g in enumerate(acceptors):
            atom_indices = g.atom_indices
            n_atoms = len(atom_indices)
            indices[i, 0] = atom_indices[0]
            n_r = n_atoms - 1
            if n_r > 0:
                r_indices = atom_indices[1:]
                for j in range(1, max_atoms):
                    indices[i, j] = r_indices[(j - 1) % n_r]

        return indices

    @staticmethod
    def _vec_angle(v1, v2):
        """计算两组向量的夹角（度）。"""
        cos_angle = np.sum(v1 * v2, axis=1) / (
            np.linalg.norm(v1, axis=1) * np.linalg.norm(v2, axis=1)
        )
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        return np.degrees(np.arccos(cos_angle))
