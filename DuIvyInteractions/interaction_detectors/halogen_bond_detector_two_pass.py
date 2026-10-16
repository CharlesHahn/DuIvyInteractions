# -*- coding: utf-8 -*-
"""卤键检测器（两轮遍历 + 稀疏存储）。

判据：X-A ≤ 4.0 Å，C-X···A 在 165°±30°，X···A-R 在 120°±30°。
参考：PLIP, Auffinger et al.

PBC：距离 calc_bonds(box)、角度 calc_angles(box)（acc_angle 经 R-padding 展开-
还原 + argmin 择优），见 doc/pbc_design.md §6C。
"""

import numpy as np
from typing import List, Tuple, Dict

from ..core.interfaces import InteractionDetectorTwoPass
from ..core.datas import Group, Interaction, InteractionSparse
from ..core.pbc import get_box


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

    def compute_pair_metrics(self, group_tuples, all_positions, box=None):
        """向量化计算距离和角度（使用缓存索引）。

        box: 盒尺寸（MIC）；None = 欧氏距离（无 PBC）。
        """
        if not hasattr(self, '_cached_c_idx'):
            raise RuntimeError(
                "compute_pair_metrics requires cached indices. "
                "Call initialize_candidates or _build_indices_from_sparse first.")

        from MDAnalysis.lib.distances import calc_bonds, calc_angles

        c = all_positions[self._cached_c_idx]
        x = all_positions[self._cached_x_idx]
        a = all_positions[self._cached_a_idx]
        r = all_positions[self._cached_r_padding]  # (n_pairs, max_r, 3)

        # X-A 距离（MIC）
        distance = calc_bonds(x, a, box=box)

        # C-X···A 角度（顶点 x；calc_angles 返回弧度）
        don_angle = np.degrees(calc_angles(c, x, a, box=box))

        # X···A-R 角度（跳过 padding 第 0 列，那是 A 本身）：展开 → 批量 MIC 角 → argmin
        r_only = r[:, 1:, :]                       # (n_pairs, max_r, 3)
        n_pairs = len(x)
        max_r = r_only.shape[1]
        x_rep = np.repeat(x, max_r, axis=0)        # (n_pairs*max_r, 3)
        a_rep = np.repeat(a, max_r, axis=0)
        r_flat = r_only.reshape(-1, 3)
        all_acc = np.degrees(calc_angles(x_rep, a_rep, r_flat, box=box))  # 顶点 a
        all_acc_angles = all_acc.reshape(n_pairs, max_r)
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
        """执行 Pass2：补全全帧数据（逐帧 PBC 距离/角度）。

        自写循环（基类 run_pass2 不暴露帧/box）：逐帧取 get_box(ts) 传入
        compute_pair_metrics，使 PBC 距离/角度逐帧正确（NPT 系综盒逐帧不同）。
        """
        if not sparse.data:
            return []

        self._build_indices_from_sparse(sparse)
        group_tuples = [entry["groups"] for entry in sparse.data.values()]
        n_groups = len(group_tuples)
        n_frames = trajectory.n_frames

        # 预分配稠密数组
        existence = np.zeros((n_groups, n_frames), dtype=bool)
        metrics = {name: np.full((n_groups, n_frames), np.nan)
                   for name in self.metric_names}

        times = np.empty(n_frames)
        for f, ts in enumerate(trajectory):
            times[f] = ts.time
            frame_metrics = self.compute_pair_metrics(
                group_tuples, ts.positions, box=get_box(ts))
            existence[:, f] = self.apply_threshold(frame_metrics)
            for name in self.metric_names:
                metrics[name][:, f] = frame_metrics[name]

        # 构建 results 列表（与基类一致：不过滤从未存在的 pair）
        results = [
            (entry["groups"], existence[i],
             {k: v[i] for k, v in metrics.items()})
            for i, entry in enumerate(sparse.data.values())
        ]

        results = self._post_process(results)
        return self._build_interaction(results, times=times)

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
