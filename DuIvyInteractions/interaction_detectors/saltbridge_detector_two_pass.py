# -*- coding: utf-8 -*-
"""盐桥检测器（两轮遍历 + 稀疏存储）。

判据：两个相反电荷中心距离 ≤ SALTBRIDGE_DIST_MAX。
参考：PLIP, Barlow and Thornton 1983（+1.5 Å 扩展）。

PBC：电荷中心聚合先以基团参考原子做 MIC 回绕（minimize_vectors 展平-回绕-
还原），中心距离 calc_bonds(box)，见 doc/pbc_design.md §6D。
"""

import numpy as np
from typing import List, Tuple, Dict

from ..core.interfaces import InteractionDetectorTwoPass
from ..core.datas import Group, Interaction, InteractionSparse
from ..core.pbc import get_box


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

    def compute_pair_metrics(self, group_tuples, all_positions, box=None):
        """向量化计算电荷中心距离（MIC）。

        box: 盒尺寸（MIC）；None = 欧氏距离（无 PBC）。
        """
        from MDAnalysis.lib.distances import calc_bonds
        pos_c = self._charge_centers(
            all_positions, self._pos_idx, self._pos_q, self._pos_valid, box=box)
        neg_c = self._charge_centers(
            all_positions, self._neg_idx, self._neg_q, self._neg_valid, box=box)
        return {"distance": calc_bonds(
            pos_c[self._pos_flat], neg_c[self._neg_flat], box=box)}

    def apply_threshold(self, metrics):
        """距离 ≤ SALTBRIDGE_DIST_MAX。"""
        return metrics["distance"] <= SALTBRIDGE_DIST_MAX

    def run_pass2(self, sparse: InteractionSparse,
                  trajectory) -> List[Interaction]:
        """执行 Pass2：补全全帧数据（逐帧 PBC 电荷中心距离）。

        自写循环（基类 run_pass2 不暴露帧/box）：逐帧取 get_box(ts) 传入
        compute_pair_metrics，使 PBC 距离逐帧正确（NPT 系综盒逐帧不同）。
        """
        if not sparse.data:
            return []

        # 从 InteractionSparse 重新构建 padding
        self._build_padding_from_sparse(sparse)
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
    def _charge_centers(positions, group_idx, group_q, group_valid, box=None):
        """padding 向量化计算全部 group 的电荷中心（PBC 感知）。

        box: 盒尺寸（MIC）；None = 欧氏（无 PBC）。
        有 box 时先以每基团第一个原子为参考做最小镜像回绕（保基团连续），
        再加权平均——防御性保证残基跨盒时中心正确。
        """
        coords = positions[group_idx]                    # (n, max_atoms, 3)
        if box is not None:
            from MDAnalysis.lib.distances import minimize_vectors
            n, ma = coords.shape[:2]
            ref = coords[:, 0:1, :]
            rel = coords - ref
            rel_mic = minimize_vectors(rel.reshape(-1, 3), box).reshape(n, ma, 3)
            coords = ref + rel_mic
        weighted = coords * group_q[:, :, None] * group_valid[:, :, None]
        q_sum = np.sum(group_q * group_valid, axis=1)
        centers = np.sum(weighted, axis=1) / q_sum[:, None]
        return centers
