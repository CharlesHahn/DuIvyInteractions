# -*- coding: utf-8 -*-
"""盐桥检测器（两轮遍历 + 稀疏存储）。

判据：两个相反电荷中心距离 ≤ SALTBRIDGE_DIST_MAX。
参考：PLIP, Barlow and Thornton 1983（+1.5 Å 扩展）。

与 PerFrame 版本的区别：
- 无第一帧 prefilter，每帧计算全部 pair → 零遗漏
- 稀疏存储（只记录 active 的 pair），内存低
"""

import numpy as np
from typing import List, Tuple, Dict

from ..core.interfaces import InteractionDetectorTwoPass
from ..core.datas import Group, Interaction


# PLIP 阈值：5.5 Å
SALTBRIDGE_DIST_MAX = 5.5  # Å


class SaltBridgeDetectorTwoPass(InteractionDetectorTwoPass):
    """盐桥检测器（两轮遍历 + 稀疏存储）。"""

    def __init__(self):
        self._pos_groups = None
        self._neg_groups = None
        self._pos_idx = None
        self._neg_idx = None
        self._pos_q = None
        self._neg_q = None
        self._pos_valid = None
        self._neg_valid = None
        self._pos_flat = None
        self._neg_flat = None

    @property
    def name(self) -> str:
        return "salt_bridge"

    @property
    def required_group_types(self) -> List[str]:
        return ["charged_positive", "charged_negative"]

    @property
    def metric_names(self) -> List[str]:
        return ["distance"]

    # ==================== detect_frame ====================

    def detect_frame(self, items: list, all_positions: np.ndarray,
                     frame: int) -> Tuple[List[int], Dict[str, List[float]]]:
        """计算全部 pair 的电荷中心距离，返回 active 的 pair。

        Args:
            items: 未使用（groups 在 __init__ 已分组）
            all_positions: (n_atoms_total, 3)
            frame: 帧号

        Returns:
            (active_indices, {"distance": [...]})
        """
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

    # ==================== 重写 detect ====================

    def detect(self, groups, trajectory=None, n_workers=1,
               topology_path=None, trajectory_path=None,
               tuple_filter=None) -> List[Interaction]:
        """盐桥检测主流程。"""
        if trajectory is None:
            raise ValueError("trajectory is required")

        # 1. 分组
        self._pos_groups = [g for g in groups if g.group_type == "charged_positive"]
        self._neg_groups = [g for g in groups if g.group_type == "charged_negative"]

        if not self._pos_groups or not self._neg_groups:
            return []

        # 2. 构建 padding
        self._pos_idx, self._pos_q, self._pos_valid = self._build_padding(
            self._pos_groups)
        self._neg_idx, self._neg_q, self._neg_valid = self._build_padding(
            self._neg_groups)

        # 3. 笛卡尔积
        n_pos = len(self._pos_groups)
        n_neg = len(self._neg_groups)
        pos_grid, neg_grid = np.meshgrid(
            np.arange(n_pos), np.arange(n_neg), indexing='ij')
        self._pos_flat = pos_grid.ravel()
        self._neg_flat = neg_grid.ravel()

        # 4. tuple_filter
        if tuple_filter is not None:
            pair_mask = np.array([
                tuple_filter((self._pos_groups[pi], self._neg_groups[ni]))
                for pi, ni in zip(self._pos_flat, self._neg_flat)
            ])
            self._pos_flat = self._pos_flat[pair_mask]
            self._neg_flat = self._neg_flat[pair_mask]

        if len(self._pos_flat) == 0:
            return []

        n_frames = trajectory.n_frames

        # 5. Pass1: 逐帧检测，稀疏存储
        sparse_results: Dict[int, dict] = {}

        for f, ts in enumerate(trajectory):
            active_indices, frame_metrics = self.detect_frame(
                None, ts.positions, f)

            for idx, pos in enumerate(active_indices):
                if pos not in sparse_results:
                    pi = self._pos_flat[pos]
                    ni = self._neg_flat[pos]
                    sparse_results[pos] = {
                        "groups": (self._pos_groups[pi], self._neg_groups[ni]),
                        "frames": [],
                        "metrics": {name: [] for name in self.metric_names},
                    }
                sparse_results[pos]["frames"].append(f)
                for name in self.metric_names:
                    sparse_results[pos]["metrics"][name].append(
                        frame_metrics[name][idx])

        if not sparse_results:
            return []

        # 6. Pass2: 稠密化
        existence, metrics = self.densify(sparse_results, n_frames)

        # 7. 构建 results
        indices = sorted(sparse_results.keys())
        results = [
            (sparse_results[i]["groups"], existence[row],
             {k: v[row] for k, v in metrics.items()})
            for row, i in enumerate(indices)
        ]

        results = self._post_process(results)
        return self._build_interaction(results)

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
