# -*- coding: utf-8 -*-
"""疏水相互作用检测器（两轮遍历 + 稀疏存储）。

判据：两个疏水碳原子之间距离在 0.5~4.0 Å。
参考：PLIP。

Pass1：KDTree 自查询发现近邻对 + 距离过滤
Pass2：从 InteractionSparse 重建索引 + 向量化计算
"""

import numpy as np
from typing import List, Tuple, Dict

from ..core.interfaces import InteractionDetectorTwoPass
from ..core.datas import Group, Interaction, InteractionSparse


# PLIP 阈值
HYDROPH_DIST_MAX = 4.0  # Å，疏水原子间最大距离
HYDROPH_MIN_DIST = 0.5  # Å，最小距离


class HydrophobicDetectorTwoPass(InteractionDetectorTwoPass):
    """疏水相互作用检测器（两轮遍历 + 稀疏存储）。"""

    @property
    def name(self) -> str:
        return "hydrophobic"

    @property
    def required_group_types(self) -> List[str]:
        return ["hydrophobic"]

    @property
    def metric_names(self) -> List[str]:
        return ["distance"]

    # ==================== Pass1：KDTree 发现 ====================

    def run_pass1(self, groups: List[Group], trajectory,
                  tuple_filter=None) -> InteractionSparse:
        """执行 Pass1：KDTree 自查询发现近邻对。"""
        hydro_groups = [g for g in groups if g.group_type == "hydrophobic"]

        if not hydro_groups:
            return InteractionSparse(interaction_type=self.name, data={})

        # 构建原子索引数组
        self._atom_indices = np.array(
            [g.atoms[0].atom_global_idx for g in hydro_groups])

        sparse_data: Dict[Tuple[int, ...], dict] = {}

        for f, ts in enumerate(trajectory):
            pairs, distances = self._kdtree_discover(ts.positions)

            for i, j, dist in zip(pairs[0], pairs[1], distances):
                g1 = hydro_groups[i]
                g2 = hydro_groups[j]

                if tuple_filter is not None and not tuple_filter((g1, g2)):
                    continue

                group_ids = (g1.group_id, g2.group_id)

                if group_ids not in sparse_data:
                    sparse_data[group_ids] = {
                        "groups": (g1, g2),
                        "frames": [],
                        "metrics": {"distance": []}
                    }

                sparse_data[group_ids]["frames"].append(f)
                sparse_data[group_ids]["metrics"]["distance"].append(dist)

        # TODO: 后处理去重 — 同一 (group_id, residue_id) 对只保留平均距离最近的 pair
        # 物理原因：残基 LEU100 有 3 个疏水碳 C1/C2/C3，配体有 CA，
        # CA-C1、CA-C2、CA-C3 都是 "CA-LEU100" 的疏水接触，应只保留最近的。
        # 实现：计算每对在 active 帧的平均距离，每 (group_id, residue_id) 只保留最优。

        return InteractionSparse(interaction_type=self.name, data=sparse_data)

    # ==================== Pass2：从 InteractionSparse 重建 ====================

    def run_pass2(self, sparse: InteractionSparse,
                  trajectory) -> List[Interaction]:
        """执行 Pass2：从 InteractionSparse 重建索引。"""
        if not sparse.data:
            return []

        self._build_indices_from_sparse(sparse)
        return super().run_pass2(sparse, trajectory)

    # ==================== compute_pair_metrics ====================

    def compute_pair_metrics(self, group_tuples, all_positions):
        """向量化计算距离。"""
        idx1 = np.array([gt[0].atoms[0].atom_global_idx for gt in group_tuples])
        idx2 = np.array([gt[1].atoms[0].atom_global_idx for gt in group_tuples])
        pos1 = all_positions[idx1]
        pos2 = all_positions[idx2]
        return {"distance": np.linalg.norm(pos1 - pos2, axis=1)}

    def apply_threshold(self, metrics):
        """0.5 < 距离 < 4.0。"""
        return (metrics["distance"] > HYDROPH_MIN_DIST) & \
               (metrics["distance"] < HYDROPH_DIST_MAX)

    # ==================== 内部辅助方法 ====================

    def _kdtree_discover(self, all_positions):
        """KDTree 自查询发现近邻对。返回 (pairs, distances)。"""
        from scipy.spatial import cKDTree

        coords = all_positions[self._atom_indices]
        tree = cKDTree(coords)
        sparse_dist = tree.sparse_distance_matrix(tree, HYDROPH_DIST_MAX)

        coo = sparse_dist.tocoo()
        # 去重 (i<j) + 排除自身 (i==j)
        mask = coo.row < coo.col
        pairs = (coo.row[mask], coo.col[mask])
        distances = coo.data[mask]

        # 最小距离过滤
        dist_mask = distances > HYDROPH_MIN_DIST
        return (pairs[0][dist_mask], pairs[1][dist_mask]), distances[dist_mask]

    def _build_indices_from_sparse(self, sparse: InteractionSparse):
        """从 InteractionSparse 重建原子索引数组（Pass2 用）。"""
        group_tuples = [entry["groups"] for entry in sparse.data.values()]
        self._pair_idx1 = np.array(
            [gt[0].atoms[0].atom_global_idx for gt in group_tuples])
        self._pair_idx2 = np.array(
            [gt[1].atoms[0].atom_global_idx for gt in group_tuples])
