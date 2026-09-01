# -*- coding: utf-8 -*-
"""氢键检测器（两轮遍历 + 稀疏存储）。

判据：D-A 距离 ≤ 4.1 Å 且 D-H···A 角度 ≥ 100°。
参考：PLIP, Hubbard & Haider 2001（距离 +0.6 Å，角度 +10°）。

Pass1：KDTree 筛选 D-A 距离 + 精确计算角度
Pass2：从 InteractionSparse 重新构建索引 + 向量化计算
"""

import numpy as np
from typing import List, Tuple, Dict

from ..core.interfaces import InteractionDetectorTwoPass
from ..core.datas import Group, Interaction, InteractionSparse


# PLIP 阈值
HBOND_DIST_MAX = 4.1         # Å，D-A 最大距离
HBOND_DON_ANGLE_MIN = 100.0  # °，D-H···A 最小角度


class HydrogenBondDetectorTwoPass(InteractionDetectorTwoPass):
    """氢键检测器（两轮遍历 + 稀疏存储）。

    Pass1：覆盖 run_pass1，KDTree 筛距离 + 精确计算角度。
    Pass2：从 InteractionSparse 重新构建索引 + 向量化计算。
    """

    @property
    def name(self) -> str:
        return "hydrogen_bond"

    @property
    def required_group_types(self) -> List[str]:
        return ["H_donor", "H_acceptor"]

    @property
    def metric_names(self) -> List[str]:
        return ["distance", "angle"]

    # ==================== Pass1：KDTree 发现 ====================

    def run_pass1(self, groups: List[Group], trajectory,
                  tuple_filter=None) -> InteractionSparse:
        """执行 Pass1：KDTree 筛距离 + 精确计算角度。"""
        # 1. 分组
        donors = [g for g in groups if g.group_type == "H_donor"]
        acceptors = [g for g in groups if g.group_type == "H_acceptor"]

        if not donors or not acceptors:
            return InteractionSparse(interaction_type=self.name, data={})

        # 2. 构建索引数组
        self._build_indices(donors, acceptors)

        # 3. 逐帧 KDTree 发现
        sparse_data: Dict[Tuple[int, ...], dict] = {}

        for f, ts in enumerate(trajectory):
            active_donors, active_acceptors, distances, angles = \
                self._kdtree_discover(ts.positions)

            for d_idx, a_idx, dist, angle in zip(
                    active_donors, active_acceptors, distances, angles):
                donor = donors[d_idx]
                acceptor = acceptors[a_idx]

                # tuple_filter
                if tuple_filter is not None and not tuple_filter((donor, acceptor)):
                    continue

                group_ids = (donor.group_id, acceptor.group_id)

                if group_ids not in sparse_data:
                    sparse_data[group_ids] = {
                        "groups": (donor, acceptor),
                        "frames": [],
                        "metrics": {"distance": [], "angle": []}
                    }

                sparse_data[group_ids]["frames"].append(f)
                sparse_data[group_ids]["metrics"]["distance"].append(dist)
                sparse_data[group_ids]["metrics"]["angle"].append(angle)

        return InteractionSparse(interaction_type=self.name, data=sparse_data)

    # ==================== Pass2：从 InteractionSparse 重建 ====================

    def run_pass2(self, sparse: InteractionSparse,
                  trajectory) -> List[Interaction]:
        """执行 Pass2：从 InteractionSparse 重新构建索引 + 向量化计算。"""
        if not sparse.data:
            return []

        # 从 InteractionSparse 重新构建索引数组
        self._build_indices_from_sparse(sparse)

        # 调用父类的 run_pass2
        return super().run_pass2(sparse, trajectory)

    # ==================== compute_pair_metrics ====================

    def compute_pair_metrics(self, group_tuples, all_positions):
        """向量化计算距离和角度。"""
        # 从 group_tuples 提取原子索引
        d_indices = np.array([gt[0].atoms[0].atom_global_idx for gt in group_tuples])
        h_indices = np.array([gt[0].atoms[1].atom_global_idx for gt in group_tuples])
        a_indices = np.array([gt[1].atoms[0].atom_global_idx for gt in group_tuples])

        d = all_positions[d_indices]
        h = all_positions[h_indices]
        a = all_positions[a_indices]

        # D-A 距离
        da_vec = d - a
        dist = np.linalg.norm(da_vec, axis=1)

        # D-H···A 角度
        hd_vec = d - h
        ha_vec = a - h
        cos_angle = np.sum(hd_vec * ha_vec, axis=1) / (
            np.linalg.norm(hd_vec, axis=1) * np.linalg.norm(ha_vec, axis=1)
        )
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        ang = np.degrees(np.arccos(cos_angle))

        return {"distance": dist, "angle": ang}

    def apply_threshold(self, metrics):
        """距离 ≤ 4.1Å 且角度 ≥ 100°。"""
        return (metrics["distance"] <= HBOND_DIST_MAX) & \
               (metrics["angle"] >= HBOND_DON_ANGLE_MIN)

    # ==================== 内部辅助方法 ====================

    def _build_indices(self, donors, acceptors):
        """构建索引数组（Pass1 用）。"""
        self._donor_d_idx = np.array([g.atoms[0].atom_global_idx for g in donors])
        self._donor_h_idx = np.array([g.atoms[1].atom_global_idx for g in donors])
        self._acceptor_a_idx = np.array([g.atoms[0].atom_global_idx for g in acceptors])

    def _build_indices_from_sparse(self, sparse: InteractionSparse):
        """从 InteractionSparse 重新构建索引数组（Pass2 用）。"""
        group_tuples = [entry["groups"] for entry in sparse.data.values()]
        donors = [gt[0] for gt in group_tuples]
        acceptors = [gt[1] for gt in group_tuples]

        self._donor_d_idx = np.array([g.atoms[0].atom_global_idx for g in donors])
        self._donor_h_idx = np.array([g.atoms[1].atom_global_idx for g in donors])
        self._acceptor_a_idx = np.array([g.atoms[0].atom_global_idx for g in acceptors])

    def _kdtree_discover(self, all_positions):
        """KDTree 发现 D-A 距离 ≤ 4.1Å 的 pair + 精确计算角度。"""
        from scipy.spatial import cKDTree

        # 获取 D 和 A 坐标
        d_coords = all_positions[self._donor_d_idx]
        a_coords = all_positions[self._acceptor_a_idx]

        # 构建两棵树
        tree_d = cKDTree(d_coords)
        tree_a = cKDTree(a_coords)

        # 一步获取 (d_idx, a_idx) 对 + 精确距离
        sparse_dist = tree_d.sparse_distance_matrix(tree_a, HBOND_DIST_MAX)

        if sparse_dist.nnz == 0:
            return np.array([]), np.array([]), np.array([]), np.array([])

        # 转换为 coo_matrix 并提取索引和距离
        coo_dist = sparse_dist.tocoo()
        active_donors = coo_dist.row
        active_acceptors = coo_dist.col
        distances = coo_dist.data

        # 计算角度（向量化）
        d = d_coords[active_donors]
        h = all_positions[self._donor_h_idx[active_donors]]
        a = a_coords[active_acceptors]

        hd_vec = d - h
        ha_vec = a - h
        cos_angle = np.sum(hd_vec * ha_vec, axis=1) / (
            np.linalg.norm(hd_vec, axis=1) * np.linalg.norm(ha_vec, axis=1))
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        angles = np.degrees(np.arccos(cos_angle))

        # 过滤角度 ≥ 100°
        angle_mask = angles >= HBOND_DON_ANGLE_MIN

        return (active_donors[angle_mask], active_acceptors[angle_mask],
                distances[angle_mask], angles[angle_mask])
