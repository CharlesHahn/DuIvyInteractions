# -*- coding: utf-8 -*-
"""水桥检测器（两轮遍历 + 稀疏存储）。

判据：水分子位于供体和受体之间，距离和角度满足 PLIP 定义。
参考：PLIP, Jiang et al. 2005。

Pass1：KDTree 发现三元组 + 内联计算角度
Pass2：从 InteractionSparse 重建索引 + 向量化计算
"""

import numpy as np
from typing import List, Tuple, Dict
from collections import defaultdict

from ..core.interfaces import InteractionDetectorTwoPass
from ..core.datas import Group, Interaction, InteractionSparse


# PLIP 阈值
WATER_BRIDGE_MAXDIST = 4.1   # Å，水 O 到极性原子最大距离
WATER_BRIDGE_OMEGA_MIN = 71  # °，受体-水O-供体H 最小角度
WATER_BRIDGE_OMEGA_MAX = 140 # °，受体-水O-供体H 最大角度
WATER_BRIDGE_THETA_MIN = 100 # °，水O-供体H-供体D 最小角度


class WaterBridgeDetectorTwoPass(InteractionDetectorTwoPass):
    """水桥检测器（两轮遍历 + 稀疏存储）。"""

    @property
    def name(self) -> str:
        return "water_bridge"

    @property
    def required_group_types(self) -> List[str]:
        return ["H_donor", "water", "H_acceptor"]

    @property
    def metric_names(self) -> List[str]:
        return ["dist_dw", "dist_wa", "theta", "omega"]

    # ==================== Pass1：KDTree 发现 ====================

    def run_pass1(self, groups: List[Group], trajectory,
                  tuple_filter=None) -> InteractionSparse:
        """执行 Pass1：KDTree 发现三元组 + 内联计算角度。"""
        from ..group_identifiers.amber_ff_identifier import WATER_RESIDUES

        donors = [g for g in groups
                  if g.group_type == "H_donor" and g.residue_name not in WATER_RESIDUES]
        waters = [g for g in groups if g.group_type == "water"]
        acceptors = [g for g in groups
                     if g.group_type == "H_acceptor" and g.residue_name not in WATER_RESIDUES]

        if not donors or not waters or not acceptors:
            return InteractionSparse(interaction_type=self.name, data={})

        # 构建索引
        self._donor_d_idx = np.array([g.atoms[0].atom_global_idx for g in donors])
        self._donor_h_idx = np.array([g.atoms[1].atom_global_idx for g in donors])
        self._water_ow_idx = np.array([g.atoms[0].atom_global_idx for g in waters])
        self._acceptor_a_idx = np.array([g.atoms[0].atom_global_idx for g in acceptors])

        sparse_data: Dict[Tuple[int, ...], dict] = {}

        for f, ts in enumerate(trajectory):
            positions = ts.positions
            triples = self._kdtree_discover(positions)

            for idx in range(len(triples[0])):
                di = triples[0][idx]
                wi = triples[1][idx]
                ai = triples[2][idx]
                donor = donors[di]
                water = waters[wi]
                acceptor = acceptors[ai]

                if tuple_filter is not None and not tuple_filter((donor, water, acceptor)):
                    continue

                group_ids = (donor.group_id, water.group_id, acceptor.group_id)

                if group_ids not in sparse_data:
                    sparse_data[group_ids] = {
                        "groups": (donor, water, acceptor),
                        "frames": [],
                        "metrics": {name: [] for name in self.metric_names}
                    }

                sparse_data[group_ids]["frames"].append(f)
                sparse_data[group_ids]["metrics"]["dist_dw"].append(float(triples[3][idx]))
                sparse_data[group_ids]["metrics"]["dist_wa"].append(float(triples[4][idx]))
                sparse_data[group_ids]["metrics"]["theta"].append(float(triples[5][idx]))
                sparse_data[group_ids]["metrics"]["omega"].append(float(triples[6][idx]))

        return InteractionSparse(interaction_type=self.name, data=sparse_data)

    # ==================== Pass2：从 InteractionSparse 重建 ====================

    def run_pass2(self, sparse: InteractionSparse,
                  trajectory) -> List[Interaction]:
        """执行 Pass2：从 InteractionSparse 重建缓存索引。"""
        if not sparse.data:
            return []

        self._build_indices_from_sparse(sparse)
        return super().run_pass2(sparse, trajectory)

    # ==================== compute_pair_metrics ====================

    def compute_pair_metrics(self, group_tuples, all_positions):
        """向量化计算 4 个指标（使用缓存索引）。"""
        if not hasattr(self, '_cached_d_idx'):
            raise RuntimeError(
                "compute_pair_metrics requires cached indices. "
                "Call _build_indices_from_sparse first via run_pass2.")

        d = all_positions[self._cached_d_idx]
        h = all_positions[self._cached_h_idx]
        ow = all_positions[self._cached_ow_idx]
        a = all_positions[self._cached_a_idx]

        # 距离
        dist_dw = np.linalg.norm(d - ow, axis=1)
        dist_wa = np.linalg.norm(ow - a, axis=1)

        # theta: 水O-供体H-供体D
        theta = self._vec_angle(d - h, ow - h)

        # omega: 受体-水O-供体H
        omega = self._vec_angle(a - ow, h - ow)

        return {"dist_dw": dist_dw, "dist_wa": dist_wa,
                "theta": theta, "omega": omega}

    def apply_threshold(self, metrics):
        """距离 < 4.1Å & theta ≥ 100° & 71° ≤ omega ≤ 140°。"""
        dist_ok = (metrics["dist_dw"] < WATER_BRIDGE_MAXDIST) & \
                  (metrics["dist_wa"] < WATER_BRIDGE_MAXDIST)
        theta_ok = metrics["theta"] >= WATER_BRIDGE_THETA_MIN
        omega_ok = (metrics["omega"] >= WATER_BRIDGE_OMEGA_MIN) & \
                   (metrics["omega"] <= WATER_BRIDGE_OMEGA_MAX)
        return dist_ok & theta_ok & omega_ok

    # ==================== 内部辅助方法 ====================

    def _kdtree_discover(self, positions):
        """KDTree 发现三元组 + 向量化计算角度。

        Returns:
            (di_arr, wi_arr, ai_arr, dist_dw_arr, dist_wa_arr, theta_arr, omega_arr)
            每个都是 numpy 数组，长度 = 满足阈值的三元组数量。
        """
        from scipy.spatial import cKDTree

        d_coords = positions[self._donor_d_idx]
        w_coords = positions[self._water_ow_idx]
        a_coords = positions[self._acceptor_a_idx]
        h_coords = positions[self._donor_h_idx]

        tree_d = cKDTree(d_coords)
        tree_a = cKDTree(a_coords)
        tree_w = cKDTree(w_coords)

        # KDTree 获取距离对
        dw_sparse = tree_w.sparse_distance_matrix(tree_d, WATER_BRIDGE_MAXDIST)
        wa_sparse = tree_w.sparse_distance_matrix(tree_a, WATER_BRIDGE_MAXDIST)

        if dw_sparse.nnz == 0 or wa_sparse.nnz == 0:
            return tuple(np.array([]) for _ in range(7))

        dw_coo = dw_sparse.tocoo()
        wa_coo = wa_sparse.tocoo()

        # 按 wi 分组
        dw_by_w = defaultdict(list)
        for idx in range(len(dw_coo.row)):
            dw_by_w[dw_coo.row[idx]].append(
                (dw_coo.col[idx], dw_coo.data[idx]))

        wa_by_w = defaultdict(list)
        for idx in range(len(wa_coo.row)):
            wa_by_w[wa_coo.row[idx]].append(
                (wa_coo.col[idx], wa_coo.data[idx]))

        # 向量化组合三元组（np.repeat + np.tile）
        common_w = set(dw_by_w) & set(wa_by_w)
        if not common_w:
            return tuple(np.array([]) for _ in range(7))

        all_di, all_wi, all_ai = [], [], []
        all_dist_dw, all_dist_wa = [], []

        for wi in common_w:
            dw_list = dw_by_w[wi]
            wa_list = wa_by_w[wi]
            n_dw = len(dw_list)
            n_wa = len(wa_list)
            n_triples = n_dw * n_wa

            dw_di = np.array([x[0] for x in dw_list])
            dw_dist = np.array([x[1] for x in dw_list])
            wa_ai = np.array([x[0] for x in wa_list])
            wa_dist = np.array([x[1] for x in wa_list])

            all_di.append(np.repeat(dw_di, n_wa))
            all_wi.append(np.full(n_triples, wi))
            all_ai.append(np.tile(wa_ai, n_dw))
            all_dist_dw.append(np.repeat(dw_dist, n_wa))
            all_dist_wa.append(np.tile(wa_dist, n_dw))

        di_arr = np.concatenate(all_di)
        wi_arr = np.concatenate(all_wi)
        ai_arr = np.concatenate(all_ai)
        dist_dw_arr = np.concatenate(all_dist_dw)
        dist_wa_arr = np.concatenate(all_dist_wa)

        # 向量化取坐标
        d = d_coords[di_arr]
        h = h_coords[di_arr]
        ow = w_coords[wi_arr]
        a = a_coords[ai_arr]

        # 向量化计算角度
        theta = self._vec_angle(d - h, ow - h)
        omega = self._vec_angle(a - ow, h - ow)

        # 向量化阈值筛选
        mask = (theta >= WATER_BRIDGE_THETA_MIN) & \
               (omega >= WATER_BRIDGE_OMEGA_MIN) & \
               (omega <= WATER_BRIDGE_OMEGA_MAX)

        return (di_arr[mask], wi_arr[mask], ai_arr[mask],
                dist_dw_arr[mask], dist_wa_arr[mask],
                theta[mask], omega[mask])

    def _build_indices_from_sparse(self, sparse: InteractionSparse):
        """从 InteractionSparse 重建缓存索引（Pass2 用）。"""
        group_tuples = [entry["groups"] for entry in sparse.data.values()]

        self._cached_d_idx = np.array(
            [gt[0].atoms[0].atom_global_idx for gt in group_tuples])
        self._cached_h_idx = np.array(
            [gt[0].atoms[1].atom_global_idx for gt in group_tuples])
        self._cached_ow_idx = np.array(
            [gt[1].atoms[0].atom_global_idx for gt in group_tuples])
        self._cached_a_idx = np.array(
            [gt[2].atoms[0].atom_global_idx for gt in group_tuples])

    @staticmethod
    def _vec_angle(v1, v2):
        """计算两组向量的夹角（度）。v1, v2: (n, 3) → (n,)"""
        cos_angle = np.sum(v1 * v2, axis=1) / (
            np.linalg.norm(v1, axis=1) * np.linalg.norm(v2, axis=1)
        )
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        return np.degrees(np.arccos(cos_angle))
