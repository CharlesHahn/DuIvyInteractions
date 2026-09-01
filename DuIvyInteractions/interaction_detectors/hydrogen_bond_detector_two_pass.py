# -*- coding: utf-8 -*-
"""氢键检测器（两轮遍历 + 稀疏存储）。

判据：D-A 距离 ≤ 4.1 Å 且 D-H···A 角度 ≥ 100°。
参考：PLIP, Hubbard & Haider 2001（距离 +0.6 Å，角度 +10°）。
"""

import numpy as np
from typing import List, Tuple, Dict
from scipy.spatial import KDTree

from ..core.interfaces import InteractionDetectorTwoPass
from ..core.datas import Group, Interaction


# PLIP 阈值
HBOND_DIST_MAX = 4.1         # Å，D-A 最大距离
HBOND_DON_ANGLE_MIN = 100.0  # °，D-H···A 最小角度


class HydrogenBondDetectorTwoPass(InteractionDetectorTwoPass):
    """氢键检测器（两轮遍历 + 稀疏存储）。"""

    @property
    def name(self) -> str:
        return "hydrogen_bond"

    @property
    def required_group_types(self) -> List[str]:
        return ["H_donor", "H_acceptor"]

    @property
    def metric_names(self) -> List[str]:
        return ["distance", "angle"]

    # ==================== 子类必须实现 ====================

    def initialize_candidates(self, groups, trajectory, tuple_filter=None):
        """构建索引数组，生成候选 items。"""
        donors = [g for g in groups if g.group_type == "H_donor"]
        acceptors = [g for g in groups if g.group_type == "H_acceptor"]

        if not donors or not acceptors:
            return []

        # 构建索引数组
        self._donor_d_idx = np.array([g.atoms[0].atom_global_idx for g in donors])
        self._donor_h_idx = np.array([g.atoms[1].atom_global_idx for g in donors])
        self._acceptor_a_idx = np.array([g.atoms[0].atom_global_idx for g in acceptors])

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

        self._pair_donor = pair_donor
        self._pair_acceptor = pair_acceptor

        items = [(donors[pair_donor[i]], acceptors[pair_acceptor[i]])
                 for i in range(len(pair_donor))]
        return items

    def discover_active_pairs(self, pair_indices, all_positions, frame):
        """Pass1：KDTree 距离判定 → 角度判定 → active pairs。

        第一层：KDTree(D, radius=4.1Å) → 满足距离准则的 pair
        第二层：对这些 pair 算 D-H···A 角度 ≥ 100°
        """
        d_pos = all_positions[self._donor_d_idx]
        h_pos = all_positions[self._donor_h_idx]
        a_pos = all_positions[self._acceptor_a_idx]

        # 建 donor → pair index 的映射
        # pair_donor[i] = di 表示第 i 个 pair 的 donor 是 donors[di]
        pair_donor = self._pair_donor
        pair_acceptor = self._pair_acceptor

        # 第一层：KDTree 距离判定
        tree = KDTree(d_pos)
        pairs = tree.query_ball_point(a_pos, HBOND_DIST_MAX)

        # 映射回全局 pair index
        donor_to_pair_indices = {}
        for i, di in enumerate(pair_donor):
            donor_to_pair_indices.setdefault(int(di), []).append(i)

        candidate_pair_idx = []
        for ai, nearby_d in enumerate(pairs):
            for di in nearby_d:
                for pi in donor_to_pair_indices.get(int(di), []):
                    if pair_acceptor[pi] == ai:
                        candidate_pair_idx.append(pi)

        if not candidate_pair_idx:
            return [], {}

        candidate_pair_idx = np.array(candidate_pair_idx)
        di = pair_donor[candidate_pair_idx]
        ai = pair_acceptor[candidate_pair_idx]

        # 第二层：角度判定
        d = d_pos[di]
        h = h_pos[di]
        a = a_pos[ai]

        dist = np.linalg.norm(d - a, axis=1)

        hd_vec = d - h
        ha_vec = a - h
        cos_angle = np.sum(hd_vec * ha_vec, axis=1) / (
            np.linalg.norm(hd_vec, axis=1) * np.linalg.norm(ha_vec, axis=1)
        )
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        ang = np.degrees(np.arccos(cos_angle))

        # 阈值过滤
        mask = (dist <= HBOND_DIST_MAX) & (ang >= HBOND_DON_ANGLE_MIN)

        active_indices = candidate_pair_idx[mask].tolist()
        active_metrics = {
            "distance": dist[mask].tolist(),
            "angle": ang[mask].tolist(),
        }

        return active_indices, active_metrics

    def compute_pair_metrics(self, pair_indices, all_positions):
        """Pass2：对 discovered pairs 算 distance + angle。"""
        d = all_positions[self._donor_d_idx[self._pair_donor[pair_indices]]]
        h = all_positions[self._donor_h_idx[self._pair_donor[pair_indices]]]
        a = all_positions[self._acceptor_a_idx[self._pair_acceptor[pair_indices]]]

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
