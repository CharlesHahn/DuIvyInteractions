# -*- coding: utf-8 -*-
"""氢键检测器（逐帧策略，KDTree 预筛选 + 向量化计算）。

判据：D-A 距离 ≤ 4.1 Å 且 D-H···A 角度 ≥ 100°。
参考：PLIP, Hubbard & Haider 2001（距离 +0.6 Å，角度 +10°）。
"""

import numpy as np
from typing import List, Tuple, Dict
from scipy.spatial import KDTree

from ..core.interfaces import InteractionDetectorPerFrame
from ..core.datas import Group, Interaction


# PLIP 阈值
HBOND_DIST_MAX = 4.1         # Å，D-A 最大距离
HBOND_DON_ANGLE_MIN = 100.0  # °，D-H···A 最小角度
PREFILTER_CUTOFF = HBOND_DIST_MAX * 3  # 12.3 Å


class HydrogenBondDetectorPerFrame(InteractionDetectorPerFrame):
    """氢键检测器（逐帧策略，KDTree 预筛选 + 向量化计算）。"""

    def __init__(self):
        self._pair_donor_idx = None
        self._pair_acceptor_idx = None
        self._donor_d_idx = None
        self._donor_h_idx = None
        self._acceptor_a_idx = None

    @property
    def name(self) -> str:
        return "hydrogen_bond"

    @property
    def required_group_types(self) -> List[str]:
        return ["H_donor", "H_acceptor"]

    @property
    def metric_names(self) -> List[str]:
        return ["distance", "angle"]

    # ==================== 抽象方法（基类要求，此处不使用） ====================

    def get_candidate_tuples(self, groups, coordinates=None):
        raise NotImplementedError("detect() is overridden, this method is not used")

    def compute_metrics_for_frame(self, tuples, all_positions, frame):
        raise NotImplementedError("detect() is overridden, this method is not used")

    def apply_threshold(self, metrics):
        return (metrics["distance"] <= HBOND_DIST_MAX) & \
               (metrics["angle"] >= HBOND_DON_ANGLE_MIN)

    # ==================== 重写 detect ====================

    def detect(self, groups, trajectory=None,
               n_workers: int = 1,
               topology_path: str = None,
               trajectory_path: str = None,
               tuple_filter=None) -> List[Interaction]:
        """氢键检测主流程。"""
        if trajectory is None:
            raise ValueError("trajectory is required")
        if n_workers > 1:
            raise NotImplementedError("PerFrame 检测器暂不支持并行")

        # 1. 分组
        from ..group_identifiers.amber_ff_identifier import WATER_RESIDUES
        donors = [g for g in groups
                  if g.group_type == "H_donor" and g.residue_name not in WATER_RESIDUES]
        acceptors = [g for g in groups
                     if g.group_type == "H_acceptor" and g.residue_name not in WATER_RESIDUES]

        if not donors or not acceptors:
            return []

        # 2. 构建索引数组
        self._donor_d_idx = np.array([g.atoms[0].atom_global_idx for g in donors])
        self._donor_h_idx = np.array([g.atoms[1].atom_global_idx for g in donors])
        self._acceptor_a_idx = np.array([g.atoms[0].atom_global_idx for g in acceptors])

        # 3. 第一帧 KDTree 预筛选
        first_pos = trajectory[0].positions
        d_pos = first_pos[self._donor_d_idx]     # (n_donors, 3)
        a_pos = first_pos[self._acceptor_a_idx]  # (n_acceptors, 3)

        tree = KDTree(d_pos)
        pairs = tree.query_ball_point(a_pos, PREFILTER_CUTOFF)

        pair_donor = []
        pair_acceptor = []
        for ai, nearby_donors in enumerate(pairs):
            for di in nearby_donors:
                pair_donor.append(di)
                pair_acceptor.append(ai)

        if not pair_donor:
            return []

        pair_donor = np.array(pair_donor)
        pair_acceptor = np.array(pair_acceptor)

        # 4. tuple_filter
        if tuple_filter is not None:
            mask = np.array([
                tuple_filter((donors[di], acceptors[ai]))
                for di, ai in zip(pair_donor, pair_acceptor)
            ])
            pair_donor = pair_donor[mask]
            pair_acceptor = pair_acceptor[mask]

        if len(pair_donor) == 0:
            return []

        self._pair_donor_idx = pair_donor
        self._pair_acceptor_idx = pair_acceptor

        n_pairs = len(pair_donor)
        n_frames = trajectory.n_frames

        # 5. 预分配结果数组
        existence = np.zeros((n_pairs, n_frames), dtype=bool)
        distance = np.zeros((n_pairs, n_frames))
        angle = np.zeros((n_pairs, n_frames))

        # 6. 逐帧计算
        for f, ts in enumerate(trajectory):
            positions = ts.positions
            dist, ang = self._compute_metrics(positions)
            distance[:, f] = dist
            angle[:, f] = ang
            existence[:, f] = (dist <= HBOND_DIST_MAX) & (ang >= HBOND_DON_ANGLE_MIN)

        # 7. 过滤从未存在的 pair
        has_any = np.any(existence, axis=1)
        if not np.any(has_any):
            return []

        tuples = [(donors[pair_donor[i]], acceptors[pair_acceptor[i]])
                  for i in range(n_pairs) if has_any[i]]
        existence = existence[has_any]
        distance = distance[has_any]
        angle = angle[has_any]

        return [Interaction(
            interaction_type=self.name,
            groups=tuples,
            existence=existence,
            metrics={"distance": distance, "angle": angle}
        )]

    # ==================== 内部辅助方法 ====================

    def _compute_metrics(self, positions):
        """向量化计算全部候选 pair 的 D-A 距离和 D-H-A 角度。

        Returns:
            dist: (n_pairs,) Å
            ang: (n_pairs,) °
        """
        # 取坐标：(n_pairs, 3)
        d_pos = positions[self._donor_d_idx[self._pair_donor_idx]]
        h_pos = positions[self._donor_h_idx[self._pair_donor_idx]]
        a_pos = positions[self._acceptor_a_idx[self._pair_acceptor_idx]]

        # D-A 距离
        da_vec = d_pos - a_pos
        dist = np.linalg.norm(da_vec, axis=1)

        # D-H···A 角度
        hd_vec = d_pos - h_pos    # H→D
        ha_vec = a_pos - h_pos    # H→A
        cos_angle = np.sum(hd_vec * ha_vec, axis=1) / (
            np.linalg.norm(hd_vec, axis=1) * np.linalg.norm(ha_vec, axis=1)
        )
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        ang = np.degrees(np.arccos(cos_angle))

        return dist, ang
