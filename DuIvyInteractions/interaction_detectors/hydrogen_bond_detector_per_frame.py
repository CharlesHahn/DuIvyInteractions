# -*- coding: utf-8 -*-
"""氢键检测器（逐帧策略，KDTree 预筛选 + 向量化计算）。

判据：D-A 距离 ≤ 4.1 Å 且 D-H···A 角度 ≥ 100°。
参考：PLIP, Hubbard & Haider 2001（距离 +0.6 Å，角度 +10°）。

PBC：候选预筛 capped_distance(box)、逐帧距离 calc_bonds(box)、角度 calc_angles(box)，
见 doc/pbc_design.md §6A。
"""

import numpy as np
from typing import List, Tuple, Dict

from ..core.interfaces import InteractionDetectorPerFrame
from ..core.datas import Group, Interaction
from ..core.pbc import get_box


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
            raise NotImplementedError("PerFrame detector does not support parallelism yet")

        # 1. 分组
        donors = [g for g in groups if g.group_type == "H_donor"]
        acceptors = [g for g in groups if g.group_type == "H_acceptor"]

        if not donors or not acceptors:
            return []

        # 2. 构建索引数组
        self._donor_d_idx = np.array([g.atoms[0].atom_global_idx for g in donors])
        self._donor_h_idx = np.array([g.atoms[1].atom_global_idx for g in donors])
        self._acceptor_a_idx = np.array([g.atoms[0].atom_global_idx for g in acceptors])

        # 3. 第一帧 KDTree 预筛选（PBC 感知）
        from MDAnalysis.lib.distances import capped_distance
        first_pos = trajectory[0].positions
        d_pos = first_pos[self._donor_d_idx]     # (n_donors, 3)
        a_pos = first_pos[self._acceptor_a_idx]  # (n_acceptors, 3)

        # ref=a_pos, conf=d_pos → pairs[:,0]=acceptor, pairs[:,1]=donor
        pairs = capped_distance(
            a_pos, d_pos, PREFILTER_CUTOFF,
            box=get_box(trajectory[0]), method='pkdtree', return_distances=False)

        if len(pairs) == 0:
            return []

        pair_acceptor = pairs[:, 0]
        pair_donor = pairs[:, 1]

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
        times = np.empty(n_frames)

        # 6. 逐帧计算（PBC 感知）
        for f, ts in enumerate(trajectory):
            times[f] = ts.time
            positions = ts.positions
            dist, ang = self._compute_metrics(positions, get_box(ts))
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
            metrics={"distance": distance, "angle": angle},
            times=times,
        )]

    # ==================== 内部辅助方法 ====================

    def _compute_metrics(self, positions, box=None):
        """向量化计算全部候选 pair 的 D-A 距离和 D-H-A 角度（MIC）。

        box: 盒尺寸（MIC）；None = 欧氏距离（无 PBC）。

        Returns:
            dist: (n_pairs,) Å
            ang: (n_pairs,) °
        """
        from MDAnalysis.lib.distances import calc_bonds, calc_angles

        # 取坐标：(n_pairs, 3)
        d_pos = positions[self._donor_d_idx[self._pair_donor_idx]]
        h_pos = positions[self._donor_h_idx[self._pair_donor_idx]]
        a_pos = positions[self._acceptor_a_idx[self._pair_acceptor_idx]]

        # D-A 距离（MIC）
        dist = calc_bonds(d_pos, a_pos, box=box)

        # D-H···A 角度（MIC，顶点 h；calc_angles 返回弧度）
        ang = np.degrees(calc_angles(d_pos, h_pos, a_pos, box=box))

        return dist, ang
