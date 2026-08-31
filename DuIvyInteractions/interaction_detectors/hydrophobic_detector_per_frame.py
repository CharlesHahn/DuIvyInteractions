# -*- coding: utf-8 -*-
"""疏水相互作用检测器（逐帧策略，KDTree 预筛选 + 向量化计算）。

判据：两个疏水碳原子之间距离在 0.5~4.0 Å。
参考：PLIP。
"""

import numpy as np
from typing import List, Tuple, Dict
from scipy.spatial import KDTree

from ..core.interfaces import InteractionDetectorPerFrame
from ..core.datas import Group, Interaction


# PLIP 阈值
HYDROPH_DIST_MAX = 4.0  # Å，疏水原子间最大距离
HYDROPH_MIN_DIST = 0.5  # Å，最小距离
PREFILTER_CUTOFF = HYDROPH_DIST_MAX * 2  # 8.0 Å


class HydrophobicDetectorPerFrame(InteractionDetectorPerFrame):
    """疏水相互作用检测器（逐帧策略，KDTree 预筛选 + 向量化计算）。"""

    def __init__(self):
        self._hydro_idx = None
        self._pair_a = None
        self._pair_b = None

    @property
    def name(self) -> str:
        return "hydrophobic"

    @property
    def required_group_types(self) -> List[str]:
        return ["hydrophobic"]

    @property
    def metric_names(self) -> List[str]:
        return ["distance"]

    # ==================== 抽象方法（基类要求，此处不使用） ====================

    def get_candidate_tuples(self, groups, coordinates=None):
        raise NotImplementedError("detect() is overridden, this method is not used")

    def compute_metrics_for_frame(self, tuples, all_positions, frame):
        raise NotImplementedError("detect() is overridden, this method is not used")

    def apply_threshold(self, metrics):
        return (metrics["distance"] > HYDROPH_MIN_DIST) & \
               (metrics["distance"] < HYDROPH_DIST_MAX)

    # ==================== 重写 detect ====================

    def detect(self, groups, trajectory=None,
               n_workers: int = 1,
               topology_path: str = None,
               trajectory_path: str = None,
               tuple_filter=None) -> List[Interaction]:
        """疏水检测主流程。"""
        if trajectory is None:
            raise ValueError("trajectory is required")
        if n_workers > 1:
            raise NotImplementedError("PerFrame 检测器暂不支持并行")

        # 1. 筛选疏水原子
        hydro = [g for g in groups if g.group_type == "hydrophobic"]
        if len(hydro) < 2:
            return []

        # 2. 构建索引数组（单原子 group）
        self._hydro_idx = np.array([g.atoms[0].atom_global_idx for g in hydro])

        # 3. 第一帧 KDTree 预筛选
        first_pos = trajectory[0].positions
        h_pos = first_pos[self._hydro_idx]
        tree = KDTree(h_pos)
        pairs = tree.query_pairs(r=PREFILTER_CUTOFF, output_type='ndarray')

        if len(pairs) == 0:
            return []

        pair_a = pairs[:, 0]
        pair_b = pairs[:, 1]

        # 4. tuple_filter
        if tuple_filter is not None:
            mask = np.array([
                tuple_filter((hydro[ai], hydro[bi]))
                for ai, bi in zip(pair_a, pair_b)
            ])
            pair_a = pair_a[mask]
            pair_b = pair_b[mask]

        if len(pair_a) == 0:
            return []

        self._pair_a = pair_a
        self._pair_b = pair_b

        n_pairs = len(pair_a)
        n_frames = trajectory.n_frames

        # 5. 预分配结果数组
        existence = np.zeros((n_pairs, n_frames), dtype=bool)
        distance = np.zeros((n_pairs, n_frames))

        # 6. 逐帧计算
        for f, ts in enumerate(trajectory):
            positions = ts.positions
            pos_a = positions[self._hydro_idx[self._pair_a]]
            pos_b = positions[self._hydro_idx[self._pair_b]]
            dist = np.linalg.norm(pos_a - pos_b, axis=1)
            distance[:, f] = dist
            existence[:, f] = (dist > HYDROPH_MIN_DIST) & (dist < HYDROPH_DIST_MAX)

        # 7. 过滤从未存在的 pair
        has_any = np.any(existence, axis=1)
        if not np.any(has_any):
            return []

        tuples = [(hydro[pair_a[i]], hydro[pair_b[i]])
                  for i in range(n_pairs) if has_any[i]]
        existence = existence[has_any]
        distance = distance[has_any]

        # 8. 去重
        tuples, existence, distance = self._deduplicate(tuples, existence, distance)

        if not tuples:
            return []

        return [Interaction(
            interaction_type=self.name,
            groups=tuples,
            existence=existence,
            metrics={"distance": distance}
        )]

    # ==================== 去重 ====================

    @staticmethod
    def _deduplicate(tuples, existence, distance):
        """去重：同一原子与同残基多个原子的接触，只保留最近的。

        对于每个 (group_id, residue_id) 对，只保留平均距离最近的那个接触。
        一个 pair 必须是它所有 key 的最优解才会被保留。
        """
        n_pairs = len(tuples)
        if n_pairs == 0:
            return tuples, existence, distance

        # 第一步：找出每个 key 的最优 pair index
        best = {}  # (group_id, residue_id) → (index, avg_dist)

        for i, (g1, g2) in enumerate(tuples):
            active = existence[i]
            if not np.any(active):
                continue
            avg_dist = float(np.mean(distance[i][active]))

            key1 = (g1.group_id, g2.residue_id)
            key2 = (g2.group_id, g1.residue_id)

            for key in [key1, key2]:
                if key not in best or avg_dist < best[key][1]:
                    best[key] = (i, avg_dist)

        # 第二步：只保留对所有 key 都是最优的 pair
        keep = []
        for i, (g1, g2) in enumerate(tuples):
            key1 = (g1.group_id, g2.residue_id)
            key2 = (g2.group_id, g1.residue_id)
            if key1 in best and key2 in best:
                if best[key1][0] == i and best[key2][0] == i:
                    keep.append(i)

        if not keep:
            return [], np.array([]), np.array([])

        keep = np.array(keep)
        return [tuples[i] for i in keep], existence[keep], distance[keep]
