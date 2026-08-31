# -*- coding: utf-8 -*-
"""金属配位检测器（逐帧策略，向量化）。

判据：金属离子到配位原子距离 < 3.0 Å。
参考：PLIP, Harding 2001。
"""

import numpy as np
from typing import List, Tuple, Dict

from ..core.interfaces import InteractionDetectorPerFrame
from ..core.datas import Group, Interaction


# PLIP 阈值
METAL_DIST_MAX = 3.0  # Å，金属到配位原子最大距离
PREFILTER_CUTOFF = METAL_DIST_MAX * 2  # 6.0 Å


class MetalCoordinationDetectorPerFrame(InteractionDetectorPerFrame):
    """金属配位检测器（逐帧策略，向量化）。"""

    def __init__(self):
        self._metal_idx = None
        self._binding_idx = None
        self._pair_metal = None
        self._pair_binding = None

    @property
    def name(self) -> str:
        return "metal_coordination"

    @property
    def required_group_types(self) -> List[str]:
        return ["metal", "metal_binding"]

    @property
    def metric_names(self) -> List[str]:
        return ["distance"]

    # ==================== 抽象方法（基类要求，此处不使用） ====================

    def get_candidate_tuples(self, groups, coordinates=None):
        raise NotImplementedError("detect() is overridden, this method is not used")

    def compute_metrics_for_frame(self, tuples, all_positions, frame):
        raise NotImplementedError("detect() is overridden, this method is not used")

    def apply_threshold(self, metrics):
        return metrics["distance"] < METAL_DIST_MAX

    # ==================== 重写 detect ====================

    def detect(self, groups, trajectory=None,
               n_workers: int = 1,
               topology_path: str = None,
               trajectory_path: str = None,
               tuple_filter=None) -> List[Interaction]:
        """金属配位检测主流程。"""
        if trajectory is None:
            raise ValueError("trajectory is required")
        if n_workers > 1:
            raise NotImplementedError("PerFrame 检测器暂不支持并行")

        # 1. 分组
        metals = [g for g in groups if g.group_type == "metal"]
        bindings = [g for g in groups if g.group_type == "metal_binding"]

        if not metals or not bindings:
            return []

        # 2. 构建索引数组（单原子 group，无 padding）
        self._metal_idx = np.array([g.atoms[0].atom_global_idx for g in metals])
        self._binding_idx = np.array([g.atoms[0].atom_global_idx for g in bindings])

        # 3. 笛卡尔积
        n_m = len(metals)
        n_b = len(bindings)
        m_grid, b_grid = np.meshgrid(np.arange(n_m), np.arange(n_b), indexing='ij')
        pair_metal = m_grid.ravel()
        pair_binding = b_grid.ravel()

        # 4. tuple_filter
        if tuple_filter is not None:
            mask = np.array([
                tuple_filter((metals[mi], bindings[bi]))
                for mi, bi in zip(pair_metal, pair_binding)
            ])
            pair_metal = pair_metal[mask]
            pair_binding = pair_binding[mask]

        if len(pair_metal) == 0:
            return []

        # 5. 第一帧：距离预过滤
        first_pos = trajectory[0].positions
        m_pos = first_pos[self._metal_idx[pair_metal]]
        b_pos = first_pos[self._binding_idx[pair_binding]]
        dist = np.linalg.norm(m_pos - b_pos, axis=1)
        mask = dist < PREFILTER_CUTOFF

        pair_metal = pair_metal[mask]
        pair_binding = pair_binding[mask]

        if len(pair_metal) == 0:
            return []

        self._pair_metal = pair_metal
        self._pair_binding = pair_binding

        n_pairs = len(pair_metal)
        n_frames = trajectory.n_frames

        # 6. 预分配结果数组
        existence = np.zeros((n_pairs, n_frames), dtype=bool)
        distance = np.zeros((n_pairs, n_frames))

        # 7. 逐帧计算
        for f, ts in enumerate(trajectory):
            positions = ts.positions
            m = positions[self._metal_idx[self._pair_metal]]
            b = positions[self._binding_idx[self._pair_binding]]
            dist = np.linalg.norm(m - b, axis=1)
            distance[:, f] = dist
            existence[:, f] = dist < METAL_DIST_MAX

        # 8. 过滤从未存在的 pair
        has_any = np.any(existence, axis=1)
        if not np.any(has_any):
            return []

        tuples = [(metals[pair_metal[i]], bindings[pair_binding[i]])
                  for i in range(n_pairs) if has_any[i]]

        return [Interaction(
            interaction_type=self.name,
            groups=tuples,
            existence=existence[has_any],
            metrics={"distance": distance[has_any]}
        )]
