# -*- coding: utf-8 -*-
"""疏水相互作用检测器。

判据：两个疏水碳原子之间距离在 0.5~4.0 Å。
参考：PLIP。
"""

import numpy as np
from typing import List, Tuple, Dict
from itertools import combinations

from ..core.interfaces import InteractionDetectorPerTuple
from ..core.datas import Group


# PLIP 阈值
HYDROPH_DIST_MAX = 4.0  # Å，疏水原子间最大距离
HYDROPH_MIN_DIST = 0.5  # Å，最小距离


class HydrophobicDetectorPerTuple(InteractionDetectorPerTuple):
    """疏水相互作用检测器。"""

    # 预过滤 cutoff = 距离阈值 × 2。设为 None 禁用。
    PREFILTER_CUTOFF = HYDROPH_DIST_MAX * 2  # 8.0 Å

    @property
    def name(self) -> str:
        return "hydrophobic"

    @property
    def required_group_types(self) -> List[str]:
        return ["hydrophobic"]

    @property
    def metric_names(self) -> List[str]:
        return ["distance"]

    def get_candidate_tuples(self, groups: List[Group],
                             coordinates: np.ndarray = None) -> List[Tuple[Group, Group]]:
        """生成所有疏水原子对（不重复）。"""
        hydro = [g for g in groups if g.group_type == "hydrophobic"]
        return list(combinations(hydro, 2))

    def filter_candidate_tuples(self, tuples: List[Tuple[Group, Group]],
                                coordinates: np.ndarray) -> List[Tuple[Group, Group]]:
        """用第一帧坐标预过滤。"""
        if self.PREFILTER_CUTOFF is None:
            return tuples

        # 每个 hydrophobic 基团只有 1 个原子，直接取坐标
        centers = {}
        for gt in tuples:
            for g in gt:
                if id(g) not in centers:
                    centers[id(g)] = coordinates[g.atoms[0].atom_global_idx]

        return [(a, b) for a, b in tuples
                if np.linalg.norm(centers[id(a)] - centers[id(b)]) < self.PREFILTER_CUTOFF]

    def compute_metrics(self, group_tuple: Tuple[Group, Group],
                        coords: np.ndarray) -> Dict[str, np.ndarray]:
        """计算疏水相互作用指标。

        Args:
            group_tuple: (hydro1, hydro2)
            coords: (F, 2, 3)

        Returns:
            {"distance": (F,)}
        """
        pos1 = coords[:, 0, :]  # (F, 3)
        pos2 = coords[:, 1, :]  # (F, 3)
        distance = np.linalg.norm(pos1 - pos2, axis=1)  # (F,)
        return {"distance": distance}

    def apply_threshold(self, metrics: Dict[str, np.ndarray]) -> np.ndarray:
        """距离在 0.5~4.0 Å 之间。"""
        return (metrics["distance"] > HYDROPH_MIN_DIST) & \
               (metrics["distance"] < HYDROPH_DIST_MAX)
