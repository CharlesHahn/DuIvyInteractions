# -*- coding: utf-8 -*-
"""金属配位检测器。

判据：金属离子到配位原子距离 < 3.0 Å。
参考：PLIP, Harding 2001。
"""

import numpy as np
from typing import List, Tuple, Dict
from itertools import product

from ..core.interfaces import InteractionDetectorPerTuple
from ..core.datas import Group


# PLIP 阈值
METAL_DIST_MAX = 3.0  # Å，金属到配位原子最大距离


class MetalCoordinationDetectorPerTuple(InteractionDetectorPerTuple):
    """金属配位检测器。"""

    # 预过滤 cutoff = 距离阈值 × 2。设为 None 禁用。
    PREFILTER_CUTOFF = METAL_DIST_MAX * 2  # 6.0 Å

    @property
    def name(self) -> str:
        return "metal_coordination"

    @property
    def required_group_types(self) -> List[str]:
        return ["metal", "metal_binding"]

    @property
    def metric_names(self) -> List[str]:
        return ["distance"]

    def get_candidate_tuples(self, groups: List[Group],
                             coordinates: np.ndarray = None) -> List[Tuple[Group, Group]]:
        """生成所有金属-配位原子对。"""
        metals = [g for g in groups if g.group_type == "metal"]
        binding = [g for g in groups if g.group_type == "metal_binding"]
        return list(product(metals, binding))

    def filter_candidate_tuples(self, tuples: List[Tuple[Group, Group]],
                                coordinates: np.ndarray) -> List[Tuple[Group, Group]]:
        """用第一帧坐标预过滤。"""
        if self.PREFILTER_CUTOFF is None:
            return tuples

        centers = {}
        for gt in tuples:
            for g in gt:
                if id(g) not in centers:
                    centers[id(g)] = coordinates[g.atoms[0].atom_global_idx]

        return [(m, a) for m, a in tuples
                if np.linalg.norm(centers[id(m)] - centers[id(a)]) < self.PREFILTER_CUTOFF]

    def compute_metrics(self, group_tuple: Tuple[Group, Group],
                        coords: np.ndarray) -> Dict[str, np.ndarray]:
        """计算金属到配位原子的距离。

        Args:
            group_tuple: (metal, binding_atom)
            coords: (F, 2, 3)

        Returns:
            {"distance": (F,)}
        """
        metal_pos = coords[:, 0, :]   # (F, 3)
        atom_pos = coords[:, 1, :]    # (F, 3)
        distance = np.linalg.norm(metal_pos - atom_pos, axis=1)  # (F,)
        return {"distance": distance}

    def apply_threshold(self, metrics: Dict[str, np.ndarray]) -> np.ndarray:
        """距离 < 3.0 Å。"""
        return metrics["distance"] < METAL_DIST_MAX
