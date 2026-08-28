# -*- coding: utf-8 -*-
"""盐桥检测器。

判据：两个相反电荷中心距离 ≤ SALTBRIDGE_DIST_MAX。
参考：PLIP, Barlow and Thornton 1983（+1.5 Å 扩展）。
"""

import numpy as np
from typing import List, Tuple, Dict

from ..core.interfaces import InteractionDetector
from ..core.datas import Group


# PLIP 阈值：5.5 Å = 0.55 nm
SALTBRIDGE_DIST_MAX = 0.55  # nm


class SaltBridgeDetector(InteractionDetector):
    """盐桥检测器。"""

    @property
    def name(self) -> str:
        return "salt_bridge"

    @property
    def required_group_types(self) -> List[str]:
        return ["charged_positive", "charged_negative"]

    @property
    def metric_names(self) -> List[str]:
        return ["distance"]

    def get_candidate_tuples(self, groups: List[Group]) -> List[Tuple[Group, Group]]:
        """生成所有正电-负电基团组。"""
        pos = [g for g in groups if g.group_type == "charged_positive"]
        neg = [g for g in groups if g.group_type == "charged_negative"]
        return [(p, n) for p in pos for n in neg]

    def compute_metrics(self, group_tuple: Tuple[Group, Group],
                        coords: np.ndarray) -> Dict[str, np.ndarray]:
        """计算电荷中心距离。

        Args:
            group_tuple: (正电基团, 负电基团)
            coords: (F, n_atoms, 3) 基团原子在全部帧的坐标（nm）

        Returns:
            {"distance": (F,)} 电荷中心距离（nm）
        """
        pos_group, neg_group = group_tuple
        n_pos = len(pos_group.atoms)

        # 拆分坐标
        pos_coords = coords[:, :n_pos, :]   # (F, n_pos, 3)
        neg_coords = coords[:, n_pos:, :]   # (F, n_neg, 3)

        # 电荷向量
        pos_q = np.array([a.atom_charge for a in pos_group.atoms])
        neg_q = np.array([a.atom_charge for a in neg_group.atoms])

        # 电荷加权中心：Σ(q_i × r_i) / Σ(q_i)
        pos_center = np.sum(pos_coords * pos_q[None, :, None], axis=1) / np.sum(pos_q)
        neg_center = np.sum(neg_coords * neg_q[None, :, None], axis=1) / np.sum(neg_q)

        # 距离
        distance = np.linalg.norm(pos_center - neg_center, axis=1)

        return {"distance": distance}

    def apply_threshold(self, metrics: Dict[str, np.ndarray]) -> np.ndarray:
        """距离 ≤ SALTBRIDGE_DIST_MAX。"""
        return metrics["distance"] < SALTBRIDGE_DIST_MAX
