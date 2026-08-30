# -*- coding: utf-8 -*-
"""氢键检测器。

判据：D-A 距离 ≤ 4.1 Å 且 D-H···A 角度 ≥ 100°。
参考：PLIP, Hubbard & Haider 2001（距离 +0.6 Å，角度 +10°）。
"""

import numpy as np
from typing import List, Tuple, Dict

from ..core.interfaces import InteractionDetector
from ..core.datas import Group


# PLIP 阈值
HBOND_DIST_MAX = 4.1    # Å，D-A 最大距离
HBOND_DON_ANGLE_MIN = 100.0  # °，D-H···A 最小角度


class HydrogenBondDetector(InteractionDetector):
    """氢键检测器。"""

    # 预过滤 cutoff = 距离阈值 × 3。设为 None 禁用。
    PREFILTER_CUTOFF = HBOND_DIST_MAX * 3  # 12.3 Å

    @property
    def name(self) -> str:
        return "hydrogen_bond"

    @property
    def required_group_types(self) -> List[str]:
        return ["H_donor", "H_acceptor"]

    @property
    def metric_names(self) -> List[str]:
        return ["distance", "angle"]

    def get_candidate_tuples(self, groups: List[Group],
                             coordinates: np.ndarray = None) -> List[Tuple[Group, Group]]:
        """生成所有供体-受体基团组。"""
        donors = [g for g in groups if g.group_type == "H_donor"]
        acceptors = [g for g in groups if g.group_type == "H_acceptor"]
        return [(d, a) for d in donors for a in acceptors]

    def filter_candidate_tuples(self, tuples: List[Tuple[Group, Group]],
                                coordinates: np.ndarray) -> List[Tuple[Group, Group]]:
        """用第一帧坐标预过滤：D-A 距离 > cutoff 的对直接排除。"""
        if self.PREFILTER_CUTOFF is None:
            return tuples

        # H_donor atoms=[D, H]，D 是 atoms[0]
        # H_acceptor atoms=[A]，A 是 atoms[0]
        # _get_atom_indices 拼接顺序：[D, H, A]
        return [(d, a) for d, a in tuples
                if self._da_distance_first_frame(d, a, coordinates) < self.PREFILTER_CUTOFF]

    def _da_distance_first_frame(self, donor: Group, acceptor: Group,
                                 coordinates: np.ndarray) -> float:
        """计算第一帧的 D-A 距离。"""
        d_pos = coordinates[donor.atoms[0].atom_global_idx]   # D
        a_pos = coordinates[acceptor.atoms[0].atom_global_idx]  # A
        return float(np.linalg.norm(d_pos - a_pos))

    def compute_metrics(self, group_tuple: Tuple[Group, Group],
                        coords: np.ndarray) -> Dict[str, np.ndarray]:
        """计算 D-A 距离和 D-H···A 角度。

        Args:
            group_tuple: (H_donor, H_acceptor)
            coords: (F, n_atoms, 3) 基团原子在全部帧的坐标（Å）
                原子顺序：[D, H, A]

        Returns:
            {"distance": (F,), "angle": (F,)}
        """
        # H_donor atoms=[D, H]，H_acceptor atoms=[A]
        # coords 拼接顺序：[D, H, A]
        d_pos = coords[:, 0, :]   # (F, 3)
        h_pos = coords[:, 1, :]   # (F, 3)
        a_pos = coords[:, 2, :]   # (F, 3)

        # D-A 距离
        da_vec = d_pos - a_pos                        # (F, 3)
        distance = np.linalg.norm(da_vec, axis=1)     # (F,)

        # D-H···A 角度
        hd_vec = d_pos - h_pos                        # (F, 3) H→D
        ha_vec = a_pos - h_pos                        # (F, 3) H→A

        cos_angle = np.sum(hd_vec * ha_vec, axis=1) / (
            np.linalg.norm(hd_vec, axis=1) * np.linalg.norm(ha_vec, axis=1)
        )
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        angle = np.degrees(np.arccos(cos_angle))      # (F,)

        return {"distance": distance, "angle": angle}

    def apply_threshold(self, metrics: Dict[str, np.ndarray]) -> np.ndarray:
        """距离 ≤ 4.1 Å 且角度 ≥ 100°。"""
        return (metrics["distance"] <= HBOND_DIST_MAX) & \
               (metrics["angle"] >= HBOND_DON_ANGLE_MIN)
