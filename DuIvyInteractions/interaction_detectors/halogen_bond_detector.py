# -*- coding: utf-8 -*-
"""卤键检测器。

判据：X-A ≤ 4.0 Å，C-X···A 在 165°±30°，X···A-R 在 120°±30°。
参考：PLIP, Auffinger et al.
"""

import numpy as np
from typing import List, Tuple, Dict

from ..core.interfaces import InteractionDetector
from ..core.datas import Group


# PLIP 阈值
HALOGEN_DIST_MAX = 4.0        # Å
HALOGEN_DON_ANGLE = 165.0     # °，C-X···A 最优角度
HALOGEN_ACC_ANGLE = 120.0     # °，X···A-R 最优角度
HALOGEN_ANGLE_DEV = 30.0      # °，角度偏差上限


class HalogenBondDetector(InteractionDetector):
    """卤键检测器。"""

    # 预过滤 cutoff = 距离阈值 × 3。设为 None 禁用。
    PREFILTER_CUTOFF = HALOGEN_DIST_MAX * 3  # 12.0 Å

    @property
    def name(self) -> str:
        return "halogen_bond"

    @property
    def required_group_types(self) -> List[str]:
        return ["halogen_donor", "halogen_acceptor"]

    @property
    def metric_names(self) -> List[str]:
        return ["distance", "don_angle", "acc_angle"]

    def get_candidate_tuples(self, groups: List[Group],
                             coordinates: np.ndarray = None) -> List[Tuple[Group, Group]]:
        """生成所有卤键供体-受体基团组。"""
        donors = [g for g in groups if g.group_type == "halogen_donor"]
        acceptors = [g for g in groups if g.group_type == "halogen_acceptor"]
        return [(d, a) for d in donors for a in acceptors]

    def filter_candidate_tuples(self, tuples: List[Tuple[Group, Group]],
                                coordinates: np.ndarray) -> List[Tuple[Group, Group]]:
        """用第一帧坐标预过滤：X-A 距离 > cutoff 的对直接排除。"""
        if self.PREFILTER_CUTOFF is None:
            return tuples

        # halogen_donor atoms=[C, X]，X 是 atoms[1]
        # halogen_acceptor atoms=[A, R1, R2, ...]，A 是 atoms[0]
        return [(d, a) for d, a in tuples
                if self._xa_distance_first_frame(d, a, coordinates) < self.PREFILTER_CUTOFF]

    def _xa_distance_first_frame(self, donor: Group, acceptor: Group,
                                 coordinates: np.ndarray) -> float:
        """计算第一帧的 X-A 距离。"""
        x_pos = coordinates[donor.atoms[1].atom_global_idx]    # X
        a_pos = coordinates[acceptor.atoms[0].atom_global_idx]  # A
        return float(np.linalg.norm(x_pos - a_pos))

    def compute_metrics(self, group_tuple: Tuple[Group, Group],
                        coords: np.ndarray) -> Dict[str, np.ndarray]:
        """计算 X-A 距离、C-X···A 角度、X···A-R 角度。

        Args:
            group_tuple: (halogen_donor, halogen_acceptor)
            coords: (F, n_atoms, 3) 基团原子在全部帧的坐标（Å）
                原子顺序：[C, X, A, R1, R2, ...]

        Returns:
            {"distance": (F,), "don_angle": (F,), "acc_angle": (F,)}
        """
        n_donor_atoms = len(group_tuple[0].atoms)  # 2

        # 取坐标
        c_pos = coords[:, 0, :]                       # C (F, 3)
        x_pos = coords[:, 1, :]                       # X (F, 3)
        a_pos = coords[:, 2, :]                       # A (F, 3)
        r_positions = coords[:, n_donor_atoms + 1:, :]  # R 原子 (F, n_r, 3)

        # X-A 距离
        xa_vec = x_pos - a_pos
        distance = np.linalg.norm(xa_vec, axis=1)

        # C-X···A 角度（don_angle）
        xc_vec = c_pos - x_pos    # X→C
        xa_vec2 = a_pos - x_pos   # X→A
        don_angle = self._vec_angle(xc_vec, xa_vec2)

        # X···A-R 角度（acc_angle）：对每个 R 计算，取最接近 120° 的
        ax_vec = x_pos - a_pos    # A→X
        n_r = r_positions.shape[1]
        all_acc_angles = np.zeros((n_r, coords.shape[0]))

        for i in range(n_r):
            ar_vec = r_positions[:, i, :] - a_pos  # A→R_i
            all_acc_angles[i] = self._vec_angle(ax_vec, ar_vec)

        # 每帧取最接近 HALOGEN_ACC_ANGLE 的 R
        diffs = np.abs(all_acc_angles - HALOGEN_ACC_ANGLE)
        best_r_idx = np.argmin(diffs, axis=0)
        frame_idx = np.arange(coords.shape[0])
        acc_angle = all_acc_angles[best_r_idx, frame_idx]

        return {"distance": distance, "don_angle": don_angle, "acc_angle": acc_angle}

    def apply_threshold(self, metrics: Dict[str, np.ndarray]) -> np.ndarray:
        """距离 ≤ 4.0 Å，don_angle 在 165°±30°，acc_angle 在 120°±30°。"""
        dist_ok = metrics["distance"] <= HALOGEN_DIST_MAX
        don_ok = (metrics["don_angle"] >= HALOGEN_DON_ANGLE - HALOGEN_ANGLE_DEV) & \
                 (metrics["don_angle"] <= HALOGEN_DON_ANGLE + HALOGEN_ANGLE_DEV)
        acc_ok = (metrics["acc_angle"] >= HALOGEN_ACC_ANGLE - HALOGEN_ANGLE_DEV) & \
                 (metrics["acc_angle"] <= HALOGEN_ACC_ANGLE + HALOGEN_ANGLE_DEV)
        return dist_ok & don_ok & acc_ok

    @staticmethod
    def _vec_angle(v1: np.ndarray, v2: np.ndarray) -> np.ndarray:
        """计算两组向量的夹角（度）。v1, v2: (F, 3) → (F,)"""
        cos_angle = np.sum(v1 * v2, axis=1) / (
            np.linalg.norm(v1, axis=1) * np.linalg.norm(v2, axis=1)
        )
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        return np.degrees(np.arccos(cos_angle))
