# -*- coding: utf-8 -*-
"""π-阳离子相互作用检测器。

判据：电荷中心到环心距离 < 6.0 Å，投影偏移 < 2.0 Å。
参考：PLIP, Gallivan and Dougherty 1999。
"""

import numpy as np
from typing import List, Tuple, Dict
from itertools import product

from ..core.interfaces import InteractionDetector
from ..core.datas import Group


# PLIP 阈值
PICATION_DIST_MAX = 6.0    # Å，电荷中心到环心最大距离
PICATION_OFFSET_MAX = 2.0  # Å，最大投影偏移
PICATION_MIN_DIST = 0.5    # Å，最小距离


class PiCationDetector(InteractionDetector):
    """π-阳离子相互作用检测器。"""

    # 预过滤 cutoff = 距离阈值 × 3。设为 None 禁用。
    PREFILTER_CUTOFF = PICATION_DIST_MAX * 3  # 18.0 Å

    @property
    def name(self) -> str:
        return "pi_cation"

    @property
    def required_group_types(self) -> List[str]:
        return ["aromatic_ring", "charged_positive"]

    @property
    def metric_names(self) -> List[str]:
        return ["distance", "offset"]

    def get_candidate_tuples(self, groups: List[Group],
                             coordinates: np.ndarray = None) -> List[Tuple[Group, Group]]:
        """生成所有环-正电荷对（笛卡尔积）。"""
        rings = [g for g in groups if g.group_type == "aromatic_ring"]
        pos = [g for g in groups if g.group_type == "charged_positive"]
        return list(product(rings, pos))

    def filter_candidate_tuples(self, tuples: List[Tuple[Group, Group]],
                                coordinates: np.ndarray) -> List[Tuple[Group, Group]]:
        """用第一帧坐标预过滤。"""
        if self.PREFILTER_CUTOFF is None:
            return tuples

        # 计算每个基团的中心（按 id 缓存）
        centers = {}
        for gt in tuples:
            for g in gt:
                if id(g) not in centers:
                    if g.group_type == "aromatic_ring":
                        idx = np.array(g.atom_indices)
                        centers[id(g)] = np.mean(coordinates[idx], axis=0)
                    elif g.group_type == "charged_positive":
                        centers[id(g)] = self._charge_center_first_frame(g, coordinates)

        return [(r, p) for r, p in tuples
                if np.linalg.norm(centers[id(r)] - centers[id(p)]) < self.PREFILTER_CUTOFF]

    def _charge_center_first_frame(self, group: Group, coordinates: np.ndarray) -> np.ndarray:
        """计算第一帧的电荷中心。"""
        idx = np.array(group.atom_indices)
        coords = coordinates[idx]  # (n, 3)
        q = np.array([a.atom_charge for a in group.atoms])
        return np.sum(coords * q[:, None], axis=0) / np.sum(q)

    def compute_metrics(self, group_tuple: Tuple[Group, Group],
                        coords: np.ndarray) -> Dict[str, np.ndarray]:
        """计算 π-阳离子指标。

        Args:
            group_tuple: (aromatic_ring, charged_positive)
            coords: (F, n_atoms, 3)
                原子顺序：[ring_atoms..., charged_atoms...]

        Returns:
            {"distance": (F,), "offset": (F,)}
        """
        ring, pos = group_tuple
        n_ring = len(ring.atoms)

        ring_coords = coords[:, :n_ring, :]    # (F, n_ring, 3)
        pos_coords = coords[:, n_ring:, :]     # (F, n_pos, 3)

        # 环心
        ring_center = np.mean(ring_coords, axis=1)  # (F, 3)

        # 电荷中心
        pos_q = np.array([a.atom_charge for a in pos.atoms])
        charge_center = np.sum(pos_coords * pos_q[None, :, None], axis=1) / np.sum(pos_q)  # (F, 3)

        # 距离
        distance = np.linalg.norm(charge_center - ring_center, axis=1)  # (F,)

        # 投影偏移
        ring_normal = self._ring_normal(ring_coords)  # (F, 3)
        offset = self._projection_distance(ring_normal, ring_center, charge_center)  # (F,)

        return {"distance": distance, "offset": offset}

    def apply_threshold(self, metrics: Dict[str, np.ndarray]) -> np.ndarray:
        """逐帧判定是否存在 π-阳离子相互作用。"""
        return (metrics["distance"] > PICATION_MIN_DIST) & \
               (metrics["distance"] < PICATION_DIST_MAX) & \
               (metrics["offset"] < PICATION_OFFSET_MAX)

    # ==================== 内部辅助方法（复制自 PiStackingDetector） ====================

    def _ring_normal(self, ring_coords: np.ndarray) -> np.ndarray:
        """计算环法向量。严格复现 PLIP：每个原子与两个环内邻居的叉积，再平均，再归一化。

        Args:
            ring_coords: (F, n, 3)
        Returns:
            normal: (F, 3) 单位向量
        """
        prev = np.roll(ring_coords, 1, axis=1)    # 前一个邻居
        next_ = np.roll(ring_coords, -1, axis=1)  # 后一个邻居
        vec1 = ring_coords - prev                  # (F, n, 3)
        vec2 = ring_coords - next_                 # (F, n, 3)
        normals = np.cross(vec1, vec2)             # (F, n, 3)
        normal = np.mean(normals, axis=1)          # (F, 3)
        norm = np.linalg.norm(normal, axis=1, keepdims=True)
        return normal / norm

    def _projection_distance(self, normal: np.ndarray,
                             plane_point: np.ndarray,
                             target_point: np.ndarray) -> np.ndarray:
        """将 target_point 投影到 (normal, plane_point) 平面，返回投影点到 plane_point 的距离。

        Args:
            normal: (F, 3) 平面法向量
            plane_point: (F, 3) 平面上一点
            target_point: (F, 3) 待投影的点
        Returns:
            dist: (F,)
        """
        # 选法向量方向使朝向 target_point
        d1 = np.linalg.norm(target_point - (plane_point + normal), axis=1)
        d2 = np.linalg.norm(target_point - (plane_point - normal), axis=1)
        sign = np.where(d1 < d2, 1.0, -1.0)[:, None]  # (F, 1)
        oriented_normal = normal * sign                  # (F, 3)

        # 正交投影
        t = target_point - plane_point                   # (F, 3)
        proj_dist = np.sum(t * oriented_normal, axis=1)  # (F,)
        proj_point = target_point - proj_dist[:, None] * oriented_normal  # (F, 3)

        return np.linalg.norm(proj_point - plane_point, axis=1)  # (F,)
