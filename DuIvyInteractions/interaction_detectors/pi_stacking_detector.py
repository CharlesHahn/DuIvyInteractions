# -*- coding: utf-8 -*-
"""π-π 堆积检测器。

判据：环心距离 ≤ 5.5 Å，法向量夹角满足 P 型或 T 型，投影偏移 ≤ 2.0 Å。
参考：PLIP, McGaughey 1998。
"""

import numpy as np
from typing import List, Tuple, Dict
from itertools import combinations

from ..core.interfaces import InteractionDetector
from ..core.datas import Group


# PLIP 阈值
PISTACK_DIST_MAX = 5.5       # Å，环心最大距离
PISTACK_ANG_DEV = 30.0       # °，角度偏差上限
PISTACK_OFFSET_MAX = 2.0     # Å，最大投影偏移
PISTACK_PLANARITY = 5.0      # °，平面性检验阈值
PISTACK_MIN_DIST = 0.5       # Å，最小距离（排除自身）


class PiStackingDetector(InteractionDetector):
    """π-π 堆积检测器。"""

    # 预过滤 cutoff = 距离阈值 × 3。设为 None 禁用。
    PREFILTER_CUTOFF = PISTACK_DIST_MAX * 3  # 16.5 Å

    def __init__(self, check_planarity: bool = False):
        """初始化。

        Args:
            check_planarity: 是否逐帧检验环平面性。默认 False。
        """
        self.check_planarity = check_planarity

    @property
    def name(self) -> str:
        return "pi_stacking"

    @property
    def required_group_types(self) -> List[str]:
        return ["aromatic_ring"]

    @property
    def metric_names(self) -> List[str]:
        names = ["distance", "angle", "offset", "pistacking_type"]
        if self.check_planarity:
            names.extend(["planarity_ring1", "planarity_ring2"])
        return names

    def get_candidate_tuples(self, groups: List[Group],
                             coordinates: np.ndarray = None) -> List[Tuple[Group, Group]]:
        """生成所有环对（不重复）。"""
        rings = [g for g in groups if g.group_type == "aromatic_ring"]
        return list(combinations(rings, 2))

    def filter_candidate_tuples(self, tuples: List[Tuple[Group, Group]],
                                coordinates: np.ndarray) -> List[Tuple[Group, Group]]:
        """用第一帧坐标预过滤。"""
        if self.PREFILTER_CUTOFF is None and not self.check_planarity:
            return tuples

        # 计算每个环的几何质心（按 id 缓存）
        centers = {}
        for gt in tuples:
            for g in gt:
                if id(g) not in centers:
                    idx = np.array(g.atom_indices)
                    centers[id(g)] = np.mean(coordinates[idx], axis=0)

        # 距离预过滤
        if self.PREFILTER_CUTOFF is not None:
            tuples = [(r1, r2) for r1, r2 in tuples
                      if np.linalg.norm(centers[id(r1)] - centers[id(r2)]) < self.PREFILTER_CUTOFF]

        # 平面性预过滤（第一帧）
        if self.check_planarity:
            tuples = [(r1, r2) for r1, r2 in tuples
                      if self._is_planar_first_frame(r1, coordinates)
                      and self._is_planar_first_frame(r2, coordinates)]

        return tuples

    def _is_planar_first_frame(self, ring: Group, coordinates: np.ndarray) -> bool:
        """检查第一帧环是否满足平面性。"""
        idx = np.array(ring.atom_indices)
        coords = coordinates[idx]  # (n, 3)
        normals = self._compute_normals_single_frame(coords)  # (n, 3)
        max_angle = self._max_normal_angle_single_frame(normals)
        return max_angle <= PISTACK_PLANARITY

    def compute_metrics(self, group_tuple: Tuple[Group, Group],
                        coords: np.ndarray) -> Dict[str, np.ndarray]:
        """计算 π-π 堆积指标。

        Args:
            group_tuple: (ring1, ring2)
            coords: (F, n_atoms, 3)
                原子顺序：[ring1_atoms..., ring2_atoms...]

        Returns:
            指标字典
        """
        ring1, ring2 = group_tuple
        n1 = len(ring1.atoms)

        ring1_coords = coords[:, :n1, :]    # (F, n1, 3)
        ring2_coords = coords[:, n1:, :]    # (F, n2, 3)

        # 环心
        center1 = np.mean(ring1_coords, axis=1)  # (F, 3)
        center2 = np.mean(ring2_coords, axis=1)  # (F, 3)

        # 环心距离
        distance = np.linalg.norm(center1 - center2, axis=1)  # (F,)

        # 法向量
        normal1 = self._ring_normal(ring1_coords)  # (F, 3)
        normal2 = self._ring_normal(ring2_coords)  # (F, 3)

        # 法向量夹角（取最小角，处理方向不确定性）
        cos_angle = np.sum(normal1 * normal2, axis=1)
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        raw_angle = np.degrees(np.arccos(cos_angle))  # (F,)
        angle = np.minimum(raw_angle, 180.0 - raw_angle)  # (F,)

        # 投影偏移量（双向取较小值）
        offset1 = self._projection_distance(normal2, center2, center1)
        offset2 = self._projection_distance(normal1, center1, center2)
        offset = np.minimum(offset1, offset2)  # (F,)

        # 堆积分类
        pistacking_type = self._classify(angle, offset)

        result = {
            "distance": distance,
            "angle": angle,
            "offset": offset,
            "pistacking_type": pistacking_type,
        }

        # 平面性（逐帧）
        if self.check_planarity:
            normals1 = self._normals_per_frame(ring1_coords)  # (F, n1, 3)
            normals2 = self._normals_per_frame(ring2_coords)  # (F, n2, 3)
            result["planarity_ring1"] = self._max_normal_angle(normals1)  # (F,)
            result["planarity_ring2"] = self._max_normal_angle(normals2)  # (F,)

        return result

    def apply_threshold(self, metrics: Dict[str, np.ndarray]) -> np.ndarray:
        """逐帧判定是否存在 π-π 堆积。"""
        base = (metrics["distance"] > PISTACK_MIN_DIST) & \
               (metrics["distance"] <= PISTACK_DIST_MAX) & \
               (metrics["pistacking_type"] != 'N')

        if self.check_planarity:
            planar = (metrics["planarity_ring1"] <= PISTACK_PLANARITY) & \
                     (metrics["planarity_ring2"] <= PISTACK_PLANARITY)
            return base & planar

        return base

    # ==================== 内部辅助方法 ====================

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

    def _classify(self, angle: np.ndarray, offset: np.ndarray) -> np.ndarray:
        """分类 P 型 / T 型 / N 型。

        Args:
            angle: (F,) 法向量夹角（°）
            offset: (F,) 投影偏移量（Å）
        Returns:
            pistacking_type: (F,) dtype='U1'，'N'/'P'/'T'
        """
        pistacking_type = np.full(angle.shape[0], 'N', dtype='U1')

        # P 型：angle ≤ 30°（值域 [0°, 90°]，含完全平行）
        p_mask = (angle <= PISTACK_ANG_DEV) & (offset < PISTACK_OFFSET_MAX)
        # T 型：angle ≥ 60°（值域 [0°, 90°]，含完全垂直）
        t_mask = (angle >= 90 - PISTACK_ANG_DEV) & (offset < PISTACK_OFFSET_MAX)

        pistacking_type[p_mask] = 'P'
        pistacking_type[t_mask] = 'T'

        return pistacking_type

    def _normals_per_frame(self, ring_coords: np.ndarray) -> np.ndarray:
        """计算每帧每个原子的法向量。

        Args:
            ring_coords: (F, n, 3)
        Returns:
            normals: (F, n, 3)
        """
        prev = np.roll(ring_coords, 1, axis=1)
        next_ = np.roll(ring_coords, -1, axis=1)
        vec1 = ring_coords - prev
        vec2 = ring_coords - next_
        return np.cross(vec1, vec2)  # (F, n, 3)

    def _max_normal_angle(self, normals: np.ndarray) -> np.ndarray:
        """计算每帧环内法向量两两夹角的最大值。

        Args:
            normals: (F, n, 3)
        Returns:
            max_angle: (F,)
        """
        n = normals.shape[1]
        max_angle = np.zeros(normals.shape[0])

        for i in range(n):
            for j in range(i + 1, n):
                cos = np.sum(normals[:, i, :] * normals[:, j, :], axis=1)
                norm_i = np.linalg.norm(normals[:, i, :], axis=1)
                norm_j = np.linalg.norm(normals[:, j, :], axis=1)
                cos = cos / (norm_i * norm_j)
                cos = np.clip(cos, -1.0, 1.0)
                angle = np.degrees(np.arccos(cos))
                max_angle = np.maximum(max_angle, angle)

        return max_angle

    def _compute_normals_single_frame(self, coords: np.ndarray) -> np.ndarray:
        """计算单帧环内每个原子的法向量。

        Args:
            coords: (n, 3)
        Returns:
            normals: (n, 3)
        """
        n = len(coords)
        normals = np.zeros((n, 3))
        for i in range(n):
            prev_idx = (i - 1) % n
            next_idx = (i + 1) % n
            vec1 = coords[i] - coords[prev_idx]
            vec2 = coords[i] - coords[next_idx]
            normals[i] = np.cross(vec1, vec2)
        return normals

    def _max_normal_angle_single_frame(self, normals: np.ndarray) -> float:
        """计算单帧环内法向量两两夹角的最大值。

        Args:
            normals: (n, 3)
        Returns:
            max_angle: float
        """
        n = len(normals)
        max_angle = 0.0
        for i in range(n):
            for j in range(i + 1, n):
                cos = np.dot(normals[i], normals[j])
                norm_i = np.linalg.norm(normals[i])
                norm_j = np.linalg.norm(normals[j])
                cos = cos / (norm_i * norm_j)
                cos = np.clip(cos, -1.0, 1.0)
                angle = np.degrees(np.arccos(cos))
                max_angle = max(max_angle, angle)
        return max_angle
