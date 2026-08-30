# -*- coding: utf-8 -*-
"""水桥检测器。

判据：水分子位于供体和受体之间，距离和角度满足 PLIP 定义。
参考：PLIP, Jiang et al. 2005。
"""

import numpy as np
from typing import List, Tuple, Dict
from scipy.spatial import KDTree

from ..core.interfaces import InteractionDetector
from ..core.datas import Group


# PLIP 阈值
WATER_BRIDGE_MINDIST = 2.5   # Å，水 O 到极性原子最小距离
WATER_BRIDGE_MAXDIST = 4.1   # Å，水 O 到极性原子最大距离
WATER_BRIDGE_OMEGA_MIN = 71  # °，受体-水O-供体H 最小角度
WATER_BRIDGE_OMEGA_MAX = 140 # °，受体-水O-供体H 最大角度
WATER_BRIDGE_THETA_MIN = 100 # °，水O-供体H-供体D 最小角度
PREFILTER_CUTOFF = WATER_BRIDGE_MAXDIST * 2  # 8.2 Å


class WaterBridgeDetector(InteractionDetector):
    """水桥检测器。"""

    @property
    def name(self) -> str:
        return "water_bridge"

    @property
    def required_group_types(self) -> List[str]:
        return ["H_donor", "water", "H_acceptor"]

    @property
    def metric_names(self) -> List[str]:
        return ["dist_dw", "dist_wa", "theta", "omega"]

    def get_candidate_tuples(self, groups: List[Group],
                             coordinates: np.ndarray = None) -> List[Tuple[Group, Group, Group]]:
        """生成 (供体, 水, 受体) 三元组，用 KDTree 预筛选。

        预筛选规则：
        1. 供体和受体排除水分子
        2. 同一个水必须同时靠近一个供体和一个受体（< 8.2 Å）
        3. 供体和受体距离 > 2.5 Å（太近则直接氢键，不需要水桥）
        """
        from ..group_identifiers.amber_ff_identifier import WATER_RESIDUES

        donors = [g for g in groups
                  if g.group_type == "H_donor" and g.residue_name not in WATER_RESIDUES]
        waters = [g for g in groups if g.group_type == "water"]
        acceptors = [g for g in groups
                     if g.group_type == "H_acceptor" and g.residue_name not in WATER_RESIDUES]

        if not donors or not waters or not acceptors:
            return []

        # 建 KDTree
        donor_pos = np.array([coordinates[d.atoms[0].atom_global_idx] for d in donors])
        acceptor_pos = np.array([coordinates[a.atoms[0].atom_global_idx] for a in acceptors])
        donor_tree = KDTree(donor_pos)
        acceptor_tree = KDTree(acceptor_pos)

        # 每个水查询附近供体和受体
        tuples = []
        for w in waters:
            w_pos = coordinates[w.atoms[0].atom_global_idx]

            nearby_d_indices = donor_tree.query_ball_point(w_pos, PREFILTER_CUTOFF)
            nearby_a_indices = acceptor_tree.query_ball_point(w_pos, PREFILTER_CUTOFF)

            if not nearby_d_indices or not nearby_a_indices:
                continue

            for di in nearby_d_indices:
                d_pos = donor_pos[di]
                for ai in nearby_a_indices:
                    a_pos = acceptor_pos[ai]
                    # D-A 距离过滤：太近则直接氢键，不需要水桥
                    if np.linalg.norm(d_pos - a_pos) < WATER_BRIDGE_MINDIST:
                        continue
                    tuples.append((donors[di], w, acceptors[ai]))

        return tuples

    def filter_candidate_tuples(self, tuples, coordinates):
        """不需要额外过滤，已在 get_candidate_tuples 中完成。"""
        return tuples

    def compute_metrics(self, group_tuple: Tuple[Group, Group, Group],
                        coords: np.ndarray) -> Dict[str, np.ndarray]:
        """计算水桥指标。

        Args:
            group_tuple: (H_donor, water, H_acceptor)
            coords: (F, 6, 3) 原子顺序：[D, H, Ow, Hw1, Hw2, A]

        Returns:
            {"dist_dw": (F,), "dist_wa": (F,), "theta": (F,), "omega": (F,)}
        """
        d_pos = coords[:, 0, :]    # D    (F, 3)
        h_pos = coords[:, 1, :]    # H    (F, 3)
        ow_pos = coords[:, 2, :]   # Ow   (F, 3)
        a_pos = coords[:, 5, :]    # A    (F, 3)

        # 距离
        dist_dw = np.linalg.norm(d_pos - ow_pos, axis=1)   # (F,)
        dist_wa = np.linalg.norm(ow_pos - a_pos, axis=1)   # (F,)

        # theta 角：水O-供体H-供体D
        vec_hd = d_pos - h_pos       # H→D
        vec_ho = ow_pos - h_pos      # H→Ow
        theta = self._vec_angle(vec_hd, vec_ho)

        # omega 角：受体-水O-供体H
        vec_oa = a_pos - ow_pos      # Ow→A
        vec_oh = h_pos - ow_pos      # Ow→H
        omega = self._vec_angle(vec_oa, vec_oh)

        return {
            "dist_dw": dist_dw,
            "dist_wa": dist_wa,
            "theta": theta,
            "omega": omega,
        }

    def apply_threshold(self, metrics: Dict[str, np.ndarray]) -> np.ndarray:
        """距离和角度阈值判定。"""
        dw_ok = (metrics["dist_dw"] > WATER_BRIDGE_MINDIST) & \
                (metrics["dist_dw"] < WATER_BRIDGE_MAXDIST)
        wa_ok = (metrics["dist_wa"] > WATER_BRIDGE_MINDIST) & \
                (metrics["dist_wa"] < WATER_BRIDGE_MAXDIST)
        theta_ok = metrics["theta"] >= WATER_BRIDGE_THETA_MIN
        omega_ok = (metrics["omega"] >= WATER_BRIDGE_OMEGA_MIN) & \
                   (metrics["omega"] <= WATER_BRIDGE_OMEGA_MAX)

        return dw_ok & wa_ok & theta_ok & omega_ok

    @staticmethod
    def _vec_angle(v1: np.ndarray, v2: np.ndarray) -> np.ndarray:
        """计算两组向量的夹角（度）。v1, v2: (F, 3) → (F,)"""
        cos_angle = np.sum(v1 * v2, axis=1) / (
            np.linalg.norm(v1, axis=1) * np.linalg.norm(v2, axis=1)
        )
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        return np.degrees(np.arccos(cos_angle))
