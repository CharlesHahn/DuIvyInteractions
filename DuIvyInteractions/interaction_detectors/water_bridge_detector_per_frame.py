# -*- coding: utf-8 -*-
"""水桥检测器（逐帧策略，KDTree 预筛选 + 向量化计算）。

判据：水分子位于供体和受体之间，距离和角度满足 PLIP 定义。
参考：PLIP, Jiang et al. 2005。
"""

import numpy as np
from typing import List, Tuple, Dict
from scipy.spatial import KDTree

from ..core.interfaces import InteractionDetectorPerFrame
from ..core.datas import Group, Interaction


# PLIP 阈值
WATER_BRIDGE_MINDIST = 2.5   # Å，水 O 到极性原子最小距离
WATER_BRIDGE_MAXDIST = 4.1   # Å，水 O 到极性原子最大距离
WATER_BRIDGE_OMEGA_MIN = 71  # °，受体-水O-供体H 最小角度
WATER_BRIDGE_OMEGA_MAX = 140 # °，受体-水O-供体H 最大角度
WATER_BRIDGE_THETA_MIN = 100 # °，水O-供体H-供体D 最小角度
PREFILTER_CUTOFF = WATER_BRIDGE_MAXDIST * 2  # 8.2 Å


class WaterBridgeDetectorPerFrame(InteractionDetectorPerFrame):
    """水桥检测器（逐帧策略，KDTree 预筛选 + 向量化计算）。"""

    def __init__(self):
        self._donor_d_idx = None
        self._donor_h_idx = None
        self._water_ow_idx = None
        self._acceptor_a_idx = None
        self._triple_d = None
        self._triple_w = None
        self._triple_a = None

    @property
    def name(self) -> str:
        return "water_bridge"

    @property
    def required_group_types(self) -> List[str]:
        return ["H_donor", "water", "H_acceptor"]

    @property
    def metric_names(self) -> List[str]:
        return ["dist_dw", "dist_wa", "theta", "omega"]

    # ==================== 抽象方法（基类要求，此处不使用） ====================

    def get_candidate_tuples(self, groups, coordinates=None):
        raise NotImplementedError("detect() is overridden, this method is not used")

    def compute_metrics_for_frame(self, tuples, all_positions, frame):
        raise NotImplementedError("detect() is overridden, this method is not used")

    def apply_threshold(self, metrics):
        dw_ok = (metrics["dist_dw"] > WATER_BRIDGE_MINDIST) & \
                (metrics["dist_dw"] < WATER_BRIDGE_MAXDIST)
        wa_ok = (metrics["dist_wa"] > WATER_BRIDGE_MINDIST) & \
                (metrics["dist_wa"] < WATER_BRIDGE_MAXDIST)
        theta_ok = metrics["theta"] >= WATER_BRIDGE_THETA_MIN
        omega_ok = (metrics["omega"] >= WATER_BRIDGE_OMEGA_MIN) & \
                   (metrics["omega"] <= WATER_BRIDGE_OMEGA_MAX)
        return dw_ok & wa_ok & theta_ok & omega_ok

    # ==================== 重写 detect ====================

    def detect(self, groups, trajectory=None,
               n_workers: int = 1,
               topology_path: str = None,
               trajectory_path: str = None,
               tuple_filter=None) -> List[Interaction]:
        """水桥检测主流程。"""
        if trajectory is None:
            raise ValueError("trajectory is required")
        if n_workers > 1:
            raise NotImplementedError("PerFrame 检测器暂不支持并行")

        # 1. 分组
        from ..group_identifiers.amber_ff_identifier import WATER_RESIDUES
        donors = [g for g in groups
                  if g.group_type == "H_donor" and g.residue_name not in WATER_RESIDUES]
        waters = [g for g in groups if g.group_type == "water"]
        acceptors = [g for g in groups
                     if g.group_type == "H_acceptor" and g.residue_name not in WATER_RESIDUES]

        if not donors or not waters or not acceptors:
            return []

        # 2. 构建索引数组
        self._donor_d_idx = np.array([g.atoms[0].atom_global_idx for g in donors])
        self._donor_h_idx = np.array([g.atoms[1].atom_global_idx for g in donors])
        self._water_ow_idx = np.array([g.atoms[0].atom_global_idx for g in waters])
        self._acceptor_a_idx = np.array([g.atoms[0].atom_global_idx for g in acceptors])

        # 3. 第一帧 KDTree 预筛选
        first_pos = trajectory[0].positions
        d_pos = first_pos[self._donor_d_idx]
        w_pos = first_pos[self._water_ow_idx]
        a_pos = first_pos[self._acceptor_a_idx]

        donor_tree = KDTree(d_pos)
        acceptor_tree = KDTree(a_pos)

        triple_d, triple_w, triple_a = [], [], []
        for wi in range(len(waters)):
            wp = w_pos[wi]
            nearby_d = donor_tree.query_ball_point(wp, PREFILTER_CUTOFF)
            nearby_a = acceptor_tree.query_ball_point(wp, PREFILTER_CUTOFF)
            if not nearby_d or not nearby_a:
                continue
            for di in nearby_d:
                for ai in nearby_a:
                    if np.linalg.norm(d_pos[di] - a_pos[ai]) < WATER_BRIDGE_MINDIST:
                        continue
                    triple_d.append(di)
                    triple_w.append(wi)
                    triple_a.append(ai)

        if not triple_d:
            return []

        triple_d = np.array(triple_d)
        triple_w = np.array(triple_w)
        triple_a = np.array(triple_a)

        # 4. tuple_filter
        if tuple_filter is not None:
            mask = np.array([
                tuple_filter((donors[d], waters[w], acceptors[a]))
                for d, w, a in zip(triple_d, triple_w, triple_a)
            ])
            triple_d = triple_d[mask]
            triple_w = triple_w[mask]
            triple_a = triple_a[mask]

        if len(triple_d) == 0:
            return []

        self._triple_d = triple_d
        self._triple_w = triple_w
        self._triple_a = triple_a

        n_triples = len(triple_d)
        n_frames = trajectory.n_frames

        # 5. 预分配结果数组
        existence = np.zeros((n_triples, n_frames), dtype=bool)
        dist_dw = np.zeros((n_triples, n_frames))
        dist_wa = np.zeros((n_triples, n_frames))
        theta = np.zeros((n_triples, n_frames))
        omega = np.zeros((n_triples, n_frames))

        # 6. 逐帧计算
        for f, ts in enumerate(trajectory):
            positions = ts.positions
            metrics = self._compute_metrics(positions)
            dist_dw[:, f] = metrics["dist_dw"]
            dist_wa[:, f] = metrics["dist_wa"]
            theta[:, f] = metrics["theta"]
            omega[:, f] = metrics["omega"]
            existence[:, f] = self.apply_threshold(metrics)

        # 7. 过滤从未存在的三元组
        has_any = np.any(existence, axis=1)
        if not np.any(has_any):
            return []

        tuples = [(donors[triple_d[i]], waters[triple_w[i]], acceptors[triple_a[i]])
                  for i in range(n_triples) if has_any[i]]

        return [Interaction(
            interaction_type=self.name,
            groups=tuples,
            existence=existence[has_any],
            metrics={
                "dist_dw": dist_dw[has_any],
                "dist_wa": dist_wa[has_any],
                "theta": theta[has_any],
                "omega": omega[has_any],
            }
        )]

    # ==================== 向量化计算 ====================

    def _compute_metrics(self, positions):
        """向量化计算全部候选三元组的指标。"""
        td = self._triple_d
        tw = self._triple_w
        ta = self._triple_a

        # 取坐标：(n_triples, 3)
        d = positions[self._donor_d_idx[td]]
        h = positions[self._donor_h_idx[td]]
        ow = positions[self._water_ow_idx[tw]]
        a = positions[self._acceptor_a_idx[ta]]

        # 距离
        dist_dw = np.linalg.norm(d - ow, axis=1)
        dist_wa = np.linalg.norm(ow - a, axis=1)

        # theta 角：水O-供体H-供体D
        vec_hd = d - h     # H→D
        vec_ho = ow - h    # H→Ow
        theta = self._vec_angle(vec_hd, vec_ho)

        # omega 角：受体-水O-供体H
        vec_oa = a - ow    # Ow→A
        vec_oh = h - ow    # Ow→H
        omega = self._vec_angle(vec_oa, vec_oh)

        return {
            "dist_dw": dist_dw,
            "dist_wa": dist_wa,
            "theta": theta,
            "omega": omega,
        }

    @staticmethod
    def _vec_angle(v1: np.ndarray, v2: np.ndarray) -> np.ndarray:
        """计算两组向量的夹角（度）。"""
        cos_angle = np.sum(v1 * v2, axis=1) / (
            np.linalg.norm(v1, axis=1) * np.linalg.norm(v2, axis=1)
        )
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        return np.degrees(np.arccos(cos_angle))
