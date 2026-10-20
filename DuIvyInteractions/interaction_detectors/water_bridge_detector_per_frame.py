# -*- coding: utf-8 -*-
"""水桥检测器（逐帧策略，KDTree 预筛选 + 向量化计算）。

判据：水分子位于供体和受体之间，距离和角度满足 PLIP 定义。
参考：PLIP, Jiang et al. 2005。

PBC：候选预筛 capped_distance(box)、逐帧距离 calc_bonds(box)、角度 calc_angles(box)，
见 doc/pbc_design.md §6B。
"""

import numpy as np
from typing import List, Tuple, Dict
from collections import defaultdict

from ..core.interfaces import InteractionDetectorPerFrame
from ..core.datas import Group, Interaction
from ..core.pbc import get_box


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
            raise NotImplementedError("PerFrame detector does not support parallelism yet")

        # 1. 分组
        # 分组：供体/受体排除水分子自身（按 water 组原子集合，跨力场通用，不依赖残基名）
        waters = [g for g in groups if g.group_type == "water"]
        water_atoms = set()
        for w in waters:
            water_atoms.update(w.atom_indices)
        donors = [g for g in groups
                  if g.group_type == "H_donor" and not (set(g.atom_indices) & water_atoms)]
        acceptors = [g for g in groups
                     if g.group_type == "H_acceptor" and not (set(g.atom_indices) & water_atoms)]

        if not donors or not waters or not acceptors:
            return []

        # 2. 构建索引数组
        self._donor_d_idx = np.array([g.atoms[0].atom_global_idx for g in donors])
        self._donor_h_idx = np.array([g.atoms[1].atom_global_idx for g in donors])
        self._water_ow_idx = np.array([g.atoms[0].atom_global_idx for g in waters])
        self._acceptor_a_idx = np.array([g.atoms[0].atom_global_idx for g in acceptors])

        # 3. 第一帧 KDTree 预筛选（PBC 感知；全对 + 按 w 分组组合三元组）
        from MDAnalysis.lib.distances import capped_distance, calc_bonds
        first_pos = trajectory[0].positions
        first_box = get_box(trajectory[0])
        d_pos = first_pos[self._donor_d_idx]
        w_pos = first_pos[self._water_ow_idx]
        a_pos = first_pos[self._acceptor_a_idx]

        # 全体 w↔d、w↔a 候选对（pairs[:,0]=w，两组对齐后按 w 分组）
        dw_pairs = capped_distance(
            w_pos, d_pos, PREFILTER_CUTOFF,
            box=first_box, method='pkdtree', return_distances=False)
        wa_pairs = capped_distance(
            w_pos, a_pos, PREFILTER_CUTOFF,
            box=first_box, method='pkdtree', return_distances=False)

        if len(dw_pairs) == 0 or len(wa_pairs) == 0:
            return []

        dw_by_w = defaultdict(list)
        for idx in range(len(dw_pairs)):
            dw_by_w[dw_pairs[idx, 0]].append(dw_pairs[idx, 1])
        wa_by_w = defaultdict(list)
        for idx in range(len(wa_pairs)):
            wa_by_w[wa_pairs[idx, 0]].append(wa_pairs[idx, 1])

        common_w = set(dw_by_w) & set(wa_by_w)
        if not common_w:
            return []

        triple_d, triple_w, triple_a = [], [], []
        for wi in common_w:
            for di in dw_by_w[wi]:
                for ai in wa_by_w[wi]:
                    # D-A 最短距排除（MIC）
                    if calc_bonds(d_pos[di:di+1], a_pos[ai:ai+1],
                                  box=first_box)[0] < WATER_BRIDGE_MINDIST:
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
        times = np.empty(n_frames)

        # 6. 逐帧计算（PBC 感知）
        for f, ts in enumerate(trajectory):
            times[f] = ts.time
            positions = ts.positions
            metrics = self._compute_metrics(positions, get_box(ts))
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
            },
            times=times,
        )]

    # ==================== 向量化计算 ====================

    def _compute_metrics(self, positions, box=None):
        """向量化计算全部候选三元组的指标（MIC）。

        box: 盒尺寸（MIC）；None = 欧氏距离（无 PBC）。
        """
        from MDAnalysis.lib.distances import calc_bonds, calc_angles

        td = self._triple_d
        tw = self._triple_w
        ta = self._triple_a

        # 取坐标：(n_triples, 3)
        d = positions[self._donor_d_idx[td]]
        h = positions[self._donor_h_idx[td]]
        ow = positions[self._water_ow_idx[tw]]
        a = positions[self._acceptor_a_idx[ta]]

        # 距离（MIC）
        dist_dw = calc_bonds(d, ow, box=box)
        dist_wa = calc_bonds(ow, a, box=box)

        # theta 角：水O-供体H-供体D（顶点 h；calc_angles 返回弧度）
        theta = np.degrees(calc_angles(ow, h, d, box=box))

        # omega 角：受体-水O-供体H（顶点 ow）
        omega = np.degrees(calc_angles(a, ow, h, box=box))

        return {
            "dist_dw": dist_dw,
            "dist_wa": dist_wa,
            "theta": theta,
            "omega": omega,
        }
