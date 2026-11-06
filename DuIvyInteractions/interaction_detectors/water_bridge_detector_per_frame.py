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


# 氢键几何判据（Baker & Hubbard 1984；设计见 doc/water_bridge_design.md）
WATER_BRIDGE_HB_DIST_MAX = 2.5      # Å，H···A 距离上限
WATER_BRIDGE_HB_ANGLE_MIN = 120.0   # °，D-H···A 角度下限
WATER_BRIDGE_MIN_DIST = 0.5         # Å，防原子重叠下界
WATER_BRIDGE_MAX_DONORS = 2         # 每水每帧最多供体数（2 对孤对）
WATER_BRIDGE_MAX_BRIDGES = 2        # 每水每帧最多桥数（2 个 H）
WATER_BRIDGE_PREFILTER = 3.0        # Å，候选预筛半径（H···O / H_w···A）
WATER_BRIDGE_OMEGA_MIN = 71.0       # °，A-Ow···H(供体) 最小角（"水在中间"约束）
WATER_BRIDGE_OMEGA_MAX = 140.0      # °，A-Ow···H(供体) 最大角
PREFILTER_CUTOFF = WATER_BRIDGE_PREFILTER  # 兼容旧名


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
        return ["dist_h_ow", "angle_dh_ow", "dist_d_ow",
                "dist_hw_a", "angle_ow_hw_a", "dist_ow_a", "angle_h_ow_a"]

    # ==================== 抽象方法（基类要求，此处不使用） ====================

    def get_candidate_tuples(self, groups, coordinates=None):
        raise NotImplementedError("detect() is overridden, this method is not used")

    def compute_metrics_for_frame(self, tuples, all_positions, frame):
        raise NotImplementedError("detect() is overridden, this method is not used")

    def apply_threshold(self, metrics):
        """H···A ≤ 2.5 Å 且 D-H···A ≥ 120°（两侧各一，见设计文档 §2.2）。"""
        dw_ok = (metrics["dist_h_ow"] > WATER_BRIDGE_MIN_DIST) & \
                (metrics["dist_h_ow"] <= WATER_BRIDGE_HB_DIST_MAX)
        wa_ok = (metrics["dist_hw_a"] > WATER_BRIDGE_MIN_DIST) & \
                (metrics["dist_hw_a"] <= WATER_BRIDGE_HB_DIST_MAX)
        theta_ok = metrics["angle_dh_ow"] >= WATER_BRIDGE_HB_ANGLE_MIN
        omega_ok = metrics["angle_ow_hw_a"] >= WATER_BRIDGE_HB_ANGLE_MIN
        # "水在中间"约束（PLIP omega 等价）：供体 H 与受体 A 在水 O 两侧
        mid_ok = (metrics["angle_h_ow_a"] >= WATER_BRIDGE_OMEGA_MIN) & \
                 (metrics["angle_h_ow_a"] <= WATER_BRIDGE_OMEGA_MAX)
        return dw_ok & wa_ok & theta_ok & omega_ok & mid_ok

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

        # 退化检查：判据依赖水的显式 H（见设计文档 §6）
        if any(sum(1 for a in w.atoms if a.atom_element == "H") < 2 for w in waters):
            print("[WARN] water_bridge: water has no explicit H, skipped")
            return []

        # 2. 构建索引数组
        self._donor_d_idx = np.array([g.atoms[0].atom_global_idx for g in donors])
        self._donor_h_idx = np.array([g.atoms[1].atom_global_idx for g in donors])
        self._water_ow_idx = np.array([g.atoms[0].atom_global_idx for g in waters])
        self._acceptor_a_idx = np.array([g.atoms[0].atom_global_idx for g in acceptors])
        # 水的两个 H（每水 2 个，形状 (n_water, 2)）
        self._water_hw_idx = np.array(
            [[w.atoms[1].atom_global_idx, w.atoms[2].atom_global_idx] for w in waters])

        # 3. 第一帧 KDTree 预筛选（PBC 感知；基于显式 H 的距离）
        #    设计原则：PerFrame 策略的候选集只基于第一帧（后续帧复用该候选集）
        from MDAnalysis.lib.distances import capped_distance
        first_pos = trajectory[0].positions
        first_box = get_box(trajectory[0])
        h_pos = first_pos[self._donor_h_idx]
        w_pos = first_pos[self._water_ow_idx]
        hw_pos = first_pos[self._water_hw_idx.ravel()]
        a_pos = first_pos[self._acceptor_a_idx]

        ho_pairs = capped_distance(
            h_pos, w_pos, PREFILTER_CUTOFF,
            box=first_box, method='pkdtree', return_distances=False)
        hwa_pairs = capped_distance(
            hw_pos, a_pos, PREFILTER_CUTOFF,
            box=first_box, method='pkdtree', return_distances=False)

        if len(ho_pairs) == 0 or len(hwa_pairs) == 0:
            return []

        dw_by_w = defaultdict(set)
        for k in range(len(ho_pairs)):
            dw_by_w[ho_pairs[k, 1]].add(ho_pairs[k, 0])
        wa_by_w = defaultdict(set)
        for k in range(len(hwa_pairs)):
            wa_by_w[hwa_pairs[k, 0] // 2].add(hwa_pairs[k, 1])

        common_w = set(dw_by_w) & set(wa_by_w)
        if not common_w:
            return []

        triple_d, triple_w, triple_a = [], [], []
        for wi in common_w:
            for di in dw_by_w[wi]:
                for ai in wa_by_w[wi]:
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
        metrics_out = {name: np.zeros((n_triples, n_frames))
                       for name in self.metric_names}
        times = np.empty(n_frames)

        # 每对的水/供体身份（配额用）
        water_uid = np.array([waters[wi].group_id for wi in triple_w])
        donor_uid = np.array([donors[di].group_id for di in triple_d])

        # 6. 逐帧计算（PBC 感知）+ 配额
        for f, ts in enumerate(trajectory):
            times[f] = ts.time
            positions = ts.positions
            frame_metrics = self._compute_metrics(positions, get_box(ts))
            ok = self.apply_threshold(frame_metrics)
            ok = self._apply_water_quota(
                ok, frame_metrics, water_uid, donor_uid)
            existence[:, f] = ok
            for name in self.metric_names:
                metrics_out[name][:, f] = frame_metrics[name]

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
            metrics={name: arr[has_any] for name, arr in metrics_out.items()},
            times=times,
        )]

    # ==================== 配额（去重） ====================

    @staticmethod
    def _apply_water_quota(ok, frame_metrics, water_uid, donor_uid):
        """三层配额（逐帧）：每水H≤1受体、每水O≤2供体、每水≤2桥。

        设计见 doc/water_bridge_design.md §3.1。

        第 3 层的候选为「受体侧保留的 (H_w, A)」×「供体侧保留的 D」组合
        （k×m，最多 4），按两侧角度之和择优取前 2——这样"同一供体 + 水的两个 H
        指向两个受体"可保留为 2 个桥。
        """
        idx = np.where(ok)[0]
        if len(idx) == 0:
            return ok
        ang = frame_metrics["angle_ow_hw_a"]
        ang_d = frame_metrics["angle_dh_ow"]
        dh = frame_metrics["dist_h_ow"]
        hw = frame_metrics["hw_id"]
        keep = np.zeros_like(ok)
        for w in np.unique(water_uid[idx]):
            cand = idx[water_uid[idx] == w]
            # 第 1 层：每 (水, H_w) 保留角度最优的 1 个受体
            by_h = {}
            for i in cand:
                h = hw[i]
                if h not in by_h or ang[i] > ang[by_h[h]]:
                    by_h[h] = i
            # 第 2 层：每水 O 保留 dist_h_ow 最小的 2 个供体（按供体去重）
            by_d = {}
            for i in cand:
                d = donor_uid[i]
                if d not in by_d or dh[i] < dh[by_d[d]]:
                    by_d[d] = i
            donors_kept = sorted(by_d.values(), key=lambda i: dh[i])[:WATER_BRIDGE_MAX_DONORS]
            # 第 3 层：受体侧候选 × 供体侧候选 → 按两侧角度和择优取前 2
            combos = [(i, j) for i in by_h.values() for j in donors_kept]
            combos.sort(key=lambda p: -(ang[p[0]] + ang_d[p[1]]))
            for i, _ in combos[:WATER_BRIDGE_MAX_BRIDGES]:
                keep[i] = True
        return keep

    # ==================== 向量化计算 ====================

    def _compute_metrics(self, positions, box=None):
        """向量化计算全部候选三元组的指标（MIC）。

        受体侧对水的 HW1/HW2 各算一次，取 angle_ow_hw_a 更优者并记录 hw_id。
        box: 盒尺寸（MIC）；None = 欧氏距离（无 PBC）。
        """
        from MDAnalysis.lib.distances import calc_bonds, calc_angles

        td = self._triple_d
        tw = self._triple_w
        ta = self._triple_a

        d = positions[self._donor_d_idx[td]]
        h = positions[self._donor_h_idx[td]]
        ow = positions[self._water_ow_idx[tw]]
        a = positions[self._acceptor_a_idx[ta]]

        # 供体侧
        dist_h_ow = calc_bonds(h, ow, box=box)
        dist_d_ow = calc_bonds(d, ow, box=box)
        angle_dh_ow = np.degrees(calc_angles(d, h, ow, box=box))

        # 受体侧：HW1 / HW2 各算一次，取角度更优者
        hw1 = positions[self._water_hw_idx[tw, 0]]
        hw2 = positions[self._water_hw_idx[tw, 1]]
        dist_hw1 = calc_bonds(hw1, a, box=box)
        dist_hw2 = calc_bonds(hw2, a, box=box)
        ang1 = np.degrees(calc_angles(ow, hw1, a, box=box))
        ang2 = np.degrees(calc_angles(ow, hw2, a, box=box))
        use2 = ang2 > ang1
        dist_hw_a = np.where(use2, dist_hw2, dist_hw1)
        angle_ow_hw_a = np.where(use2, ang2, ang1)
        dist_ow_a = calc_bonds(ow, a, box=box)

        return {
            "dist_h_ow": dist_h_ow,
            "angle_dh_ow": angle_dh_ow,
            "dist_d_ow": dist_d_ow,
            "dist_hw_a": dist_hw_a,
            "angle_ow_hw_a": angle_ow_hw_a,
            "dist_ow_a": dist_ow_a,
            "hw_id": use2.astype(float),
            "angle_h_ow_a": np.degrees(calc_angles(h, ow, a, box=box)),
        }
