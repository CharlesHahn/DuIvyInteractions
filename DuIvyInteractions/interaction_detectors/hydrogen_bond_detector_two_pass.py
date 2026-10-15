# -*- coding: utf-8 -*-
"""氢键检测器（两轮遍历 + 稀疏存储）。

判据：D-A 距离 ≤ 4.1 Å 且 D-H···A 角度 ≥ 100°。
参考：PLIP, Hubbard & Haider 2001（距离 +0.6 Å，角度 +10°）。

Pass1：KDTree 筛选 D-A 距离 + 精确计算角度（PBC 感知）
Pass2：从 InteractionSparse 重新构建索引 + 向量化计算（PBC 感知）

PBC：距离 calc_bonds(box)、角度 calc_angles(box)、候选发现 capped_distance(box)，
见 doc/pbc_design.md §6A。
"""

import numpy as np
from typing import List, Tuple, Dict

from ..core.interfaces import InteractionDetectorTwoPass
from ..core.datas import Group, Interaction, InteractionSparse
from ..core.pbc import get_box


# PLIP 阈值
HBOND_DIST_MAX = 4.1         # Å，D-A 最大距离
HBOND_DON_ANGLE_MIN = 100.0  # °，D-H···A 最小角度


class HydrogenBondDetectorTwoPass(InteractionDetectorTwoPass):
    """氢键检测器（两轮遍历 + 稀疏存储）。

    Pass1：覆盖 run_pass1，KDTree 筛距离 + 精确计算角度。
    Pass2：从 InteractionSparse 重新构建索引 + 向量化计算。
    """

    @property
    def name(self) -> str:
        return "hydrogen_bond"

    @property
    def required_group_types(self) -> List[str]:
        return ["H_donor", "H_acceptor"]

    @property
    def metric_names(self) -> List[str]:
        return ["distance", "angle"]

    # ==================== Pass1：KDTree 发现 ====================

    def run_pass1(self, groups: List[Group], trajectory,
                  tuple_filter=None) -> InteractionSparse:
        """执行 Pass1：KDTree 筛距离 + 精确计算角度。"""
        # 1. 分组
        donors = [g for g in groups if g.group_type == "H_donor"]
        acceptors = [g for g in groups if g.group_type == "H_acceptor"]

        if not donors or not acceptors:
            return InteractionSparse(interaction_type=self.name, data={})

        # 2. 构建索引数组
        self._build_indices(donors, acceptors)

        # 3. 逐帧 KDTree 发现
        sparse_data: Dict[Tuple[int, ...], dict] = {}

        for f, ts in enumerate(trajectory):
            active_donors, active_acceptors, distances, angles = \
                self._kdtree_discover(ts.positions, get_box(ts))

            for d_idx, a_idx, dist, angle in zip(
                    active_donors, active_acceptors, distances, angles):
                donor = donors[d_idx]
                acceptor = acceptors[a_idx]

                # tuple_filter
                if tuple_filter is not None and not tuple_filter((donor, acceptor)):
                    continue

                group_ids = (donor.group_id, acceptor.group_id)

                if group_ids not in sparse_data:
                    sparse_data[group_ids] = {
                        "groups": (donor, acceptor),
                        "frames": [],
                        "metrics": {"distance": [], "angle": []}
                    }

                sparse_data[group_ids]["frames"].append(f)
                sparse_data[group_ids]["metrics"]["distance"].append(dist)
                sparse_data[group_ids]["metrics"]["angle"].append(angle)

        return InteractionSparse(interaction_type=self.name, data=sparse_data)

    # ==================== Pass2：从 InteractionSparse 重建 ====================

    def run_pass2(self, sparse: InteractionSparse,
                  trajectory) -> List[Interaction]:
        """执行 Pass2：补全全帧数据（逐帧 PBC 距离/角度）。

        自写循环（基类 run_pass2 不暴露帧/box）：逐帧取 get_box(ts) 传入
        compute_pair_metrics，使 PBC 距离/角度逐帧正确（NPT 系综盒逐帧不同）。
        """
        if not sparse.data:
            return []

        # 从 InteractionSparse 重新构建索引数组
        self._build_indices_from_sparse(sparse)

        group_tuples = [entry["groups"] for entry in sparse.data.values()]
        n_groups = len(group_tuples)
        n_frames = trajectory.n_frames

        # 预分配稠密数组
        existence = np.zeros((n_groups, n_frames), dtype=bool)
        metrics = {name: np.full((n_groups, n_frames), np.nan)
                   for name in self.metric_names}

        times = np.empty(n_frames)
        for f, ts in enumerate(trajectory):
            times[f] = ts.time
            frame_metrics = self.compute_pair_metrics(
                group_tuples, ts.positions, box=get_box(ts))
            existence[:, f] = self.apply_threshold(frame_metrics)
            for name in self.metric_names:
                metrics[name][:, f] = frame_metrics[name]

        # 构建 results 列表（与基类一致：不过滤从未存在的 pair）
        results = [
            (entry["groups"], existence[i],
             {k: v[i] for k, v in metrics.items()})
            for i, entry in enumerate(sparse.data.values())
        ]

        results = self._post_process(results)
        return self._build_interaction(results, times=times)

    # ==================== compute_pair_metrics ====================

    def compute_pair_metrics(self, group_tuples, all_positions, box=None):
        """向量化计算距离和角度（使用 _build_indices_from_sparse 缓存的索引）。

        box: 盒尺寸（MIC）；None = 欧氏距离（无 PBC）。
        """
        if not hasattr(self, '_donor_d_idx'):
            raise RuntimeError(
                "compute_pair_metrics requires cached indices. "
                "Call _build_indices_from_sparse first via run_pass2.")

        from MDAnalysis.lib.distances import calc_bonds, calc_angles

        d = all_positions[self._donor_d_idx]
        h = all_positions[self._donor_h_idx]
        a = all_positions[self._acceptor_a_idx]

        # D-A 距离（MIC）
        dist = calc_bonds(d, a, box=box)

        # D-H···A 角度（MIC，顶点 h；calc_angles 返回弧度）
        ang = np.degrees(calc_angles(d, h, a, box=box))

        return {"distance": dist, "angle": ang}

    def apply_threshold(self, metrics):
        """距离 ≤ 4.1Å 且角度 ≥ 100°。"""
        return (metrics["distance"] <= HBOND_DIST_MAX) & \
               (metrics["angle"] >= HBOND_DON_ANGLE_MIN)

    # ==================== 内部辅助方法 ====================

    def _build_indices(self, donors, acceptors):
        """构建索引数组（Pass1 用）。"""
        self._donor_d_idx = np.array([g.atoms[0].atom_global_idx for g in donors])
        self._donor_h_idx = np.array([g.atoms[1].atom_global_idx for g in donors])
        self._acceptor_a_idx = np.array([g.atoms[0].atom_global_idx for g in acceptors])

    def _build_indices_from_sparse(self, sparse: InteractionSparse):
        """从 InteractionSparse 重新构建索引数组（Pass2 用）。"""
        group_tuples = [entry["groups"] for entry in sparse.data.values()]
        donors = [gt[0] for gt in group_tuples]
        acceptors = [gt[1] for gt in group_tuples]

        self._donor_d_idx = np.array([g.atoms[0].atom_global_idx for g in donors])
        self._donor_h_idx = np.array([g.atoms[1].atom_global_idx for g in donors])
        self._acceptor_a_idx = np.array([g.atoms[0].atom_global_idx for g in acceptors])

    def _kdtree_discover(self, all_positions, box=None):
        """PBC 感知发现 D-A 距离 ≤ 4.1Å 的 pair + 精确计算角度。

        box: 盒尺寸（MIC）；None = 欧氏距离（无 PBC）。
        """
        from MDAnalysis.lib.distances import capped_distance, calc_angles

        # 获取 D 和 A 坐标
        d_coords = all_positions[self._donor_d_idx]
        a_coords = all_positions[self._acceptor_a_idx]

        # 交叉查询（MIC；capped_distance 返回 (pairs[:,0]=d, pairs[:,1]=a, dists)）
        pairs, distances = capped_distance(
            d_coords, a_coords, HBOND_DIST_MAX,
            box=box, method='pkdtree', return_distances=True)

        if len(pairs) == 0:
            return np.array([]), np.array([]), np.array([]), np.array([])

        active_donors = pairs[:, 0]
        active_acceptors = pairs[:, 1]

        # 计算角度（MIC，顶点 h；calc_angles 返回弧度）
        angles = np.degrees(calc_angles(
            d_coords[active_donors],
            all_positions[self._donor_h_idx[active_donors]],
            a_coords[active_acceptors], box=box))

        # 过滤角度 ≥ 100°
        angle_mask = angles >= HBOND_DON_ANGLE_MIN

        return (active_donors[angle_mask], active_acceptors[angle_mask],
                distances[angle_mask], angles[angle_mask])
