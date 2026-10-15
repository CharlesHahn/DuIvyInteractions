# -*- coding: utf-8 -*-
"""金属配位检测器（两轮遍历 + 稀疏存储）。

判据：金属离子到配位原子距离 < 3.0 Å。
参考：PLIP, Harding 2001。

Pass1：KDTree 交叉查询发现近邻对（PBC 感知）
Pass2：从 InteractionSparse 重建索引 + 向量化计算（PBC 感知）

PBC：距离计算传 box（最小镜像约定），见 doc/pbc_design.md。
"""

import numpy as np
from typing import List, Tuple, Dict

from ..core.interfaces import InteractionDetectorTwoPass
from ..core.datas import Group, Interaction, InteractionSparse
from ..core.pbc import get_box


# PLIP 阈值
METAL_DIST_MAX = 3.0  # Å


class MetalCoordinationDetectorTwoPass(InteractionDetectorTwoPass):
    """金属配位检测器（两轮遍历 + 稀疏存储）。"""

    @property
    def name(self) -> str:
        return "metal_coordination"

    @property
    def required_group_types(self) -> List[str]:
        return ["metal", "metal_binding"]

    @property
    def metric_names(self) -> List[str]:
        return ["distance"]

    # ==================== Pass1：KDTree 发现 ====================

    def run_pass1(self, groups: List[Group], trajectory,
                  tuple_filter=None) -> InteractionSparse:
        """执行 Pass1：KDTree 交叉查询发现金属-配位近邻对。"""
        metals = [g for g in groups if g.group_type == "metal"]
        bindings = [g for g in groups if g.group_type == "metal_binding"]

        if not metals or not bindings:
            return InteractionSparse(interaction_type=self.name, data={})

        # 构建原子索引数组
        self._metal_indices = np.array(
            [g.atoms[0].atom_global_idx for g in metals])
        self._binding_indices = np.array(
            [g.atoms[0].atom_global_idx for g in bindings])

        sparse_data: Dict[Tuple[int, ...], dict] = {}

        for f, ts in enumerate(trajectory):
            active_metals, active_bindings, distances = \
                self._kdtree_discover(ts.positions, get_box(ts))

            for m_idx, b_idx, dist in zip(
                    active_metals, active_bindings, distances):
                metal = metals[m_idx]
                binding = bindings[b_idx]

                if tuple_filter is not None and not tuple_filter((metal, binding)):
                    continue

                group_ids = (metal.group_id, binding.group_id)

                if group_ids not in sparse_data:
                    sparse_data[group_ids] = {
                        "groups": (metal, binding),
                        "frames": [],
                        "metrics": {"distance": []}
                    }

                sparse_data[group_ids]["frames"].append(f)
                sparse_data[group_ids]["metrics"]["distance"].append(dist)

        return InteractionSparse(interaction_type=self.name, data=sparse_data)

    # ==================== Pass2：从 InteractionSparse 重建 ====================

    def run_pass2(self, sparse: InteractionSparse,
                  trajectory) -> List[Interaction]:
        """执行 Pass2：补全全帧数据（逐帧 PBC 距离）。

        自写循环（基类 run_pass2 不暴露帧/box）：逐帧取 get_box(ts) 传入
        compute_pair_metrics，使 PBC 距离逐帧正确（NPT 系综盒逐帧不同）。
        """
        if not sparse.data:
            return []

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
        """向量化计算距离（使用 _build_indices_from_sparse 缓存的索引）。

        box: 盒尺寸（MIC）；None = 欧氏距离（无 PBC）。
        """
        if not hasattr(self, '_pass2_metal_idx'):
            raise RuntimeError(
                "compute_pair_metrics requires cached indices. "
                "Call _build_indices_from_sparse first via run_pass2.")
        from MDAnalysis.lib.distances import calc_bonds
        m_pos = all_positions[self._pass2_metal_idx]
        b_pos = all_positions[self._pass2_binding_idx]
        return {"distance": calc_bonds(m_pos, b_pos, box=box)}

    def apply_threshold(self, metrics):
        """距离 < METAL_DIST_MAX。"""
        return metrics["distance"] < METAL_DIST_MAX

    # ==================== 内部辅助方法 ====================

    def _kdtree_discover(self, all_positions, box=None):
        """PBC 感知的交叉查询发现金属-配位近邻对。

        box: 盒尺寸（MIC）；None = 欧氏距离（无 PBC）。
        """
        from MDAnalysis.lib.distances import capped_distance

        m_coords = all_positions[self._metal_indices]
        b_coords = all_positions[self._binding_indices]

        pairs, distances = capped_distance(
            m_coords, b_coords, METAL_DIST_MAX,
            box=box, method='pkdtree', return_distances=True)

        # 空对：capped_distance 返回 pairs.shape=(0,2)（实测）
        if len(pairs) == 0:
            return np.array([]), np.array([]), np.array([])
        return pairs[:, 0], pairs[:, 1], distances

    def _build_indices_from_sparse(self, sparse: InteractionSparse):
        """从 InteractionSparse 重建原子索引数组（Pass2 用）。"""
        group_tuples = [entry["groups"] for entry in sparse.data.values()]
        self._pass2_metal_idx = np.array(
            [gt[0].atoms[0].atom_global_idx for gt in group_tuples])
        self._pass2_binding_idx = np.array(
            [gt[1].atoms[0].atom_global_idx for gt in group_tuples])
