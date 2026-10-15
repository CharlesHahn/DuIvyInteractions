# -*- coding: utf-8 -*-
"""疏水相互作用检测器（两轮遍历 + 稀疏存储）。

判据：两个疏水碳原子之间距离在 0.5~4.0 Å。
参考：PLIP。

Pass1：KDTree 自查询发现近邻对 + 距离过滤（PBC 感知）
Pass2：从 InteractionSparse 重建索引 + 向量化计算（PBC 感知）

PBC：距离计算传 box（最小镜像约定），见 doc/pbc_design.md。
"""

import numpy as np
from typing import List, Tuple, Dict

from ..core.interfaces import InteractionDetectorTwoPass
from ..core.datas import Group, Interaction, InteractionSparse
from ..core.pbc import get_box


# PLIP 阈值
HYDROPH_DIST_MAX = 4.0  # Å，疏水原子间最大距离
HYDROPH_MIN_DIST = 0.5  # Å，最小距离


class HydrophobicDetectorTwoPass(InteractionDetectorTwoPass):
    """疏水相互作用检测器（两轮遍历 + 稀疏存储）。"""

    @property
    def name(self) -> str:
        return "hydrophobic"

    @property
    def required_group_types(self) -> List[str]:
        return ["hydrophobic"]

    @property
    def metric_names(self) -> List[str]:
        return ["distance"]

    # ==================== Pass1：KDTree 发现 ====================

    def run_pass1(self, groups: List[Group], trajectory,
                  tuple_filter=None) -> InteractionSparse:
        """执行 Pass1：KDTree 自查询发现近邻对。"""
        hydro_groups = [g for g in groups if g.group_type == "hydrophobic"]

        if not hydro_groups:
            return InteractionSparse(interaction_type=self.name, data={})

        # 构建原子索引数组
        self._atom_indices = np.array(
            [g.atoms[0].atom_global_idx for g in hydro_groups])

        sparse_data: Dict[Tuple[int, ...], dict] = {}

        for f, ts in enumerate(trajectory):
            pairs, distances = self._kdtree_discover(ts.positions,
                                                     get_box(ts))

            for i, j, dist in zip(pairs[0], pairs[1], distances):
                g1 = hydro_groups[i]
                g2 = hydro_groups[j]

                if tuple_filter is not None and not tuple_filter((g1, g2)):
                    continue

                group_ids = (g1.group_id, g2.group_id)

                if group_ids not in sparse_data:
                    sparse_data[group_ids] = {
                        "groups": (g1, g2),
                        "frames": [],
                        "metrics": {"distance": []}
                    }

                sparse_data[group_ids]["frames"].append(f)
                sparse_data[group_ids]["metrics"]["distance"].append(dist)

        # TODO: 后处理去重 — 同一 (group_id, residue_id) 对只保留平均距离最近的 pair
        # 物理原因：残基 LEU100 有 3 个疏水碳 C1/C2/C3，配体有 CA，
        # CA-C1、CA-C2、CA-C3 都是 "CA-LEU100" 的疏水接触，应只保留最近的。
        # 实现：计算每对在 active 帧的平均距离，每 (group_id, residue_id) 只保留最优。

        return InteractionSparse(interaction_type=self.name, data=sparse_data)

    # ==================== Pass2：从 InteractionSparse 重建 ====================

    def run_pass2(self, sparse: InteractionSparse,
                  trajectory) -> List[Interaction]:
        """执行 Pass2：补全全帧数据（逐帧 PBC 距离）。

        自写循环（基类 run_pass2 不暴露帧/box）：逐帧取 get_box(ts) 传入
        compute_pair_metrics，使 PBC 距离逐帧正确（NPT 系综盒逐帧不同）。
        结果构建与基类 run_pass2 一致（保留全部 pair、_post_process + _build_interaction）。
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
        if not hasattr(self, '_pair_idx1'):
            raise RuntimeError(
                "compute_pair_metrics requires cached indices. "
                "Call _build_indices_from_sparse first via run_pass2.")
        from MDAnalysis.lib.distances import calc_bonds
        pos1 = all_positions[self._pair_idx1]
        pos2 = all_positions[self._pair_idx2]
        return {"distance": calc_bonds(pos1, pos2, box=box)}

    def apply_threshold(self, metrics):
        """0.5 < 距离 < 4.0。"""
        return (metrics["distance"] > HYDROPH_MIN_DIST) & \
               (metrics["distance"] < HYDROPH_DIST_MAX)

    # ==================== 内部辅助方法 ====================

    def _kdtree_discover(self, all_positions, box=None):
        """PBC 感知的自查询发现近邻对。返回 (pairs, distances)。

        box: 盒尺寸（MIC）；None = 欧氏距离（无 PBC）。
        """
        from MDAnalysis.lib.distances import self_capped_distance

        coords = all_positions[self._atom_indices]
        pairs, distances = self_capped_distance(
            coords, HYDROPH_DIST_MAX, min_cutoff=HYDROPH_MIN_DIST,
            box=box, method='pkdtree', return_distances=True)
        if len(pairs) == 0:
            return (np.array([]), np.array([])), np.array([])
        return (pairs[:, 0], pairs[:, 1]), distances

    def _build_indices_from_sparse(self, sparse: InteractionSparse):
        """从 InteractionSparse 重建原子索引数组（Pass2 用）。"""
        group_tuples = [entry["groups"] for entry in sparse.data.values()]
        self._pair_idx1 = np.array(
            [gt[0].atoms[0].atom_global_idx for gt in group_tuples])
        self._pair_idx2 = np.array(
            [gt[1].atoms[0].atom_global_idx for gt in group_tuples])
