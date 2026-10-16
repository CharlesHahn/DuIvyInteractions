# -*- coding: utf-8 -*-
"""π-阳离子相互作用检测器（两轮遍历 + 稀疏存储）。

判据：电荷中心到环心距离 < 6.0 Å，投影偏移 < 2.0 Å。
参考：PLIP, Gallivan and Dougherty 1999。

PBC：环原子回绕（minimize_vectors 展平-回绕-还原）→ 环心/法向量/电荷中心；
距离 calc_bonds(box)；投影偏移 t 向量做 MIC，见 doc/pbc_design.md §6E。
"""

import numpy as np
from typing import List, Tuple, Dict

from ..core.interfaces import InteractionDetectorTwoPass
from ..core.datas import Group, Interaction, InteractionSparse
from ..core.pbc import get_box


# PLIP 阈值
PICATION_DIST_MAX = 6.0    # Å
PICATION_OFFSET_MAX = 2.0  # Å
PICATION_MIN_DIST = 0.5    # Å


class PiCationDetectorTwoPass(InteractionDetectorTwoPass):
    """π-阳离子相互作用检测器（两轮遍历 + 稀疏存储）。"""

    @property
    def name(self) -> str:
        return "pi_cation"

    @property
    def required_group_types(self) -> List[str]:
        return ["aromatic_ring", "charged_positive"]

    @property
    def metric_names(self) -> List[str]:
        return ["distance", "offset"]

    # ==================== initialize_candidates ====================

    def initialize_candidates(self, groups, trajectory, tuple_filter=None):
        """分组 → padding → 笛卡尔积 → tuple_filter → 设置缓存索引。"""
        rings = [g for g in groups if g.group_type == "aromatic_ring"]
        cations = [g for g in groups if g.group_type == "charged_positive"]

        if not rings or not cations:
            return []

        # 构建索引
        ring_idx, ring_prev, ring_next, ring_valid, _ = \
            self._build_circular_padding(rings)
        cation_idx, cation_q, cation_valid = self._build_padding(cations)

        # 笛卡尔积
        n_r = len(rings)
        n_c = len(cations)
        r_grid, c_grid = np.meshgrid(np.arange(n_r), np.arange(n_c), indexing='ij')
        pair_ring = r_grid.ravel()
        pair_cation = c_grid.ravel()

        # tuple_filter
        if tuple_filter is not None:
            mask = np.array([
                tuple_filter((rings[ri], cations[ci]))
                for ri, ci in zip(pair_ring, pair_cation)
            ])
            pair_ring = pair_ring[mask]
            pair_cation = pair_cation[mask]

        if len(pair_ring) == 0:
            return []

        # 设置缓存索引
        self._cached_ring_idx = ring_idx[pair_ring]
        self._cached_ring_prev = ring_prev[pair_ring]
        self._cached_ring_next = ring_next[pair_ring]
        self._cached_ring_valid = ring_valid[pair_ring]
        self._cached_cation_idx = cation_idx[pair_cation]
        self._cached_cation_q = cation_q[pair_cation]
        self._cached_cation_valid = cation_valid[pair_cation]

        return [(rings[pair_ring[i]], cations[pair_cation[i]])
                for i in range(len(pair_ring))]

    # ==================== compute_pair_metrics ====================

    def compute_pair_metrics(self, group_tuples, all_positions, box=None):
        """向量化计算距离和投影偏移（使用缓存索引，PBC 感知）。

        box: 盒尺寸（MIC）；None = 欧氏距离（无 PBC）。
        """
        if not hasattr(self, '_cached_ring_idx'):
            raise RuntimeError(
                "compute_pair_metrics requires cached indices. "
                "Call initialize_candidates or _build_indices_from_sparse first.")

        from MDAnalysis.lib.distances import minimize_vectors, calc_bonds

        # 环坐标（逐对已按 pair 展开）
        r_coords = all_positions[self._cached_ring_idx]
        r_prev = all_positions[self._cached_ring_prev]
        r_next = all_positions[self._cached_ring_next]
        valid = self._cached_ring_valid

        if box is not None:
            n_pairs, ma = r_coords.shape[:2]
            # ① 回绕环原子（以每环第一个原子为参考）
            ref = r_coords[:, 0:1, :]
            rel = r_coords - ref
            rel_mic = minimize_vectors(rel.reshape(-1, 3), box).reshape(n_pairs, ma, 3)
            r_coords_w = ref + rel_mic
            # prev/next 是同一环的原子（环形邻居），用同一参考回绕
            rel_p = r_prev - ref
            rel_n = r_next - ref
            r_prev_w = ref + minimize_vectors(rel_p.reshape(-1, 3), box).reshape(n_pairs, ma, 3)
            r_next_w = ref + minimize_vectors(rel_n.reshape(-1, 3), box).reshape(n_pairs, ma, 3)
            r_coords, r_prev, r_next = r_coords_w, r_prev_w, r_next_w

        # ② 环心（回绕后原子平均）
        r_masked = r_coords * valid[:, :, None]
        r_center = np.sum(r_masked, axis=1) / np.sum(valid, axis=1)[:, None]

        # ③ 环法向量（回绕后的相对向量 cross）
        atom_normals = np.cross(r_coords - r_prev, r_coords - r_next)
        n_masked = atom_normals * valid[:, :, None]
        normal = np.sum(n_masked, axis=1) / np.sum(valid, axis=1)[:, None]
        norm = np.linalg.norm(normal, axis=1, keepdims=True)
        normal = normal / norm

        # ④ 电荷中心（回绕）
        c_coords = all_positions[self._cached_cation_idx]
        if box is not None:
            n_pairs, ma_c = c_coords.shape[:2]
            ref_c = c_coords[:, 0:1, :]
            rel_c = c_coords - ref_c
            rel_c_mic = minimize_vectors(
                rel_c.reshape(-1, 3), box).reshape(n_pairs, ma_c, 3)
            c_coords = ref_c + rel_c_mic
        c_weighted = c_coords * self._cached_cation_q[:, :, None] * \
                     self._cached_cation_valid[:, :, None]
        q_sum = np.sum(self._cached_cation_q * self._cached_cation_valid, axis=1)
        c_center = np.sum(c_weighted, axis=1) / q_sum[:, None]

        # ⑤ 距离
        from MDAnalysis.lib.distances import calc_bonds
        distance = calc_bonds(r_center, c_center, box=box)

        # ⑥ 投影偏移（t 向量做 MIC）
        offset = self._projection_distance(normal, r_center, c_center, box=box)

        return {"distance": distance, "offset": offset}

    def apply_threshold(self, metrics):
        """0.5 < distance < 6.0 & offset < 2.0。"""
        return (metrics["distance"] > PICATION_MIN_DIST) & \
               (metrics["distance"] < PICATION_DIST_MAX) & \
               (metrics["offset"] < PICATION_OFFSET_MAX)

    # ==================== Pass2 ====================

    def run_pass2(self, sparse: InteractionSparse,
                  trajectory) -> List[Interaction]:
        """执行 Pass2：补全全帧数据（逐帧 PBC 距离/投影）。

        自写循环（基类 run_pass2 不暴露帧/box）：逐帧取 get_box(ts) 传入
        compute_pair_metrics，使 PBC 距离/投影逐帧正确（NPT 系综盒逐帧不同）。
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

    # ==================== 内部辅助方法 ====================

    def _build_indices_from_sparse(self, sparse: InteractionSparse):
        """从 InteractionSparse 重建缓存索引（Pass2 用）。"""
        group_tuples = [entry["groups"] for entry in sparse.data.values()]
        rings = [gt[0] for gt in group_tuples]
        cations = [gt[1] for gt in group_tuples]

        ring_idx, ring_prev, ring_next, ring_valid, _ = \
            self._build_circular_padding(rings)
        cation_idx, cation_q, cation_valid = self._build_padding(cations)

        self._cached_ring_idx = ring_idx
        self._cached_ring_prev = ring_prev
        self._cached_ring_next = ring_next
        self._cached_ring_valid = ring_valid
        self._cached_cation_idx = cation_idx
        self._cached_cation_q = cation_q
        self._cached_cation_valid = cation_valid

    @staticmethod
    def _build_circular_padding(rings: List[Group]):
        """构建环形 padding + 显式邻居索引。"""
        n = len(rings)
        max_atoms = max(len(r.atoms) for r in rings)
        indices = np.zeros((n, max_atoms), dtype=int)
        prev_idx = np.zeros((n, max_atoms), dtype=int)
        next_idx = np.zeros((n, max_atoms), dtype=int)
        valid = np.zeros((n, max_atoms), dtype=bool)
        n_atoms = np.zeros(n, dtype=int)

        for i, ring in enumerate(rings):
            na = len(ring.atoms)
            n_atoms[i] = na
            atom_indices = ring.atom_indices
            for j in range(max_atoms):
                indices[i, j] = atom_indices[j % na]
            for j in range(max_atoms):
                prev_idx[i, j] = atom_indices[(j - 1) % na]
                next_idx[i, j] = atom_indices[(j + 1) % na]
            valid[i, :na] = True

        return indices, prev_idx, next_idx, valid, n_atoms

    @staticmethod
    def _build_padding(groups: List[Group]):
        """构建 padding 索引矩阵。"""
        n = len(groups)
        max_atoms = max(len(g.atoms) for g in groups)
        indices = np.zeros((n, max_atoms), dtype=int)
        charges = np.zeros((n, max_atoms))
        valid = np.zeros((n, max_atoms), dtype=bool)
        for i, g in enumerate(groups):
            na = len(g.atoms)
            indices[i, :na] = g.atom_indices
            charges[i, :na] = [a.atom_charge for a in g.atoms]
            valid[i, :na] = True
        return indices, charges, valid

    @staticmethod
    def _projection_distance(normal, plane_point, target_point, box=None):
        """投影距离（PBC 感知：t 向量做 MIC）。

        box: 盒尺寸（MIC）；None = 欧氏（无 PBC）。
        """
        from MDAnalysis.lib.distances import minimize_vectors
        d1 = np.linalg.norm(target_point - (plane_point + normal), axis=1)
        d2 = np.linalg.norm(target_point - (plane_point - normal), axis=1)
        sign = np.where(d1 < d2, 1.0, -1.0)[:, None]
        oriented_normal = normal * sign
        t = target_point - plane_point
        if box is not None:
            t = minimize_vectors(t.reshape(-1, 3), box).reshape(-1, 3)
        proj_dist = np.sum(t * oriented_normal, axis=1)
        proj_point = plane_point + (t - proj_dist[:, None] * oriented_normal)
        return np.linalg.norm(proj_point - plane_point, axis=1)
