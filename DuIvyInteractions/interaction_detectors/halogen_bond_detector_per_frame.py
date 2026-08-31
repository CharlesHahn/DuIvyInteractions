# -*- coding: utf-8 -*-
"""卤键检测器（逐帧策略，向量化）。

判据：X-A ≤ 4.0 Å，C-X···A 在 165°±30°，X···A-R 在 120°±30°。
参考：PLIP, Auffinger et al.
"""

import numpy as np
from typing import List, Tuple, Dict

from ..core.interfaces import InteractionDetectorPerFrame
from ..core.datas import Group, Interaction


# PLIP 阈值
HALOGEN_DIST_MAX = 4.0        # Å
HALOGEN_DON_ANGLE = 165.0     # °，C-X···A 最优角度
HALOGEN_ACC_ANGLE = 120.0     # °，X···A-R 最优角度
HALOGEN_ANGLE_DEV = 30.0      # °，角度偏差上限
PREFILTER_CUTOFF = HALOGEN_DIST_MAX * 3  # 12.0 Å


class HalogenBondDetectorPerFrame(InteractionDetectorPerFrame):
    """卤键检测器（逐帧策略，向量化）。"""

    def __init__(self):
        self._donor_c_idx = None
        self._donor_x_idx = None
        self._acceptor_indices = None
        self._pair_donor = None
        self._pair_acceptor = None

    @property
    def name(self) -> str:
        return "halogen_bond"

    @property
    def required_group_types(self) -> List[str]:
        return ["halogen_donor", "halogen_acceptor"]

    @property
    def metric_names(self) -> List[str]:
        return ["distance", "don_angle", "acc_angle"]

    # ==================== 抽象方法（基类要求，此处不使用） ====================

    def get_candidate_tuples(self, groups, coordinates=None):
        raise NotImplementedError("detect() is overridden, this method is not used")

    def compute_metrics_for_frame(self, tuples, all_positions, frame):
        raise NotImplementedError("detect() is overridden, this method is not used")

    def apply_threshold(self, metrics):
        dist_ok = metrics["distance"] <= HALOGEN_DIST_MAX
        don_ok = (metrics["don_angle"] >= HALOGEN_DON_ANGLE - HALOGEN_ANGLE_DEV) & \
                 (metrics["don_angle"] <= HALOGEN_DON_ANGLE + HALOGEN_ANGLE_DEV)
        acc_ok = (metrics["acc_angle"] >= HALOGEN_ACC_ANGLE - HALOGEN_ANGLE_DEV) & \
                 (metrics["acc_angle"] <= HALOGEN_ACC_ANGLE + HALOGEN_ANGLE_DEV)
        return dist_ok & don_ok & acc_ok

    # ==================== 重写 detect ====================

    def detect(self, groups, trajectory=None,
               n_workers: int = 1,
               topology_path: str = None,
               trajectory_path: str = None,
               tuple_filter=None) -> List[Interaction]:
        """卤键检测主流程。"""
        if trajectory is None:
            raise ValueError("trajectory is required")
        if n_workers > 1:
            raise NotImplementedError("PerFrame 检测器暂不支持并行")

        # 1. 分组
        donors = [g for g in groups if g.group_type == "halogen_donor"]
        acceptors = [g for g in groups if g.group_type == "halogen_acceptor"]

        if not donors or not acceptors:
            return []

        # 2. 构建数据结构
        # donor: 固定 2 原子 [C, X]
        self._donor_c_idx = np.array([g.atoms[0].atom_global_idx for g in donors])
        self._donor_x_idx = np.array([g.atoms[1].atom_global_idx for g in donors])

        # acceptor: [A, R1, R2, ...] 环形 R padding
        self._acceptor_indices = self._build_acceptor_padding(acceptors)

        # 3. 笛卡尔积
        n_d = len(donors)
        n_a = len(acceptors)
        d_grid, a_grid = np.meshgrid(np.arange(n_d), np.arange(n_a), indexing='ij')
        pair_donor = d_grid.ravel()
        pair_acceptor = a_grid.ravel()

        # 4. tuple_filter
        if tuple_filter is not None:
            mask = np.array([
                tuple_filter((donors[di], acceptors[ai]))
                for di, ai in zip(pair_donor, pair_acceptor)
            ])
            pair_donor = pair_donor[mask]
            pair_acceptor = pair_acceptor[mask]

        if len(pair_donor) == 0:
            return []

        # 5. 第一帧：X-A 距离预过滤
        first_pos = trajectory[0].positions
        x_pos = first_pos[self._donor_x_idx[pair_donor]]
        a_pos = first_pos[self._acceptor_indices[pair_acceptor, 0]]
        dist = np.linalg.norm(x_pos - a_pos, axis=1)
        mask = dist < PREFILTER_CUTOFF

        pair_donor = pair_donor[mask]
        pair_acceptor = pair_acceptor[mask]

        if len(pair_donor) == 0:
            return []

        self._pair_donor = pair_donor
        self._pair_acceptor = pair_acceptor

        n_pairs = len(pair_donor)
        n_frames = trajectory.n_frames

        # 6. 预分配结果数组
        existence = np.zeros((n_pairs, n_frames), dtype=bool)
        distance = np.zeros((n_pairs, n_frames))
        don_angle = np.zeros((n_pairs, n_frames))
        acc_angle = np.zeros((n_pairs, n_frames))

        # 7. 逐帧计算
        for f, ts in enumerate(trajectory):
            positions = ts.positions
            metrics = self._compute_metrics(positions)
            distance[:, f] = metrics["distance"]
            don_angle[:, f] = metrics["don_angle"]
            acc_angle[:, f] = metrics["acc_angle"]
            existence[:, f] = self.apply_threshold(metrics)

        # 8. 过滤从未存在的 pair
        has_any = np.any(existence, axis=1)
        if not np.any(has_any):
            return []

        tuples = [(donors[pair_donor[i]], acceptors[pair_acceptor[i]])
                  for i in range(n_pairs) if has_any[i]]

        return [Interaction(
            interaction_type=self.name,
            groups=tuples,
            existence=existence[has_any],
            metrics={
                "distance": distance[has_any],
                "don_angle": don_angle[has_any],
                "acc_angle": acc_angle[has_any],
            }
        )]

    # ==================== 数据结构构建 ====================

    @staticmethod
    def _build_acceptor_padding(acceptors: List[Group]):
        """构建 acceptor padding 索引矩阵。

        存储格式：[A, R1, R2, R3, ...]，R 部分环形填充。
        position 0 是 A，position 1+ 是 R（可能重复）。

        Returns:
            indices: (n_acceptors, max_atoms) int
        """
        n = len(acceptors)
        max_atoms = max(len(g.atoms) for g in acceptors)
        indices = np.zeros((n, max_atoms), dtype=int)

        for i, g in enumerate(acceptors):
            atom_indices = g.atom_indices
            n_atoms = len(atom_indices)
            # position 0 = A
            indices[i, 0] = atom_indices[0]
            # position 1+ = R 部分环形填充
            n_r = n_atoms - 1
            if n_r > 0:
                r_indices = atom_indices[1:]
                for j in range(1, max_atoms):
                    indices[i, j] = r_indices[(j - 1) % n_r]

        return indices

    # ==================== 向量化计算 ====================

    def _compute_metrics(self, positions):
        """向量化计算全部候选对的指标。"""
        pd = self._pair_donor
        pa = self._pair_acceptor

        # donor 坐标
        c = positions[self._donor_c_idx[pd]]   # (n_pairs, 3)
        x = positions[self._donor_x_idx[pd]]   # (n_pairs, 3)

        # acceptor 坐标：A 和全部 R（R 部分环形填充）
        a = positions[self._acceptor_indices[pa, 0]]    # (n_pairs, 3)
        r = positions[self._acceptor_indices[pa, 1:]]   # (n_pairs, max_r, 3)

        # X-A 距离
        xa = x - a
        distance = np.linalg.norm(xa, axis=1)

        # C-X···A 角度（don_angle）
        xc = c - x
        xa2 = a - x
        don_angle = self._vec_angle(xc, xa2)

        # X···A-R 角度（acc_angle）：对每个 R 计算，取最接近 120° 的
        ax = x - a                              # (n_pairs, 3)
        ar = r - a[:, None, :]                  # (n_pairs, max_r, 3)

        # 向量化角度计算
        cos = np.sum(ax[:, None, :] * ar, axis=2)  # (n_pairs, max_r)
        norm_ax = np.linalg.norm(ax, axis=1, keepdims=True)  # (n_pairs, 1)
        norm_ar = np.linalg.norm(ar, axis=2)                  # (n_pairs, max_r)
        cos = cos / (norm_ax * norm_ar)
        cos = np.clip(cos, -1.0, 1.0)
        all_acc_angles = np.degrees(np.arccos(cos))  # (n_pairs, max_r)

        # 取最接近 120° 的 R
        diffs = np.abs(all_acc_angles - HALOGEN_ACC_ANGLE)
        best_r = np.argmin(diffs, axis=1)  # (n_pairs,)
        acc_angle = all_acc_angles[np.arange(len(pd)), best_r]

        return {
            "distance": distance,
            "don_angle": don_angle,
            "acc_angle": acc_angle,
        }

    @staticmethod
    def _vec_angle(v1, v2):
        """计算两组向量的夹角（度）。"""
        cos_angle = np.sum(v1 * v2, axis=1) / (
            np.linalg.norm(v1, axis=1) * np.linalg.norm(v2, axis=1)
        )
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        return np.degrees(np.arccos(cos_angle))
