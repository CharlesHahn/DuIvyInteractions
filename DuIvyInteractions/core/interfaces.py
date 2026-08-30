# -*- coding: utf-8 -*-
"""核心接口定义：Reader（读取器）、GroupIdentifier（识别器）、InteractionDetector（检测器）。"""

from abc import ABC, abstractmethod
from typing import List, Tuple, Dict
from functools import partial
import numpy as np

from .datas import Group, Interaction, SystemData


class Reader(ABC):
    """数据读取器接口。"""

    @property
    @abstractmethod
    def name(self) -> str:
        """读取器名称。"""
        ...

    @abstractmethod
    def read(self, source: str) -> SystemData:
        """从文件读取数据。

        Args:
            source: 文件路径

        Returns:
            SystemData 实例
        """
        ...


class GroupIdentifier(ABC):
    """基团识别器接口。"""

    @property
    @abstractmethod
    def name(self) -> str:
        """识别器名称。"""
        ...

    @abstractmethod
    def identify(self, system_data: SystemData) -> List[Group]:
        """从 SystemData 识别基团。

        Args:
            system_data: 体系数据

        Returns:
            识别到的基团列表
        """
        ...


class InteractionDetector(ABC):
    """相互作用检测器基类（模板方法模式）。

    子类只需实现：name, required_group_types, metric_names,
    get_candidate_tuples, compute_metrics, apply_threshold。

    detect 方法由基类固化算法骨架：
    遍历基团组 → 加载坐标 → 计算指标 → 过滤 → 构建 Interaction。
    """

    # ==================== 子类必须实现 ====================

    @property
    @abstractmethod
    def name(self) -> str:
        """检测器名称（如 "salt_bridge"）。"""
        ...

    @property
    @abstractmethod
    def required_group_types(self) -> List[str]:
        """需要的基团类型列表，Pipeline 据此过滤传入的基团。"""
        ...

    @property
    @abstractmethod
    def metric_names(self) -> List[str]:
        """指标名称列表（如 ["distance", "angle"]）。"""
        ...

    @abstractmethod
    def get_candidate_tuples(self, groups: List[Group],
                             coordinates: np.ndarray = None) -> List[Tuple[Group, ...]]:
        """生成候选基团组。子类实现具体的组合逻辑。

        Args:
            groups: 基团列表
            coordinates: 第一帧坐标数组，可用于边生成边预筛选。默认 None。
        """
        ...

    @abstractmethod
    def compute_metrics(self, group_tuple: Tuple[Group, ...],
                        coords: np.ndarray) -> Dict[str, np.ndarray]:
        """计算单个基团组在全部帧的指标。

        Args:
            group_tuple: 基团组
            coords: (F, n_atoms, 3) 这个基团组的原子在全部帧的坐标

        Returns:
            {name: (F,)} 每帧的指标值
        """
        ...

    @abstractmethod
    def apply_threshold(self, metrics: Dict[str, np.ndarray]) -> np.ndarray:
        """根据指标判断每帧是否存在。

        Args:
            metrics: {name: (F,)} 每帧的指标值

        Returns:
            (F,) bool 每帧是否存在
        """
        ...

    # ==================== 可选覆盖 ====================

    def filter_candidate_tuples(self, tuples: List[Tuple[Group, ...]],
                                coordinates: np.ndarray) -> List[Tuple[Group, ...]]:
        """用一帧坐标过滤候选基团组。默认不过滤。"""
        return tuples

    # ==================== 基类固化（模板方法） ====================

    def detect(self, groups: List[Group], trajectory=None,
               n_workers: int = 1,
               topology_path: str = None,
               trajectory_path: str = None,
               tuple_filter=None) -> List[Interaction]:
        """检测全部帧的相互作用。

        Args:
            groups: 基团列表（已由 Pipeline 按 required_group_types 过滤）
            trajectory: MDAnalysis 轨迹对象（串行时必填，与 trajectory_path 二选一）
            n_workers: 并行 worker 数（1=串行）
            topology_path: 拓扑文件路径（并行时必填）
            trajectory_path: 轨迹文件路径（并行时必填）
            tuple_filter: 可选的基团组过滤函数 (Tuple[Group,...]) -> bool

        Returns:
            检测到的相互作用列表
        """
        # 确定轨迹来源
        if n_workers > 1:
            if topology_path is None or trajectory_path is None:
                raise ValueError(
                    "Parallel execution requires topology_path and trajectory_path")
            import MDAnalysis as mda
            traj_for_filter = mda.Universe(topology_path, trajectory_path).trajectory
        else:
            if trajectory is None:
                raise ValueError("Serial execution requires trajectory")
            traj_for_filter = trajectory

        first_frame = traj_for_filter[0].positions
        tuples = self.get_candidate_tuples(groups, first_frame)

        # 用户自定义过滤（如只检测不同蛋白之间的相互作用）
        if tuple_filter is not None:
            tuples = [t for t in tuples if tuple_filter(t)]

        tuples = self.filter_candidate_tuples(tuples, first_frame)

        if n_workers <= 1:
            results = [r for gt in tuples
                       if (r := self._process_tuple(gt, trajectory)) is not None]
        else:
            from concurrent.futures import ProcessPoolExecutor
            with ProcessPoolExecutor(
                    max_workers=n_workers,
                    initializer=_worker_init,
                    initargs=(topology_path, trajectory_path)) as pool:
                worker = partial(_worker_process_tuple, detector=self)
                results = list(pool.map(worker, tuples))
            results = [r for r in results if r is not None]

        results = self._post_process(results)
        return self._build_interaction(results)

    # ==================== 内部辅助方法 ====================

    def _process_tuple(self, gt: Tuple[Group, ...],
                       trajectory) -> tuple:
        """处理单个基团组。返回 (gt, existence, metrics) 或 None。"""
        n_frames = trajectory.n_frames
        indices = self._get_atom_indices(gt)
        coords = np.zeros((n_frames, len(indices), 3))
        for f, ts in enumerate(trajectory):
            coords[f] = ts.positions[indices]

        metrics = self.compute_metrics(gt, coords)
        existence = self.apply_threshold(metrics)

        if np.any(existence):
            return (gt, existence, metrics)
        return None

    def _get_atom_indices(self, gt: Tuple[Group, ...]) -> np.ndarray:
        """提取基团组中所有原子的全局索引。"""
        return np.array([idx for g in gt for idx in g.atom_indices])

    def _post_process(self, results: list) -> list:
        """后处理钩子：在检测完成后、构建 Interaction 之前调用。

        子类可覆盖此方法实现跨对逻辑（如去重）。
        默认不做任何处理。

        Args:
            results: [(gt, existence, metrics), ...] 列表
        Returns:
            处理后的 results 列表
        """
        return results

    def _build_interaction(self, results: list) -> List[Interaction]:
        """将结果列表构建为 Interaction 对象。"""
        if not results:
            return []
        groups = [r[0] for r in results]
        existence = np.array([r[1] for r in results])
        metrics = {k: np.array([r[2][k] for r in results])
                   for k in results[0][2]}
        return [Interaction(
            interaction_type=self.name,
            groups=groups,
            existence=existence,
            metrics=metrics
        )]


# ==================== 并行 worker（模块级函数，可 pickle） ====================

_worker_trajectory = None  # 每个 worker 进程独立的轨迹对象


def _worker_init(topology_path: str, trajectory_path: str):
    """每个 worker 进程启动时执行一次，加载轨迹。"""
    global _worker_trajectory
    import MDAnalysis as mda
    u = mda.Universe(topology_path, trajectory_path)
    _worker_trajectory = u.trajectory


def _worker_process_tuple(gt: Tuple[Group, ...],
                          detector: InteractionDetector):
    """并行 worker：用已加载的轨迹处理单个基团组。"""
    return detector._process_tuple(gt, _worker_trajectory)
