# -*- coding: utf-8 -*-
"""核心接口定义：Reader、GroupIdentifier、InteractionDetectorPerTuple、InteractionDetectorPerFrame、InteractionDetectorTwoPass。"""

from abc import ABC, abstractmethod
from typing import List, Tuple, Dict
from functools import partial
import numpy as np

from .datas import Group, Interaction, InteractionSparse, SystemData


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


class InteractionDetectorPerTuple(ABC):
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


class InteractionDetectorPerFrame(ABC):
    """相互作用检测器基类（逐帧策略）。

    与 PerTuple 的区别：
    - PerTuple: for each tuple → load ALL frames → vectorized over frames
    - PerFrame: for each frame → process ALL tuples → vectorized over tuples

    适用场景：候选 tuple 数量极大（如水桥），逐 tuple 遍历轨迹不可接受。
    """

    # ==================== 子类必须实现 ====================

    @property
    @abstractmethod
    def name(self) -> str:
        """检测器名称（如 "water_bridge"）。"""
        ...

    @property
    @abstractmethod
    def required_group_types(self) -> List[str]:
        """需要的基团类型列表。"""
        ...

    @property
    @abstractmethod
    def metric_names(self) -> List[str]:
        """指标名称列表。"""
        ...

    @abstractmethod
    def get_candidate_tuples(self, groups: List[Group],
                             coordinates: np.ndarray = None
                             ) -> List[Tuple[Group, ...]]:
        """生成候选基团组。"""
        ...

    @abstractmethod
    def compute_metrics_for_frame(
        self,
        tuples: List[Tuple[Group, ...]],
        all_positions: np.ndarray,
        frame: int,
    ) -> Dict[str, np.ndarray]:
        """对全部 tuple 计算单帧指标。

        Args:
            tuples: 全部候选基团组
            all_positions: (n_atoms_total, 3) 当前帧全部原子坐标（Å）
            frame: 帧号

        Returns:
            {name: (n_tuples,)} 每个 tuple 在当前帧的指标值
        """
        ...

    @abstractmethod
    def apply_threshold(self, metrics: Dict[str, np.ndarray]) -> np.ndarray:
        """根据指标判断是否存在。shape: (n_tuples,) → (n_tuples,) bool。"""
        ...

    # ==================== 可选覆盖 ====================

    def filter_candidate_tuples(self, tuples: List[Tuple[Group, ...]],
                                coordinates: np.ndarray) -> List[Tuple[Group, ...]]:
        """用一帧坐标过滤候选基团组。默认不过滤。"""
        return tuples

    def _post_process(self, results: list) -> list:
        """后处理钩子。默认不做任何处理。"""
        return results

    # ==================== 基类固化（模板方法） ====================

    def detect(self, groups: List[Group], trajectory=None,
               n_workers: int = 1,
               topology_path: str = None,
               trajectory_path: str = None,
               tuple_filter=None) -> List[Interaction]:
        """检测全部帧的相互作用。

        接口与 PerTuple.detect 完全一致，便于替换。
        """
        if trajectory is None:
            raise ValueError("trajectory is required")

        # 1. 生成候选 tuple
        first_frame = trajectory[0].positions
        tuples = self.get_candidate_tuples(groups, first_frame)

        if tuple_filter is not None:
            tuples = [t for t in tuples if tuple_filter(t)]

        tuples = self.filter_candidate_tuples(tuples, first_frame)

        if not tuples:
            return []

        n_tuples = len(tuples)
        n_frames = trajectory.n_frames

        # 2. 预分配结果数组
        existence = np.zeros((n_tuples, n_frames), dtype=bool)
        metrics = {name: np.zeros((n_tuples, n_frames))
                   for name in self.metric_names}

        # 3. 逐帧处理
        for f, ts in enumerate(trajectory):
            frame_metrics = self.compute_metrics_for_frame(
                tuples, ts.positions, f)
            existence[:, f] = self.apply_threshold(frame_metrics)
            for name in self.metric_names:
                metrics[name][:, f] = frame_metrics[name]

        # 4. 过滤从未存在的 tuple
        has_any = np.any(existence, axis=1)
        if not np.any(has_any):
            return []

        results = [(tuples[i], existence[i],
                    {k: v[i] for k, v in metrics.items()})
                   for i in range(n_tuples) if has_any[i]]

        # 5. 后处理
        results = self._post_process(results)

        # 6. 构建 Interaction
        return self._build_interaction(results)

    # ==================== 内部辅助方法 ====================

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


class InteractionDetectorTwoPass(ABC):
    """相互作用检测器基类（两轮遍历 + 稀疏存储）。

    与 PerTuple / PerFrame 的区别：
    - PerTuple: for each tuple → load ALL frames → vectorized over frames
    - PerFrame: for each frame → compute ALL tuples → dense storage
    - TwoPass:  Pass1 discover active tuples → Pass2 fill full metrics

    两阶段职责：
    - Pass1（发现）：逐帧精确检测哪些基团组存在相互作用，结果稀疏存储。
    - Pass2（补全）：对 Pass1 发现的基团组遍历全帧，补全完整时间序列数据。

    公共接口：
    - run_pass1(groups, trajectory, tuple_filter) → InteractionSparse
    - run_pass2(sparse, trajectory) → List[Interaction]
    - detect(groups, trajectory, ...) → List[Interaction]（便捷方法 = Pass1 + Pass2）

    Pass1 和 Pass2 只通过 InteractionSparse 传递数据，无耦合。

    子类必须实现：
    - name, required_group_types, metric_names: 元信息
    - initialize_candidates: 初始化候选基团组
    - compute_pair_metrics: 计算指标
    - apply_threshold: 判定 existence

    子类可覆盖：
    - run_pass1: 直接覆盖实现自定义逻辑（如 KDTree）
    - run_pass2: 直接覆盖实现自定义逻辑
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
        """需要的基团类型列表。"""
        ...

    @property
    @abstractmethod
    def metric_names(self) -> List[str]:
        """指标名称列表（如 ["distance", "angle"]）。"""
        ...

    def initialize_candidates(self, groups: List[Group], trajectory,
                              tuple_filter=None) -> list:
        """初始化候选基团组。子类实现。

        Args:
            groups: 基团列表
            trajectory: MDAnalysis 轨迹对象
            tuple_filter: 可选的基团组过滤函数

        Returns:
            items: 候选基团组列表
        """
        return []

    def compute_pair_metrics(self, group_tuples: List[Tuple[Group, ...]],
                             all_positions: np.ndarray) -> Dict[str, np.ndarray]:
        """对给定的基团组列表计算指标。子类实现。

        Args:
            group_tuples: 基团组列表，每组是一个 tuple
            all_positions: (n_atoms_total, 3) 当前帧全部原子坐标（Å）

        Returns:
            {name: (n_groups,)} 每个基团组的指标值
        """
        return {}

    def apply_threshold(self, metrics: Dict[str, np.ndarray]) -> np.ndarray:
        """根据指标判定 existence。子类实现。

        Args:
            metrics: {name: (n_groups,)} 每个基团组的指标值

        Returns:
            (n_groups,) bool，每组是否存在相互作用
        """
        return np.array([])

    # ==================== 子类可覆盖（有默认实现） ====================

    def run_pass1(self, groups: List[Group], trajectory,
                  tuple_filter=None) -> InteractionSparse:
        """执行 Pass1：逐帧发现 active pairs，返回稀疏结果。

        默认实现：调用 initialize_candidates + compute_pair_metrics + apply_threshold。
        子类可覆盖实现自定义逻辑（如 KDTree、级联筛选）。

        Args:
            groups: 基团列表
            trajectory: MDAnalysis 轨迹对象
            tuple_filter: 可选的基团组过滤函数

        Returns:
            InteractionSparse 对象
        """
        group_tuples = self.initialize_candidates(groups, trajectory, tuple_filter)
        if not group_tuples:
            return InteractionSparse(interaction_type=self.name, data={})

        sparse_data: Dict[Tuple[int, ...], dict] = {}

        for f, ts in enumerate(trajectory):
            all_metrics = self.compute_pair_metrics(group_tuples, ts.positions)
            mask = self.apply_threshold(all_metrics)

            for idx in np.where(mask)[0]:
                group_tuple = group_tuples[idx]
                group_ids = tuple(g.group_id for g in group_tuple)

                if group_ids not in sparse_data:
                    sparse_data[group_ids] = {
                        "groups": group_tuple,
                        "frames": [],
                        "metrics": {name: [] for name in self.metric_names}
                    }

                sparse_data[group_ids]["frames"].append(f)
                for name in self.metric_names:
                    sparse_data[group_ids]["metrics"][name].append(
                        all_metrics[name][idx])

        return InteractionSparse(interaction_type=self.name, data=sparse_data)

    def run_pass2(self, sparse: InteractionSparse,
                  trajectory) -> List[Interaction]:
        """执行 Pass2：补全全帧数据。

        默认实现：调用 compute_pair_metrics + apply_threshold。
        子类可覆盖实现自定义逻辑。

        Args:
            sparse: run_pass1 返回的稀疏结果
            trajectory: MDAnalysis 轨迹对象

        Returns:
            List[Interaction]，完整的相互作用结果
        """
        if not sparse.data:
            return []

        group_tuples = [entry["groups"] for entry in sparse.data.values()]
        n_groups = len(group_tuples)
        n_frames = trajectory.n_frames

        # 预分配稠密数组
        existence = np.zeros((n_groups, n_frames), dtype=bool)
        metrics = {name: np.full((n_groups, n_frames), np.nan)
                   for name in self.metric_names}

        # 逐帧计算全部 discovered groups 的 metric
        for f, ts in enumerate(trajectory):
            frame_metrics = self.compute_pair_metrics(group_tuples, ts.positions)
            existence[:, f] = self.apply_threshold(frame_metrics)
            for name in self.metric_names:
                metrics[name][:, f] = frame_metrics[name]

        # 构建 results 列表
        results = [
            (entry["groups"], existence[i],
             {k: v[i] for k, v in metrics.items()})
            for i, entry in enumerate(sparse.data.values())
        ]

        results = self._post_process(results)
        return self._build_interaction(results)

    # ==================== 公共接口 ====================

    def detect(self, groups: List[Group], trajectory=None,
               n_workers: int = 1,
               topology_path: str = None,
               trajectory_path: str = None,
               tuple_filter=None) -> List[Interaction]:
        """便捷方法：Pass1 + Pass2。

        接口与 PerTuple / PerFrame.detect 完全一致，便于替换。
        """
        if trajectory is None:
            raise ValueError("trajectory is required")

        sparse = self.run_pass1(groups, trajectory, tuple_filter)
        return self.run_pass2(sparse, trajectory)

    # ==================== 内部辅助方法 ====================

    def _post_process(self, results: list) -> list:
        """后处理钩子。默认不做任何处理。"""
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
                          detector: InteractionDetectorPerTuple):
    """并行 worker：用已加载的轨迹处理单个基团组。"""
    return detector._process_tuple(gt, _worker_trajectory)
