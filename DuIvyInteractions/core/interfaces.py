# -*- coding: utf-8 -*-
"""核心接口定义：Reader（读取器）、GroupIdentifier（识别器）、InteractionDetector（检测器）。"""

from abc import ABC, abstractmethod
from typing import List, Tuple, Dict
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
    def get_candidate_tuples(self, groups: List[Group]) -> List[Tuple[Group, ...]]:
        """生成候选基团组。子类实现具体的组合逻辑。"""
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

    def detect(self, groups: List[Group], trajectory) -> List[Interaction]:
        """检测全部帧的相互作用。

        Args:
            groups: 基团列表（已由 Pipeline 按 required_group_types 过滤）
            trajectory: MDAnalysis 轨迹对象

        Returns:
            检测到的相互作用列表
        """
        tuples = self.get_candidate_tuples(groups)

        # 可选预过滤：加载第一帧坐标
        tuples = self.filter_candidate_tuples(tuples, trajectory[0].positions)

        n_frames = trajectory.n_frames
        results = []

        for gt in tuples:
            indices = self._get_atom_indices(gt)
            coords = np.zeros((n_frames, len(indices), 3))
            for f, ts in enumerate(trajectory):
                coords[f] = ts.positions[indices]

            metrics = self.compute_metrics(gt, coords)
            existence = self.apply_threshold(metrics)
            if np.any(existence):
                results.append((gt, existence, metrics))

        return self._build_interaction(results)

    # ==================== 内部辅助方法 ====================

    def _get_atom_indices(self, gt: Tuple[Group, ...]) -> np.ndarray:
        """提取基团组中所有原子的全局索引。"""
        return np.array([idx for g in gt for idx in g.atom_indices])

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
