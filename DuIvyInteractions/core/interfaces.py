# -*- coding: utf-8 -*-
"""核心接口定义：Reader（读取器）、GroupIdentifier（识别器）、InteractionDetector（检测器）。"""

from abc import ABC, abstractmethod
from typing import List

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
    """相互作用检测器接口。"""

    @property
    @abstractmethod
    def name(self) -> str:
        """检测器名称。"""
        ...

    @property
    @abstractmethod
    def required_group_types(self) -> List[str]:
        """需要的基团类型列表。"""
        ...

    @abstractmethod
    def detect(self, groups: List[Group], trajectory) -> List[Interaction]:
        """对全部帧检测相互作用。

        内部遍历轨迹，逐帧向量化计算，累积为 Interaction 对象。

        Args:
            groups: 所有基团列表
            trajectory: MDAnalysis 轨迹对象（支持迭代和随机访问）

        Returns:
            检测到的相互作用列表，每个 Interaction 包含全帧数据
        """
        ...
