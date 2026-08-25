# -*- coding: utf-8 -*-
"""核心接口定义：GroupIdentifier（识别器）和 InteractionDetector（检测器）。"""

from abc import ABC, abstractmethod
from typing import List
import numpy as np

from .data import Group, Interaction


class GroupIdentifier(ABC):
    """基团识别器接口。"""

    @property
    @abstractmethod
    def name(self) -> str:
        """识别器名称。"""
        ...

    @abstractmethod
    def identify(self, source: str) -> List[Group]:
        """从数据源识别基团。

        Args:
            source: 数据源（文件路径、SMILES 等）

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
    def detect_frame(self, groups: List[Group], coordinates: np.ndarray,
                     frame: int, time_ps: float) -> List[Interaction]:
        """检测单帧的相互作用。

        Args:
            groups: 所有基团列表
            coordinates: 坐标数组
            frame: 帧号
            time_ps: 时间（皮秒）

        Returns:
            检测到的相互作用列表
        """
        ...
