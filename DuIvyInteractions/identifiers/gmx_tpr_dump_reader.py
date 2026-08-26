# -*- coding: utf-8 -*-
"""从 gmx dump 文本读取数据。"""

from typing import List
from ..core.interfaces import Reader
from ..core.data import SystemData


class GmxTprDumpReader(Reader):
    """从 gmx dump 文本读取数据。"""

    @property
    def name(self) -> str:
        return "gmx_tpr_dump"

    def read(self, source: str) -> SystemData:
        """从 gmx dump 文本读取数据。

        Args:
            source: gmx dump 输出的文本文件路径

        Returns:
            SystemData 实例
        """
        # 1. 解析文本
        moltypes = self._parse_dump(source)
        
        # 2. 分配全局索引
        self._assign_global_indices(moltypes)
        
        # 3. 构建 SystemData
        return self._build_system_data(moltypes, source)

    def _parse_dump(self, source: str) -> list:
        """解析 gmx dump 文本。

        Args:
            source: 文件路径

        Returns:
            解析后的分子类型列表
        """
        # TODO: 复用 parse_tpr_dump.py 的逻辑
        ...

    def _assign_global_indices(self, moltypes: list) -> None:
        """分配全局索引。

        Args:
            moltypes: 分子类型列表
        """
        # TODO: 为原子和残基分配全局索引
        ...

    def _build_system_data(self, moltypes: list, source: str) -> SystemData:
        """构建 SystemData。

        Args:
            moltypes: 分子类型列表
            source: 文件路径（用于提取 system_name）

        Returns:
            SystemData 实例
        """
        # TODO: 检查 Bond 和 Constraint 的关系
        # TODO: 构建 SystemData
        ...
