# -*- coding: utf-8 -*-
"""Interaction Exporter 基类。

将 Interaction 数据导出为 DuIvyTools 支持的 xvg 和 xpm 格式。
"""

from abc import ABC, abstractmethod
from typing import List, Dict, Optional

import numpy as np

from DuIvyTools.DuIvyTools.FileParser.xvgParser import XVG
from DuIvyTools.DuIvyTools.FileParser.xpmParser import XPM

from ..core.datas import Interaction


class InteractionExporter(ABC):
    """Interaction 导出器基类。

    将 Interaction 数据导出为 DuIvyTools 支持的 xvg 和 xpm 格式。
    所有相互作用类型的导出逻辑相同，子类只需定义名称和标签。
    """

    # ==================== 子类必须实现 ====================

    @property
    @abstractmethod
    def name(self) -> str:
        """相互作用名称，如 'Hydrogen Bond'。"""
        ...

    @property
    @abstractmethod
    def metric_labels(self) -> Dict[str, str]:
        """metric 名称到标签的映射。

        示例：{"distance": "D-A Distance (Å)", "angle": "D-H···A Angle (°)"}
        """
        ...

    # ==================== 通用导出方法 ====================

    def to_xvg_metric(
        self,
        interaction: Interaction,
        metric_name: str,
        pair_indices: Optional[List[int]] = None,
        title: Optional[str] = None,
        xlabel: str = "Frame",
        ylabel: Optional[str] = None,
    ) -> XVG:
        """将某个 metric 转换为 XVG 对象。

        Args:
            interaction: Interaction 数据
            metric_name: metric 名称，如 "distance", "angle"
            pair_indices: 要导出的 pair 索引，None 表示全部
            title: 图表标题，None 则自动生成
            xlabel: X 轴标签，默认为 "Frame"
            ylabel: Y 轴标签，None 则从 metric_labels 获取

        Returns:
            XVG 对象

        Raises:
            ValueError: metric_name 不存在或 pair_indices 超出范围
        """
        # 验证 interaction 非空
        if not interaction.metrics:
            raise ValueError("No metrics found in interaction")

        # 验证 metric_name
        if metric_name not in interaction.metrics:
            raise ValueError(
                f"Metric '{metric_name}' not found. "
                f"Available: {list(interaction.metrics.keys())}"
            )

        # 确定并验证 pair 索引
        pair_indices = self._validate_pair_indices(interaction, pair_indices)

        # 获取数据
        data = interaction.metrics[metric_name]  # (n_pairs, n_frames)
        n_frames = interaction.n_frames

        # 构建 XVG 对象
        xvg = XVG("", is_file=False, new_file=True)

        # 设置标题和标签
        xvg.title = title or f"{self.name} - {self.metric_labels.get(metric_name, metric_name)}"
        xvg.xlabel = xlabel
        xvg.ylabel = ylabel or self.metric_labels.get(metric_name, metric_name)

        # 设置图例
        xvg.legends = self.get_pair_legends(interaction, pair_indices)

        # 构建数据列
        # 第 0 列：帧号（X 轴）
        xvg.data_columns = [list(range(n_frames))]

        # 后续列：每个 pair 的数据
        for idx in pair_indices:
            xvg.data_columns.append(data[idx].tolist())

        # 设置维度
        xvg.column_num = len(xvg.data_columns)
        xvg.row_num = n_frames

        # 设置数据头
        xvg.data_heads = [xlabel] + xvg.legends

        return xvg

    def to_xpm_existence(
        self,
        interaction: Interaction,
        pair_indices: Optional[List[int]] = None,
        title: Optional[str] = None,
        xlabel: str = "Frame",
        ylabel: str = "Pair",
    ) -> XPM:
        """将 existence 转换为 XPM 对象（热力图）。

        Args:
            interaction: Interaction 数据
            pair_indices: 要导出的 pair 索引，None 表示全部
            title: 图表标题，None 则自动生成
            xlabel: X 轴标签
            ylabel: Y 轴标签

        Returns:
            XPM 对象

        Raises:
            ValueError: interaction 为空或 pair_indices 超出范围
        """
        # 验证 interaction 非空
        if interaction.n_pairs == 0:
            raise ValueError("No pairs found in interaction")

        # 确定并验证 pair 索引
        pair_indices = self._validate_pair_indices(interaction, pair_indices)

        # 获取数据
        existence = interaction.existence[pair_indices]  # (n_selected, n_frames)
        n_frames = interaction.n_frames
        n_pairs = len(pair_indices)

        # 构建 XPM 对象
        xpm = XPM("", is_file=False, new_file=True)

        # 设置基本属性
        xpm.title = title or f"{self.name} Existence"
        xpm.xlabel = xlabel
        xpm.ylabel = ylabel
        xpm.legend = "Existence"
        xpm.type = "Discrete"

        # 设置维度
        xpm.width = n_frames
        xpm.height = n_pairs

        # 设置坐标轴
        xpm.xaxis = list(range(n_frames))
        xpm.yaxis = list(range(n_pairs))

        # 构建值矩阵（0 或 1）
        xpm.value_matrix = existence.astype(float).tolist()

        # 刷新颜色和字符
        xpm.refresh_by_value_matrix(is_Continuous=False)

        return xpm

    # ==================== 保存方法 ====================

    def save_xvg(
        self,
        interaction: Interaction,
        metric_name: str,
        path: str,
        **kwargs
    ) -> None:
        """保存某个 metric 为 xvg 文件。

        Args:
            interaction: Interaction 数据
            metric_name: metric 名称
            path: 输出文件路径
            **kwargs: 传递给 to_xvg_metric 的参数
        """
        xvg = self.to_xvg_metric(interaction, metric_name, **kwargs)
        xvg.save(path)

    def save_xpm(
        self,
        interaction: Interaction,
        path: str,
        **kwargs
    ) -> None:
        """保存 existence 为 xpm 文件。

        Args:
            interaction: Interaction 数据
            path: 输出文件路径
            **kwargs: 传递给 to_xpm_existence 的参数
        """
        xpm = self.to_xpm_existence(interaction, **kwargs)
        xpm.save(path)

    # ==================== 辅助方法 ====================

    def get_pair_label(self, interaction: Interaction, pair_idx: int) -> str:
        """生成 pair 的标签，如 'ARG73-D927'。

        Args:
            interaction: Interaction 数据
            pair_idx: pair 索引

        Returns:
            标签字符串
        """
        g_tuple = interaction.groups[pair_idx]
        labels = []
        for g in g_tuple:
            label = f"{g.residue_name}{g.residue_id}"
            labels.append(label)
        return "-".join(labels)

    def get_pair_legends(self, interaction: Interaction, pair_indices: List[int]) -> List[str]:
        """生成多个 pair 的图例列表。

        Args:
            interaction: Interaction 数据
            pair_indices: pair 索引列表

        Returns:
            图例列表
        """
        return [self.get_pair_label(interaction, i) for i in pair_indices]

    # ==================== 内部辅助方法 ====================

    def _validate_pair_indices(
        self, interaction: Interaction, pair_indices: Optional[List[int]]
    ) -> List[int]:
        """验证并返回 pair 索引列表。

        Args:
            interaction: Interaction 数据
            pair_indices: pair 索引列表，None 表示全部

        Returns:
            验证后的 pair 索引列表

        Raises:
            ValueError: pair_indices 超出范围
        """
        max_idx = interaction.n_pairs - 1

        if pair_indices is None:
            return list(range(interaction.n_pairs))

        # 验证每个索引
        for idx in pair_indices:
            if not isinstance(idx, int):
                raise ValueError(f"pair_index must be int, got {type(idx)}")
            if idx < 0 or idx > max_idx:
                raise ValueError(f"pair_index {idx} out of range [0, {max_idx}]")

        return pair_indices
