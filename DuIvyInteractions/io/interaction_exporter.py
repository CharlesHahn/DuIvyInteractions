# -*- coding: utf-8 -*-
"""Interaction Exporter 基类。

将 Interaction 数据导出为 DuIvyTools 支持的 xvg 和 xpm 格式。
"""

import csv
import string
from abc import ABC, abstractmethod
from typing import List, Dict, Optional

import numpy as np

from DuIvyTools.DuIvyTools.FileParser.xvgParser import XVG
from DuIvyTools.DuIvyTools.FileParser.xpmParser import XPM

from ..core.datas import Interaction

# DIT.mplstyle 颜色循环（DuIvyTools 默认配色）
DIT_COLORS = [
    '#38A7D0', '#F67088', '#66C2A5', '#FC8D62', '#8DA0CB',
    '#E78AC3', '#A6D854', '#FFD92F', '#E5C494', '#B3B3B3',
]

# XPM 字符表：与 DuIvyTools refresh 的 letters 方案一致
XPM_LETTERS = string.ascii_letters + "0123456789!@#$%^&*()-_=+{}|;"


def _validate_xpm_indices(value_matrix: np.ndarray, n_colors: int) -> None:
    """校验 value_matrix 是 [0, n_colors) 内的整数索引。"""
    arr = np.asarray(value_matrix)
    if arr.dtype.kind not in "iu":
        raise TypeError(f"value_matrix must be integer, got {arr.dtype}")
    if arr.size and (arr.min() < 0 or arr.max() >= n_colors):
        raise ValueError(f"value_matrix index out of range [0, {n_colors}): "
                         f"min={arr.min()}, max={arr.max()}")


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
        xlabel: str = "Time (ps)",
        ylabel: Optional[str] = None,
    ) -> XVG:
        """将某个 metric 转换为 XVG 对象。

        Args:
            interaction: Interaction 数据
            metric_name: metric 名称，如 "distance", "angle"
            pair_indices: 要导出的 pair 索引，None 表示全部
            title: 图表标题，None 则自动生成
            xlabel: X 轴标签，默认为 "Time (ps)"
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
        # 第 0 列：时间（X 轴）
        xvg.data_columns = [interaction.times.tolist()]

        # 后续列：每个 pair 的数据
        for idx in pair_indices:
            xvg.data_columns.append(data[idx].tolist())

        # 设置维度
        xvg.column_num = len(xvg.data_columns)
        xvg.row_num = n_frames

        # 设置数据头
        xvg.data_heads = [xlabel] + xvg.legends

        return xvg

    def to_xvg_count(
        self,
        interaction: Interaction,
        title: Optional[str] = None,
        xlabel: str = "Time (ps)",
        ylabel: Optional[str] = None,
    ) -> XVG:
        """将每帧的活跃相互作用数量转换为 XVG 对象。

        Args:
            interaction: Interaction 数据
            title: 图表标题，None 则自动生成
            xlabel: X 轴标签
            ylabel: Y 轴标签

        Returns:
            XVG 对象
        """
        count = np.sum(interaction.existence, axis=0).astype(int)

        xvg = XVG("", is_file=False, new_file=True)
        xvg.title = title or f"{self.name} Count"
        xvg.xlabel = xlabel
        xvg.ylabel = ylabel or "Count"
        xvg.data_columns = [interaction.times.tolist(), count.tolist()]
        xvg.column_num = 2
        xvg.row_num = interaction.n_frames
        xvg.data_heads = [xlabel, ylabel or "Count"]

        return xvg

    def to_xpm_existence(
        self,
        interaction: Interaction,
        pair_indices: Optional[List[int]] = None,
        title: Optional[str] = None,
        xlabel: str = "Time (ps)",
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

        # 构建 XPM（value_matrix 即颜色索引，colors/notes 按索引对齐）
        pair_legends = self.get_pair_legends(interaction, pair_indices)
        return self._build_discrete_xpm(
            value_matrix=interaction.existence[pair_indices].astype(int),
            colors=["#FFFFFF", DIT_COLORS[0]],
            notes=["No", "Yes"],
            title=title or f"{self.name} Existence",
            legend=" ".join(f"{i}:{label}" for i, label in enumerate(pair_legends)),
            xlabel=xlabel,
            ylabel=ylabel,
            times=interaction.times,
        )

    def _build_discrete_xpm(
        self,
        value_matrix: np.ndarray,   # (n_pairs, n_frames) 颜色索引 0..k-1
        colors: List[str],          # 长度 k，与索引一一对应
        notes: List[str],           # 长度 k
        title: str,
        legend: str,
        xlabel: str,
        ylabel: str,
        times: np.ndarray,
    ) -> XPM:
        """构建 Discrete 型 XPM：value_matrix 即颜色索引，colors/notes 按索引对齐。"""
        n_pairs, n_frames = value_matrix.shape
        if len(colors) != len(notes):
            raise ValueError(f"colors({len(colors)}) != notes({len(notes)})")
        if len(colors) > len(XPM_LETTERS):
            raise ValueError(f"too many colors ({len(colors)}), max {len(XPM_LETTERS)}")
        _validate_xpm_indices(value_matrix, len(colors))
        xpm = XPM("", is_file=False, new_file=True)
        xpm.title, xpm.legend = title, legend
        xpm.xlabel, xpm.ylabel = xlabel, ylabel
        xpm.type = "Discrete"
        xpm.width, xpm.height = n_frames, n_pairs
        xpm.color_num, xpm.char_per_pixel = len(colors), 1
        xpm.chars, xpm.colors, xpm.notes = (
            list(XPM_LETTERS[:len(colors)]), list(colors), list(notes))
        xpm.xaxis, xpm.yaxis = times.tolist(), list(range(n_pairs))
        xpm.value_matrix = value_matrix.astype(int).tolist()
        xpm.datalines = ["".join(XPM_LETTERS[v] for v in row)
                         for row in value_matrix.astype(int)]
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

    def save_xvg_count(
        self,
        interaction: Interaction,
        path: str,
        **kwargs
    ) -> None:
        """保存每帧活跃相互作用数量为 xvg 文件。

        Args:
            interaction: Interaction 数据
            path: 输出文件路径
            **kwargs: 传递给 to_xvg_count 的参数
        """
        xvg = self.to_xvg_count(interaction, **kwargs)
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

    def to_csv_summary(self, interaction: Interaction, path: str) -> None:
        """每对的汇总统计（仅活跃帧）保存为 CSV。

        列：pair_label, occupancy, avg_{metric}, std_{metric}, ...

        Args:
            interaction: Interaction 数据
            path: 输出文件路径
        """
        occ = interaction.occupancy()

        # 表头：只输出数值指标，跳过字符串指标
        headers = ["pair_label", "occupancy"]
        numeric_metrics = []
        for name, arr in interaction.metrics.items():
            is_str = hasattr(arr, 'dtype') and arr.dtype.kind in ('U', 'S', 'O')
            if not is_str:
                numeric_metrics.append(name)
                headers.append(f"avg_{name}")
                headers.append(f"std_{name}")

        # 逐行
        rows = []
        for i in range(interaction.n_pairs):
            label = self.get_pair_label(interaction, i)
            row = [label, f"{occ[i]:.4f}"]
            for name in numeric_metrics:
                active = interaction.metrics[name][i][interaction.existence[i]]
                if active.size > 0:
                    row.append(f"{np.nanmean(active):.4f}")
                    row.append(f"{np.nanstd(active):.4f}")
                else:
                    row.append("")
                    row.append("")
            rows.append(row)

        # 写文件
        with open(path, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(headers)
            writer.writerows(rows)

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
