# -*- coding: utf-8 -*-
"""π-堆积导出器。"""

from typing import Dict, Optional, List

import numpy as np

from .interaction_exporter import InteractionExporter, DIT_COLORS
from DuIvyTools.DuIvyTools.FileParser.xpmParser import XPM
from ..core.datas import Interaction


class PiStackingExporter(InteractionExporter):
    """π-堆积导出器。"""

    @property
    def name(self) -> str:
        return "Pi Stacking"

    @property
    def metric_labels(self) -> Dict[str, str]:
        return {
            "distance": "Ring Center Distance (Å)",
            "angle": "Normal Vector Angle (°)",
            "offset": "Offset (Å)",
            "pistacking_type": "Type (P/T/N)",
        }

    def get_pair_label(self, interaction: Interaction, pair_idx: int) -> str:
        """生成π-堆积标签。

        格式：'{残基名}{全局残基号}({起始原子号}-{结束原子号})···{残基名}{全局残基号}({起始原子号}-{结束原子号})'
        """
        ring1, ring2 = interaction.groups[pair_idx]
        ring1_label = self._ring_label(ring1)
        ring2_label = self._ring_label(ring2)
        return f"{ring1_label}···{ring2_label}"

    @staticmethod
    def _ring_label(ring) -> str:
        """生成单个芳香环的标签。"""
        indices = [a.atom_global_idx for a in ring.atoms]
        start, end = min(indices), max(indices)
        return f"{ring.residue_name}{ring.residue_id}({start}-{end})"

    def to_xpm_stacking_type(
        self,
        interaction: Interaction,
        pair_indices: Optional[List[int]] = None,
        title: Optional[str] = None,
        xlabel: str = "Time (ps)",
        ylabel: str = "Pair",
    ) -> XPM:
        """输出堆积类型的 XPM（0=无, 1=T 型, 2=P 型）。

        Args:
            interaction: Interaction 数据
            pair_indices: 要导出的 pair 索引，None 表示全部
            title: 图表标题
            xlabel: X 轴标签
            ylabel: Y 轴标签

        Returns:
            XPM 对象
        """
        if interaction.n_pairs == 0:
            raise ValueError("No pairs found in interaction")

        pair_indices = self._validate_pair_indices(interaction, pair_indices)

        existence = interaction.existence[pair_indices]
        pstype = interaction.metrics["pistacking_type"][pair_indices]

        # 构建三值矩阵：0=无, 1=T, 2=P
        vm = np.zeros(existence.shape, dtype=int)
        vm[(existence) & (pstype == 'T')] = 1
        vm[(existence) & (pstype == 'P')] = 2

        # 构建 XPM（value_matrix 即颜色索引，colors/notes 按索引对齐）
        pair_legends = self.get_pair_legends(interaction, pair_indices)
        return self._build_discrete_xpm(
            value_matrix=vm,
            colors=["#FFFFFF", DIT_COLORS[1], DIT_COLORS[0]],
            notes=["None", "T-shaped", "Parallel"],
            title=title or f"{self.name} Type",
            legend=" ".join(f"{i}:{label}" for i, label in enumerate(pair_legends)),
            xlabel=xlabel,
            ylabel=ylabel,
            times=interaction.times,
        )

    def save_xpm_stacking_type(
        self,
        interaction: Interaction,
        path: str,
        **kwargs
    ) -> None:
        """保存堆积类型 XPM 文件。"""
        xpm = self.to_xpm_stacking_type(interaction, **kwargs)
        xpm.save(path)
