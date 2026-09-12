# -*- coding: utf-8 -*-
"""π-堆积导出器。"""

from typing import Dict

from .interaction_exporter import InteractionExporter
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
