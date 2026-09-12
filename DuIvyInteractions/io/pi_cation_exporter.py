# -*- coding: utf-8 -*-
"""π-阳离子导出器。"""

from typing import Dict

from .interaction_exporter import InteractionExporter
from ..core.datas import Interaction


class PiCationExporter(InteractionExporter):
    """π-阳离子导出器。"""

    @property
    def name(self) -> str:
        return "Pi Cation"

    @property
    def metric_labels(self) -> Dict[str, str]:
        return {
            "distance": "Ring-Charge Distance (Å)",
            "offset": "Offset (Å)",
        }

    def get_pair_label(self, interaction: Interaction, pair_idx: int) -> str:
        """生成π-阳离子标签。

        格式：'{残基名}{全局残基号}({起始原子号}-{结束原子号})···{残基名}{全局残基号}({起始原子号}-{结束原子号})'
        """
        ring, cation = interaction.groups[pair_idx]
        ring_label = self._group_label(ring)
        cation_label = self._group_label(cation)
        return f"{ring_label}···{cation_label}"

    @staticmethod
    def _group_label(group) -> str:
        """生成单个基团的标签。"""
        indices = [a.atom_global_idx for a in group.atoms]
        start, end = min(indices), max(indices)
        return f"{group.residue_name}{group.residue_id}({start}-{end})"
