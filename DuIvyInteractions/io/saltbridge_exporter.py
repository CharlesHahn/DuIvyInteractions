# -*- coding: utf-8 -*-
"""盐桥导出器。"""

from typing import Dict

from .interaction_exporter import InteractionExporter
from ..core.datas import Interaction


class SaltBridgeExporter(InteractionExporter):
    """盐桥导出器。"""

    @property
    def name(self) -> str:
        return "Salt Bridge"

    @property
    def metric_labels(self) -> Dict[str, str]:
        return {
            "distance": "Charge Center Distance (Å)",
        }

    def get_pair_label(self, interaction: Interaction, pair_idx: int) -> str:
        """生成盐桥标签。

        格式：'{残基名}{全局残基号}({起始原子号}-{结束原子号})···{残基名}{全局残基号}({起始原子号}-{结束原子号})'
        """
        pos, neg = interaction.groups[pair_idx]
        pos_label = self._group_label(pos)
        neg_label = self._group_label(neg)
        return f"{pos_label}···{neg_label}"

    @staticmethod
    def _group_label(group) -> str:
        """生成单个带电基团的标签。"""
        indices = [a.atom_global_idx for a in group.atoms]
        start, end = min(indices), max(indices)
        return f"{group.residue_name}{group.residue_id}({start}-{end})"
