# -*- coding: utf-8 -*-
"""疏水导出器。"""

from typing import Dict

from .interaction_exporter import InteractionExporter
from ..core.datas import Interaction


class HydrophobicExporter(InteractionExporter):
    """疏水导出器。"""

    @property
    def name(self) -> str:
        return "Hydrophobic"

    @property
    def metric_labels(self) -> Dict[str, str]:
        return {
            "distance": "Distance (Å)",
        }

    def get_pair_label(self, interaction: Interaction, pair_idx: int) -> str:
        """生成疏水标签。

        格式：'{残基名}{全局残基号}:{原子名}({全局原子号})···{残基名}{全局残基号}:{原子名}({全局原子号})'
        """
        g1, g2 = interaction.groups[pair_idx]
        a1, a2 = g1.atoms[0], g2.atoms[0]
        label1 = (f"{g1.residue_name}{g1.residue_id}:"
                  f"{a1.atom_name}({a1.atom_global_idx})")
        label2 = (f"{g2.residue_name}{g2.residue_id}:"
                  f"{a2.atom_name}({a2.atom_global_idx})")
        return f"{label1}···{label2}"
