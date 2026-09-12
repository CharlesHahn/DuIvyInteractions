# -*- coding: utf-8 -*-
"""金属配位导出器。"""

from typing import Dict

from .interaction_exporter import InteractionExporter
from ..core.datas import Interaction


class MetalCoordinationExporter(InteractionExporter):
    """金属配位导出器。"""

    @property
    def name(self) -> str:
        return "Metal Coordination"

    @property
    def metric_labels(self) -> Dict[str, str]:
        return {
            "distance": "Metal-Ligand Distance (Å)",
        }

    def get_pair_label(self, interaction: Interaction, pair_idx: int) -> str:
        """生成金属配位标签。

        格式：'{残基名}{全局残基号}:{原子名}({全局原子号})···{残基名}{全局残基号}:{原子名}({全局原子号})'
        """
        metal, binding = interaction.groups[pair_idx]
        m, b = metal.atoms[0], binding.atoms[0]
        m_label = (f"{metal.residue_name}{metal.residue_id}:"
                   f"{m.atom_name}({m.atom_global_idx})")
        b_label = (f"{binding.residue_name}{binding.residue_id}:"
                   f"{b.atom_name}({b.atom_global_idx})")
        return f"{m_label}···{b_label}"
