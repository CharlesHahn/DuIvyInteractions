# -*- coding: utf-8 -*-
"""氢键导出器。"""

from typing import Dict

from .interaction_exporter import InteractionExporter
from ..core.datas import Interaction


class HydrogenBondExporter(InteractionExporter):
    """氢键导出器。"""

    @property
    def name(self) -> str:
        return "Hydrogen Bond"

    @property
    def metric_labels(self) -> Dict[str, str]:
        return {
            "distance": "D-A Distance (Å)",
            "angle": "D-H···A Angle (°)",
        }

    def get_pair_label(self, interaction: Interaction, pair_idx: int) -> str:
        """生成氢键标签。

        格式：'{残基名}{全局残基号}:{原子名}({全局原子号})-{原子名}({全局原子号})···{残基名}{全局残基号}:{原子名}({全局原子号})'
        """
        donor, acceptor = interaction.groups[pair_idx]
        d, h = donor.atoms[0], donor.atoms[1]
        a = acceptor.atoms[0]
        d_label = (f"{donor.residue_name}{donor.residue_id}:"
                   f"{d.atom_name}({d.atom_global_idx})-"
                   f"{h.atom_name}({h.atom_global_idx})")
        a_label = (f"{acceptor.residue_name}{acceptor.residue_id}:"
                   f"{a.atom_name}({a.atom_global_idx})")
        return f"{d_label}···{a_label}"
