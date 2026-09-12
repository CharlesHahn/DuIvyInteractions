# -*- coding: utf-8 -*-
"""卤键导出器。"""

from typing import Dict

from .interaction_exporter import InteractionExporter
from ..core.datas import Interaction


class HalogenBondExporter(InteractionExporter):
    """卤键导出器。"""

    @property
    def name(self) -> str:
        return "Halogen Bond"

    @property
    def metric_labels(self) -> Dict[str, str]:
        return {
            "distance": "X···A Distance (Å)",
            "don_angle": "C-X···A Angle (°)",
            "acc_angle": "X···A-R Angle (°)",
        }

    def get_pair_label(self, interaction: Interaction, pair_idx: int) -> str:
        """生成卤键标签。

        格式：'{残基名}{全局残基号}:{原子名}({全局原子号})-{原子名}({全局原子号})···{残基名}{全局残基号}:{原子名}({全局原子号})'
        供体展示 C-X，受体展示 A（R 逐帧动态选择，无法确定）。
        """
        donor, acceptor = interaction.groups[pair_idx]
        c, x = donor.atoms[0], donor.atoms[1]
        a = acceptor.atoms[0]
        d_label = (f"{donor.residue_name}{donor.residue_id}:"
                   f"{c.atom_name}({c.atom_global_idx})-"
                   f"{x.atom_name}({x.atom_global_idx})")
        a_label = (f"{acceptor.residue_name}{acceptor.residue_id}:"
                   f"{a.atom_name}({a.atom_global_idx})")
        return f"{d_label}···{a_label}"
