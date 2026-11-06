# -*- coding: utf-8 -*-
"""水桥导出器。"""

from typing import Dict

from .interaction_exporter import InteractionExporter
from ..core.datas import Interaction


class WaterBridgeExporter(InteractionExporter):
    """水桥导出器。"""

    @property
    def name(self) -> str:
        return "Water Bridge"

    @property
    def metric_labels(self) -> Dict[str, str]:
        return {
            "dist_h_ow": "H···Ow Distance (Å)",
            "angle_dh_ow": "D-H···Ow Angle (°)",
            "dist_d_ow": "D···Ow Distance (Å)",
            "dist_hw_a": "Hw···A Distance (Å)",
            "angle_ow_hw_a": "Ow-Hw···A Angle (°)",
            "dist_ow_a": "Ow···A Distance (Å)",
            "angle_h_ow_a": "H(donor)-Ow···A Angle (°)",
        }

    def get_pair_label(self, interaction: Interaction, pair_idx: int) -> str:
        """生成水桥标签。

        格式：供体(D-H)···水(O)···受体(A)
        """
        donor, water, acceptor = interaction.groups[pair_idx]
        d, h = donor.atoms[0], donor.atoms[1]
        ow = water.atoms[0]
        a = acceptor.atoms[0]
        d_label = (f"{donor.residue_name}{donor.residue_id}:"
                   f"{d.atom_name}({d.atom_global_idx})-"
                   f"{h.atom_name}({h.atom_global_idx})")
        w_label = (f"{water.residue_name}{water.residue_id}:"
                   f"{ow.atom_name}({ow.atom_global_idx})")
        a_label = (f"{acceptor.residue_name}{acceptor.residue_id}:"
                   f"{a.atom_name}({a.atom_global_idx})")
        return f"{d_label}···{w_label}···{a_label}"
