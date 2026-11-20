# -*- coding: utf-8 -*-
"""io - 结果文件读写（Interaction 数据的序列化/反序列化）。"""

from .h5 import save_interactions, save_interactions_atomic, load_interactions
from .interaction_exporter import InteractionExporter
from .hydrogen_bond_exporter import HydrogenBondExporter
from .saltbridge_exporter import SaltBridgeExporter
from .pi_stacking_exporter import PiStackingExporter
from .pi_cation_exporter import PiCationExporter
from .halogen_bond_exporter import HalogenBondExporter
from .hydrophobic_exporter import HydrophobicExporter
from .metal_coordination_exporter import MetalCoordinationExporter
from .water_bridge_exporter import WaterBridgeExporter

# interaction_type（检测器名）→ Exporter 子类（dii export 与 web 共用的单一来源）
EXPORTER_CLASSES = {
    "hydrogen_bond": HydrogenBondExporter,
    "pi_stacking": PiStackingExporter,
    "salt_bridge": SaltBridgeExporter,
    "hydrophobic": HydrophobicExporter,
    "halogen_bond": HalogenBondExporter,
    "metal_coordination": MetalCoordinationExporter,
    "water_bridge": WaterBridgeExporter,
    "pi_cation": PiCationExporter,
}

__all__ = [
    "save_interactions", "save_interactions_atomic", "load_interactions",
    "InteractionExporter", "HydrogenBondExporter", "SaltBridgeExporter",
    "PiStackingExporter", "PiCationExporter", "HalogenBondExporter",
    "HydrophobicExporter", "MetalCoordinationExporter", "WaterBridgeExporter",
    "EXPORTER_CLASSES",
]
