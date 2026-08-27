# -*- coding: utf-8 -*-
"""core - 核心基础：数据类 + 接口定义 + 常量。"""

from .constants import GROUP_TYPES, BOND_TYPES
from .datas import (
    Group, AtomData, BondData, ResidueData, InterResidueBond,
    SystemData, Interaction, HydrogenBond, PiStacking, WaterBridge
)
from .interfaces import Reader, GroupIdentifier, InteractionDetector

__all__ = [
    # 常量
    "GROUP_TYPES", "BOND_TYPES",
    # 数据类
    "Group", "AtomData", "BondData", "ResidueData", "InterResidueBond",
    "SystemData", "Interaction", "HydrogenBond", "PiStacking", "WaterBridge",
    # 接口
    "Reader", "GroupIdentifier", "InteractionDetector",
]
