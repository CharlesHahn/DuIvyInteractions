# -*- coding: utf-8 -*-
"""core - 核心基础：数据类 + 接口定义 + 常量。"""

from .constants import GROUP_TYPES, BOND_TYPES
from .datas import (
    Group, AtomData, BondData, ResidueData, InterResidueBond,
    SystemData, Interaction
)
from .interfaces import Reader, GroupIdentifier, InteractionDetectorPerTuple

__all__ = [
    # 常量
    "GROUP_TYPES", "BOND_TYPES",
    # 数据类
    "Group", "AtomData", "BondData", "ResidueData", "InterResidueBond",
    "SystemData", "Interaction",
    # 接口
    "Reader", "GroupIdentifier", "InteractionDetectorPerTuple",
]
