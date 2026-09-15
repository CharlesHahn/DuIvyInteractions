# -*- coding: utf-8 -*-
"""Pipeline：串联 Reader → Identifier → Detector → h5 保存。"""

import os
from typing import List

import MDAnalysis as mda

from .system_readers import GmxTprReader
from .group_identifiers import IDENTIFIER_CLASSES
from .group_identifiers.amber_ff_identifier import WATER_RESIDUES
from .interaction_detectors import (
    HydrogenBondDetectorTwoPass, HydrogenBondDetectorPerFrame, HydrogenBondDetectorPerTuple,
    PiStackingDetectorTwoPass, PiStackingDetectorPerFrame, PiStackingDetectorPerTuple,
    SaltBridgeDetectorTwoPass, SaltBridgeDetectorPerFrame, SaltBridgeDetectorPerTuple,
    HydrophobicDetectorTwoPass, HydrophobicDetectorPerFrame, HydrophobicDetectorPerTuple,
    HalogenBondDetectorTwoPass, HalogenBondDetectorPerFrame, HalogenBondDetectorPerTuple,
    MetalCoordinationDetectorTwoPass, MetalCoordinationDetectorPerFrame, MetalCoordinationDetectorPerTuple,
    WaterBridgeDetectorTwoPass, WaterBridgeDetectorPerFrame, WaterBridgeDetectorPerTuple,
    PiCationDetectorTwoPass, PiCationDetectorPerFrame, PiCationDetectorPerTuple,
)
from .io.h5 import save_interactions


ALL_INTERACTIONS = (
    "hydrogen_bond", "pi_stacking", "salt_bridge", "hydrophobic",
    "halogen_bond", "metal_coordination", "water_bridge", "pi_cation",
)

# 类型名 → (TwoPass, PerFrame, PerTuple)
DETECTOR_CLASSES = {
    "hydrogen_bond": (HydrogenBondDetectorTwoPass,
                      HydrogenBondDetectorPerFrame,
                      HydrogenBondDetectorPerTuple),
    "pi_stacking": (PiStackingDetectorTwoPass,
                    PiStackingDetectorPerFrame,
                    PiStackingDetectorPerTuple),
    "salt_bridge": (SaltBridgeDetectorTwoPass,
                    SaltBridgeDetectorPerFrame,
                    SaltBridgeDetectorPerTuple),
    "hydrophobic": (HydrophobicDetectorTwoPass,
                    HydrophobicDetectorPerFrame,
                    HydrophobicDetectorPerTuple),
    "halogen_bond": (HalogenBondDetectorTwoPass,
                     HalogenBondDetectorPerFrame,
                     HalogenBondDetectorPerTuple),
    "metal_coordination": (MetalCoordinationDetectorTwoPass,
                           MetalCoordinationDetectorPerFrame,
                           MetalCoordinationDetectorPerTuple),
    "water_bridge": (WaterBridgeDetectorTwoPass,
                     WaterBridgeDetectorPerFrame,
                     WaterBridgeDetectorPerTuple),
    "pi_cation": (PiCationDetectorTwoPass,
                  PiCationDetectorPerFrame,
                  PiCationDetectorPerTuple),
}

STRATEGY_INDEX = {"two_pass": 0, "per_frame": 1, "per_tuple": 2}


class Pipeline:
    """串联 Reader → Identifier → Detector → h5 保存。"""

    def __init__(self, ff: str, strategy: str = "two_pass"):
        self.ff = ff
        self.strategy = strategy

    def run(self, tpr: str, xtc: str, output: str,
            interactions: List[str] = None) -> None:
        """运行相互作用检测并保存 h5。

        Args:
            tpr: 拓扑文件路径
            xtc: 轨迹文件路径
            output: 输出目录
            interactions: 要检测的相互作用类型列表，None 表示全部
        """
        # 1. 读取 + 识别（只做一次）
        sd = GmxTprReader().read(tpr)
        identifier = self._make_identifier()
        groups = identifier.identify(sd)
        # 2. 加载轨迹（只做一次）
        u = mda.Universe(tpr, xtc)
        os.makedirs(output, exist_ok=True)
        # 3. 逐类型检测 + 保存（单个失败不中断其余）
        names = ALL_INTERACTIONS if interactions is None else interactions
        for name in names:
            try:
                detector = self._make_detector(name)
                filtered = self._filter_groups(groups, detector)
                results = detector.detect(filtered, trajectory=u.trajectory)
                save_interactions(results, os.path.join(output, f"{name}.h5"))
            except Exception as e:
                print(f"[WARN] {name} 检测失败: {e}")

    def _make_identifier(self):
        """按力场构造基团识别器。"""
        if self.ff not in IDENTIFIER_CLASSES:
            raise ValueError(
                f"未知力场: '{self.ff}'。可用: {', '.join(IDENTIFIER_CLASSES)}")
        return IDENTIFIER_CLASSES[self.ff]()

    def _make_detector(self, name: str):
        """按类型和策略构造检测器。"""
        cls = DETECTOR_CLASSES[name][STRATEGY_INDEX[self.strategy]]
        return cls()

    @staticmethod
    def _filter_groups(groups, detector) -> List:
        """按 required_group_types 过滤；除水桥外排除水分子。"""
        needs_water = "water" in detector.required_group_types
        return [g for g in groups
                if g.group_type in detector.required_group_types
                and (needs_water or g.residue_name not in WATER_RESIDUES)]
