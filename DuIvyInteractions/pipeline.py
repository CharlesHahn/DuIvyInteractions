# -*- coding: utf-8 -*-
"""Pipeline：串联 Reader → Identifier → Detector → h5 保存。"""

import os
from typing import List, Optional, Tuple

import MDAnalysis as mda

from .core.datas import Group
from .system_readers import GmxTprReader
from .group_identifiers import IDENTIFIER_CLASSES
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
from .selection import (
    CROSS_MIN_SIDES, Assignment,
    resolve_selection, validate_selections, assign_groups,
    filter_groups, build_tuple_filter,
)


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
            interactions: List[str] = None,
            select: List[str] = None, pairing: str = "any") -> None:
        """运行相互作用检测并保存 h5（--select/--pairing 限定参与池与配对）。"""
        sd = GmxTprReader().read(tpr)
        identifier = self._make_identifier()
        groups = identifier.identify(sd)
        water_residues = identifier.WATER_RESIDUES
        u = mda.Universe(tpr, xtc)
        os.makedirs(output, exist_ok=True)
        groups, assignment = self._apply_selection(
            groups, u, water_residues, select, pairing)
        names = ALL_INTERACTIONS if interactions is None else interactions
        for name in names:
            self._detect_and_save(name, groups, u, water_residues,
                                  assignment, pairing, output)

    def _detect_and_save(self, name: str, groups: List[Group], u: mda.Universe,
                         water_residues: frozenset,
                         assignment: Optional[Assignment], pairing: str,
                         output: str) -> None:
        """单类型检测并保存 h5；失败只告警，不中断其余类型。"""
        try:
            detector = self._make_detector(name)
            filtered = self._filter_groups(groups, detector, water_residues)
            tuple_filter = None
            if assignment is not None:
                tuple_filter = build_tuple_filter(pairing, assignment, name)
            results = detector.detect(filtered, trajectory=u.trajectory,
                                      tuple_filter=tuple_filter)
            save_interactions(results, os.path.join(output, f"{name}.h5"))
        except Exception as e:
            print(f"[WARN] {name} detection failed: {e}")

    @staticmethod
    def _validate_pairing(select: Optional[List[str]], pairing: str) -> None:
        """between 需要至少 2 个 select。"""
        if pairing == "between" and len(select or []) < CROSS_MIN_SIDES:
            raise ValueError("--pairing between 需要至少 2 个 --select")

    def _apply_selection(self, groups: List[Group], u: mda.Universe,
                         water_residues: frozenset,
                         select: Optional[List[str]], pairing: str
                         ) -> Tuple[List[Group], Optional[Assignment]]:
        """应用 --select/--pairing 过滤；select 为空则返回 (原基团, None)。"""
        self._validate_pairing(select, pairing)
        if not select:
            return groups, None
        u.trajectory[0]  # 几何选择按第一帧静态求值
        sets = [resolve_selection(u, s) for s in select]
        validate_selections(sets)
        assignment = assign_groups(groups, sets, water_residues)
        kept = filter_groups(groups, assignment)
        if not kept:
            raise ValueError(
                "参与池内无基团（0 个基团保留）。若使用了 around/sphzone "
                "几何选择，基团可能因部分原子在圈外被排除，请改用 byres around …")
        return kept, assignment

    def _make_identifier(self):
        """按力场构造基团识别器。"""
        if self.ff not in IDENTIFIER_CLASSES:
            raise ValueError(
                f"unknown force field: '{self.ff}'. Available: {', '.join(IDENTIFIER_CLASSES)}")
        return IDENTIFIER_CLASSES[self.ff]()

    def _make_detector(self, name: str):
        """按类型和策略构造检测器。"""
        cls = DETECTOR_CLASSES[name][STRATEGY_INDEX[self.strategy]]
        return cls()

    @staticmethod
    def _filter_groups(groups, detector, water_residues) -> List:
        """按 required_group_types 过滤；除水桥外排除水分子。

        water_residues: 当前力场的水残基名集合（来自 identifier.WATER_RESIDUES），
        跨力场正确排除（Amber SOL / CHARMM TIP3 / OPLS HO4/HO5）。
        """
        needs_water = "water" in detector.required_group_types
        return [g for g in groups
                if g.group_type in detector.required_group_types
                and (needs_water or g.residue_name not in water_residues)]
