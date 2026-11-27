# -*- coding: utf-8 -*-
"""Pipeline：串联 Reader → Identifier → Detector → h5 保存。"""

import os
from datetime import datetime
from typing import List, Optional, Set, Tuple

import MDAnalysis as mda

import DuIvyInteractions

from .core.datas import Group, Interaction
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
from .io.h5 import (
    RunMeta, save_interactions, save_interactions_atomic, load_interactions,
)
from .redundancy import FLAG_PI_OVERLAP, pi_overlap_flags
from .selection import (
    CROSS_MIN_SIDES, WATER_SIDE, Assignment,
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
            select: List[str] = None, pairing: str = "any",
            mark_redundant: bool = False) -> None:
        """运行相互作用检测并保存 h5（--select/--pairing 限定参与池与配对；
        mark_redundant 标记疏水与 π-π 的重叠）。"""
        names = _prepare_run_names(interactions, mark_redundant, self.strategy)
        run_meta = build_run_meta(
            tpr, xtc, output, self.ff, self.strategy, interactions, select,
            pairing, mark_redundant)
        sd = GmxTprReader().read(tpr)
        identifier = self._make_identifier()
        groups = identifier.identify(sd)
        water_residues = identifier.WATER_RESIDUES
        u = mda.Universe(tpr, xtc)
        os.makedirs(output, exist_ok=True)
        groups, assignment = self._apply_selection(groups, u, select, pairing)
        succeeded = {
            name for name in names
            if self._detect_and_save(name, groups, u, water_residues,
                                     assignment, pairing, output, run_meta)
        }
        if mark_redundant:
            _finish_marking(output, succeeded, run_meta)

    def _detect_and_save(self, name: str, groups: List[Group], u: mda.Universe,
                         water_residues: frozenset,
                         assignment: Optional[Assignment], pairing: str,
                         output: str, run: Optional[RunMeta] = None) -> bool:
        """单类型检测并保存 h5；失败只告警，返回是否成功。"""
        try:
            detector = self._make_detector(name)
            filtered = self._filter_groups(groups, detector, water_residues,
                                           assignment)
            tuple_filter = None
            if assignment is not None:
                tuple_filter = build_tuple_filter(pairing, assignment, name)
            results = detector.detect(filtered, trajectory=u.trajectory,
                                      tuple_filter=tuple_filter)
            save_interactions(results, os.path.join(output, f"{name}.h5"),
                              run=run)
            return True
        except Exception as e:
            print(f"[WARN] {name} detection failed: {e}")
            return False

    @staticmethod
    def _validate_pairing(select: Optional[List[str]], pairing: str) -> None:
        """between 需要至少 2 个 select。"""
        if pairing == "between" and len(select or []) < CROSS_MIN_SIDES:
            raise ValueError("--pairing between 需要至少 2 个 --select")

    def _apply_selection(self, groups: List[Group], u: mda.Universe,
                         select: Optional[List[str]], pairing: str
                         ) -> Tuple[List[Group], Optional[Assignment]]:
        """应用 --select/--pairing 过滤；select 为空则返回 (原基团, None)。"""
        self._validate_pairing(select, pairing)
        if not select:
            return groups, None
        u.trajectory[0]  # 几何选择按第一帧静态求值
        sets = [resolve_selection(u, s) for s in select]
        validate_selections(sets)
        assignment = assign_groups(groups, sets)
        kept = filter_groups(groups, assignment)
        # 空池检查基于"非水基团"：水基团恒真保留（桥接水），
        # 若用户选择未命中任何真实化学基团（如 name CA 只选到裸原子），
        # 应报错而非静默产出空结果。
        kept_non_water = [g for g in kept if assignment[g.group_id] != WATER_SIDE]
        if not kept_non_water:
            raise ValueError(
                "参与池内无可用基团（非水基团 0 个）。若使用了 around/sphzone "
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
    def _filter_groups(groups, detector, water_residues,
                       assignment: Optional[Assignment] = None) -> List:
        """按 required_group_types 过滤；除水桥外排除水残基基团。

        用户显式选择（assignment 非 None）时，**显式归属到具体 select 的
        水残基基团豁免水排除**（如 --select "resname SOL" 让水参与氢键）；
        无 select（assignment=None）时行为与既有完全一致。
        """
        needs_water = "water" in detector.required_group_types
        result = []
        for g in groups:
            if g.group_type not in detector.required_group_types:
                continue
            if needs_water or g.residue_name not in water_residues:
                result.append(g)
                continue
            # 水残基基团：仅显式归属到具体 select 索引（int）时豁免；
            # WATER_SIDE（str）与未归属（None）均不豁免（显式优于隐式）
            side = assignment.get(g.group_id) if assignment is not None else None
            if isinstance(side, int):
                result.append(g)
        return result


def build_run_meta(tpr: str, xtc: str, output: str, ff: str, strategy: str,
                   interactions: Optional[List[str]],
                   select: Optional[List[str]], pairing: str,
                   mark_redundant: bool,
                   argv: Optional[List[str]] = None,
                   cwd: Optional[str] = None,
                   created_at: Optional[str] = None) -> RunMeta:
    """构建一次 dii run 的运行参数（argv/cwd/created_at 缺省时自动生成）。"""
    types = None if interactions is None else _unique(interactions)
    selects = _unique(select)
    tpr, xtc, output = (os.path.abspath(p) for p in (tpr, xtc, output))
    if argv is None:
        argv = _canonical_argv(tpr, xtc, output, ff, strategy, types, selects,
                               pairing, mark_redundant)
    return RunMeta(
        topology=tpr, trajectory=xtc, output=output, ff=ff, strategy=strategy,
        interactions=types or (), select=selects, pairing=pairing,
        mark_redundant=bool(mark_redundant), argv=tuple(argv),
        cwd=cwd if cwd is not None else os.getcwd(),
        created_at=created_at if created_at is not None else _now_iso(),
        dii_version=_dii_version())


def _canonical_argv(tpr: str, xtc: str, output: str, ff: str, strategy: str,
                    interactions: Optional[Tuple[str, ...]],
                    select: Tuple[str, ...], pairing: str,
                    mark_redundant: bool) -> List[str]:
    """生成与 dii run 参数顺序一致的规范命令行。"""
    types = "all" if interactions is None else ",".join(interactions)
    argv = ["dii", "run", "-t", tpr, "-f", xtc, "-o", output, "--ff", ff,
            "--strategy", strategy, "--interactions", types]
    for value in select:
        argv += ["--select", value]
    if pairing != "any":
        argv += ["--pairing", pairing]
    if mark_redundant:
        argv.append("--mark-redundant")
    return argv


def _unique(values: Optional[List[str]]) -> Tuple[str, ...]:
    """去重并保持传入顺序。"""
    return tuple(dict.fromkeys(values or ()))


def _now_iso() -> str:
    """当前本地时间（含时区、秒精度）的 ISO 字符串。"""
    return datetime.now().astimezone().isoformat(timespec="seconds")


def _dii_version() -> str:
    """DuIvyInteractions 包版本号（未定义时为 unknown）。"""
    return str(getattr(DuIvyInteractions, "__version__", "unknown"))


def _prepare_run_names(interactions: Optional[List[str]], mark_redundant: bool,
                       strategy: str) -> List[str]:
    """构建实际运行类型列表；mark 时校验并自动补跑 pi_stacking。"""
    names = list(ALL_INTERACTIONS if interactions is None else interactions)
    if not mark_redundant:
        return names
    if "hydrophobic" not in names:
        raise ValueError("--mark-redundant requires hydrophobic")
    if strategy == "per_tuple":
        raise ValueError("--mark-redundant does not support per_tuple")
    if "pi_stacking" not in names:
        names.append("pi_stacking")
        print("[INFO] auto-running pi_stacking for --mark-redundant")
    return names


def _finish_marking(output: str, succeeded: Set[str],
                    run: Optional[RunMeta] = None) -> None:
    """确认依赖类型成功并写回 flags；依赖失败则报错。"""
    missing = {"hydrophobic", "pi_stacking"} - succeeded
    if missing:
        raise ValueError(
            "mark-redundant skipped: "
            + ", ".join(f"{name} detection failed" for name in sorted(missing)))
    n_pairs, n_frames = _mark_hydrophobic_redundancy(output, run)
    print(f"[INFO] marked {n_pairs} hydrophobic pairs ({n_frames} pair-frames) "
          "overlapping pi-stacking (kept, not removed)")


def _mark_hydrophobic_redundancy(output: str,
                                 run: Optional[RunMeta] = None) -> Tuple[int, int]:
    """读取 hydrophobic/pi_stacking h5，计算并原子写回 flags。"""
    hydro_path = os.path.join(output, "hydrophobic.h5")
    pi_path = os.path.join(output, "pi_stacking.h5")
    hydrophobic, pi_interaction = _load_marking_inputs(hydro_path, pi_path)
    if not hydrophobic:
        return (0, 0)
    counts = _apply_overlap_flags(hydrophobic, pi_interaction)
    _save_marked(hydrophobic, hydro_path, run)
    return counts


def _load_marking_inputs(hydro_path: str, pi_path: str
                         ) -> Tuple[List[Interaction], Optional[Interaction]]:
    """读取两个 h5，并取出 pi_stacking Interaction（无则 None）。"""
    try:
        hydrophobic = load_interactions(hydro_path)
        pi_stacking = load_interactions(pi_path)
    except (OSError, ValueError) as e:
        raise ValueError(f"cannot read marking inputs: {e}") from e
    pi_matches = [
        it for it in pi_stacking if it.interaction_type == "pi_stacking"]
    if pi_stacking and not pi_matches:
        raise ValueError("no pi_stacking interaction in pi_stacking.h5")
    if len(pi_matches) > 1:
        raise ValueError("multiple pi_stacking interactions in pi_stacking.h5")
    return hydrophobic, pi_matches[0] if pi_matches else None


def _apply_overlap_flags(hydrophobic: List[Interaction],
                         pi_interaction: Optional[Interaction]
                         ) -> Tuple[int, int]:
    """为每个 hydrophobic Interaction 写入 flags，返回 (对数, 对-帧数)。"""
    total_pairs = 0
    total_frames = 0
    for interaction in hydrophobic:
        if interaction.interaction_type != "hydrophobic":
            continue
        try:
            flags = pi_overlap_flags(interaction, pi_interaction)
        except ValueError as e:
            raise ValueError(f"mark-redundant failed: {e}") from e
        interaction.flags[FLAG_PI_OVERLAP] = flags
        total_pairs += int(flags.any(axis=1).sum())
        total_frames += int(flags.sum())
    return (total_pairs, total_frames)


def _save_marked(hydrophobic: List[Interaction], hydro_path: str,
                 run: Optional[RunMeta] = None) -> None:
    """原子写回带 flags 的 hydrophobic h5（保留 run/ 溯源组）。"""
    try:
        save_interactions_atomic(hydrophobic, hydro_path, run=run)
    except OSError as e:
        raise ValueError(f"cannot write {hydro_path}: {e}") from e
