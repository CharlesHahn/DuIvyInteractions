# -*- coding: utf-8 -*-
"""疏水-芳香重叠标记：判定疏水对-帧是否落在 π-π 环对内。"""

from typing import Dict, List, Optional, Set, Tuple

import numpy as np

from .core.datas import Interaction


FLAG_PI_OVERLAP = "redundant_pi_stacking"
HYDROPHOBIC = "hydrophobic"
PI_STACKING = "pi_stacking"


def pi_overlap_flags(hydrophobic: Interaction,
                     pi_stacking: Optional[Interaction]) -> np.ndarray:
    """返回疏水对-帧与 π-π 环对重叠的 (n_pairs, n_frames) bool 矩阵。"""
    _validate_inputs(hydrophobic, pi_stacking)
    flags = np.zeros_like(hydrophobic.existence, dtype=bool)
    if pi_stacking is None or pi_stacking.n_pairs == 0:
        return flags

    atom_rings = _build_atom_rings(pi_stacking)
    pair_rings = _hydrophobic_ring_sets(hydrophobic, atom_rings)
    for frame in range(hydrophobic.n_frames):
        active = _active_ring_pairs(pi_stacking, frame)
        if not active:
            continue
        for pair_idx in np.flatnonzero(hydrophobic.existence[:, frame]):
            if _pair_hits(active, pair_rings[pair_idx]):
                flags[pair_idx, frame] = True
    return flags & hydrophobic.existence


def _validate_inputs(hydrophobic: Interaction,
                     pi_stacking: Optional[Interaction]) -> None:
    """校验类型、成对结构与帧/时间轴一致性。"""
    if hydrophobic.interaction_type != HYDROPHOBIC:
        raise ValueError(
            f"expected hydrophobic interaction, got '{hydrophobic.interaction_type}'")
    if any(len(pair) != 2 for pair in hydrophobic.groups):
        raise ValueError("hydrophobic groups must contain 2 groups per pair")
    if pi_stacking is None:
        return
    if pi_stacking.interaction_type != PI_STACKING:
        raise ValueError(
            f"expected pi_stacking interaction, got '{pi_stacking.interaction_type}'")
    if any(len(pair) != 2 for pair in pi_stacking.groups):
        raise ValueError("pi_stacking groups must contain 2 groups per pair")
    if hydrophobic.n_frames != pi_stacking.n_frames:
        raise ValueError(
            f"n_frames mismatch: {hydrophobic.n_frames} != {pi_stacking.n_frames}")
    if not np.allclose(hydrophobic.times, pi_stacking.times):
        raise ValueError("times mismatch between hydrophobic and pi_stacking")


def _build_atom_rings(pi_stacking: Interaction) -> Dict[int, Set[int]]:
    """建立 原子全局索引 → 所属芳香环 group_id 集合。"""
    atom_rings: Dict[int, Set[int]] = {}
    seen_rings: Set[int] = set()
    for pair in pi_stacking.groups:
        for ring in pair:
            if ring.group_id in seen_rings:
                continue
            seen_rings.add(ring.group_id)
            for atom_idx in ring.atom_indices:
                atom_rings.setdefault(atom_idx, set()).add(ring.group_id)
    return atom_rings


def _hydrophobic_ring_sets(hydrophobic: Interaction,
                           atom_rings: Dict[int, Set[int]]
                           ) -> List[Tuple[Set[int], Set[int]]]:
    """为每个疏水对记录两端原子所属的芳香环集合。"""
    result: List[Tuple[Set[int], Set[int]]] = []
    for group1, group2 in hydrophobic.groups:
        atom1 = group1.atoms[0].atom_global_idx
        atom2 = group2.atoms[0].atom_global_idx
        result.append((atom_rings.get(atom1, set()),
                       atom_rings.get(atom2, set())))
    return result


def _active_ring_pairs(pi_stacking: Interaction,
                       frame: int) -> Set[Tuple[int, int]]:
    """返回该帧活跃 π-π 环对的无序 group_id 键集合。"""
    active: Set[Tuple[int, int]] = set()
    for pair_idx in np.flatnonzero(pi_stacking.existence[:, frame]):
        ring1, ring2 = pi_stacking.groups[pair_idx]
        active.add(_ring_pair_key(ring1.group_id, ring2.group_id))
    return active


def _pair_hits(active: Set[Tuple[int, int]],
               rings: Tuple[Set[int], Set[int]]) -> bool:
    """判断疏水对两端是否分属某个活跃 π-π 环对。"""
    rings1, rings2 = rings
    for ring1 in rings1:
        for ring2 in rings2:
            if ring1 == ring2:
                continue
            if _ring_pair_key(ring1, ring2) in active:
                return True
    return False


def _ring_pair_key(ring1: int, ring2: int) -> Tuple[int, int]:
    """将环对规范化为 (较小 group_id, 较大 group_id)。"""
    return (ring1, ring2) if ring1 < ring2 else (ring2, ring1)
