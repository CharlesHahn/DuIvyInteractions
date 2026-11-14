# -*- coding: utf-8 -*-
"""疏水-芳香重叠标记的单元测试与真实数据回归。"""

from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pytest

from DuIvyInteractions.core.datas import AtomData, Group, Interaction
from DuIvyInteractions.group_identifiers import AmberFFGroupIdentifier
from DuIvyInteractions.interaction_detectors import (
    HydrophobicDetectorPerFrame, HydrophobicDetectorTwoPass,
    PiStackingDetectorPerFrame, PiStackingDetectorTwoPass,
)
from DuIvyInteractions.redundancy import pi_overlap_flags
from DuIvyInteractions.system_readers import GmxTprReader


TPR_FILE = Path(__file__).parent.parent / "test_MD_case_amber" / "md.tpr"
XTC_FILE = Path(__file__).parent.parent / "test_MD_case_amber" / "md1ns.xtc"


def _atom(atom_idx: int) -> AtomData:
    return AtomData(
        atom_global_idx=atom_idx, atom_idx_in_residue=0,
        atom_name="C", atom_type="ca", atom_element="C",
        atom_charge=0.0, atom_mass=12.0,
    )


def _ring(group_id: int, atom_indices: list) -> Group:
    return Group(
        group_id=group_id, group_type="aromatic_ring", molecule="mol",
        residue_name="PHE", residue_id=group_id,
        atoms=[_atom(i) for i in atom_indices],
    )


def _hydrophobic(group_id: int, atom_idx: int) -> Group:
    return Group(
        group_id=group_id, group_type="hydrophobic", molecule="mol",
        residue_name="PHE", residue_id=group_id,
        atoms=[_atom(atom_idx)],
    )


def _interaction(interaction_type: str, pairs: list,
                 existence: np.ndarray) -> Interaction:
    n_frames = existence.shape[1]
    return Interaction(
        interaction_type=interaction_type,
        groups=pairs,
        existence=existence,
        metrics={},
        times=np.arange(n_frames, dtype=float) * 10.0,
    )


@pytest.fixture
def synthetic():
    """三环体系：ring1-ring2 在帧 0/2 形成 π-π，ring3 不参与。"""
    ring1 = _ring(10, [1, 2, 3])
    ring2 = _ring(11, [4, 5, 6])
    ring3 = _ring(12, [7, 8])
    # 故意反序存储，验证环对无序规范化
    pi_stacking = _interaction(
        "pi_stacking", [(ring2, ring1)],
        np.array([[True, False, True]]),
    )
    pairs = [
        (_hydrophobic(20, 1), _hydrophobic(21, 4)),   # 环1-环2：应标记
        (_hydrophobic(22, 1), _hydrophobic(23, 99)),  # 环1-非环：不标记
        (_hydrophobic(24, 1), _hydrophobic(25, 2)),   # 环1内部：不标记
        (_hydrophobic(26, 7), _hydrophobic(27, 4)),   # 环3-环2：不标记
        (_hydrophobic(28, 5), _hydrophobic(29, 2)),   # 环2-环1：应标记
    ]
    hydrophobic = _interaction(
        "hydrophobic", pairs, np.ones((len(pairs), 3), dtype=bool))
    return hydrophobic, pi_stacking


class TestPiOverlapFlagsSynthetic:
    """合成体系上的判定语义。"""

    def test_flags_ring_pair_and_frames(self, synthetic):
        hydrophobic, pi_stacking = synthetic
        flags = pi_overlap_flags(hydrophobic, pi_stacking)
        assert flags.shape == (5, 3)
        assert flags.dtype == bool
        np.testing.assert_array_equal(flags[0], [True, False, True])
        np.testing.assert_array_equal(flags[1], [False, False, False])
        np.testing.assert_array_equal(flags[2], [False, False, False])
        np.testing.assert_array_equal(flags[3], [False, False, False])
        np.testing.assert_array_equal(flags[4], [True, False, True])

    def test_flags_subset_existence(self, synthetic):
        hydrophobic, pi_stacking = synthetic
        hydrophobic.existence[0, 1] = False
        flags = pi_overlap_flags(hydrophobic, pi_stacking)
        assert not np.any(flags & ~hydrophobic.existence)
        assert flags[0].tolist() == [True, False, True]

    def test_none_pi_returns_all_false(self, synthetic):
        hydrophobic, _ = synthetic
        flags = pi_overlap_flags(hydrophobic, None)
        assert flags.shape == hydrophobic.existence.shape
        assert not flags.any()

    def test_empty_pi_returns_all_false(self, synthetic):
        hydrophobic, _ = synthetic
        empty = _interaction(
            "pi_stacking", [], np.zeros((0, hydrophobic.n_frames), dtype=bool))
        flags = pi_overlap_flags(hydrophobic, empty)
        assert not flags.any()

    def test_empty_hydrophobic(self):
        hydrophobic = _interaction(
            "hydrophobic", [], np.zeros((0, 2), dtype=bool))
        flags = pi_overlap_flags(hydrophobic, None)
        assert flags.shape == (0, 2)

    def test_wrong_interaction_type_raises(self, synthetic):
        hydrophobic, _ = synthetic
        with pytest.raises(ValueError, match="expected hydrophobic"):
            pi_overlap_flags(
                _interaction("salt_bridge", hydrophobic.groups,
                             hydrophobic.existence),
                None)

    def test_n_frames_mismatch_raises(self, synthetic):
        hydrophobic, pi_stacking = synthetic
        with pytest.raises(ValueError, match="n_frames mismatch"):
            pi_overlap_flags(hydrophobic,
                             _interaction("pi_stacking", pi_stacking.groups,
                                          np.array([[True, False]])))

    def test_times_mismatch_raises(self, synthetic):
        hydrophobic, pi_stacking = synthetic
        pi_stacking.times = pi_stacking.times + 1.0
        with pytest.raises(ValueError, match="times mismatch"):
            pi_overlap_flags(hydrophobic, pi_stacking)


class TestFusedRingSemantics:
    """稠合环共享原子：只匹配 π-π 结果中实际配对的环。"""

    def test_shared_atom_hits_paired_ring(self):
        ring_a = _ring(20, [1, 2, 3])
        ring_b = _ring(21, [2, 3, 4])
        ring_c = _ring(22, [5, 6])
        ring_d = _ring(23, [9, 10])
        # ring_b 只出现在未活跃的 π-π 对中，用于注册共享原子 2/3 的第二个环
        pi_stacking = _interaction(
            "pi_stacking", [(ring_a, ring_c), (ring_b, ring_d)],
            np.array([[True], [False]]))
        pairs = [
            (_hydrophobic(30, 2), _hydrophobic(31, 5)),  # 共享原子 2 → 命中 ring_a
            (_hydrophobic(32, 1), _hydrophobic(33, 5)),  # ring_a 独有 → 命中
            (_hydrophobic(34, 4), _hydrophobic(35, 5)),  # ring_b 独有 → 不命中
        ]
        hydrophobic = _interaction(
            "hydrophobic", pairs, np.ones((3, 1), dtype=bool))
        flags = pi_overlap_flags(hydrophobic, pi_stacking)
        assert flags[:, 0].tolist() == [True, True, False]


@pytest.fixture(scope="module")
def two_pass_result():
    system_data = GmxTprReader().read(str(TPR_FILE))
    groups = AmberFFGroupIdentifier().identify(system_data)
    hydrophobic = [g for g in groups if g.group_type == "hydrophobic"]
    rings = [g for g in groups if g.group_type == "aromatic_ring"]
    universe = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    hydro = HydrophobicDetectorTwoPass().detect(
        hydrophobic, trajectory=universe.trajectory)[0]
    universe = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    pi = PiStackingDetectorTwoPass().detect(
        rings, trajectory=universe.trajectory)[0]
    return hydro, pi


@pytest.fixture(scope="module")
def per_frame_result():
    system_data = GmxTprReader().read(str(TPR_FILE))
    groups = AmberFFGroupIdentifier().identify(system_data)
    hydrophobic = [g for g in groups if g.group_type == "hydrophobic"]
    rings = [g for g in groups if g.group_type == "aromatic_ring"]
    universe = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    hydro = HydrophobicDetectorPerFrame().detect(
        hydrophobic, trajectory=universe.trajectory)[0]
    universe = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    pi = PiStackingDetectorPerFrame().detect(
        rings, trajectory=universe.trajectory)[0]
    return hydro, pi


class TestRealDataBaselines:
    """D927 真实数据基线（策略相关，见 doc/redundancy_marking_design.md）。"""

    def test_two_pass_counts(self, two_pass_result):
        hydro, pi = two_pass_result
        flags = pi_overlap_flags(hydro, pi)
        assert hydro.n_pairs == 2366
        assert pi.n_pairs == 10
        assert int(flags.any(axis=1).sum()) == 67
        assert int(flags.sum()) == 749

    def test_per_frame_counts(self, per_frame_result):
        """per_frame 的疏水同残基去重会改变基线，不能套 two_pass 数字。"""
        hydro, pi = per_frame_result
        flags = pi_overlap_flags(hydro, pi)
        assert hydro.n_pairs == 705
        assert int(flags.any(axis=1).sum()) == 9
        assert int(flags.sum()) == 143

    def test_flags_only_on_ring_ring_pairs(self, two_pass_result):
        hydro, pi = two_pass_result
        flags = pi_overlap_flags(hydro, pi)
        ring_atoms = {
            atom_idx
            for pair in pi.groups
            for ring in pair
            for atom_idx in ring.atom_indices
        }
        for pair_idx in np.flatnonzero(flags.any(axis=1)):
            atom1 = hydro.groups[pair_idx][0].atoms[0].atom_global_idx
            atom2 = hydro.groups[pair_idx][1].atoms[0].atom_global_idx
            assert atom1 in ring_atoms
            assert atom2 in ring_atoms
