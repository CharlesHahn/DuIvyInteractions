# -*- coding: utf-8 -*-
"""π-阳离子 PerFrame 检测器集成测试。

验证 PiCationDetectorPerFrame 的结果正确性。
"""

import pytest
from pathlib import Path

import numpy as np
import MDAnalysis as mda

from DuIvyInteractions.system_readers import GmxTprReader
from DuIvyInteractions.group_identifiers import AmberFFGroupIdentifier
from DuIvyInteractions.interaction_detectors import PiCationDetectorPerFrame


TPR_FILE = Path(__file__).parent.parent / "test_MD_case" / "md.tpr"
XTC_FILE = Path(__file__).parent.parent / "test_MD_case" / "md1ns.xtc"

MOL_D927 = "seg_1_D927"


@pytest.fixture(scope="module")
def system_data():
    reader = GmxTprReader()
    return reader.read(str(TPR_FILE))


@pytest.fixture(scope="module")
def groups(system_data):
    return AmberFFGroupIdentifier().identify(system_data)


@pytest.fixture(scope="module")
def relevant_groups(groups):
    return [g for g in groups
            if g.group_type in ("aromatic_ring", "charged_positive")]


@pytest.fixture(scope="module")
def pi_cations(relevant_groups):
    u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    detector = PiCationDetectorPerFrame()
    return detector.detect(relevant_groups, trajectory=u.trajectory)


# ============================================================
# 基本功能
# ============================================================

class TestPiCationPerFrameBasic:

    def test_returns_list(self, pi_cations):
        assert isinstance(pi_cations, list)

    def test_interaction_type(self, pi_cations):
        for it in pi_cations:
            assert it.interaction_type == "pi_cation"

    def test_metrics_keys(self, pi_cations):
        for it in pi_cations:
            assert "distance" in it.metrics
            assert "offset" in it.metrics


# ============================================================
# 指标值范围
# ============================================================

class TestPiCationPerFrameMetrics:

    def test_distance_positive(self, pi_cations):
        for it in pi_cations:
            assert np.all(it.metrics["distance"] >= 0)

    def test_offset_positive(self, pi_cations):
        for it in pi_cations:
            assert np.all(it.metrics["offset"] >= 0)

    def test_metrics_shape(self, pi_cations):
        for it in pi_cations:
            assert it.metrics["distance"].shape == (it.n_pairs, it.n_frames)
            assert it.metrics["offset"].shape == (it.n_pairs, it.n_frames)


# ============================================================
# 活跃帧值范围
# ============================================================

class TestPiCationPerFrameActiveFrames:

    def test_active_distance_range(self, pi_cations):
        for it in pi_cations:
            for i in range(it.n_pairs):
                active = it.metrics["distance"][i][it.existence[i]]
                if len(active) > 0:
                    assert np.all(active < 6.0)

    def test_active_offset_range(self, pi_cations):
        for it in pi_cations:
            for i in range(it.n_pairs):
                active = it.metrics["offset"][i][it.existence[i]]
                if len(active) > 0:
                    assert np.all(active < 2.0)


# ============================================================
# 候选对
# ============================================================

class TestPiCationPerFramePairs:

    def test_pairs_have_correct_types(self, pi_cations):
        for it in pi_cations:
            for ring, pos in it.groups:
                assert ring.group_type == "aromatic_ring"
                assert pos.group_type == "charged_positive"


# ============================================================
# 结果数据
# ============================================================

class TestPiCationPerFrameResults:

    def test_has_results(self, pi_cations):
        assert len(pi_cations) > 0

    def test_n_pairs(self, pi_cations):
        assert pi_cations[0].n_pairs == 8

    def test_top_pair(self, pi_cations):
        it = pi_cations[0]
        occ = it.occupancy()
        top = np.argmax(occ)
        r, p = it.groups[top]
        assert r.residue_name == "TRP" and r.residue_id == 38
        assert p.residue_name == "LYS" and p.residue_id == 47
        assert occ[top] > 0.6

    def test_d927_participation(self, pi_cations):
        it = pi_cations[0]
        d927 = sum(1 for r, p in it.groups
                   if r.molecule == MOL_D927 or p.molecule == MOL_D927)
        assert d927 >= 1

    def test_lys_dominant(self, pi_cations):
        it = pi_cations[0]
        lys = sum(1 for _, p in it.groups if p.residue_name == "LYS")
        arg = sum(1 for _, p in it.groups if p.residue_name == "ARG")
        assert lys > arg

    def test_occupancy_range(self, pi_cations):
        occ = pi_cations[0].occupancy()
        assert np.all(occ >= 0) and np.all(occ <= 1)
