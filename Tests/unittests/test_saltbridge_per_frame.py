# -*- coding: utf-8 -*-
"""盐桥 PerFrame 检测器集成测试。

验证 SaltBridgeDetectorPerFrame 的结果正确性。
"""

import pytest
from pathlib import Path

import numpy as np
import MDAnalysis as mda

from DuIvyInteractions.system_readers import GmxTprReader
from DuIvyInteractions.group_identifiers import AmberFFGroupIdentifier
from DuIvyInteractions.interaction_detectors import (
    SaltBridgeDetectorPerTuple, SaltBridgeDetectorPerFrame
)


TPR_FILE = Path(__file__).parent.parent / "test_MD_case" / "md.tpr"
XTC_FILE = Path(__file__).parent.parent / "test_MD_case" / "md1ns.xtc"


@pytest.fixture(scope="module")
def system_data():
    reader = GmxTprReader()
    return reader.read(str(TPR_FILE))


@pytest.fixture(scope="module")
def groups(system_data):
    return AmberFFGroupIdentifier().identify(system_data)


@pytest.fixture(scope="module")
def charged_groups(groups):
    return [g for g in groups
            if g.group_type in ("charged_positive", "charged_negative")]


@pytest.fixture(scope="module")
def salt_bridges(charged_groups):
    u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    detector = SaltBridgeDetectorPerFrame()
    return detector.detect(charged_groups, trajectory=u.trajectory)


@pytest.fixture(scope="module")
def salt_bridges_per_tuple(charged_groups):
    u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    detector = SaltBridgeDetectorPerTuple()
    return detector.detect(charged_groups, trajectory=u.trajectory)


# ============================================================
# 基本功能
# ============================================================

class TestSaltBridgePerFrameBasic:

    def test_has_results(self, salt_bridges):
        assert len(salt_bridges) > 0

    def test_interaction_type(self, salt_bridges):
        assert salt_bridges[0].interaction_type == "salt_bridge"

    def test_n_pairs(self, salt_bridges):
        assert salt_bridges[0].n_pairs == 47

    def test_all_opposite_charge(self, salt_bridges):
        it = salt_bridges[0]
        for g1, g2 in it.groups:
            assert g1.group_type == "charged_positive"
            assert g2.group_type == "charged_negative"

    def test_occupancy_range(self, salt_bridges):
        occ = salt_bridges[0].occupancy()
        assert np.all(occ >= 0) and np.all(occ <= 1)

    def test_top_pair_occupancy(self, salt_bridges):
        it = salt_bridges[0]
        occ = it.occupancy()
        assert np.max(occ) == 1.0

    def test_metrics_shape(self, salt_bridges):
        it = salt_bridges[0]
        assert it.metrics["distance"].shape == (it.n_pairs, it.n_frames)

    def test_distance_range(self, salt_bridges):
        it = salt_bridges[0]
        for i in range(it.n_pairs):
            active = it.metrics["distance"][i][it.existence[i]]
            if len(active) > 0:
                assert np.all(active <= 5.5)

    def test_residue_coverage(self, salt_bridges):
        names = set()
        for g1, g2 in salt_bridges[0].groups:
            names.update([g1.residue_name, g2.residue_name])
        assert {"ARG", "LYS", "ASP", "GLU"}.issubset(names)


# ============================================================
# 与 PerTuple 版本交叉验证
# ============================================================

class TestSaltBridgePerFrameVsPerTuple:

    def test_n_pairs_match(self, salt_bridges, salt_bridges_per_tuple):
        assert salt_bridges[0].n_pairs == salt_bridges_per_tuple[0].n_pairs

    def test_pairs_match(self, salt_bridges, salt_bridges_per_tuple):
        fp = {(g1.group_id, g2.group_id) for g1, g2 in salt_bridges[0].groups}
        tp = {(g1.group_id, g2.group_id) for g1, g2 in salt_bridges_per_tuple[0].groups}
        assert fp == tp

    def test_existence_match(self, salt_bridges, salt_bridges_per_tuple):
        iff, it = salt_bridges[0], salt_bridges_per_tuple[0]
        fo = {(g1.group_id, g2.group_id): i for i, (g1, g2) in enumerate(iff.groups)}
        to = {(g1.group_id, g2.group_id): i for i, (g1, g2) in enumerate(it.groups)}
        for k in fo:
            assert np.array_equal(iff.existence[fo[k]], it.existence[to[k]])

    def test_distance_match(self, salt_bridges, salt_bridges_per_tuple):
        iff, it = salt_bridges[0], salt_bridges_per_tuple[0]
        fo = {(g1.group_id, g2.group_id): i for i, (g1, g2) in enumerate(iff.groups)}
        to = {(g1.group_id, g2.group_id): i for i, (g1, g2) in enumerate(it.groups)}
        for k in fo:
            assert np.allclose(iff.metrics["distance"][fo[k]],
                               it.metrics["distance"][to[k]], atol=1e-6)
