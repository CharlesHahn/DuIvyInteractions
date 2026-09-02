# -*- coding: utf-8 -*-
"""π-阳离子 TwoPass 检测器集成测试。

从真实 tpr + xtc 输入，验证检测结果。
与 PerFrame 版本对比，验证结果一致性。
"""

import pytest
from pathlib import Path

import numpy as np

from DuIvyInteractions.input_readers import GmxTprReader
from DuIvyInteractions.group_identifiers import AmberFFGroupIdentifier
from DuIvyInteractions.interaction_detectors import (
    PiCationDetectorTwoPass, PiCationDetectorPerFrame)
import MDAnalysis as mda


TPR_FILE = Path(__file__).parent.parent / "test_MD_case" / "md.tpr"
XTC_FILE = Path(__file__).parent.parent / "test_MD_case" / "md1ns.xtc"


@pytest.fixture(scope="module")
def system_data():
    reader = GmxTprReader()
    return reader.read(str(TPR_FILE))


@pytest.fixture(scope="module")
def groups(system_data):
    identifier = AmberFFGroupIdentifier()
    return identifier.identify(system_data)


@pytest.fixture(scope="module")
def relevant_groups(groups):
    return [g for g in groups
            if g.group_type in ("aromatic_ring", "charged_positive")]


@pytest.fixture(scope="module")
def pi_cations(relevant_groups):
    u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    detector = PiCationDetectorTwoPass()
    return detector.detect(relevant_groups, trajectory=u.trajectory)


@pytest.fixture(scope="module")
def pi_cations_per_frame(relevant_groups):
    u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    detector = PiCationDetectorPerFrame()
    return detector.detect(relevant_groups, trajectory=u.trajectory)


# ============================================================
# 基本功能测试
# ============================================================

class TestPiCationTwoPassBasic:

    def test_has_results(self, pi_cations):
        assert len(pi_cations) > 0

    def test_interaction_type(self, pi_cations):
        assert pi_cations[0].interaction_type == "pi_cation"

    def test_n_pairs(self, pi_cations):
        assert pi_cations[0].n_pairs > 0

    def test_metrics_keys(self, pi_cations):
        for interaction in pi_cations:
            assert "distance" in interaction.metrics
            assert "offset" in interaction.metrics


# ============================================================
# 输出格式测试
# ============================================================

class TestPiCationTwoPassFormat:

    def test_groups_are_ring_cation(self, pi_cations):
        it = pi_cations[0]
        for g1, g2 in it.groups:
            assert g1.group_type == "aromatic_ring"
            assert g2.group_type == "charged_positive"

    def test_existence_is_2d(self, pi_cations):
        for interaction in pi_cations:
            assert interaction.existence.ndim == 2

    def test_metrics_shape(self, pi_cations):
        for interaction in pi_cations:
            n_frames = interaction.n_frames
            for key in ["distance", "offset"]:
                assert interaction.metrics[key].shape == (
                    interaction.n_pairs, n_frames)


# ============================================================
# 阈值测试
# ============================================================

class TestPiCationTwoPassThreshold:

    def test_active_distance_range(self, pi_cations):
        it = pi_cations[0]
        for i in range(it.n_pairs):
            active = it.metrics["distance"][i][it.existence[i]]
            if len(active) > 0:
                assert np.all(active > 0.5)
                assert np.all(active < 6.0)

    def test_active_offset_range(self, pi_cations):
        it = pi_cations[0]
        for i in range(it.n_pairs):
            active = it.metrics["offset"][i][it.existence[i]]
            if len(active) > 0:
                assert np.all(active < 2.0)


# ============================================================
# 与 PerFrame 对比测试
# ============================================================

class TestPiCationTwoPassVsPerFrame:

    def test_two_pass_finds_at_least_as_many(self, pi_cations, pi_cations_per_frame):
        assert pi_cations[0].n_pairs >= pi_cations_per_frame[0].n_pairs

    def test_per_frame_pairs_are_subset(self, pi_cations, pi_cations_per_frame):
        two_pairs = {(g1.group_id, g2.group_id)
                     for g1, g2 in pi_cations[0].groups}
        per_pairs = {(g1.group_id, g2.group_id)
                     for g1, g2 in pi_cations_per_frame[0].groups}
        assert per_pairs.issubset(two_pairs)


# ============================================================
# 检测器元信息测试
# ============================================================

class TestPiCationTwoPassDetectorMeta:

    def test_name(self):
        assert PiCationDetectorTwoPass().name == "pi_cation"

    def test_required_group_types(self):
        det = PiCationDetectorTwoPass()
        assert "aromatic_ring" in det.required_group_types
        assert "charged_positive" in det.required_group_types

    def test_metric_names(self):
        det = PiCationDetectorTwoPass()
        assert "distance" in det.metric_names
        assert "offset" in det.metric_names
