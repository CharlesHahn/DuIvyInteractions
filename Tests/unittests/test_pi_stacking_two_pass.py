# -*- coding: utf-8 -*-
"""π-π 堆积 TwoPass 检测器集成测试。

从真实 tpr + xtc 输入，验证检测结果。
与 PerFrame 版本对比，验证结果一致性。
"""

import pytest
from pathlib import Path

import numpy as np

from DuIvyInteractions.system_readers import GmxTprReader
from DuIvyInteractions.group_identifiers import AmberFFGroupIdentifier
from DuIvyInteractions.interaction_detectors import (
    PiStackingDetectorTwoPass, PiStackingDetectorPerFrame)
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
def aromatic_groups(groups):
    return [g for g in groups if g.group_type == "aromatic_ring"]


@pytest.fixture(scope="module")
def pi_stackings(aromatic_groups):
    u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    detector = PiStackingDetectorTwoPass(check_planarity=False)
    return detector.detect(aromatic_groups, trajectory=u.trajectory)


@pytest.fixture(scope="module")
def pi_stackings_per_frame(aromatic_groups):
    u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    detector = PiStackingDetectorPerFrame(check_planarity=False)
    return detector.detect(aromatic_groups, trajectory=u.trajectory)


@pytest.fixture(scope="module")
def pi_stackings_planar(aromatic_groups):
    u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    detector = PiStackingDetectorTwoPass(check_planarity=True)
    return detector.detect(aromatic_groups, trajectory=u.trajectory)


# ============================================================
# 基本功能测试
# ============================================================

class TestPiStackingTwoPassBasic:

    def test_has_results(self, pi_stackings):
        assert len(pi_stackings) > 0

    def test_interaction_type(self, pi_stackings):
        assert pi_stackings[0].interaction_type == "pi_stacking"

    def test_n_pairs(self, pi_stackings):
        assert pi_stackings[0].n_pairs > 0

    def test_metrics_keys(self, pi_stackings):
        for interaction in pi_stackings:
            assert "distance" in interaction.metrics
            assert "angle" in interaction.metrics
            assert "offset" in interaction.metrics
            assert "pistacking_type" in interaction.metrics


# ============================================================
# 输出格式测试
# ============================================================

class TestPiStackingTwoPassFormat:

    def test_groups_are_rings(self, pi_stackings):
        it = pi_stackings[0]
        for g1, g2 in it.groups:
            assert g1.group_type == "aromatic_ring"
            assert g2.group_type == "aromatic_ring"

    def test_existence_is_2d(self, pi_stackings):
        for interaction in pi_stackings:
            assert interaction.existence.ndim == 2

    def test_metrics_shape(self, pi_stackings):
        for interaction in pi_stackings:
            n_frames = interaction.n_frames
            for key in ["distance", "angle", "offset"]:
                assert interaction.metrics[key].shape == (
                    interaction.n_pairs, n_frames)

    def test_pistacking_type_values(self, pi_stackings):
        """pistacking_type 应只含 P, T, N。"""
        it = pi_stackings[0]
        for i in range(it.n_pairs):
            types = set(it.metrics["pistacking_type"][i])
            assert types.issubset({'P', 'T', 'N'})


# ============================================================
# 阈值测试
# ============================================================

class TestPiStackingTwoPassThreshold:

    def test_active_distance_range(self, pi_stackings):
        it = pi_stackings[0]
        for i in range(it.n_pairs):
            active = it.metrics["distance"][i][it.existence[i]]
            if len(active) > 0:
                assert np.all(active > 0.5)
                assert np.all(active <= 5.5)

    def test_active_type_not_N(self, pi_stackings):
        """活跃帧的 pistacking_type 不应为 N。"""
        it = pi_stackings[0]
        for i in range(it.n_pairs):
            active_types = it.metrics["pistacking_type"][i][it.existence[i]]
            for t in active_types:
                assert t in ('P', 'T')


# ============================================================
# 平面性测试
# ============================================================

class TestPiStackingTwoPassPlanarity:

    def test_planarity_metrics_exist(self, pi_stackings_planar):
        it = pi_stackings_planar[0]
        assert "planarity_ring1" in it.metrics
        assert "planarity_ring2" in it.metrics

    def test_planarity_pairs_le_no_planarity(self, pi_stackings, pi_stackings_planar):
        """有平面性检查时，pair 数应 ≤ 无平面性检查。"""
        n_no = pi_stackings[0].n_pairs if pi_stackings else 0
        n_yes = pi_stackings_planar[0].n_pairs if pi_stackings_planar else 0
        assert n_yes <= n_no


# ============================================================
# 与 PerFrame 对比测试
# ============================================================

class TestPiStackingTwoPassVsPerFrame:

    def test_two_pass_finds_at_least_as_many(self, pi_stackings, pi_stackings_per_frame):
        n_two = pi_stackings[0].n_pairs if pi_stackings else 0
        n_per = pi_stackings_per_frame[0].n_pairs if pi_stackings_per_frame else 0
        assert n_two >= n_per

    def test_per_frame_pairs_are_subset(self, pi_stackings, pi_stackings_per_frame):
        if not pi_stackings or not pi_stackings_per_frame:
            pytest.skip("No results")
        two_pairs = {(g1.group_id, g2.group_id)
                     for g1, g2 in pi_stackings[0].groups}
        per_pairs = {(g1.group_id, g2.group_id)
                     for g1, g2 in pi_stackings_per_frame[0].groups}
        assert per_pairs.issubset(two_pairs)


# ============================================================
# 检测器元信息测试
# ============================================================

class TestPiStackingTwoPassDetectorMeta:

    def test_name(self):
        assert PiStackingDetectorTwoPass().name == "pi_stacking"

    def test_required_group_types(self):
        assert "aromatic_ring" in PiStackingDetectorTwoPass().required_group_types

    def test_metric_names_without_planarity(self):
        det = PiStackingDetectorTwoPass(check_planarity=False)
        assert "planarity_ring1" not in det.metric_names
        assert "planarity_ring2" not in det.metric_names

    def test_metric_names_with_planarity(self):
        det = PiStackingDetectorTwoPass(check_planarity=True)
        assert "planarity_ring1" in det.metric_names
        assert "planarity_ring2" in det.metric_names
