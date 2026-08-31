# -*- coding: utf-8 -*-
"""π-π 堆积 PerFrame 检测器集成测试。

验证 PiStackingDetectorPerFrame 的结果正确性。
"""

import pytest
from pathlib import Path

import numpy as np
import MDAnalysis as mda

from DuIvyInteractions.input_readers import GmxTprReader
from DuIvyInteractions.group_identifiers import AmberFFGroupIdentifier
from DuIvyInteractions.interaction_detectors import PiStackingDetectorPerFrame


TPR_FILE = Path(__file__).parent.parent / "test_MD_case" / "md.tpr"
XTC_FILE = Path(__file__).parent.parent / "test_MD_case" / "md1ns.xtc"

MOL_RBD = "seg_0_RBD_pro"
MOL_D927 = "seg_1_D927"
MOL_KRAS = "seg_2_KRAS_pro"


@pytest.fixture(scope="module")
def system_data():
    reader = GmxTprReader()
    return reader.read(str(TPR_FILE))


@pytest.fixture(scope="module")
def groups(system_data):
    return AmberFFGroupIdentifier().identify(system_data)


@pytest.fixture(scope="module")
def aromatic_groups(groups):
    return [g for g in groups if g.group_type == "aromatic_ring"]


@pytest.fixture(scope="module")
def pi_stackings(aromatic_groups):
    u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    detector = PiStackingDetectorPerFrame(check_planarity=False)
    return detector.detect(aromatic_groups, trajectory=u.trajectory)


@pytest.fixture(scope="module")
def pi_stackings_with_planarity(aromatic_groups):
    u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    detector = PiStackingDetectorPerFrame(check_planarity=True)
    return detector.detect(aromatic_groups, trajectory=u.trajectory)


# ============================================================
# 基本功能
# ============================================================

class TestPiStackingPerFrameBasic:

    def test_returns_list(self, pi_stackings):
        assert isinstance(pi_stackings, list)

    def test_interaction_type(self, pi_stackings):
        for it in pi_stackings:
            assert it.interaction_type == "pi_stacking"

    def test_metrics_keys(self, pi_stackings):
        for it in pi_stackings:
            for key in ["distance", "angle", "offset", "pistacking_type"]:
                assert key in it.metrics

    def test_no_planarity_by_default(self, pi_stackings):
        for it in pi_stackings:
            assert "planarity_ring1" not in it.metrics
            assert "planarity_ring2" not in it.metrics


# ============================================================
# 平面性参数
# ============================================================

class TestPiStackingPerFramePlanarity:

    def test_planarity_metrics_present(self, pi_stackings_with_planarity):
        """check_planarity=True 时应包含平面性指标。"""
        if len(pi_stackings_with_planarity) > 0:
            for it in pi_stackings_with_planarity:
                assert "planarity_ring1" in it.metrics
                assert "planarity_ring2" in it.metrics

    def test_planarity_shape(self, pi_stackings_with_planarity):
        """平面性指标形状应为 (n_pairs, n_frames)。"""
        if len(pi_stackings_with_planarity) > 0:
            for it in pi_stackings_with_planarity:
                assert it.metrics["planarity_ring1"].shape == (it.n_pairs, it.n_frames)
                assert it.metrics["planarity_ring2"].shape == (it.n_pairs, it.n_frames)


# ============================================================
# 指标值范围
# ============================================================

class TestPiStackingPerFrameMetrics:

    def test_distance_positive(self, pi_stackings):
        for it in pi_stackings:
            assert np.all(it.metrics["distance"] >= 0)

    def test_angle_range(self, pi_stackings):
        for it in pi_stackings:
            assert np.all(it.metrics["angle"] >= 0)
            assert np.all(it.metrics["angle"] <= 90)

    def test_offset_positive(self, pi_stackings):
        for it in pi_stackings:
            assert np.all(it.metrics["offset"] >= 0)

    def test_pistacking_type_values(self, pi_stackings):
        for it in pi_stackings:
            types = set(it.metrics["pistacking_type"].flatten())
            assert types.issubset({'N', 'P', 'T'})

    def test_existence_matches_type(self, pi_stackings):
        for it in pi_stackings:
            for i in range(it.n_pairs):
                active = it.existence[i]
                active_types = it.metrics["pistacking_type"][i][active]
                assert np.all(active_types != 'N')


# ============================================================
# P/T 分类
# ============================================================

class TestPiStackingPerFrameClassification:

    def test_p_type_angle(self, pi_stackings):
        for it in pi_stackings:
            p_mask = it.metrics["pistacking_type"] == 'P'
            if np.any(p_mask):
                assert np.all(it.metrics["angle"][p_mask] <= 30.0)

    def test_t_type_angle(self, pi_stackings):
        for it in pi_stackings:
            t_mask = it.metrics["pistacking_type"] == 'T'
            if np.any(t_mask):
                assert np.all(it.metrics["angle"][t_mask] >= 60.0)

    def test_active_offset(self, pi_stackings):
        for it in pi_stackings:
            for i in range(it.n_pairs):
                active = it.existence[i]
                if np.any(active):
                    assert np.all(it.metrics["offset"][i][active] < 2.0)


# ============================================================
# 结果数据
# ============================================================

class TestPiStackingPerFrameResults:

    def test_has_results(self, pi_stackings):
        assert len(pi_stackings) > 0

    def test_n_pairs(self, pi_stackings):
        assert pi_stackings[0].n_pairs == 10

    def test_aromatic_ring_count(self, aromatic_groups):
        assert len(aromatic_groups) == 38

    def test_ring_distribution(self, aromatic_groups):
        from collections import Counter
        mol_counts = Counter(g.molecule for g in aromatic_groups)
        assert mol_counts[MOL_RBD] == 17
        assert mol_counts[MOL_D927] == 3
        assert mol_counts[MOL_KRAS] == 18

    def test_top_pair(self, pi_stackings):
        it = pi_stackings[0]
        occ = it.occupancy()
        top = np.argmax(occ)
        r1, r2 = it.groups[top]
        assert r1.residue_name == "HIS" and r1.residue_id == 237
        assert r2.residue_name == "TYR" and r2.residue_id == 238
        assert occ[top] > 0.5

    def test_top_pair_is_t_type(self, pi_stackings):
        it = pi_stackings[0]
        occ = it.occupancy()
        top = np.argmax(occ)
        pt = it.metrics["pistacking_type"][top]
        assert np.sum(pt == 'T') > np.sum(pt == 'P')

    def test_d927_participation(self, pi_stackings):
        it = pi_stackings[0]
        d927 = sum(1 for r1, r2 in it.groups
                   if r1.molecule == MOL_D927 or r2.molecule == MOL_D927)
        assert d927 >= 3

    def test_t_type_dominant(self, pi_stackings):
        it = pi_stackings[0]
        assert int(np.sum(it.metrics["pistacking_type"] == 'T')) > \
               int(np.sum(it.metrics["pistacking_type"] == 'P'))

    def test_active_distance_range(self, pi_stackings):
        it = pi_stackings[0]
        for i in range(it.n_pairs):
            active = it.metrics["distance"][i][it.existence[i]]
            if len(active) > 0:
                assert np.all(active <= 5.5)

    def test_active_angle_range(self, pi_stackings):
        it = pi_stackings[0]
        for i in range(it.n_pairs):
            active = it.metrics["angle"][i][it.existence[i]]
            if len(active) > 0:
                assert np.all((active <= 30.0) | (active >= 60.0))
