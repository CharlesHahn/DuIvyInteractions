# -*- coding: utf-8 -*-
"""卤键 PerFrame 检测器集成测试。

验证 HalogenBondDetectorPerFrame 的结果正确性。
"""

import pytest
from pathlib import Path

import numpy as np
import MDAnalysis as mda

from DuIvyInteractions.input_readers import GmxTprReader
from DuIvyInteractions.group_identifiers import AmberFFGroupIdentifier
from DuIvyInteractions.interaction_detectors import HalogenBondDetectorPerFrame


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
def halogen_bonds(groups):
    filtered = [g for g in groups
                if g.group_type in ("halogen_donor", "halogen_acceptor")]
    u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    detector = HalogenBondDetectorPerFrame()
    return detector.detect(filtered, trajectory=u.trajectory)


# ============================================================
# 基本功能
# ============================================================

class TestHalogenBondPerFrameBasic:

    def test_has_results(self, halogen_bonds):
        assert len(halogen_bonds) > 0

    def test_interaction_type(self, halogen_bonds):
        assert halogen_bonds[0].interaction_type == "halogen_bond"

    def test_metrics_keys(self, halogen_bonds):
        assert set(halogen_bonds[0].metrics.keys()) == {"distance", "don_angle", "acc_angle"}


# ============================================================
# 结果数据
# ============================================================

class TestHalogenBondPerFrameResults:

    def test_n_pairs(self, halogen_bonds):
        assert halogen_bonds[0].n_pairs == 1

    def test_active_pair_is_d927(self, halogen_bonds):
        it = halogen_bonds[0]
        occ = it.occupancy()
        active = np.where(occ > 0)[0]
        assert len(active) == 1
        g1, g2 = it.groups[active[0]]
        assert g1.molecule == MOL_D927

    def test_occupancy_range(self, halogen_bonds):
        occ = halogen_bonds[0].occupancy()
        assert np.all(occ >= 0) and np.all(occ <= 1)

    def test_metrics_shape(self, halogen_bonds):
        it = halogen_bonds[0]
        for key in ["distance", "don_angle", "acc_angle"]:
            assert it.metrics[key].shape == (it.n_pairs, it.n_frames)

    def test_active_distance_range(self, halogen_bonds):
        it = halogen_bonds[0]
        for i in range(it.n_pairs):
            active = it.metrics["distance"][i][it.existence[i]]
            if len(active) > 0:
                assert np.all(active <= 4.0)

    def test_active_don_angle_range(self, halogen_bonds):
        it = halogen_bonds[0]
        for i in range(it.n_pairs):
            active = it.metrics["don_angle"][i][it.existence[i]]
            if len(active) > 0:
                assert np.all(active >= 135.0) and np.all(active <= 195.0)

    def test_active_acc_angle_range(self, halogen_bonds):
        it = halogen_bonds[0]
        for i in range(it.n_pairs):
            active = it.metrics["acc_angle"][i][it.existence[i]]
            if len(active) > 0:
                assert np.all(active >= 90.0) and np.all(active <= 150.0)
