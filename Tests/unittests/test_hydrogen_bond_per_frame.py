# -*- coding: utf-8 -*-
"""氢键 PerFrame 检测器集成测试。

验证 HydrogenBondDetectorPerFrame 的结果正确性。
"""

import pytest
from pathlib import Path

import numpy as np
import MDAnalysis as mda

from DuIvyInteractions.input_readers import GmxTprReader
from DuIvyInteractions.group_identifiers import AmberFFGroupIdentifier
from DuIvyInteractions.interaction_detectors import HydrogenBondDetectorPerFrame


TPR_FILE = Path(__file__).parent.parent / "test_MD_case" / "md.tpr"
XTC_FILE = Path(__file__).parent.parent / "test_MD_case" / "md1ns.xtc"

MOL_RBD = "seg_0_RBD_pro"
MOL_KRAS = "seg_2_KRAS_pro"
WATER_RESIDUES = {"SOL", "HOH", "WAT"}
PROTEIN_MOLECULES = {MOL_RBD, MOL_KRAS}


@pytest.fixture(scope="module")
def system_data():
    reader = GmxTprReader()
    return reader.read(str(TPR_FILE))


@pytest.fixture(scope="module")
def groups(system_data):
    return AmberFFGroupIdentifier().identify(system_data)


@pytest.fixture(scope="module")
def protein_groups(groups):
    return [g for g in groups
            if g.residue_name not in WATER_RESIDUES
            and g.molecule in PROTEIN_MOLECULES]


@pytest.fixture(scope="module")
def hydrogen_bonds(protein_groups):
    u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    detector = HydrogenBondDetectorPerFrame()
    return detector.detect(
        protein_groups, trajectory=u.trajectory,
        tuple_filter=lambda gt: gt[0].molecule != gt[1].molecule
    )


# ============================================================
# 基本功能
# ============================================================

class TestHydrogenBondPerFrameBasic:

    def test_has_results(self, hydrogen_bonds):
        assert len(hydrogen_bonds) > 0

    def test_interaction_type(self, hydrogen_bonds):
        assert hydrogen_bonds[0].interaction_type == "hydrogen_bond"

    def test_n_pairs(self, hydrogen_bonds):
        assert hydrogen_bonds[0].n_pairs == 119

    def test_all_inter_protein(self, hydrogen_bonds):
        it = hydrogen_bonds[0]
        for g1, g2 in it.groups:
            assert g1.molecule != g2.molecule

    def test_occupancy_range(self, hydrogen_bonds):
        occ = hydrogen_bonds[0].occupancy()
        assert np.all(occ >= 0) and np.all(occ <= 1)

    def test_top_pair(self, hydrogen_bonds):
        it = hydrogen_bonds[0]
        occ = it.occupancy()
        top = np.argmax(occ)
        g1, g2 = it.groups[top]
        assert g1.residue_name == "ARG" and g1.residue_id == 73
        assert g2.residue_name == "ASP" and g2.residue_id == 175
        assert occ[top] > 0.8

    def test_metrics_shape(self, hydrogen_bonds):
        it = hydrogen_bonds[0]
        assert it.metrics["distance"].shape == (it.n_pairs, it.n_frames)
        assert it.metrics["angle"].shape == (it.n_pairs, it.n_frames)

    def test_distance_range(self, hydrogen_bonds):
        it = hydrogen_bonds[0]
        for i in range(it.n_pairs):
            active = it.metrics["distance"][i][it.existence[i]]
            if len(active) > 0:
                assert np.all(active <= 4.1)

    def test_angle_range(self, hydrogen_bonds):
        it = hydrogen_bonds[0]
        for i in range(it.n_pairs):
            active = it.metrics["angle"][i][it.existence[i]]
            if len(active) > 0:
                assert np.all(active >= 100.0)
