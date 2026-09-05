# -*- coding: utf-8 -*-
"""氢键检测器集成测试：从真实 tpr + xtc 输入，验证完整的氢键检测结果。

使用 D927 体系的 RBD 和 KRAS 蛋白之间的氢键作为测试数据。
"""

import pytest
from collections import Counter
from pathlib import Path

import numpy as np

from DuIvyInteractions.system_readers import GmxTprReader
from DuIvyInteractions.group_identifiers import AmberFFGroupIdentifier
from DuIvyInteractions.interaction_detectors import HydrogenBondDetectorPerTuple
import MDAnalysis as mda


# 测试数据路径
TPR_FILE = Path(__file__).parent.parent / "test_MD_case" / "md.tpr"
XTC_FILE = Path(__file__).parent.parent / "test_MD_case" / "md1ns.xtc"

# GmxTprReader 的分子名前缀
MOL_RBD = "seg_0_RBD_pro"
MOL_KRAS = "seg_2_KRAS_pro"

# 水分子残基名
WATER_RESIDUES = {"SOL", "HOH", "WAT"}

# 蛋白分子
PROTEIN_MOLECULES = {MOL_RBD, MOL_KRAS}


@pytest.fixture(scope="module")
def system_data():
    """读取 tpr 文件，返回 SystemData。"""
    reader = GmxTprReader()
    return reader.read(str(TPR_FILE))


@pytest.fixture(scope="module")
def groups(system_data):
    """运行基团识别，返回 Group 列表。"""
    identifier = AmberFFGroupIdentifier()
    return identifier.identify(system_data)


@pytest.fixture(scope="module")
def protein_groups(groups):
    """排除水分子后的蛋白基团。"""
    return [g for g in groups
            if g.residue_name not in WATER_RESIDUES
            and g.molecule in PROTEIN_MOLECULES]


@pytest.fixture(scope="module")
def hydrogen_bonds(protein_groups):
    """运行氢键检测（仅蛋白间，排除水），返回 Interaction 列表。"""
    detector = HydrogenBondDetectorPerTuple()
    return detector.detect(
        protein_groups,
        n_workers=32,
        topology_path=str(TPR_FILE),
        trajectory_path=str(XTC_FILE),
        tuple_filter=lambda gt: gt[0].molecule != gt[1].molecule
    )


# ============================================================
# 氢键检测结果测试
# ============================================================

class TestHydrogenBondDetection:
    """验证 RBD ↔ KRAS 界面氢键检测。"""

    def test_has_results(self, hydrogen_bonds):
        """应检测到氢键。"""
        assert len(hydrogen_bonds) > 0

    def test_interaction_type(self, hydrogen_bonds):
        """相互作用类型应为 hydrogen_bond。"""
        assert hydrogen_bonds[0].interaction_type == "hydrogen_bond"

    def test_n_pairs(self, hydrogen_bonds):
        """应检测到 119 对氢键。"""
        assert hydrogen_bonds[0].n_pairs == 119

    def test_all_inter_protein(self, hydrogen_bonds):
        """所有氢键应来自不同蛋白。"""
        it = hydrogen_bonds[0]
        for g1, g2 in it.groups:
            assert g1.molecule != g2.molecule

    def test_occupancy_range(self, hydrogen_bonds):
        """占位率应在 [0, 1] 范围内。"""
        occupancy = hydrogen_bonds[0].occupancy()
        assert np.all(occupancy >= 0)
        assert np.all(occupancy <= 1)

    def test_top_pair_occupancy(self, hydrogen_bonds):
        """最高占位率的氢键应为 ARG73(RBD) → ASP175(KRAS)，占位率 > 0.8。"""
        it = hydrogen_bonds[0]
        occupancy = it.occupancy()
        top_idx = np.argmax(occupancy)
        g1, g2 = it.groups[top_idx]
        assert g1.residue_name == "ARG"
        assert g1.residue_id == 73
        assert g2.residue_name == "ASP"
        assert g2.residue_id == 175
        assert occupancy[top_idx] > 0.8

    def test_metrics_shape(self, hydrogen_bonds):
        """metrics 数组形状应为 (n_pairs, n_frames)。"""
        it = hydrogen_bonds[0]
        n_frames = it.n_frames
        assert it.metrics["distance"].shape == (it.n_pairs, n_frames)
        assert it.metrics["angle"].shape == (it.n_pairs, n_frames)

    def test_distance_range(self, hydrogen_bonds):
        """活跃帧的 D-A 距离应 ≤ 4.1 Å。"""
        it = hydrogen_bonds[0]
        for i in range(it.n_pairs):
            active_dist = it.metrics["distance"][i][it.existence[i]]
            if len(active_dist) > 0:
                assert np.all(active_dist <= 4.1)

    def test_angle_range(self, hydrogen_bonds):
        """活跃帧的 D-H···A 角度应 ≥ 100°。"""
        it = hydrogen_bonds[0]
        for i in range(it.n_pairs):
            active_angle = it.metrics["angle"][i][it.existence[i]]
            if len(active_angle) > 0:
                assert np.all(active_angle >= 100.0)
