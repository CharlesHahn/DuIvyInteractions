# -*- coding: utf-8 -*-
"""CHARMM36 真实体系水桥相互作用测试（含水分子的 fullatom 轨迹）。

SMO-BST 全原子体系（`SMO-BST_md_fullatom.tpr` + `..._100ps.xtc`，100ps 片段/11 帧），
含 63055 个 SOL 水——水桥检测需要水分子，故使用全原子轨迹而非删水复合物。

⚠️ 断言基线（pair 数：PerFrame 1329、TwoPass 1678）来自本策略实现输出，
且轨迹仅 11 帧（100ps 片段）——统计意义有限，【尚未经人工核验】，请谨慎引用。
"""

import pytest
from pathlib import Path

import MDAnalysis as mda

from DuIvyInteractions.system_readers import GmxTprReader
from DuIvyInteractions.group_identifiers.charmm_ff_identifier import CharmmFFGroupIdentifier
from DuIvyInteractions.interaction_detectors import (
    WaterBridgeDetectorPerFrame,
    WaterBridgeDetectorTwoPass,
)


TPR_FILE = Path(__file__).parent.parent / "test_MD_case_charmm36" / "SMO-BST_md_fullatom.tpr"
XTC_FILE = Path(__file__).parent.parent / "test_MD_case_charmm36" / "SMO-BST_md_fullatom_100ps.xtc"


@pytest.fixture(scope="module")
def groups():
    sd = GmxTprReader().read(str(TPR_FILE))
    return CharmmFFGroupIdentifier().identify(sd)


@pytest.fixture(scope="module")
def trajectory():
    u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    return u.trajectory


@pytest.fixture(scope="module")
def water_bridges_pf(groups, trajectory):
    return WaterBridgeDetectorPerFrame().detect(groups, trajectory=trajectory)


@pytest.fixture(scope="module")
def water_bridges_tp(groups, trajectory):
    return WaterBridgeDetectorTwoPass().detect(groups, trajectory=trajectory)


class TestWaterBridgePerFrame:

    def test_has_results(self, water_bridges_pf):
        assert len(water_bridges_pf) > 0

    def test_n_pairs(self, water_bridges_pf):
        assert water_bridges_pf[0].n_pairs == 1329


class TestWaterBridgeTwoPass:

    def test_has_results(self, water_bridges_tp):
        assert len(water_bridges_tp) > 0

    def test_n_pairs(self, water_bridges_tp):
        assert water_bridges_tp[0].n_pairs == 1678