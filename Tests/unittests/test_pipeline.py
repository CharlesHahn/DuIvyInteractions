# -*- coding: utf-8 -*-
"""Pipeline 单元测试：注册表完整性、构造器、基团过滤、run 全链路。

真实数据：Tests/test_MD_case/md.tpr + md1ns.xtc（KRAS-RBD D927 体系）。
"""

import shutil
from pathlib import Path

import pytest

from DuIvyInteractions.pipeline import (
    Pipeline, ALL_INTERACTIONS, DETECTOR_CLASSES, STRATEGY_INDEX,
)
from DuIvyInteractions.system_readers import GmxTprReader
from DuIvyInteractions.group_identifiers import AmberFFGroupIdentifier
from DuIvyInteractions.interaction_detectors import (
    HydrogenBondDetectorTwoPass, HydrogenBondDetectorPerFrame,
    HydrogenBondDetectorPerTuple,
)
from DuIvyInteractions.io.h5 import load_interactions
from DuIvyInteractions.core.datas import Group


TPR_FILE = Path(__file__).parent.parent / "test_MD_case" / "md.tpr"
XTC_FILE = Path(__file__).parent.parent / "test_MD_case" / "md1ns.xtc"

TMP_DIR = Path(__file__).parent.parent.parent / "test_temp" / "pipeline_test"


@pytest.fixture(scope="module", autouse=True)
def _tmp_cleanup():
    if TMP_DIR.exists():
        shutil.rmtree(TMP_DIR, ignore_errors=True)
    TMP_DIR.mkdir(parents=True)
    yield
    shutil.rmtree(TMP_DIR, ignore_errors=True)


@pytest.fixture(scope="module")
def real_groups():
    """从真实 tpr 识别基团（供过滤与全链路测试）。"""
    sd = GmxTprReader().read(str(TPR_FILE))
    return AmberFFGroupIdentifier().identify(sd)


class TestConstants:
    """注册表完整性。"""

    def test_all_interactions_count(self):
        assert len(ALL_INTERACTIONS) == 8

    def test_detector_classes_covers_all(self):
        assert set(DETECTOR_CLASSES) == set(ALL_INTERACTIONS)

    def test_each_type_has_three_strategies(self):
        for name in ALL_INTERACTIONS:
            assert len(DETECTOR_CLASSES[name]) == 3

    def test_strategy_index_mapping(self):
        assert STRATEGY_INDEX == {"two_pass": 0, "per_frame": 1, "per_tuple": 2}


class TestMakeIdentifier:
    """按力场构造识别器。"""

    def test_amber_returns_identifier(self):
        assert isinstance(Pipeline("amber")._make_identifier(),
                          AmberFFGroupIdentifier)

    def test_unknown_ff_raises(self):
        with pytest.raises(ValueError, match="unknown force field"):
            Pipeline("nope")._make_identifier()


class TestMakeDetector:
    """按类型和策略构造检测器。"""

    def test_hydrogen_bond_two_pass(self):
        d = Pipeline("amber", "two_pass")._make_detector("hydrogen_bond")
        assert isinstance(d, HydrogenBondDetectorTwoPass)

    def test_hydrogen_bond_per_frame(self):
        d = Pipeline("amber", "per_frame")._make_detector("hydrogen_bond")
        assert isinstance(d, HydrogenBondDetectorPerFrame)

    def test_hydrogen_bond_per_tuple(self):
        d = Pipeline("amber", "per_tuple")._make_detector("hydrogen_bond")
        assert isinstance(d, HydrogenBondDetectorPerTuple)

    def test_all_types_all_strategies_instantiate(self):
        for name in ALL_INTERACTIONS:
            for strategy in STRATEGY_INDEX:
                d = Pipeline("amber", strategy)._make_detector(name)
                assert d.name == name

    def test_unknown_type_raises(self):
        with pytest.raises(KeyError):
            Pipeline("amber")._make_detector("not_a_type")


class TestFilterGroups:
    """按 required_group_types 过滤；除水桥外排除水分子。"""

    def _make_group(self, gid, gtype, residue_name="LYS"):
        from DuIvyInteractions.core.datas import AtomData
        atom = AtomData(atom_global_idx=10 + gid, atom_idx_in_residue=0,
                        atom_name="X", atom_type="", atom_element="X",
                        atom_charge=0.0, atom_mass=0.0)
        return Group(group_id=gid, group_type=gtype, molecule="mol",
                     residue_name=residue_name, residue_id=1, atoms=[atom])

    def test_water_bridge_keeps_water(self):
        groups = [self._make_group(0, "H_donor"),
                  self._make_group(1, "water", "SOL")]
        det = Pipeline("amber", "two_pass")._make_detector("water_bridge")
        filtered = Pipeline._filter_groups(groups, det)
        assert {g.group_type for g in filtered} == {"H_donor", "water"}

    def test_saltbridge_excludes_water(self):
        groups = [self._make_group(0, "charged_positive"),
                  self._make_group(1, "water", "SOL"),
                  self._make_group(2, "charged_negative", "ASP")]
        det = Pipeline("amber", "two_pass")._make_detector("salt_bridge")
        filtered = Pipeline._filter_groups(groups, det)
        assert {g.group_type for g in filtered} == {
            "charged_positive", "charged_negative"}

    def test_real_groups_saltbridge_filter(self, real_groups):
        det = Pipeline("amber", "two_pass")._make_detector("salt_bridge")
        filtered = Pipeline._filter_groups(real_groups, det)
        types = {g.group_type for g in filtered}
        assert types == {"charged_positive", "charged_negative"}
        assert all(g.residue_name not in ("SOL", "HOH") for g in filtered)


class TestRun:
    """run 全链路（真实数据）：Reader → Identifier → Detector → h5。"""

    def test_run_pipeline_with_failure_tolerance(self, capsys):
        """一次真实 run 同时验证：成功类型生成 h5 且可回读、未知类型
        WARN 不中断、不生成文件。"""
        out = TMP_DIR / "run_mixed"
        Pipeline("amber", "two_pass").run(
            str(TPR_FILE), str(XTC_FILE), str(out),
            interactions=["salt_bridge", "not_a_type"])

        # 已知类型：h5 生成、可回读、有非空结果
        h5 = out / "salt_bridge.h5"
        assert h5.exists()
        its = load_interactions(str(h5))
        assert len(its) == 1
        assert its[0].interaction_type == "salt_bridge"
        assert its[0].n_pairs > 0

        # 未知类型：不生成 h5、不中断其余
        assert not (out / "not_a_type.h5").exists()
        captured = capsys.readouterr()
        assert "not_a_type detection failed" in captured.out
