# -*- coding: utf-8 -*-
"""DII 命令行单元测试：parser 结构、main 分派、export 错误分支。

不跑真实检测（run 全链路见 test_pipeline.py），用 monkeypatch 与合成 h5。
"""

import argparse
import shutil
import sys
from pathlib import Path

import numpy as np
import pytest

from DuIvyInteractions.DII import build_parser, main, _run_export
from DuIvyInteractions.io.h5 import save_interactions
from DuIvyInteractions.core.datas import Interaction, Group, AtomData


TMP_DIR = Path(__file__).parent.parent.parent / "test_temp" / "dii_test"


@pytest.fixture(scope="module", autouse=True)
def _tmp_cleanup():
    if TMP_DIR.exists():
        shutil.rmtree(TMP_DIR, ignore_errors=True)
    TMP_DIR.mkdir(parents=True)
    yield
    shutil.rmtree(TMP_DIR, ignore_errors=True)


def _make_atom(gidx: int, name: str, element: str) -> AtomData:
    return AtomData(
        atom_global_idx=gidx, atom_idx_in_residue=0,
        atom_name=name, atom_type="", atom_element=element,
        atom_charge=0.0, atom_mass=0.0,
    )


def _make_group(gid: int, gtype: str, res_name: str) -> Group:
    return Group(group_id=gid, group_type=gtype, molecule="test",
                 residue_name=res_name, residue_id=1,
                 atoms=[_make_atom(10 + gid, "X", "X")])


def _make_saltbridge(n_pairs: int = 2, n_frames: int = 5) -> Interaction:
    """构造盐桥 Interaction（参考 test_exporters）。"""
    pos = _make_group(0, "charged_positive", "LYS")
    neg = _make_group(1, "charged_negative", "ASP")
    return Interaction(
        interaction_type="salt_bridge",
        groups=[(pos, neg)] * n_pairs,
        existence=np.full((n_pairs, n_frames), True, dtype=bool),
        metrics={"distance": np.full((n_pairs, n_frames), 3.5)},
        times=np.arange(n_frames, dtype=float) * 10.0,
    )


class TestBuildParser:
    """parser 结构与参数校验。"""

    def test_prog(self):
        assert build_parser().prog == "dii"

    def test_run_subcommand_parses(self):
        ns = build_parser().parse_args(
            ["run", "-t", "t.tpr", "-f", "t.xtc", "-o", "out", "--ff", "amber"])
        assert ns.command == "run"
        assert ns.tpr == "t.tpr"
        assert ns.xtc == "t.xtc"
        assert ns.output == "out"
        assert ns.ff == "amber"
        assert ns.interactions == "all"
        assert ns.strategy == "two_pass"

    def test_export_subcommand_parses(self):
        ns = build_parser().parse_args(["export", "-i", "in.h5", "-o", "out"])
        assert ns.command == "export"
        assert ns.input == "in.h5"
        assert ns.output == "out"

    def test_strategy_choice_accepted(self):
        ns = build_parser().parse_args(
            ["run", "-t", "t", "-f", "x", "-o", "o", "--ff", "amber",
             "--strategy", "per_frame", "--interactions", "salt_bridge"])
        assert ns.strategy == "per_frame"
        assert ns.interactions == "salt_bridge"

    def test_invalid_ff_raises(self):
        with pytest.raises(SystemExit):
            build_parser().parse_args(
                ["run", "-t", "t", "-f", "x", "-o", "o", "--ff", "badff"])

    def test_invalid_strategy_raises(self):
        with pytest.raises(SystemExit):
            build_parser().parse_args(
                ["run", "-t", "t", "-f", "x", "-o", "o", "--ff", "amber",
                 "--strategy", "bad"])

    def test_run_missing_required_raises(self):
        with pytest.raises(SystemExit):
            build_parser().parse_args(["run"])

    def test_no_command_raises(self):
        with pytest.raises(SystemExit):
            build_parser().parse_args([])


class TestMainRun:
    """main run 分支：--interactions 解析与分派（monkeypatch Pipeline）。"""

    def _run_main(self, monkeypatch, argv):
        calls = {}

        class FakePipeline:
            def __init__(self, ff, strategy):
                calls["ff"] = ff
                calls["strategy"] = strategy

            def run(self, tpr, xtc, output, interactions):
                calls["tpr"] = tpr
                calls["xtc"] = xtc
                calls["output"] = output
                calls["interactions"] = interactions

        import DuIvyInteractions.pipeline as pipeline_mod
        monkeypatch.setattr(pipeline_mod, "Pipeline", FakePipeline)
        monkeypatch.setattr(sys, "argv", argv)
        main()
        return calls

    def test_all_default_passes_none(self, monkeypatch):
        calls = self._run_main(
            monkeypatch,
            ["dii", "run", "-t", "t.tpr", "-f", "t.xtc", "-o", "out", "--ff", "amber"])
        assert calls["interactions"] is None
        assert calls["ff"] == "amber"
        assert calls["strategy"] == "two_pass"
        assert calls["tpr"] == "t.tpr"
        assert calls["xtc"] == "t.xtc"
        assert calls["output"] == "out"

    def test_explicit_types_stripped_lowered(self, monkeypatch):
        calls = self._run_main(
            monkeypatch,
            ["dii", "run", "-t", "t", "-f", "x", "-o", "o", "--ff", "amber",
             "--interactions", "hydrogen_bond, salt_bridge"])
        assert calls["interactions"] == ["hydrogen_bond", "salt_bridge"]

    def test_types_case_insensitive(self, monkeypatch):
        """类型名大小写不敏感（内部统一 lower）。"""
        calls = self._run_main(
            monkeypatch,
            ["dii", "run", "-t", "t", "-f", "x", "-o", "o", "--ff", "amber",
             "--interactions", "Hydrogen_Bond, Salt_Bridge"])
        assert calls["interactions"] == ["hydrogen_bond", "salt_bridge"]

    def test_all_mixed_with_types_raises(self, monkeypatch):
        with pytest.raises(SystemExit, match="cannot be mixed"):
            self._run_main(
                monkeypatch,
                ["dii", "run", "-t", "t", "-f", "x", "-o", "o", "--ff", "amber",
                 "--interactions", "all,hydrogen_bond"])

    def test_unknown_type_raises(self, monkeypatch):
        with pytest.raises(SystemExit, match="unknown interaction type"):
            self._run_main(
                monkeypatch,
                ["dii", "run", "-t", "t", "-f", "x", "-o", "o", "--ff", "amber",
                 "--interactions", "nope"])


class TestMainExport:
    """main export 分支：分派到 _run_export。"""

    def test_dispatch_passes_args(self, monkeypatch):
        captured = {}
        import DuIvyInteractions.DII as dii_mod

        def fake_run_export(args):
            captured["input"] = args.input
            captured["output"] = args.output

        monkeypatch.setattr(dii_mod, "_run_export", fake_run_export)
        monkeypatch.setattr(sys, "argv", ["dii", "export", "-i", "in.h5", "-o", "out"])
        main()
        assert captured == {"input": "in.h5", "output": "out"}


class TestRunExportErrors:
    """_run_export 的错误分支与成功导出。"""

    def test_corrupted_h5_raises(self):
        p = TMP_DIR / "bad.h5"
        p.write_text("garbage not h5", encoding="utf-8")
        args = argparse.Namespace(input=str(p), output=str(TMP_DIR / "bad_out"))
        with pytest.raises(SystemExit, match="cannot read h5 file"):
            _run_export(args)

    def test_missing_h5_raises(self):
        args = argparse.Namespace(input=str(TMP_DIR / "nofile.h5"),
                                  output=str(TMP_DIR / "no_out"))
        with pytest.raises(SystemExit, match="cannot read h5 file"):
            _run_export(args)

    def test_empty_h5_raises(self):
        p = TMP_DIR / "empty.h5"
        save_interactions([], str(p))
        args = argparse.Namespace(input=str(p), output=str(TMP_DIR / "empty_out"))
        with pytest.raises(SystemExit, match="no interaction results"):
            _run_export(args)

    def test_unknown_type_raises(self):
        it = _make_saltbridge(n_pairs=1)
        it.interaction_type = "foo"
        p = TMP_DIR / "unknown.h5"
        save_interactions([it], str(p))
        args = argparse.Namespace(input=str(p), output=str(TMP_DIR / "unknown_out"))
        with pytest.raises(SystemExit, match="unknown interaction type"):
            _run_export(args)

    def test_success_creates_files(self):
        p = TMP_DIR / "ok.h5"
        save_interactions([_make_saltbridge(n_pairs=1)], str(p))
        out = TMP_DIR / "ok_out"
        _run_export(argparse.Namespace(input=str(p), output=str(out)))
        names = {f.name for f in out.iterdir()}
        assert any(n.startswith("salt_bridge_") and n.endswith("_count.xvg")
                   for n in names)
        assert any(n.startswith("salt_bridge_") and n.endswith("_existence.xpm")
                   for n in names)
        assert any(n.startswith("salt_bridge_") and n.endswith("_summary.csv")
                   for n in names)
