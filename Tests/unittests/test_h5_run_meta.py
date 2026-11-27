# -*- coding: utf-8 -*-
"""h5 运行溯源组 run/ 的单元测试。

覆盖 RunMeta 写入/读取往返、旧文件无 run/ 组、退回 profile、
文件指纹、原子写、dii run 命令行规范形态与 format 1.0 向后兼容。
"""

import os
import re
import shutil
from datetime import datetime
from pathlib import Path

import h5py
import numpy as np
import pytest

import DuIvyInteractions
from DuIvyInteractions.core.datas import AtomData, Group, Interaction
from DuIvyInteractions.io.h5 import (
    FORMAT_VERSION, RUN_GROUP, RunMeta, load_interactions, load_run_meta,
    save_interactions, save_interactions_atomic, save_run_group,
)
from DuIvyInteractions.pipeline import build_run_meta

try:                          # Python < 3.11 无 tomllib
    import tomllib
except ModuleNotFoundError:   # pragma: no cover - 取决于解释器版本
    tomllib = None


REPO_ROOT = Path(__file__).resolve().parents[2]
TMP_DIR = REPO_ROOT / "test_temp" / "run_meta_test"


def _pyproject_version() -> str:
    """pyproject.toml 中 [project] 的 version。"""
    text = (REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    if tomllib is None:  # Python < 3.11：按行解析
        return re.search(r'(?m)^version\s*=\s*"([^"]+)"', text).group(1)
    return tomllib.loads(text)["project"]["version"]


@pytest.fixture(scope="module", autouse=True)
def setup_tmp_dir():
    """模块级临时目录（测试后清理）。"""
    if TMP_DIR.exists():
        shutil.rmtree(TMP_DIR, ignore_errors=True)
    TMP_DIR.mkdir(parents=True)
    yield
    shutil.rmtree(TMP_DIR, ignore_errors=True)


def _make_interaction() -> Interaction:
    """最小合法 Interaction（1 对 × 3 帧）。"""
    atom1 = AtomData(atom_global_idx=1, atom_idx_in_residue=1, atom_name="CG",
                     atom_type="ca", atom_element="C", atom_charge=-0.1,
                     atom_mass=12.011)
    atom2 = AtomData(atom_global_idx=2, atom_idx_in_residue=2, atom_name="CZ",
                     atom_type="ca", atom_element="C", atom_charge=-0.1,
                     atom_mass=12.011)
    ring = Group(group_id=1, group_type="aromatic_ring", molecule="LIG",
                 residue_name="LIG", residue_id=1, atoms=[atom1],
                 metadata={"ring_size": 6})
    hydro = Group(group_id=2, group_type="hydrophobic", molecule="LIG",
                  residue_name="LIG", residue_id=2, atoms=[atom2], metadata={})
    return Interaction(
        interaction_type="pi_stacking", groups=[(ring, hydro)],
        existence=np.array([[True, False, True]]),
        metrics={"distance": np.array([[3.5, 3.6, 3.7]])},
        times=np.array([0.0, 10.0, 20.0]))


def _make_run_meta(topology: str = "/tmp/含空格 目录/sys.tpr",
                   trajectory: str = "/tmp/traj.xtc",
                   output: str = "/tmp/out dir") -> RunMeta:
    """含中文、空格与 '=' 的运行参数。"""
    return RunMeta(
        topology=topology, trajectory=trajectory, output=output,
        ff="charmm", strategy="two_pass",
        interactions=("pi_stacking", "hydrogen_bond"),
        select=("segid A and resname 中文配体", "resname SOL"),
        pairing="between", mark_redundant=True,
        argv=("dii", "run", "-t", topology, "--ff", "charmm",
              "--select", "name CA=1"),
        cwd="/tmp/work dir", created_at="2026-11-26T18:00:00+08:00",
        dii_version="0.1.0")


def _write_run_file(name: str, run: RunMeta) -> str:
    """写入一个带 run/ 的 h5，返回路径。"""
    path = str(TMP_DIR / name)
    save_interactions([_make_interaction()], path, run=run)
    return path


class TestRunMetaRoundtrip:
    """RunMeta 写入后读出的往返一致性。"""

    def test_roundtrip_all_fields(self):
        run = _make_run_meta()
        record = load_run_meta(_write_run_file("roundtrip.h5", run))
        assert record is not None
        assert record.meta == run

    def test_roundtrip_empty_lists(self):
        run = RunMeta(topology="/tmp/a.tpr", trajectory="/tmp/a.xtc",
                      output="/tmp/out", ff="amber", strategy="two_pass")
        record = load_run_meta(_write_run_file("empty_lists.h5", run))
        assert record.meta.interactions == ()
        assert record.meta.select == ()
        assert record.meta.argv == ()
        assert isinstance(record.meta.interactions, tuple)

    def test_missing_run_group_returns_none(self):
        path = str(TMP_DIR / "old_format.h5")
        save_interactions([_make_interaction()], path)
        assert load_run_meta(path) is None

    def test_load_run_meta_missing_fields_defaults(self):
        path = str(TMP_DIR / "bare_run_group.h5")
        with h5py.File(path, "w") as f:
            f.create_group(RUN_GROUP)
        record = load_run_meta(path)
        assert record.meta.topology == ""
        assert record.meta.pairing == "any"
        assert record.meta.mark_redundant is False
        assert record.meta.interactions == ()
        assert (record.topology_size, record.topology_mtime) == (-1, -1.0)
        assert (record.trajectory_size, record.trajectory_mtime) == (-1, -1.0)

    def test_save_run_group_on_direct_handle(self):
        run = _make_run_meta()
        path = str(TMP_DIR / "direct_group.h5")
        with h5py.File(path, "w") as f:
            save_run_group(f, run)
        assert load_run_meta(path).meta == run


class TestBackwardCompatibility:
    """新增可选组不影响既有读取路径与 format 1.0。"""

    def test_load_interactions_unaffected(self):
        path = _write_run_file("compat_load.h5", _make_run_meta())
        loaded = load_interactions(path)
        assert len(loaded) == 1
        assert loaded[0].interaction_type == "pi_stacking"
        assert loaded[0].n_pairs == 1
        assert loaded[0].n_frames == 3
        assert np.array_equal(loaded[0].existence, [[True, False, True]])
        assert np.allclose(loaded[0].metrics["distance"], [[3.5, 3.6, 3.7]])

    def test_format_version_unchanged(self):
        path = _write_run_file("format_check.h5", _make_run_meta())
        with h5py.File(path, "r") as f:
            assert f.attrs["format_version"] == FORMAT_VERSION == "1.0"
            assert f.attrs["n_interactions"] == 1
            assert set(f.attrs.keys()) == {"format_version", "n_interactions"}
            assert RUN_GROUP in f

    def test_wrong_format_version_still_rejected(self):
        path = _write_run_file("wrong_version.h5", _make_run_meta())
        with h5py.File(path, "r+") as f:
            f.attrs["format_version"] = "9.9"
        with pytest.raises(ValueError, match="Unsupported format version"):
            load_interactions(path)


class TestFingerprints:
    """run/inputs 的文件指纹取自 os.stat。"""

    def test_fingerprints_match_stat(self):
        tpr = TMP_DIR / "fake.tpr"
        xtc = TMP_DIR / "fake.xtc"
        tpr.write_bytes(b"tpr" * 10)
        xtc.write_bytes(b"xtc" * 7)
        run = _make_run_meta(topology=str(tpr), trajectory=str(xtc))
        record = load_run_meta(_write_run_file("fingerprint.h5", run))
        assert record.topology_size == os.stat(tpr).st_size
        assert record.topology_mtime == pytest.approx(os.stat(tpr).st_mtime)
        assert record.trajectory_size == os.stat(xtc).st_size
        assert record.trajectory_mtime == pytest.approx(os.stat(xtc).st_mtime)

    def test_fingerprints_missing_file_default(self):
        run = _make_run_meta(topology=str(TMP_DIR / "nope.tpr"),
                             trajectory=str(TMP_DIR / "nope.xtc"))
        record = load_run_meta(_write_run_file("fingerprint_missing.h5", run))
        assert record.topology_size == -1
        assert record.topology_mtime == -1.0
        assert record.trajectory_size == -1
        assert record.trajectory_mtime == -1.0


class TestAtomicSave:
    """原子写与回写都保留 run/ 组。"""

    def test_atomic_save_with_run(self):
        run = _make_run_meta()
        path = str(TMP_DIR / "atomic.h5")
        save_interactions_atomic([_make_interaction()], path, run=run)
        assert not list(TMP_DIR.glob(f".{Path(path).name}.*.tmp"))
        assert load_run_meta(path).meta == run
        assert len(load_interactions(path)) == 1

    def test_rewrite_preserves_run_group(self):
        run = _make_run_meta()
        path = _write_run_file("rewrite.h5", run)
        loaded = load_interactions(path)
        save_interactions_atomic(loaded, path, run=run)
        assert load_run_meta(path).meta == run
        assert len(load_interactions(path)) == 1


class TestBuildRunMeta:
    """dii run 参数与规范命令行的构建。"""

    def test_build_argv_minimal(self):
        meta = build_run_meta("/tmp/md.tpr", "/tmp/md.xtc", "/tmp/out",
                              "amber", "two_pass", None, None, "any", False)
        assert meta.argv == (
            "dii", "run", "-t", "/tmp/md.tpr", "-f", "/tmp/md.xtc",
            "-o", "/tmp/out", "--ff", "amber", "--strategy", "two_pass",
            "--interactions", "all")

    def test_build_argv_with_selects(self):
        meta = build_run_meta(
            "/tmp/md.tpr", "/tmp/md.xtc", "/tmp/out", "amber", "per_frame",
            ["hydrogen_bond", "salt_bridge"], ["segid A", "resname 中文"],
            "any", False)
        assert meta.argv[-4:] == (
            "--select", "segid A", "--select", "resname 中文")
        assert meta.select == ("segid A", "resname 中文")

    def test_build_argv_between_and_mark(self):
        meta = build_run_meta(
            "/tmp/md.tpr", "/tmp/md.xtc", "/tmp/out", "charmm", "two_pass",
            ["hydrophobic"], ["name CA", "name CB"], "between", True)
        assert meta.argv[-3:] == ("--pairing", "between", "--mark-redundant")
        assert meta.mark_redundant is True

    def test_build_defaults_injected(self):
        meta = build_run_meta("/tmp/md.tpr", "/tmp/md.xtc", "/tmp/out",
                              "amber", "two_pass", [], [], "any", False)
        assert meta.cwd == os.getcwd()
        assert datetime.fromisoformat(meta.created_at).tzinfo is not None
        expected = getattr(DuIvyInteractions, "__version__", "unknown")
        assert meta.dii_version == str(expected)

    def test_build_abspath(self):
        meta = build_run_meta("rel/md.tpr", "rel/md.xtc", "rel/out",
                              "amber", "two_pass", None, None, "any", False)
        assert meta.topology == os.path.abspath("rel/md.tpr")
        assert meta.trajectory == os.path.abspath("rel/md.xtc")
        assert meta.output == os.path.abspath("rel/out")
        assert meta.argv[3] == os.path.abspath("rel/md.tpr")

    def test_build_dedup_preserves_order(self):
        meta = build_run_meta("/tmp/md.tpr", "/tmp/md.xtc", "/tmp/out",
                              "amber", "two_pass", ["b", "a", "b"],
                              ["x", "x"], "any", False)
        assert meta.interactions == ("b", "a")
        assert meta.select == ("x",)

    def test_build_explicit_argv_kept(self):
        meta = build_run_meta("/tmp/md.tpr", "/tmp/md.xtc", "/tmp/out",
                              "amber", "two_pass", None, None, "any", False,
                              argv=["dii", "run", "--ff=amber"])
        assert meta.argv == ("dii", "run", "--ff=amber")

    def test_build_created_at_explicit(self):
        meta = build_run_meta("/tmp/md.tpr", "/tmp/md.xtc", "/tmp/out",
                              "amber", "two_pass", None, None, "any", False,
                              cwd="/tmp/cwd", created_at="2026-01-02T03:04:05")
        assert meta.cwd == "/tmp/cwd"
        assert meta.created_at == "2026-01-02T03:04:05"


class TestVersion:
    """包版本号、pyproject.toml 与 run/dii_version 的一致性。"""

    def test_package_version_defined(self):
        assert DuIvyInteractions.__version__
        assert DuIvyInteractions.__version__ != "unknown"

    def test_package_version_matches_pyproject(self):
        assert DuIvyInteractions.__version__ == _pyproject_version()

    def test_build_run_meta_uses_package_version(self):
        meta = build_run_meta("/tmp/md.tpr", "/tmp/md.xtc", "/tmp/out",
                              "amber", "two_pass", None, None, "any", False)
        assert meta.dii_version == DuIvyInteractions.__version__


class TestLegacyH5Layout:
    """自造旧格式 h5（无 run/ 组）的读取兼容性。"""

    def test_legacy_h5_without_run_group(self, tmp_path):
        path = str(tmp_path / "legacy.h5")
        interaction = _make_interaction()
        save_interactions([interaction], path)
        assert load_run_meta(path) is None
        loaded = load_interactions(path)
        assert len(loaded) == 1
        assert loaded[0].n_pairs == interaction.n_pairs
        assert loaded[0].n_frames == interaction.n_frames
        with h5py.File(path, "r") as f:
            assert set(f.attrs.keys()) == {"format_version", "n_interactions"}
            assert RUN_GROUP not in f
