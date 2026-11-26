# -*- coding: utf-8 -*-
"""3D 视图自动载入 h5 的 run/ 记录。

后端：POST /api/structure/auto 的 200/409/400/404 契约、指纹宽容度与
_select_structure 复用同一 _install_structure。
模板/JS：structureConfig 注入 run_meta/auto_attempt，前端自动尝试与 auto=0 保护。
只用 tmp_path 构造最小 PDB + h5，不依赖真实轨迹。
"""

import inspect
import json
import re
from pathlib import Path

import numpy as np
import pytest

from DuIvyInteractions.core.datas import AtomData, Group, Interaction
from DuIvyInteractions.io.h5 import RunMeta, load_run_meta, save_interactions

WEB_DIR = Path(__file__).resolve().parents[2] / "DuIvyInteractions" / "web"
STRUCTURE_JS = WEB_DIR / "static" / "structure.js"
PARTIAL_TEMPLATE = WEB_DIR / "templates" / "_structure_view.html"

TYPE = "hydrogen_bond"
PAIR_URL = f"/pair/{TYPE}/0"
CONFIG_PATTERN = r'<script id="structureConfig"[^>]*>(.*?)</script>'
BOX = (10.0, 10.0, 10.0)
ATOMS = [("N1", (0.0, 0.0, 0.0), "N"),
         ("H1", (0.0, 0.0, 1.0), "H"),
         ("O1", (0.0, 0.0, 2.8), "O"),
         ("O2", (0.0, 0.0, 4.0), "O")]
RUN_META_KEYS = {"topology", "trajectory", "ff", "strategy", "created_at",
                 "dii_version", "argv", "interactions"}


# ============================================================
# 最小结构 / 数据源 / 客户端
# ============================================================


def _cryst1() -> str:
    return ("CRYST1%9.3f%9.3f%9.3f  90.00  90.00  90.00 P 1           1"
            % BOX)


def _pdb_text() -> str:
    lines = [_cryst1()]
    lines += [f"ATOM  {serial:5d} {name:>4s} ALA A   1    "
              f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}"
              for serial, (name, (x, y, z), element)
              in enumerate(ATOMS, start=1)]
    return "\n".join(lines + ["END"]) + "\n"


def _write_pdb(path: Path) -> Path:
    path.write_text(_pdb_text(), encoding="utf-8")
    return path


def _atom(index: int, name: str, atype: str = "") -> AtomData:
    return AtomData(atom_global_idx=index, atom_idx_in_residue=0,
                    atom_name=name, atom_type=atype, atom_element="C",
                    atom_charge=0.0, atom_mass=0.0)


def _group(group_id: int, group_type: str, atoms) -> Group:
    return Group(group_id=group_id, group_type=group_type, molecule="test",
                 residue_name="LYS", residue_id=1, atoms=list(atoms))


def _interaction() -> Interaction:
    existence = np.ones((1, 1), dtype=bool)
    groups = [(_group(0, "H_donor", [_atom(0, "N1", "N"),
                                     _atom(1, "H1", "H")]),
               _group(1, "H_acceptor", [_atom(2, "O1", "O"),
                                        _atom(3, "O2", "O")]))]
    return Interaction(interaction_type=TYPE, groups=groups,
                       existence=existence,
                       metrics={"distance": np.full((1, 1), 3.0)},
                       times=np.zeros(1))


def _run_meta(topology, trajectory) -> RunMeta:
    return RunMeta(topology=str(topology), trajectory=str(trajectory),
                   output="", ff="amber", strategy="two_pass",
                   interactions=(TYPE,),
                   argv=("dii", "run", "-t", str(topology)),
                   cwd="", created_at="2026-11-26T00:00:00",
                   dii_version="0.1.0")


def _save_h5(directory: Path, run=None) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    save_interactions([_interaction()], str(directory / f"{TYPE}.h5"),
                      run=run)
    return directory


@pytest.fixture(autouse=True)
def isolate_recent(monkeypatch, tmp_path):
    monkeypatch.setattr("DuIvyInteractions.visualizers.data.RECENT_PATH",
                        tmp_path / "recent.json")


@pytest.fixture
def pdb_pair(tmp_path):
    return (_write_pdb(tmp_path / "top.pdb"),
            _write_pdb(tmp_path / "traj.pdb"))


@pytest.fixture
def h5_dir(tmp_path, pdb_pair):
    return _save_h5(tmp_path / "out", _run_meta(*pdb_pair))


def _client(**kwargs):
    from fastapi.testclient import TestClient
    from DuIvyInteractions.web.app import create_app
    return TestClient(create_app(**kwargs))


def _select(client, h5_dir) -> None:
    response = client.post("/select", data={"path": str(h5_dir)})
    assert response.status_code in (200, 303)


def _auto(client, **data):
    payload = {"type": TYPE, "idx": 0}
    payload.update(data)
    return client.post("/api/structure/auto", data=payload)


@pytest.fixture
def client(h5_dir):
    client = _client()
    _select(client, h5_dir)
    return client


def _config(html: str) -> dict:
    match = re.search(CONFIG_PATTERN, html, re.S)
    assert match, "structureConfig script is missing"
    return json.loads(match.group(1))


def _page_config(client) -> dict:
    return _config(client.get(PAIR_URL, params={"view": "3d"}).text)


def _structure_js() -> str:
    return STRUCTURE_JS.read_text(encoding="utf-8")


def _function_body(text: str, name: str) -> str:
    """取 structure.js 里 2 空格缩进的函数体（静态切片，不做 JS 解析）。"""
    match = re.search(r"  function %s\(.*?\n  \}\n" % re.escape(name), text,
                      re.S)
    assert match, f"missing function {name}"
    return match.group(0)


# ============================================================
# 后端：200 自动载入
# ============================================================


class TestAutoSuccess:
    def test_installs_recorded_structure(self, client, pdb_pair):
        response = _auto(client)
        assert response.status_code == 200
        body = response.json()
        assert body["ok"] is True
        assert body["source"] == "run"
        assert body["n_atoms"] == len(ATOMS)
        assert body["n_frames"] == 1
        assert body["has_box"] is True
        assert Path(body["topology"]) == pdb_pair[0]
        assert Path(body["trajectory"]) == pdb_pair[1]

    def test_structure_state_is_set(self, client):
        assert _auto(client).status_code == 200
        response = client.get("/api/structure/interactions",
                              params={"type": TYPE, "idx": 0})
        assert response.status_code == 200
        assert response.json()["n_atoms"] == len(ATOMS)
        cfg = _page_config(client)
        assert cfg["has_structure"] is True
        assert cfg["structure_version"] == 1

    def test_repeated_auto_bumps_structure_version(self, client):
        assert _auto(client).status_code == 200
        first = _page_config(client)["structure_version"]
        assert _auto(client).status_code == 200
        second = _page_config(client)["structure_version"]
        assert second == first + 1

    def test_unrecorded_size_is_not_compared(self, tmp_path):
        """h5 保存时文件还不存在（指纹 -1）：之后出现即接受，不按大小判变更。"""
        topology = tmp_path / "late_top.pdb"
        trajectory = tmp_path / "late_traj.pdb"
        h5_dir = _save_h5(tmp_path / "out", _run_meta(topology, trajectory))
        record = load_run_meta(str(h5_dir / f"{TYPE}.h5"))
        assert record.topology_size == -1 and record.trajectory_size == -1
        _write_pdb(topology)
        _write_pdb(trajectory)
        client = _client()
        _select(client, h5_dir)
        assert _auto(client).status_code == 200


# ============================================================
# 后端：409 回退分支
# ============================================================


class TestAutoNoRunMeta:
    @pytest.fixture
    def bare_client(self, tmp_path):
        client = _client()
        _select(client, _save_h5(tmp_path / "bare"))
        return client

    def test_conflict_body(self, bare_client):
        response = _auto(bare_client)
        assert response.status_code == 409
        assert response.json() == {
            "detail": "no run metadata in this h5",
            "reason": "no-run-meta"}

    def test_no_structure_installed(self, bare_client):
        _auto(bare_client)
        cfg = _page_config(bare_client)
        assert cfg["has_structure"] is False
        assert cfg["run_meta"] is None


class TestAutoMissingFile:
    @pytest.fixture
    def missing_client(self, tmp_path):
        topology = tmp_path / "gone.tpr"
        trajectory = tmp_path / "gone.xtc"
        client = _client()
        _select(client, _save_h5(tmp_path / "out",
                                 _run_meta(topology, trajectory)))
        return client, topology, trajectory

    def test_conflict_reason_and_paths(self, missing_client):
        client, topology, trajectory = missing_client
        response = _auto(client)
        assert response.status_code == 409
        body = response.json()
        assert body["reason"] == "missing"
        assert body["detail"] == "recorded file is missing"
        assert Path(body["topology"]) == topology
        assert Path(body["trajectory"]) == trajectory

    def test_missing_trajectory_alone_is_enough(self, tmp_path, pdb_pair):
        topology, trajectory = pdb_pair[0], tmp_path / "gone.xtc"
        client = _client()
        _select(client, _save_h5(tmp_path / "out",
                                 _run_meta(topology, trajectory)))
        response = _auto(client)
        assert response.status_code == 409
        assert response.json()["reason"] == "missing"


class TestAutoChangedFile:
    @pytest.fixture
    def changed_client(self, tmp_path, pdb_pair):
        topology, trajectory = pdb_pair
        client = _client()
        _select(client, _save_h5(tmp_path / "out",
                                 _run_meta(topology, trajectory)))
        trajectory.write_text(trajectory.read_text(encoding="utf-8")
                              + "REMARK appended\n", encoding="utf-8")
        return client, pdb_pair

    def test_conflict_reason(self, changed_client):
        client, _ = changed_client
        response = _auto(client)
        assert response.status_code == 409
        body = response.json()
        assert body["reason"] == "changed"
        assert body["detail"] == "recorded file changed"

    def test_response_keeps_recorded_paths(self, changed_client):
        client, (topology, trajectory) = changed_client
        body = _auto(client).json()
        assert Path(body["topology"]) == topology
        assert Path(body["trajectory"]) == trajectory


# ============================================================
# 后端：其它失败路径（404 / 400 / 复用 select 逻辑）
# ============================================================


class TestAutoFailurePaths:
    def test_unknown_type_returns_404(self, client):
        response = _auto(client, type="no_such_type")
        assert response.status_code == 404
        assert "interaction type not loaded" in response.json()["detail"]

    def test_path_outside_root_returns_400(self, tmp_path):
        root = tmp_path / "root"
        root.mkdir()
        outside = _write_pdb(tmp_path / "outside.pdb")
        h5_dir = _save_h5(root / "out", _run_meta(outside, outside))
        client = _client(root=str(root))
        _select(client, h5_dir)
        response = _auto(client)
        assert response.status_code == 400
        assert response.json()["detail"] == "path outside root"

    def test_auto_400_when_universe_cannot_open(self, tmp_path):
        """文件在、指纹一致，但内容不是可打开的拓扑 → 400（与 select 同语义）。"""
        topology = tmp_path / "bad.tpr"
        trajectory = tmp_path / "bad.xtc"
        topology.write_text("not a structure\n", encoding="utf-8")
        trajectory.write_text("not a structure\n", encoding="utf-8")
        client = _client()
        _select(client, _save_h5(tmp_path / "out",
                                 _run_meta(topology, trajectory)))
        response = _auto(client)
        assert response.status_code == 400
        assert "cannot open structure" in response.json()["detail"]


class TestSelectUnchanged:
    def test_select_payload_keys(self, client, pdb_pair):
        response = client.post("/api/structure/select",
                               data={"topology": str(pdb_pair[0]),
                                     "trajectory": str(pdb_pair[1])})
        assert response.status_code == 200
        assert set(response.json()) == {"ok", "n_atoms", "n_frames",
                                        "has_box"}

    def test_select_empty_path_still_400(self, client):
        response = client.post("/api/structure/select",
                               data={"topology": "", "trajectory": ""})
        assert response.status_code == 400
        assert response.json()["detail"] == "structure path required"

    def test_select_and_auto_share_install_helper(self):
        from DuIvyInteractions.web import app

        assert "_install_structure(" in inspect.getsource(
            app._select_structure)
        assert "_install_structure(" in inspect.getsource(app._auto_structure)
        assert "load_universe" not in inspect.getsource(app._select_structure)


# ============================================================
# 模板：structureConfig 注入 run_meta / auto_attempt
# ============================================================


class TestStructurePageAutoConfig:
    def test_page_exposes_run_meta(self, client, pdb_pair):
        cfg = _page_config(client)
        assert cfg["auto_url"] == "/api/structure/auto"
        assert cfg["auto_attempt"] is True
        meta = cfg["run_meta"]
        assert set(meta) == RUN_META_KEYS
        assert Path(meta["topology"]) == pdb_pair[0]
        assert Path(meta["trajectory"]) == pdb_pair[1]
        assert meta["ff"] == "amber"
        assert meta["strategy"] == "two_pass"
        assert meta["interactions"] == [TYPE]
        assert meta["argv"][:2] == ["dii", "run"]

    def test_page_without_run_meta_disables_auto(self, tmp_path):
        client = _client()
        _select(client, _save_h5(tmp_path / "bare"))
        cfg = _page_config(client)
        assert cfg["run_meta"] is None
        assert cfg["auto_attempt"] is False

    def test_page_with_structure_disables_auto(self, client, pdb_pair):
        response = client.post("/api/structure/select",
                               data={"topology": str(pdb_pair[0]),
                                     "trajectory": str(pdb_pair[1])})
        assert response.status_code == 200
        cfg = _page_config(client)
        assert cfg["has_structure"] is True
        assert cfg["run_meta"] is not None
        assert cfg["auto_attempt"] is False

    def test_partial_template_declares_auto_config(self):
        text = PARTIAL_TEMPLATE.read_text(encoding="utf-8")
        assert '"/api/structure/auto"' in text
        for key in ('"auto_url"', '"run_meta"', '"auto_attempt"'):
            assert key in text, key

    def test_partial_template_defaults_when_context_missing(self):
        """旧调用方不传 run_meta/auto_attempt 时仍渲染为 null/false。"""
        from DuIvyInteractions.web.app import TEMPLATES

        html = TEMPLATES.get_template("_structure_view.html").render(
            title="structure", ds_name="", datasource="", types=[],
            active_type=None, type=TYPE, idx=0, pair_label="pair 1",
            n_atoms=4, n_frames=1, has_structure=False, existence=[1],
            topology="", trajectory="", browse_mode="picker")
        cfg = _config(html)
        assert cfg["run_meta"] is None
        assert cfg["auto_attempt"] is False
        assert cfg["auto_url"] == "/api/structure/auto"


# ============================================================
# 前端：structure.js 自动尝试流程
# ============================================================


class TestStructureScriptAuto:
    def test_auto_constants(self):
        text = _structure_js()
        assert 'var AUTO_URL_KEY = "auto_url"' in text
        assert 'var AUTO_TRIED_FLAG = "auto"' in text
        assert 'var AUTO_SKIP_VALUE = "0"' in text

    def test_posts_to_auto_url_then_reloads(self):
        body = _function_body(_structure_js(), "autoLoadStructure")
        assert "cfg[AUTO_URL_KEY]" in body
        assert 'method: "POST"' in body
        assert "autoRequestBody()" in body
        assert "window.location.reload()" in body

    def test_request_body_carries_type_and_idx(self):
        body = _function_body(_structure_js(), "autoRequestBody")
        assert "encodeURIComponent(cfg.type)" in body
        assert "cfg.idx" in body

    def test_conflict_prefills_paths_and_hints(self):
        text = _structure_js()
        fallback = _function_body(text, "autoFallback")
        assert "AUTO_CONFLICT_STATUS" in fallback
        assert "prefillAutoPaths(body)" in fallback
        assert "AUTO_UNAVAILABLE_TEXT" in fallback
        assert "AUTO_CONFIRM_HINT" in fallback
        prefill = _function_body(text, "prefillAutoPaths")
        assert "topologyInput.value" in prefill
        assert "trajectoryInput.value" in prefill
        assert "confirm the paths" in text

    def test_skip_flag_is_written_to_url(self):
        text = _structure_js()
        body = _function_body(text, "markAutoSkipped")
        assert "history.replaceState" in body
        assert "searchParams.set(AUTO_TRIED_FLAG, AUTO_SKIP_VALUE)" in body
        skip = _function_body(text, "autoSkipRequested")
        assert "searchParams.get(AUTO_TRIED_FLAG) === AUTO_SKIP_VALUE" in skip
        assert "autoSkipRequested()" in _function_body(text, "shouldAutoLoad")

    def test_auto_only_when_context_allows(self):
        text = _structure_js()
        body = _function_body(text, "shouldAutoLoad")
        assert "cfg.auto_attempt" in body
        assert "cfg[AUTO_URL_KEY]" in body
        assert "autoSkipRequested()" in body

    def test_auto_runs_before_viewer_init(self):
        text = _structure_js()
        assert "} else if (shouldAutoLoad()) {" in text
        assert "autoLoadStructure();" in text
        assert (text.index("function autoLoadStructure()")
                < text.index("function initViewer()"))

    def test_source_label_uses_run_meta_basenames(self):
        text = _structure_js()
        body = _function_body(text, "runSourceLabel")
        assert "cfg.has_structure" in body
        assert "cfg.run_meta" in body
        assert "basename(meta.topology)" in body
        assert "basename(meta.trajectory)" in body
        assert "runSourceLabel()" in _function_body(text, "statusText")

    def test_manual_select_flow_is_untouched(self):
        body = _function_body(_structure_js(), "selectStructure")
        assert "SELECT_URL" in body
        assert "topologyInput.value.trim()" in body
        assert "window.location.reload()" in body
        assert "openButton.disabled = true" in body
