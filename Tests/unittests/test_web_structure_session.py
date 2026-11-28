# -*- coding: utf-8 -*-
"""3D 结构会话缓存与 CSP 回归。

会话相关响应（拓扑 mmCIF / 原始轨迹）由当前结构动态生成，必须禁止缓存：复用上一次
会话的轨迹会让 Mol* 抛 element count mismatch。前端另加会话版本号 ?v=… 作双保险
（结构每次重选自增）。Mol* 内置 WASM 经 fetch(data:…base64) 载入，3D 页 CSP 的
connect-src 需放行 data:。

只用最小 PDB + h5，不依赖真实轨迹。
"""

import json
import re
from pathlib import Path

import numpy as np
import pytest

from DuIvyInteractions.core.datas import AtomData, Group, Interaction
from DuIvyInteractions.io.h5 import save_interactions

WEB_DIR = Path(__file__).resolve().parents[2] / "DuIvyInteractions" / "web"
STRUCTURE_JS = WEB_DIR / "static" / "structure.js"
PAIR_URL = "/pair/hydrogen_bond/0"
CONFIG_PATTERN = r'<script id="structureConfig"[^>]*>(.*?)</script>'
BOX = (10.0, 10.0, 10.0)
ATOMS = [("N1", (0.0, 0.0, 0.0), "N"),
         ("H1", (0.0, 0.0, 1.0), "H"),
         ("O1", (0.0, 0.0, 2.8), "O"),
         ("O2", (0.0, 0.0, 4.0), "O")]


# ============================================================
# 最小结构 / 数据源
# ============================================================


def _cryst1() -> str:
    return ("CRYST1%9.3f%9.3f%9.3f  90.00  90.00  90.00 P 1           1"
            % BOX)


def _atom_lines() -> list:
    return [f"ATOM  {serial:5d} {name:>4s} ALA A   1    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}"
            for serial, (name, (x, y, z), element)
            in enumerate(ATOMS, start=1)]


def _pdb_text() -> str:
    return "\n".join([_cryst1()] + _atom_lines() + ["END"]) + "\n"


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
    return Interaction(interaction_type="hydrogen_bond", groups=groups,
                       existence=existence,
                       metrics={"distance": np.full((1, 1), 3.0)},
                       times=np.zeros(1))


@pytest.fixture(autouse=True)
def isolate_recent(monkeypatch, tmp_path):
    monkeypatch.setattr("DuIvyInteractions.visualizers.data.RECENT_PATH",
                        tmp_path / "recent.json")


@pytest.fixture
def pdb_file(tmp_path):
    path = tmp_path / "mini.pdb"
    path.write_text(_pdb_text(), encoding="utf-8")
    return path


@pytest.fixture
def h5_dir(tmp_path):
    directory = tmp_path / "out"
    directory.mkdir()
    save_interactions([_interaction()], str(directory / "hydrogen_bond.h5"))
    return directory


@pytest.fixture
def client(h5_dir):
    from fastapi.testclient import TestClient

    from DuIvyInteractions.web.app import create_app
    client = TestClient(create_app())
    client.post("/select", data={"path": str(h5_dir)})
    return client


def _open_structure(client, pdb_file):
    return client.post("/api/structure/select",
                       data={"topology": str(pdb_file),
                             "trajectory": str(pdb_file)})


def _config(html: str) -> dict:
    match = re.search(CONFIG_PATTERN, html, re.S)
    assert match, "structureConfig script is missing"
    return json.loads(match.group(1))


def _config_of(client) -> dict:
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
# A. 会话相关响应禁止缓存
# ============================================================


class TestSessionResponseCache:

    def test_trajectory_response_forbids_storage(self, client, pdb_file):
        _open_structure(client, pdb_file)
        response = client.get("/api/structure/trajectory")
        assert response.status_code == 200
        assert "no-store" in response.headers["cache-control"]

    def test_topology_response_forbids_storage(self, client, pdb_file):
        _open_structure(client, pdb_file)
        response = client.get("/api/structure/topology",
                              params={"type": "hydrogen_bond", "idx": 0})
        assert response.status_code == 200
        assert "no-store" in response.headers["cache-control"]

    def test_versioned_urls_are_accepted(self, client, pdb_file):
        """前端 ?v=… 只是缓存键，服务端不得因此 400/404。"""
        _open_structure(client, pdb_file)
        version = _config_of(client)["structure_version"]
        topology = client.get("/api/structure/topology",
                              params={"type": "hydrogen_bond", "idx": 0,
                                      "v": version})
        trajectory = client.get("/api/structure/trajectory",
                                params={"v": version})
        assert topology.status_code == 200
        assert trajectory.status_code == 200
        assert "no-store" in trajectory.headers["cache-control"]

    @pytest.mark.parametrize("asset", ("structure.js", "app.css"))
    def test_static_assets_still_revalidate(self, client, asset):
        response = client.get("/static/" + asset)
        assert response.headers["cache-control"] == "no-cache"

    def test_cache_control_constant_is_no_store(self):
        from DuIvyInteractions.web.app import SESSION_CACHE_CONTROL
        assert SESSION_CACHE_CONTROL == "no-store"


# ============================================================
# B. 会话版本号：后端自增 + 页面注入
# ============================================================


class TestStructureVersion:

    def test_version_increases_on_reselect(self, client, pdb_file):
        assert _open_structure(client, pdb_file).status_code == 200
        first = _config_of(client)["structure_version"]
        assert _open_structure(client, pdb_file).status_code == 200
        second = _config_of(client)["structure_version"]
        assert isinstance(first, int) and isinstance(second, int)
        assert second > first

    def test_page_config_carries_integer_version(self, client, pdb_file):
        _open_structure(client, pdb_file)
        config = _config_of(client)
        assert isinstance(config["structure_version"], int)
        assert config["structure_version"] >= 1

    def test_version_defaults_to_zero_without_structure(self, client):
        config = _config_of(client)
        assert config["has_structure"] is False
        assert config["structure_version"] == 0


# ============================================================
# B/D. 前端：版本号拼进拓扑/轨迹 URL、原子数不一致提示
# ============================================================


class TestStructureScript:

    def test_version_parameter_constant(self):
        assert 'VERSION_PARAM = "v"' in _structure_js()

    def test_urls_carry_version_parameter(self):
        text = _structure_js()
        encoded = "encodeURIComponent(cfg.structure_version || 0)"
        for name in ("topologyUrl", "trajectoryUrl"):
            body = _function_body(text, name)
            assert "VERSION_PARAM" in body, name
            assert 'VERSION_PARAM + "="' in body, name
            assert encoded in body, name
        topology = _function_body(text, "topologyUrl")
        assert "type=" in topology and "idx=" in topology

    def test_trajectory_loader_uses_versioned_url(self):
        body = _function_body(_structure_js(), "loadTrajectoryOnce")
        assert "url: trajectoryUrl()" in body
        assert "coordinates-url" in body

    def test_version_falls_back_to_zero(self):
        """旧页面/退化时缺 structure_version，URL 仍可取 v=0。"""
        assert _structure_js().count("cfg.structure_version || 0") >= 2

    def test_mismatch_hint_branch(self):
        text = _structure_js()
        assert '"element count mismatch"' in text
        assert "atom counts differ" in text
        assert "re-select the trajectory" in text
        body = _function_body(text, "loadFailureText")
        assert "MISMATCH_MARKER" in body
        assert "toLowerCase().indexOf(MISMATCH_MARKER)" in body
        assert "MISMATCH_HINT" in body


# ============================================================
# C. CSP：3D 页放行 WASM 的 data: URL，其余指令不放宽
# ============================================================


class TestStructureCsp:

    def test_3d_view_allows_data_wasm(self, client):
        csp = client.get(PAIR_URL, params={"view": "3d"}).headers[
            "content-security-policy"]
        assert "connect-src 'self' data:" in csp
        assert "connect-src *" not in csp

    def test_3d_view_keeps_other_directives(self, client):
        csp = client.get(PAIR_URL, params={"view": "3d"}).headers[
            "content-security-policy"]
        for directive in ("default-src 'self'",
                          "script-src 'self' 'unsafe-inline' 'unsafe-eval'",
                          "style-src 'self' 'unsafe-inline'",
                          "img-src 'self' data:"):
            assert directive in csp, directive

    def test_merged_pair_page_allows_data_wasm(self, client):
        """合并后的单对页恒含 Mol*：connect-src 必须放行 data:（WASM data: URL）"""
        csp = client.get(PAIR_URL).headers["content-security-policy"]
        assert "connect-src 'self' data:" in csp
        assert "connect-src *" not in csp
