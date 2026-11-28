# -*- coding: utf-8 -*-
"""单对页合并布局测试：3D 在首、summary/groups 上移、frames 分页在标题行、无视图切换。

不依赖真实轨迹：用 tmp_path 构造最小 PDB + h5（含可选 run/ 记录）。
"""

import json
import re
from pathlib import Path

import numpy as np
import pytest

from DuIvyInteractions.core.datas import AtomData, Group, Interaction
from DuIvyInteractions.io.h5 import RunMeta, save_interactions

WEB_DIR = Path(__file__).resolve().parents[2] / "DuIvyInteractions" / "web"
APP_CSS = WEB_DIR / "static" / "app.css"
TYPE = "hydrogen_bond"
PAIRS = 2
FRAMES = 3
PAIR_URL = f"/pair/{TYPE}/0"
CONFIG_PATTERN = r'<script id="structureConfig"[^>]*>(.*?)</script>'
BOX = (10.0, 10.0, 10.0)
ATOMS = [("N1", (0.0, 0.0, 0.0), "N"),
         ("H1", (0.0, 0.0, 1.0), "H"),
         ("O1", (0.0, 0.0, 2.8), "O"),
         ("O2", (0.0, 0.0, 4.0), "O")]
# 页面区块（自上而下）；顺序断言用它们在 HTML 里的起始下标
BLOCKS = ('<div class="breadcrumb">', "structure-source", "<h2>summary",
          "<h2>groups", 'id="pairChart"', 'id="framesPanel"')


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
    groups = [(_group(0, "H_donor", [_atom(0, "N1", "N"),
                                     _atom(1, "H1", "H")]),
               _group(1, "H_acceptor", [_atom(2, "O1", "O"),
                                        _atom(3, "O2", "O")])),
              (_group(2, "H_donor", [_atom(4, "N2", "N"),
                                     _atom(5, "H2", "H")]),
               _group(3, "H_acceptor", [_atom(6, "O3", "O"),
                                        _atom(7, "O4", "O")]))]
    return Interaction(interaction_type=TYPE, groups=groups,
                       existence=np.ones((PAIRS, FRAMES), dtype=bool),
                       metrics={"distance": np.full((PAIRS, FRAMES), 3.0)},
                       times=np.zeros(FRAMES))


def _run_meta(topology: Path, trajectory: Path) -> RunMeta:
    return RunMeta(topology=str(topology), trajectory=str(trajectory),
                   output="", ff="amber", strategy="two_pass",
                   interactions=(TYPE,), argv=("dii", "run"),
                   cwd="", created_at="2026-11-26T00:00:00",
                   dii_version="0.1.0")


def _save_h5(directory: Path, run=None) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    save_interactions([_interaction()], str(directory / f"{TYPE}.h5"), run=run)
    return directory


@pytest.fixture(autouse=True)
def isolate_recent(monkeypatch, tmp_path):
    monkeypatch.setattr("DuIvyInteractions.visualizers.data.RECENT_PATH",
                        tmp_path / "recent.json")


def _client(**kwargs):
    from fastapi.testclient import TestClient
    from DuIvyInteractions.web.app import create_app
    return TestClient(create_app(**kwargs))


def _select(client, directory: Path) -> None:
    client.post("/select", data={"path": str(directory)})


def _config(html: str) -> dict:
    match = re.search(CONFIG_PATTERN, html, re.S)
    assert match, "structureConfig script is missing"
    return json.loads(match.group(1))


def _details_tag(html: str, marker: str) -> str:
    match = re.search(r"<details[^>]*%s[^>]*>" % re.escape(marker), html)
    assert match, f"missing <details> with {marker}"
    return match.group(0)


@pytest.fixture
def bare_client(tmp_path):
    client = _client()
    _select(client, _save_h5(tmp_path / "bare"))
    return client


@pytest.fixture
def run_client(tmp_path):
    """结构已选且与 h5 的 run/ 记录一致 → auto_attempt False（来源面板折叠）。"""
    topology = _write_pdb(tmp_path / "md.pdb")
    trajectory = _write_pdb(tmp_path / "md.xtc.pdb")
    client = _client()
    _select(client, _save_h5(tmp_path / "run",
                             run=_run_meta(topology, trajectory)))
    response = client.post("/api/structure/select", data={
        "topology": str(topology), "trajectory": str(trajectory)})
    assert response.status_code == 200
    return client


# ============================================================
# 区块顺序与结构
# ============================================================


class TestBlockOrder:
    def test_blocks_appear_top_to_bottom(self, bare_client):
        html = bare_client.get(PAIR_URL).text
        positions = [html.find(block) for block in BLOCKS]
        assert all(position >= 0 for position in positions), positions
        assert positions == sorted(positions), positions

    def test_header_shows_frames(self, bare_client):
        html = bare_client.get(PAIR_URL).text
        assert re.search(r'<span id="pairMeta">pair 1 / 2 · 3 frames', html)

    def test_both_asset_bundles_loaded(self, bare_client):
        html = bare_client.get(PAIR_URL).text
        for asset in ("/static/vendor/molstar.js", "/static/structure.js",
                      "/static/plotly-basic.min.js", "/static/pair.js"):
            assert asset in html, asset


# ============================================================
# 无视图切换按钮
# ============================================================


class TestNoViewSwitch:
    def test_no_view_switch_left(self, bare_client):
        html = bare_client.get(PAIR_URL).text
        assert "view-switch" not in html
        assert "btn-active" not in html
        assert "?view=" not in html

    def test_view_param_renders_same_layout(self, bare_client):
        plain = bare_client.get(PAIR_URL).text
        legacy = bare_client.get(PAIR_URL, params={"view": "3d"}).text
        assert legacy == plain


# ============================================================
# 折叠面板
# ============================================================


class TestPanels:
    def test_frames_panel_is_plain_panel_with_pager_in_head(self, bare_client):
        html = bare_client.get(PAIR_URL).text
        assert '<div class="panel frames-panel" id="framesPanel">' in html
        assert '<details class="panel frames-panel"' not in html
        head = re.search(
            r'<div class="section-head compact">\s*<h2>frames</h2>\s*'
            r'<div class="actions">.*?id="framesPrev".*?id="framesInfo"'
            r'.*?id="framesNext".*?</div>\s*</div>', html, re.S)
        assert head, "pager must live in the frames title row"
        assert 'id="framesTable"' in html

    def test_source_panel_open_without_structure(self, bare_client):
        html = bare_client.get(PAIR_URL).text
        assert "open" in _details_tag(html, "structure-source")
        assert "structure source · not selected" in html
        assert "structure not selected" in html

    def test_source_panel_with_run_record(self, run_client):
        html = run_client.get(PAIR_URL).text
        assert "structure source · md.pdb + md.xtc.pdb" in html
        assert _config(html)["auto_attempt"] is False
        assert "open" not in _details_tag(html, "structure-source")


# ============================================================
# 内容不丢
# ============================================================


class TestContentPreserved:
    def test_summary_groups_curves_frames_present(self, bare_client):
        html = bare_client.get(PAIR_URL).text
        for token in ("occupancy", 'class="meta-grid"', "meta-table",
                      'id="pairChart"', 'id="pairFigure"', 'id="pairConfig"',
                      'id="framesTable"'):
            assert token in html, token

    def test_viewer_block_present(self, bare_client):
        html = bare_client.get(PAIR_URL).text
        for element_id in ("molstarViewer", "frameLabel", "existenceBar",
                           "groupLegend", "structureStatus", "structureConfig",
                           "structureForm", "topologyInput", "trajectoryInput",
                           "openStructure"):
            assert 'id="%s"' % element_id in html, element_id


# ============================================================
# structure source 表单高度对齐
# ============================================================


class TestStructureFormHeight:
    def test_input_and_button_share_height(self):
        css = APP_CSS.read_text(encoding="utf-8")
        assert (".structure-form input,\n.structure-form .btn "
                "{ height: 34px; box-sizing: border-box; }") in css
        assert ".structure-form input { padding: 0 10px; }" in css
        assert ".structure-form .btn { padding: 0 12px; }" in css
