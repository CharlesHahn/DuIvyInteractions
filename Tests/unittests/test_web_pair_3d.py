# -*- coding: utf-8 -*-
"""单对页 3D 区块测试：/structure 302 → /pair，单页始终内嵌 Mol* viewer。

布局合并后不再有 curves/3D 视图切换；?view=3d 仅为旧链接兼容而接受。
不依赖真实轨迹：用 tmp_path 构造最小 h5，覆盖"未选结构"的退化分支。
"""

import json
import re
from pathlib import Path

import numpy as np
import pytest

from DuIvyInteractions.core.datas import AtomData, Group, Interaction
from DuIvyInteractions.io.h5 import save_interactions

WEB_DIR = Path(__file__).resolve().parents[2] / "DuIvyInteractions" / "web"
TEMPLATES_DIR = WEB_DIR / "templates"
STRUCTURE_JS = WEB_DIR / "static" / "structure.js"
PAIR_TEMPLATE = TEMPLATES_DIR / "pair.html"
PARTIAL_TEMPLATE = TEMPLATES_DIR / "_structure_view.html"

FRAMES = 12
TYPE = "hydrogen_bond"
CONFIG_PATTERN = r'<script id="structureConfig"[^>]*>(.*?)</script>'
PAIR_URL = f"/pair/{TYPE}/0"


def _atom(index: int, name: str, atype: str = "") -> AtomData:
    return AtomData(atom_global_idx=index, atom_idx_in_residue=0,
                    atom_name=name, atom_type=atype, atom_element="C",
                    atom_charge=0.0, atom_mass=0.0)


def _group(group_id: int, group_type: str, atoms) -> Group:
    return Group(group_id=group_id, group_type=group_type, molecule="test",
                 residue_name="LYS", residue_id=1, atoms=list(atoms))


def _interaction() -> Interaction:
    existence = np.ones((2, FRAMES), dtype=bool)
    groups = [(_group(0, "H_donor", [_atom(0, "N1", "N"),
                                     _atom(1, "H1", "H")]),
               _group(1, "H_acceptor", [_atom(2, "O1", "O")])),
              (_group(2, "H_donor", [_atom(3, "N2", "N"),
                                     _atom(4, "H2", "H")]),
               _group(3, "H_acceptor", [_atom(5, "O2", "O")]))]
    return Interaction(interaction_type=TYPE, groups=groups,
                       existence=existence,
                       metrics={"distance": np.full((2, FRAMES), 3.0)},
                       times=np.zeros(FRAMES))


@pytest.fixture(autouse=True)
def isolate_recent(monkeypatch, tmp_path):
    monkeypatch.setattr("DuIvyInteractions.visualizers.data.RECENT_PATH",
                        tmp_path / "recent.json")


@pytest.fixture
def h5_dir(tmp_path):
    directory = tmp_path / "out"
    directory.mkdir()
    save_interactions([_interaction()], str(directory / f"{TYPE}.h5"))
    return directory


def _client(**kwargs):
    from fastapi.testclient import TestClient
    from DuIvyInteractions.web.app import create_app
    return TestClient(create_app(**kwargs))


@pytest.fixture
def client(h5_dir):
    client = _client()
    client.post("/select", data={"path": str(h5_dir)})
    return client


def _config(html: str) -> dict:
    match = re.search(CONFIG_PATTERN, html, re.S)
    assert match, "structureConfig script is missing"
    return json.loads(match.group(1))


def _structure_js() -> str:
    return STRUCTURE_JS.read_text(encoding="utf-8")


# ============================================================
# /structure 旧路由 → 单对页 3D 视图
# ============================================================


class TestStructureRedirect:
    def test_redirects_to_pair_3d_view(self, client):
        response = client.get("/structure/hydrogen_bond/0",
                              follow_redirects=False)
        assert response.status_code == 302
        assert response.headers["location"] == f"{PAIR_URL}?view=3d"

    def test_redirect_keeps_caller_frame(self, client):
        response = client.get("/structure/hydrogen_bond/0?frame=5",
                              follow_redirects=False)
        assert response.status_code == 302
        location = response.headers["location"]
        assert location.startswith(f"{PAIR_URL}?view=3d")
        assert "frame=5" in location

    def test_redirect_keeps_explicit_zero_frame(self, client):
        location = client.get("/structure/hydrogen_bond/0?frame=0",
                              follow_redirects=False).headers["location"]
        assert "frame=0" in location

    def test_redirect_target_renders(self, client):
        response = client.get("/structure/hydrogen_bond/0")
        assert response.status_code == 200
        assert 'id="molstarViewer"' in response.text


# ============================================================
# 单页合并：曲线与 3D 同页，始终加载两套资源
# ============================================================


class TestMergedPairPage:
    def test_pair_page_always_embeds_viewer(self, client):
        response = client.get(PAIR_URL)
        assert response.status_code == 200
        html = response.text
        for element_id in ("molstarViewer", "existenceBar", "groupLegend",
                           "structureStatus", "structureForm",
                           "topologyInput", "trajectoryInput",
                           "openStructure", "structureConfig",
                           "pairChart", "framesTable"):
            assert 'id="%s"' % element_id in html, element_id

    def test_pair_page_loads_both_bundles(self, client):
        html = client.get(PAIR_URL).text
        molstar_at = html.index("/static/vendor/molstar.js")
        structure_at = html.index("/static/structure.js")
        assert molstar_at < structure_at
        for asset in ("/static/vendor/molstar.css", "/static/vendor/molstar.js",
                      "/static/structure.js", "/static/plotly-basic.min.js",
                      "/static/pair.js"):
            assert asset in html, asset

    def test_view_param_no_longer_changes_rendering(self, client):
        plain = client.get(PAIR_URL)
        legacy = client.get(PAIR_URL, params={"view": "3d"})
        assert legacy.status_code == 200
        assert legacy.text == plain.text

    def test_illegal_view_renders_same_merged_page(self, client):
        response = client.get(PAIR_URL, params={"view": "xxx"})
        assert response.status_code == 200
        assert 'id="molstarViewer"' in response.text
        assert 'id="pairChart"' in response.text

    def test_merged_page_allows_molstar_eval(self, client):
        csp = client.get(PAIR_URL).headers["content-security-policy"]
        assert "'unsafe-eval'" in csp
        assert "connect-src 'self' data:" in csp


# ============================================================
# 帧号：注入、钳制、存在性条初始状态
# ============================================================


class TestFrameState:
    def test_frame_is_injected(self, client):
        cfg = _config(client.get(PAIR_URL, params={"frame": 7}).text)
        assert cfg["frame"] == 7
        assert cfg["n_frames"] == FRAMES

    def test_frame_is_clamped(self, client):
        cfg = _config(client.get(PAIR_URL, params={"frame": 99999}).text)
        assert cfg["frame"] == cfg["n_frames"] - 1 == FRAMES - 1

    def test_negative_frame_is_clamped(self, client):
        cfg = _config(client.get(PAIR_URL, params={"frame": -5}).text)
        assert cfg["frame"] == 0

    def test_frame_label_and_bar_start_at_frame(self, client):
        html = client.get(PAIR_URL, params={"frame": 7}).text
        bar = re.search(r'id="existenceBar"[^>]*', html).group(0)
        assert 'aria-valuenow="7"' in bar


# ============================================================
# prev/next pair 链接：保留帧号，不再带 view
# ============================================================


class TestPairNavigation:
    def test_next_pair_link_keeps_frame(self, client):
        html = client.get(PAIR_URL, params={"frame": 7}).text
        assert f'href="/pair/{TYPE}/1?frame=7"' in html
        assert "?view=" not in html

    def test_previous_pair_link_keeps_frame(self, client):
        html = client.get(f"/pair/{TYPE}/1", params={"frame": 7}).text
        assert f'href="/pair/{TYPE}/0?frame=7"' in html


# ============================================================
# 退化分支：未选择拓扑/轨迹
# ============================================================


class TestDegenerateNoStructure:
    def test_page_without_structure_shows_hint(self, client):
        html = client.get(PAIR_URL).text
        assert "structure not selected" in html
        cfg = _config(html)
        assert cfg["has_structure"] is False
        assert cfg["n_frames"] == FRAMES

    def test_page_without_structure_keeps_page_panels(self, client):
        html = client.get(PAIR_URL).text
        assert "summary" in html and "occupancy" in html
        assert "groups" in html
        assert 'id="pairChart"' in html

    def test_structure_script_hints_structure_selection(self):
        text = _structure_js()
        assert "select topology and trajectory to load the structure" in text


# ============================================================
# structure.js：读入帧号 + 地址栏同步 + 非 3D 页面静默退出
# ============================================================


class TestStructureScript:
    def test_reads_injected_frame(self):
        text = _structure_js()
        assert "cfg.frame" in text
        assert "var frame = clampFrame(cfg.frame || 0)" in text

    def test_initial_frame_is_applied_after_load(self):
        text = _structure_js()
        assert "return selectFrame(frame);" in text

    def test_frame_change_updates_url(self):
        text = _structure_js()
        assert "history.replaceState" in text
        assert 'var FRAME_PARAM = "frame"' in text
        assert "searchParams.set(FRAME_PARAM" in text
        for name in ("requestFrame", "pollFrame"):
            match = re.search(r"  function %s\(.*?\n  \}\n" % name, text, re.S)
            assert match, name
            assert "syncFrameUrl(" in match.group(0), name

    def test_script_returns_without_viewer_host(self):
        text = _structure_js()
        assert 'var host = document.getElementById("molstarViewer")' in text
        assert "if (!host) return;" in text


# ============================================================
# 模板结构：3D 区块 partial（单对页 include，可独立渲染回归）
# ============================================================


class TestPartialTemplate:
    def test_partial_exists(self):
        assert PARTIAL_TEMPLATE.is_file()

    def test_pair_template_includes_partial(self):
        text = PAIR_TEMPLATE.read_text(encoding="utf-8")
        assert '{% include "_structure_view.html" %}' in text, PAIR_TEMPLATE
        assert '{% if view' not in text, "include must be unconditional"

    def test_partial_keeps_viewer_dom_ids(self):
        """viewer 类 id 已迁到 base 常驻坞；partial 只保留槽位与表单。"""
        text = PARTIAL_TEMPLATE.read_text(encoding="utf-8")
        for element_id in ("structureConfig", "structureForm", "filepickModal",
                           "topologyInput", "trajectoryInput",
                           "openStructure", "structureDockSlot"):
            assert 'id="%s"' % element_id in text, element_id
        base = (TEMPLATES_DIR / "base.html").read_text(encoding="utf-8")
        for element_id in ("molstarViewer", "frameLabel", "existenceBar",
                           "groupLegend", "structureStatus"):
            assert 'id="%s"' % element_id in base, element_id

    def test_partial_template_renders_standalone(self):
        from DuIvyInteractions.web.app import TEMPLATES

        html = TEMPLATES.get_template("_structure_view.html").render(
            title="structure", ds_name="", datasource="", types=[],
            active_type=None, type=TYPE, idx=0, pair_label="pair 1",
            n_atoms=4, n_frames=1, has_structure=False, existence=[1],
            topology="", trajectory="", browse_mode="picker")
        assert 'id="structureDockSlot"' in html
        assert 'id="structureConfig"' in html
        assert "structure not selected" in html

    def test_partial_source_panel_is_details(self):
        from DuIvyInteractions.web.app import TEMPLATES

        html = TEMPLATES.get_template("_structure_view.html").render(
            title="structure", ds_name="", datasource="", types=[],
            active_type=None, type=TYPE, idx=0, pair_label="pair 1",
            n_atoms=4, n_frames=1, has_structure=True, existence=[1],
            topology="/data/md.tpr", trajectory="/data/md.xtc",
            browse_mode="picker")
        assert '<details class="panel structure-source"' in html
        assert "structure source · md.tpr + md.xtc" in html
