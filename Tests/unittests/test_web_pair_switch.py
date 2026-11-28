# -*- coding: utf-8 -*-
"""单对页页内切 pair：/panels 接口、三段 HTML 同源、前端无跳转切换。

用真实 h5（测试数据目录）渲染，不依赖浏览器；前端行为只做静态断言。
"""

import json
import re
from pathlib import Path

import pytest

H5_DATA_DIR = Path(__file__).resolve().parents[1] / "interaction_h5data"
WEB_DIR = Path(__file__).resolve().parents[2] / "DuIvyInteractions" / "web"
PAIR_JS = WEB_DIR / "static" / "pair.js"
STRUCTURE_JS = WEB_DIR / "static" / "structure.js"
PAIR_TEMPLATE = WEB_DIR / "templates" / "pair.html"
TYPE = "pi_stacking"
PAIR_URL = f"/pair/{TYPE}/1"
PANEL_KEYS = {"idx", "n_pairs", "head", "summary", "groups"}
SWAP_TARGETS = (("head", "pairHead"), ("summary", "pairSummary"),
                ("groups", "pairGroups"))
FIGURE_PATTERN = r'<script id="pairFigure"[^>]*>(.*?)</script>'


@pytest.fixture(autouse=True)
def isolate_recent(monkeypatch, tmp_path):
    monkeypatch.setattr("DuIvyInteractions.visualizers.data.RECENT_PATH",
                        tmp_path / "recent.json")


def _client():
    from fastapi.testclient import TestClient
    from DuIvyInteractions.web.app import create_app
    return TestClient(create_app())


@pytest.fixture
def client():
    client = _client()
    client.post("/select", data={"path": str(H5_DATA_DIR)})
    return client


def _panels(client, idx: int = 1, name: str = TYPE) -> dict:
    response = client.get(f"/api/pair/{name}/{idx}/panels")
    assert response.status_code == 200
    return response.json()


# ============================================================
# /api/pair/{name}/{idx}/panels
# ============================================================


class TestPanelsEndpoint:
    def test_payload_shape(self, client):
        data = _panels(client)
        assert set(data) == PANEL_KEYS
        assert data["idx"] == 1
        assert data["n_pairs"] == 10

    def test_fragments_carry_swap_targets(self, client):
        data = _panels(client)
        assert 'class="section-head"' in data["head"]
        assert 'id="pairTitle"' in data["head"]
        assert 'id="pairMeta"' in data["head"]
        assert 'id="pairSummaryBody"' in data["summary"]
        assert 'id="pairGroupsBody"' in data["groups"]

    def test_fragments_match_whole_page(self, client):
        """三段 HTML 与整页 include 的输出逐字一致（外包容器后子串相等）。"""
        page = client.get(PAIR_URL).text
        data = _panels(client)
        for key, wrapper in SWAP_TARGETS:
            block = '<div id="%s">%s</div>' % (wrapper, data[key])
            assert block in page, key

    def test_neighbour_pair_differs(self, client):
        assert _panels(client, 0)["summary"] != _panels(client, 1)["summary"]

    def test_out_of_range_returns_404(self, client):
        assert client.get(
            f"/api/pair/{TYPE}/10/panels").status_code == 404

    def test_unknown_type_returns_404(self, client):
        assert client.get("/api/pair/nope/0/panels").status_code == 404

    def test_pair_data_contract_unchanged(self, client):
        from DuIvyInteractions.web.app import PAIR_API_KEYS
        data = client.get(f"/api/pair/{TYPE}/1").json()
        assert set(data) == set(PAIR_API_KEYS)
        page = client.get(PAIR_URL).text
        match = re.search(FIGURE_PATTERN, page, re.S)
        assert match and json.loads(match.group(1)) == data["figure"]


# ============================================================
# 模板：可替换容器 + 无 JS 兜底
# ============================================================


class TestPairTemplate:
    def test_page_has_swap_targets_and_nav(self, client):
        html = client.get(PAIR_URL).text
        for _, wrapper in SWAP_TARGETS:
            assert 'id="%s"' % wrapper in html, wrapper
        assert html.count("data-pair-nav") == 2
        assert 'id="pairPrev"' in html and 'id="pairNext"' in html
        prev = re.search(r'<a[^>]*id="pairPrev"[^>]*>', html).group(0)
        assert 'href="/pair/%s/0?frame=0"' % TYPE in prev
        next_link = re.search(r'<a[^>]*id="pairNext"[^>]*>', html).group(0)
        assert 'href="/pair/%s/2?frame=0"' % TYPE in next_link

    def test_no_js_still_renders_content(self, client):
        html = client.get(PAIR_URL).text
        for token in ("summary", "occupancy", "meta-table", 'id="pairChart"',
                      'id="framesTable"', 'id="pairFigure"'):
            assert token in html, token

    def test_head_meta_wrapped_for_swap(self):
        text = PAIR_TEMPLATE.read_text(encoding="utf-8")
        assert '{% include "_pair_head.html" %}' in text
        assert '{% include "_pair_summary.html" %}' in text
        assert '{% include "_pair_groups.html" %}' in text


# ============================================================
# pair.js：页内切换（静态断言）
# ============================================================


def _pair_js() -> str:
    return PAIR_JS.read_text(encoding="utf-8")


def _function_body(text: str, name: str) -> str:
    match = re.search(r"  function %s\(.*?\n  \}\n" % re.escape(name), text,
                      re.S)
    assert match, f"missing function {name}"
    return match.group(0)


class TestPairScript:
    def test_calls_panels_and_pair_api(self):
        text = _pair_js()
        body = _function_body(text, "switchPair")
        assert "/panels" in text
        assert "panelUrl(" in body and "dataUrl(" in body
        assert "Promise.all" in body

    def test_replaces_sections_without_innerhtml(self):
        text = _pair_js()
        assert "innerHTML" not in text
        body = _function_body(text, "replaceSection")
        assert "DOMParser" in body
        assert "replaceChildren" in body

    def test_updates_chart_in_place(self):
        text = _pair_js()
        body = _function_body(text, "reactChart")
        assert "Plotly.react" in body
        assert "PLOT_CONFIG" in body

    def test_calls_structure_update_pair(self):
        text = _pair_js()
        body = _function_body(text, "updateStructure")
        assert "window.DiiStructure" in body
        assert "api.updatePair(payload)" in body
        assert "isMounted()" in body

    def test_updates_history_without_navigation(self):
        text = _pair_js()
        assert "history.replaceState" in text
        for token in ("location.href", "location.assign", "location.reload"):
            assert token not in text, token
        body = _function_body(text, "syncHistory")
        assert '"/pair/" + encodeURIComponent(cfg.type) + "/" + target' in body

    def test_switching_guard_and_frame_reset(self):
        text = _pair_js()
        body = _function_body(text, "switchPair")
        assert "if (switching" in body
        assert "switching = false;" in body
        apply = _function_body(text, "applyPairData")
        assert "frames.state.start = 0;" in apply
        assert "cfg.idx = target;" in apply

    def test_nav_binding_uses_data_attribute(self):
        text = _pair_js()
        assert 'closest("[data-pair-nav]")' in text
        assert "switchPair(cfg.idx + Number(" in text


# ============================================================
# structure.js：updatePair 重建背景层
# ============================================================


class TestStructureBackgroundRebuild:
    def test_update_pair_rebuilds_background(self):
        text = STRUCTURE_JS.read_text(encoding="utf-8")
        body = _function_body(text, "updatePair")
        assert "buildBackgroundLayers()" in body
        background = _function_body(text, "buildBackgroundLayers")
        assert "buildPolymerLayer(" in background
        assert "buildLigandLayer(" in background
        for token in ("TOPOLOGY_URL", "TRAJECTORY_URL", "loadTrajectoryOnce("):
            assert token not in body, token

    def test_clears_background_refs_on_pair_switch(self):
        text = STRUCTURE_JS.read_text(encoding="utf-8")
        body = _function_body(text, "clearResidueLayers")
        assert "backgroundRefs" in body
        assert "residueRefs = [];" in body
