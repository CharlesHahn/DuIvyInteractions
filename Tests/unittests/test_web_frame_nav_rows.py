# -*- coding: utf-8 -*-
"""类型页帧导航控件 + aggregates 行按需加载（/api/rows 与 rows.js）。

覆盖三块：
1. `/api/rows/<name>` 行窗口（offset/limit/筛选/排序）与 `/type` 聚合表一致；
2. 类型页模板的帧导航控件（frameInput/framePrev/frameNext）与行加载容器；
3. `rows.js` / `matrix.js` 静态契约与 app.css 追加样式。
"""

import json
import re
from pathlib import Path

import numpy as np
import pytest

from DuIvyInteractions.core.datas import AtomData, Group, Interaction
from DuIvyInteractions.io.h5 import save_interactions

TYPE = "hydrogen_bond"
PAIRS = 240
FRAMES = 6
FIRST_PAGE = 50
REPO_ROOT = Path(__file__).resolve().parents[2]
STATIC_DIR = REPO_ROOT / "DuIvyInteractions" / "web" / "static"
TEMPLATE_DIR = REPO_ROOT / "DuIvyInteractions" / "web" / "templates"
ROWS_CONFIG_RE = re.compile(
    r'<script id="rowsConfig" type="application/json">(.*?)</script>', re.S)
ROWS_BODY_RE = re.compile(r'<tbody id="rowsBody">(.*?)</tbody>', re.S)
AGG_TABLE_RE = re.compile(
    r'<div class="table-wrap short">\s*<table>(.*?)</table>', re.S)
TH_RE = re.compile(r"<th>(.*?)</th>", re.S)
TD_RE = re.compile(r"<td>\s*(.*?)\s*</td>", re.S)
TAG_RE = re.compile(r"<[^>]+>")
TOTAL_RE = re.compile(r"aggregates · (\d+) pairs")
FRAME_INPUT_RE = re.compile(r'<input[^>]*id="frameInput"[^>]*>')


@pytest.fixture(autouse=True)
def isolate_recent(monkeypatch, tmp_path):
    monkeypatch.setattr("DuIvyInteractions.visualizers.data.RECENT_PATH",
                        tmp_path / "recent.json")


def _atom(gidx: int, name: str, element: str) -> AtomData:
    return AtomData(atom_global_idx=gidx, atom_idx_in_residue=0, atom_name=name,
                    atom_type="", atom_element=element, atom_charge=0.0,
                    atom_mass=0.0)


def _group(gid: int, gtype: str, res: str, rid: int,
           atoms) -> Group:
    return Group(group_id=gid, group_type=gtype, molecule="test",
                 residue_name=res, residue_id=rid, atoms=atoms)


def _pair(i: int):
    """第 i 对：TYR{i} 供体（N-H）+ ASP{1000+i} 受体（O）。"""
    donor = _group(2 * i, "H_donor", "TYR", i,
                   [_atom(100000 + 2 * i, "N", "N"),
                    _atom(100000 + 2 * i + 1, "H", "H")])
    acceptor = _group(2 * i + 1, "H_acceptor", "ASP", 1000 + i,
                      [_atom(200000 + 2 * i, "O", "O")])
    return donor, acceptor


def _interaction() -> Interaction:
    """240 对 × 6 帧：occupancy 只有 4/6 与 5/6 两档，末对无活跃帧。"""
    existence = np.zeros((PAIRS, FRAMES), dtype=bool)
    for i in range(PAIRS):
        for frame in range(FRAMES):
            existence[i, frame] = (i + frame) % 4 != 0
    existence[PAIRS - 1, :] = False  # 空指标显示（—）与页面一致性用例
    metrics = {"distance": np.where(
        existence, 3.0 + np.arange(PAIRS, dtype=float)[:, None] / 1000.0,
        np.nan)}
    return Interaction(interaction_type=TYPE,
                       groups=[_pair(i) for i in range(PAIRS)],
                       existence=existence, metrics=metrics,
                       times=np.arange(FRAMES, dtype=float) * 10.0)


@pytest.fixture(scope="module")
def h5_dir(tmp_path_factory):
    directory = tmp_path_factory.mktemp("rows_out")
    save_interactions([_interaction()], str(directory / (TYPE + ".h5")))
    return directory


@pytest.fixture(scope="module")
def client(h5_dir):
    from fastapi.testclient import TestClient
    from DuIvyInteractions.web.app import create_app
    client = TestClient(create_app())
    client.post("/select", data={"path": str(h5_dir)})
    return client


def _rows(client, query: str = "") -> dict:
    response = client.get("/api/rows/" + TYPE + query)
    assert response.status_code == 200, response.text
    return response.json()


def _page(client, query: str = "") -> str:
    response = client.get("/type/" + TYPE + query)
    assert response.status_code == 200, response.text
    return response.text


def _page_total(html: str) -> int:
    match = TOTAL_RE.search(html)
    assert match, "aggregates heading missing"
    return int(match.group(1))


def _cell_texts(html: str):
    body = ROWS_BODY_RE.search(html).group(1)
    return [TAG_RE.sub("", cell).strip() for cell in TD_RE.findall(body)]


def _header_texts(html: str):
    table = AGG_TABLE_RE.search(html)
    assert table, "aggregates table missing"
    return [TAG_RE.sub("", head).strip() for head in TH_RE.findall(table.group(1))]


def _button(html: str, element_id: str) -> str:
    pattern = re.compile(r'<button[^>]*id="%s"[^>]*>' % re.escape(element_id))
    match = pattern.search(html)
    assert match, "%s button missing" % element_id
    return match.group(0)


def _all_rows(client) -> list:
    """按窗口顺序拼接全表（用于计算期望的筛选结果）。"""
    first = _rows(client, "?limit=200")["rows"]
    rest = _rows(client, "?offset=200&limit=200")["rows"]
    return first + rest


# ============================================================
# /api/rows：行窗口
# ============================================================


class TestRowsApiWindow:

    def test_default_window_is_first_page(self, client):
        data = _rows(client)
        assert data["offset"] == 0
        assert data["limit"] == FIRST_PAGE
        assert data["total"] == PAIRS
        assert len(data["rows"]) == FIRST_PAGE
        assert data["columns"][:3] == ["index", "pair_label", "occupancy"]
        assert data["units"]["distance"] == "Å"

    def test_offset_returns_next_slice(self, client):
        ordered = _all_rows(client)
        assert len(ordered) == PAIRS
        window = _rows(client, "?offset=10&limit=50")["rows"]
        assert window == ordered[10:60]

    def test_limit_capped_at_max(self, client):
        data = _rows(client, "?limit=999")
        assert data["limit"] == 200
        assert len(data["rows"]) == 200

    def test_negative_offset_and_zero_limit_clamped(self, client):
        data = _rows(client, "?offset=-5&limit=3")
        assert data["offset"] == 0
        assert [row["idx"] for row in data["rows"]] == \
            [row["idx"] for row in _rows(client, "?limit=3")["rows"]]
        single = _rows(client, "?limit=0")
        assert single["limit"] == 1
        assert len(single["rows"]) == 1

    def test_non_numeric_window_falls_back_to_defaults(self, client):
        data = _rows(client, "?offset=abc&limit=xyz")
        assert data["offset"] == 0
        assert data["limit"] == FIRST_PAGE
        assert len(data["rows"]) == FIRST_PAGE

    def test_invalid_min_occupancy_falls_back_to_zero(self, client):
        assert _rows(client, "?min_occupancy=nan")["total"] == PAIRS
        assert _rows(client, "?min_occupancy=")["total"] == PAIRS

    def test_offset_at_total_returns_empty(self, client):
        data = _rows(client, "?offset=%d" % PAIRS)
        assert data["total"] == PAIRS
        assert data["offset"] == PAIRS
        assert data["rows"] == []

    def test_unknown_type_returns_json_404(self, client):
        response = client.get("/api/rows/not_a_type")
        assert response.status_code == 404
        assert response.headers["content-type"].startswith("application/json")

    def test_entry_shape_aligns_with_columns(self, client):
        data = _rows(client, "?limit=1")
        row = data["rows"][0]
        assert set(row) == {"idx", "label", "values"}
        assert len(row["values"]) == len(data["columns"])
        assert row["values"][0] == str(row["idx"])
        assert row["values"][1] == row["label"]


# ============================================================
# /api/rows：筛选/排序与 /type 聚合表一致
# ============================================================


class TestRowsApiConsistency:

    def test_columns_and_units_match_page_header(self, client):
        data = _rows(client)
        expected = []
        for column in data["columns"]:
            base = (column[4:] if column.startswith(("avg_", "std_")) else "")
            unit = data["units"].get(base, "")
            expected.append(column + (" (%s)" % unit if base and unit else ""))
        assert expected == _header_texts(_page(client))

    def test_values_match_page_first_row(self, client):
        data = _rows(client)
        cells = _cell_texts(_page(client))[:len(data["columns"])]
        assert data["rows"][0]["values"] == cells

    def test_empty_metric_cells_match_page(self, client):
        """无活跃帧的对：占有率 0.0%，avg/std 显示 —（与页面同一格式化）。"""
        query = "?search=TYR239"
        api = _rows(client, query + "&limit=1")
        assert api["total"] == 1
        row = api["rows"][0]
        assert row["values"][2] == "0.0%"
        assert row["values"][3] == "—"
        assert row["values"][4] == "—"
        cells = _cell_texts(_page(client, query))
        assert row["values"] == cells[:len(api["columns"])]

    @pytest.mark.parametrize("query", (
        "", "?search=TYR11", "?min_occupancy=0.8",
        "?search=TYR11&min_occupancy=0.8"))
    def test_total_matches_page(self, client, query):
        assert _rows(client, query)["total"] == _page_total(_page(client, query))

    def test_search_filters_and_matches_page(self, client):
        data = _rows(client, "?search=TYR11&limit=200")
        expected = sum("tyr11" in row["label"].lower()
                       for row in _all_rows(client))
        assert data["total"] == expected == 11
        assert all("tyr11" in row["label"].lower() for row in data["rows"])

    def test_min_occupancy_filters_and_matches_page(self, client):
        data = _rows(client, "?min_occupancy=0.8&limit=200")
        assert _page_total(_page(client, "?min_occupancy=0.8")) == data["total"]
        assert 0 < data["total"] < PAIRS
        assert all(float(row["values"][2].rstrip("%")) >= 80.0
                   for row in data["rows"])

    @pytest.mark.parametrize("query", (
        "?sort_by=label&sort_dir=asc",
        "?sort_by=index&sort_dir=desc",
        "?sort_by=occupancy&sort_dir=asc"))
    def test_sort_order_matches_page(self, client, query):
        rows = _rows(client, query + "&limit=5")["rows"]
        cells = _cell_texts(_page(client, query))
        width = len(rows[0]["values"])
        assert rows[0]["values"] == cells[:width]

    def test_index_sort_desc_puts_last_pair_first(self, client):
        row = _rows(client, "?sort_by=index&sort_dir=desc&limit=1")["rows"][0]
        assert row["idx"] == PAIRS - 1


# ============================================================
# 模板：帧导航 + 行加载容器
# ============================================================


class TestTypeTemplate:

    def test_frame_nav_and_rows_controls_present(self, client):
        html = _page(client)
        for element_id in ("frameInput", "framePrev", "frameNext",
                           "rowsMore", "rowsConfig",
                           "rowsBody", "rowsStatus"):
            assert 'id="%s"' % element_id in html, element_id

    def test_first_frame_disables_prev_only(self, client):
        html = _page(client)
        assert "disabled" in _button(html, "framePrev")
        assert "disabled" not in _button(html, "frameNext")

    def test_last_frame_disables_next_only(self, client):
        html = _page(client, "?frame=%d" % (FRAMES - 1))
        assert "disabled" not in _button(html, "framePrev")
        assert "disabled" in _button(html, "frameNext")

    def test_frame_input_carries_clamped_range(self, client):
        tag = FRAME_INPUT_RE.search(_page(client, "?frame=99999")).group(0)
        assert 'type="number"' in tag
        assert 'min="0"' in tag
        assert 'max="%d"' % (FRAMES - 1) in tag
        assert 'value="%d"' % (FRAMES - 1) in tag

    def test_frame_number_anchor_removed(self, client):
        """帧号只由 frameInput 承载，页面不再保留隐藏的 frameNumber 锚点。"""
        template = (TEMPLATE_DIR / "type.html").read_text(encoding="utf-8")
        assert "frameNumber" not in template
        assert "frameNumber" not in _page(client, "?frame=5")

    def test_rows_config_matches_query(self, client):
        html = _page(client,
                     "?search=TYR11&min_occupancy=0.8"
                     "&sort_by=index&sort_dir=asc")
        cfg = json.loads(ROWS_CONFIG_RE.search(html).group(1))
        assert cfg["url"] == "/api/rows/" + TYPE
        assert cfg["type"] == TYPE
        assert cfg["search"] == "TYR11"
        assert cfg["min_occupancy"] == 0.8
        assert cfg["sort_by"] == "index" and cfg["sort_dir"] == "asc"
        assert cfg["total"] == _page_total(html)
        assert cfg["limit"] == FIRST_PAGE
        assert cfg["offset"] == 0

    def test_page_two_first_screen_offsets_config(self, client):
        """?page=2 仍被后端接受：首屏 50 行从第 50 行开始，JS 续载不重叠。"""
        html = _page(client, "?page=2")
        cfg = json.loads(ROWS_CONFIG_RE.search(html).group(1))
        assert cfg["offset"] == FIRST_PAGE
        assert "showing %d of %d" % (2 * FIRST_PAGE, PAIRS) in html

    def test_first_page_is_server_rendered(self, client):
        html = _page(client)
        assert "showing %d of %d" % (FIRST_PAGE, PAIRS) in html
        assert "load more (%d remaining)" % (PAIRS - FIRST_PAGE) in html
        assert len(_cell_texts(html)) == FIRST_PAGE * len(_rows(client)["columns"])

    def test_old_pager_removed_from_aggregates_panel(self, client):
        html = _page(client)
        start = html.index("<h2>aggregates · ")
        panel = html[start:html.index('id="rowsConfig"')]
        assert "‹ previous" not in panel
        assert "next ›" not in panel
        assert 'class="pager"' not in html
        assert "page 1 /" not in html

    def test_rows_script_loaded(self, client):
        html = _page(client)
        assert '<script src="/static/rows.js"></script>' in html
        assert html.index('src="/static/matrix.js"') < \
            html.index('src="/static/rows.js"')


# ============================================================
# 静态资源契约
# ============================================================


class TestStaticAssets:

    def test_rows_js_contract(self):
        source = (STATIC_DIR / "rows.js").read_text(encoding="utf-8")
        for token in ("/api/rows/", "offset=", "limit=", "rowsMore",
                      "rowsBody", "rowsStatus", "load more",
                      "createElement", "textContent"):
            assert token in source, token
        assert '"&offset=" + (baseOffset + rendered)' in source
        assert '"all " + total + " rows loaded"' in source
        assert "catch(function (error)" in source
        assert "failed to load rows: " in source
        assert "innerHTML" not in source

    def test_matrix_js_frame_nav_contract(self):
        source = (STATIC_DIR / "matrix.js").read_text(encoding="utf-8")
        for token in ("frameInput", "framePrev", "frameNext",
                      "stepFrame", "syncFrameNav", "commitFrameInput"):
            assert token in source, token
        assert 'frameInput.addEventListener("keydown"' in source
        assert 'frameInput.addEventListener("blur"' in source
        # 既有 URL / 帧号同步行为保持
        assert "history.replaceState" in source
        assert "frameInput.value = frame" in source
        assert "frameNumber" not in source

    def test_css_has_frame_nav_and_rows_status(self):
        css = (STATIC_DIR / "app.css").read_text(encoding="utf-8")
        for selector in (".frame-nav", ".frame-head #frameInput",
                         ".rows-more", ".rows-status"):
            assert selector in css, selector
