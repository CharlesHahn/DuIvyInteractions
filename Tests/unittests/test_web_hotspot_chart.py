# -*- coding: utf-8 -*-
"""热点残基分布图 Web 测试：数据集页 /residues 与类型页 /type/<name> 共用 partial。

合成数据覆盖分段、柱色、图例、Y 轴、hover title 与空数据；
真实数据（Tests/interaction_h5data）覆盖水排除与类型页单色图例。
"""

import re
from pathlib import Path
from urllib.parse import quote

import numpy as np
import pytest
from fastapi.testclient import TestClient

from DuIvyInteractions.core.datas import AtomData, Group, Interaction
from DuIvyInteractions.io.h5 import save_interactions
from DuIvyInteractions.visualizers import (
    hotspot_segments, parse_datasource, residue_hotspots)
from DuIvyInteractions.web import app as web_app

FRAMES = 4
NAN = float("nan")
MOL_A = "seg_0_pro"
MOL_B = "seg_1_lig"
H5_DATA_DIR = Path(__file__).parent.parent / "interaction_h5data"
STATIC_CSS = (Path(__file__).resolve().parents[2]
              / "DuIvyInteractions" / "web" / "static" / "app.css")
CHART_MARKER = '<div class="residue-chart">'
GLOBAL_SCRIPTS = frozenset({"/static/load.js", "/static/browse.js"})
SEG_RE = re.compile(r'<div class="residue-seg" style="flex-grow: (\d+); flex-basis: 0">')
SEG_LABEL_RE = re.compile(r'<div class="residue-seg-label" title="([^"]*)">([^<]*)</div>')
BAR_GROUP_RE = re.compile(r'<div class="bars big residue-bars">(.*?)</div>', re.S)
BAR_RE = re.compile(
    r'<span style="height: ([\d.]+)%; background: ([^"]*)"\s*title="([^"]*)"', re.S)
LEGEND_RE = re.compile(
    r'residue-legend-dot" style="background: var\(--c-([a-z_]+)\)"></span>([a-z_]+)</span>')
CHIP_RE = re.compile(r'<span class="residue-seg-chip">([^<]*)</span>')

HB_SPECS = [
    [("H_donor", "LYS", 10, MOL_A), ("H_acceptor", "ASP", 100, MOL_A)],
    [("H_donor", "LYS", 20, MOL_A), ("H_acceptor", "GLU", 200, MOL_A)],
]
HB_EXIST = [[True, True, False, False], [False, False, True, False]]
HB_DIST = [[2.0, 2.5, NAN, NAN], [NAN, NAN, 3.0, NAN]]
HYDRO_SPECS = [[("hydrophobic", "LEU", 50, MOL_A),
                ("hydrophobic", "LIG", 60, MOL_B)]]
HYDRO_EXIST = [[True, True, True, True]]
HYDRO_DIST = [[3.5, 3.6, 3.4, 3.5]]
WATER_SPECS = [[("H_donor", "LYS", 10, MOL_A), ("water", "SOL", 500, MOL_B),
                ("H_acceptor", "ASP", 100, MOL_A)]]
WATER_EXIST = [[True, False, False, False]]
WATER_DIST = [[3.0, NAN, NAN, NAN]]


# ============================================================
# fixture 构造
# ============================================================

def _atom(gidx: int) -> AtomData:
    """占位原子（热点聚合不读原子）。"""
    return AtomData(atom_global_idx=gidx, atom_idx_in_residue=0, atom_name="C",
                    atom_type="", atom_element="C", atom_charge=0.0,
                    atom_mass=0.0)


def _group(gid: int, gtype: str, res: str, rid: int, molecule: str) -> Group:
    """按 (group_type, residue_name, residue_id, molecule) 造基团。"""
    return Group(group_id=gid, group_type=gtype, molecule=molecule,
                 residue_name=res, residue_id=rid, atoms=[_atom(gid)])


def _interaction(itype: str, specs, existence, distances) -> Interaction:
    """specs 每对为 [(group_type, residue_name, residue_id, molecule), ...]。"""
    groups = [tuple(_group(10 * i + pos, *spec)
                    for pos, spec in enumerate(pair))
              for i, pair in enumerate(specs)]
    return Interaction(interaction_type=itype, groups=groups,
                       existence=np.asarray(existence, dtype=bool),
                       metrics={"distance": np.asarray(distances, dtype=float)},
                       times=np.arange(FRAMES, dtype=float) * 10.0)


def _save(directory: Path, itype: str, specs, existence, distances) -> None:
    """一个类型写一个 <type>.h5。"""
    save_interactions([_interaction(itype, specs, existence, distances)],
                      str(directory / f"{itype}.h5"))


@pytest.fixture(autouse=True)
def isolate_recent(monkeypatch, tmp_path):
    """隔离最近数据源文件，避免写入用户缓存。"""
    monkeypatch.setattr("DuIvyInteractions.visualizers.data.RECENT_PATH",
                        tmp_path / "recent.json")


@pytest.fixture
def dataset(tmp_path):
    """两分子三分量数据集：蛋白段 5 残基 + 配体段 1 残基。"""
    directory = tmp_path / "small"
    directory.mkdir()
    _save(directory, "hydrogen_bond", HB_SPECS, HB_EXIST, HB_DIST)
    _save(directory, "hydrophobic", HYDRO_SPECS, HYDRO_EXIST, HYDRO_DIST)
    _save(directory, "water_bridge", WATER_SPECS, WATER_EXIST, WATER_DIST)
    return directory


@pytest.fixture
def zero_dir(tmp_path):
    """两对残基但全帧无命中（occupancy 恒为 0）。"""
    directory = tmp_path / "zero"
    directory.mkdir()
    _save(directory, "hydrogen_bond", HB_SPECS,
          [[False] * FRAMES for _ in HB_SPECS], HB_DIST)
    return directory


@pytest.fixture(scope="module")
def real_page(tmp_path_factory):
    """真实数据（8 类 h5）的总热点页与类型页 HTML（模块内复用一次载入）。"""
    if not (H5_DATA_DIR / "water_bridge.h5").exists():
        pytest.skip("real h5 test data not available")
    monkeypatch = pytest.MonkeyPatch()
    monkeypatch.setattr("DuIvyInteractions.visualizers.data.RECENT_PATH",
                        tmp_path_factory.mktemp("recent") / "recent.json")
    client = TestClient(web_app.create_app())
    response = client.post("/select", data={"path": str(H5_DATA_DIR.resolve())},
                           follow_redirects=False)
    assert response.status_code == 303
    ds = quote(str(H5_DATA_DIR.resolve()), safe="/")
    pages = {"residues": client.get(f"/residues?ds={ds}").text,
             "water_bridge": client.get("/type/water_bridge").text}
    monkeypatch.undo()
    return pages


@pytest.fixture(scope="module")
def real_segments():
    """真实数据的分段（后端 API 作为渲染规则基准）。"""
    if not (H5_DATA_DIR / "water_bridge.h5").exists():
        pytest.skip("real h5 test data not available")
    return hotspot_segments(residue_hotspots(parse_datasource(str(H5_DATA_DIR))))


@pytest.fixture(scope="module")
def real_water_segments():
    """真实数据水桥类型的分段（水已被排除）。"""
    if not (H5_DATA_DIR / "water_bridge.h5").exists():
        pytest.skip("real h5 test data not available")
    rows = residue_hotspots(parse_datasource(str(H5_DATA_DIR)),
                            type_filter="water_bridge")
    return hotspot_segments(rows)


# ============================================================
# 页面解析辅助
# ============================================================

def _select_client(*paths) -> TestClient:
    """创建应用并逐个载入数据源。"""
    client = TestClient(web_app.create_app())
    for path in paths:
        response = client.post("/select", data={"path": str(path)},
                               follow_redirects=False)
        assert response.status_code == 303
    return client


def _ds(directory: Path) -> str:
    """数据集的绝对路径（与 /select 记录的批次键一致）。"""
    return str(Path(directory).resolve())


def _url(ds: str, **params) -> str:
    """热点页 URL（ds 保留斜杠，其它参数完整转义）。"""
    query = "&".join(f"{key}={quote(str(value), safe='')}"
                     for key, value in params.items() if value is not None)
    base = f"/residues?ds={quote(ds, safe='/')}"
    return f"{base}&{query}" if query else base


def _get_page(client: TestClient, ds: str, **params) -> str:
    """GET 热点页并断言 200。"""
    response = client.get(_url(ds, **params))
    assert response.status_code == 200, response.text
    return response.text


def _chart(html: str) -> str:
    """页面里的残基分布图片段（无图时返回空串）。"""
    if CHART_MARKER not in html:
        return ""
    return html.split(CHART_MARKER, 1)[1].split('<div class="table-wrap', 1)[0]


def _segments(chart: str):
    """[(flex_grow, molecule_label), ...]（按 DOM 顺序）。"""
    labels = [(match.group(1), match.group(2))
              for match in SEG_LABEL_RE.finditer(chart)]
    grows = [int(value) for value in SEG_RE.findall(chart)]
    assert len(grows) == len(labels)
    return list(zip(grows, [text for _, text in labels]))


def _bar_groups(chart: str):
    """每段的柱组 HTML（列表元素为段内 span 串）。"""
    return BAR_GROUP_RE.findall(chart)


def _bars(chart: str):
    """[(height, background, title), ...]（全部段、DOM 顺序）。"""
    return [match.groups() for group in _bar_groups(chart)
            for match in BAR_RE.finditer(group)]


def _legend(chart: str):
    """图例中的类型名（按 DOM 顺序）。"""
    return [match.group(2) for match in LEGEND_RE.finditer(chart)]


def _rules(css: str, selector: str):
    """CSS 中某个选择器的声明块文本列表（.a 与 .a > b 不互相匹配）。"""
    pattern = re.compile(re.escape(selector) + r"\s*\{([^}]*)\}")
    return pattern.findall(css)


# ============================================================
# 数据集页：存在性与分段
# ============================================================

class TestDatasetChart:

    def test_chart_rendered_before_table(self, dataset):
        html = _get_page(_select_client(dataset), _ds(dataset))
        chart = _chart(html)
        assert chart, "residue-chart missing"
        assert "residue-track" in chart
        assert html.index("residue-chart") < html.index("hotspots-table")

    def test_segment_count_matches_segments(self, dataset):
        chart = _chart(_get_page(_select_client(dataset), _ds(dataset)))
        assert len(SEG_RE.findall(chart)) == 2

    def test_bar_count_matches_residue_sum(self, dataset):
        chart = _chart(_get_page(_select_client(dataset), _ds(dataset)))
        groups = _bar_groups(chart)
        assert [len(BAR_RE.findall(group)) for group in groups] == [5, 1]
        assert len(_bars(chart)) == 6

    def test_segment_labels_use_molecule_value(self, dataset):
        chart = _chart(_get_page(_select_client(dataset), _ds(dataset)))
        segments = _segments(chart)
        assert [label for _, label in segments] == [MOL_A, MOL_B]
        assert [grow for grow, _ in segments] == [5, 1]
        for molecule in (MOL_A, MOL_B):
            assert f'title="{molecule}"' in chart

    def test_segment_chip_list_matches_segments(self, dataset):
        chart = _chart(_get_page(_select_client(dataset), _ds(dataset)))
        chips = CHIP_RE.findall(chart)
        assert len(chips) == len(SEG_RE.findall(chart)) == 2
        assert chips == [f"{MOL_A} · 5", f"{MOL_B} · 1"]
        assert '<p class="residue-chart-segments muted">' in chart


# ============================================================
# 数据集页：柱色、图例、Y 轴、hover
# ============================================================

class TestDatasetChartEncoding:

    def test_bar_colour_uses_dominant_type(self, dataset):
        chart = _chart(_get_page(_select_client(dataset), _ds(dataset)))
        colours = {title.split(" · ")[0]: background
                   for _, background, title in _bars(chart)}
        assert colours["LEU50"] == "var(--c-hydrophobic)"
        assert colours["LIG60"] == "var(--c-hydrophobic)"
        assert colours["ASP100"] == "var(--c-hydrogen_bond)"
        assert colours["LYS20"] == "var(--c-hydrogen_bond)"

    def test_legend_lists_dominant_types_sorted(self, dataset):
        chart = _chart(_get_page(_select_client(dataset), _ds(dataset)))
        assert _legend(chart) == ["hydrogen_bond", "hydrophobic"]

    def test_yaxis_has_three_ticks(self, dataset):
        chart = _chart(_get_page(_select_client(dataset), _ds(dataset)))
        assert "<span>100%</span>" in chart
        assert "<span>50%</span>" in chart
        assert "<span>0</span>" in chart

    def test_hover_title_carries_label_molecule_occupancy_pairs(self, dataset):
        chart = _chart(_get_page(_select_client(dataset), _ds(dataset)))
        titles = [title for _, _, title in _bars(chart)]
        assert ("LEU50 · seg_0_pro · 100.0% · 1 pair · 1 type") in titles
        lys = next(title for title in titles if title.startswith("LYS10 · "))
        assert lys == "LYS10 · seg_0_pro · 50.0% · 2 pairs · 2 types"

    def test_singular_pair_title(self, dataset):
        chart = _chart(_get_page(_select_client(dataset), _ds(dataset)))
        titles = [title for _, _, title in _bars(chart)]
        leu = next(title for title in titles if title.startswith("LEU50 · "))
        assert " · 1 pair · " in leu
        assert "1 pairs" not in leu

    def test_singular_type_title(self, dataset):
        chart = _chart(_get_page(_select_client(dataset), _ds(dataset)))
        titles = [title for _, _, title in _bars(chart)]
        leu = next(title for title in titles if title.startswith("LEU50 · "))
        assert leu.endswith(" · 1 type")
        assert "1 types" not in leu
        lys = next(title for title in titles if title.startswith("LYS10 · "))
        assert lys.endswith(" · 2 types")

    def test_bar_height_is_occupancy_percent(self, dataset):
        chart = _chart(_get_page(_select_client(dataset), _ds(dataset)))
        heights = {title.split(" · ")[0]: height
                   for height, _, title in _bars(chart)}
        assert heights["LEU50"] == "100.0"
        assert heights["LYS20"] == "25.0"
        assert heights["ASP100"] == "50.0"


# ============================================================
# 数据集页：空数据与水
# ============================================================

class TestDatasetChartEmpty:

    def test_no_chart_when_rows_empty(self, dataset):
        html = _get_page(_select_client(dataset), _ds(dataset), search="NOPE")
        assert "no residues match" in html
        assert CHART_MARKER not in html

    def test_zero_occupancy_residues_render_no_chart(self, zero_dir):
        html = _get_page(_select_client(zero_dir), _ds(zero_dir))
        assert "LYS10" in html          # 表格仍列出零占用残基
        assert CHART_MARKER not in html  # 图表不渲染零占用柱

    def test_water_residue_absent_from_segments(self, dataset):
        html = _get_page(_select_client(dataset), _ds(dataset))
        chart = _chart(html)
        assert "SOL" not in chart
        assert [label for _, label in _segments(chart)] == [MOL_A, MOL_B]


# ============================================================
# 类型页：同一 partial
# ============================================================

class TestTypePageChart:

    def test_chart_inside_top_residues_panel(self, dataset):
        html = _select_client(dataset).get("/type/hydrophobic").text
        chart = _chart(html)
        assert chart, "type page chart missing"
        assert 'class="residue-chart"' in html
        assert html.index("residue-chart") < html.index("hotspots-table")
        panel = html.split("top residues · this type", 1)[1]
        assert "residue-chart" in panel.split("hotspots-table", 1)[0]

    def test_type_page_segments_match_backend(self, dataset):
        html = _select_client(dataset).get("/type/hydrophobic").text
        chart = _chart(html)
        segments = _segments(chart)
        assert [label for _, label in segments] == [MOL_A, MOL_B]
        assert len(_bars(chart)) == 2

    def test_type_page_legend_single_colour(self, dataset):
        chart = _chart(_select_client(dataset).get("/type/hydrophobic").text)
        assert _legend(chart) == ["hydrophobic"]
        assert {background for _, background, _ in _bars(chart)} == {
            "var(--c-hydrophobic)"}


# ============================================================
# 无 JS 依赖
# ============================================================

class TestNoJavascript:

    def test_dataset_page_adds_no_script(self, dataset):
        """图表页不加载图表专用脚本：只允许全局外壳脚本（load.js / browse.js）。"""
        html = _get_page(_select_client(dataset), _ds(dataset))
        scripts = set(re.findall(r'<script src="([^"]+)"', html))
        assert scripts <= GLOBAL_SCRIPTS, scripts
        assert "<script" not in _chart(html)

    def test_type_page_script_set_unchanged(self, dataset):
        """语义：仅全局脚本 + matrix/rows，不混入图表专用脚本。"""
        html = _select_client(dataset).get("/type/hydrophobic").text
        scripts = set(re.findall(r'<script src="([^"]+)"', html))
        assert scripts - GLOBAL_SCRIPTS == {"/static/matrix.js",
                                            "/static/rows.js"}
        assert "<script" not in _chart(html)


# ============================================================
# CSS
# ============================================================

class TestChartCss:

    def test_bar_width_not_capped(self):
        rules = _rules(STATIC_CSS.read_text(encoding="utf-8"),
                       ".residue-bars span")
        assert any("max-width: none" in rule for rule in rules)
        assert any("flex: 1 1 0" in rule for rule in rules)

    def test_segment_track_and_label_defined(self):
        css = STATIC_CSS.read_text(encoding="utf-8")
        assert _rules(css, ".residue-chart")
        assert _rules(css, ".residue-seg")
        assert _rules(css, ".residue-track")
        assert _rules(css, ".residue-seg-label")

    def test_plot_grid_reserves_yaxis_column(self):
        rules = _rules(STATIC_CSS.read_text(encoding="utf-8"),
                       ".residue-plot")
        assert any("grid-template-columns" in rule for rule in rules)
        assert any("96px" in rule for rule in _rules(
            STATIC_CSS.read_text(encoding="utf-8"), ".bars.residue-bars"))

    def test_segment_chip_rules_defined(self):
        css = STATIC_CSS.read_text(encoding="utf-8")
        assert any("flex-wrap: wrap" in rule
                   for rule in _rules(css, ".residue-chart-segments"))
        assert any("white-space: nowrap" in rule
                   for rule in _rules(css, ".residue-seg-chip"))


# ============================================================
# 真实数据
# ============================================================

class TestRealDataChart:

    def test_segment_and_bar_counts(self, real_page, real_segments):
        chart = _chart(real_page["residues"])
        assert len(SEG_RE.findall(chart)) == len(real_segments) == 6
        assert len(_bars(chart)) == sum(seg["n_residues"]
                                        for seg in real_segments) == 316
        assert [label for _, label in _segments(chart)] == [
            seg["molecule"] for seg in real_segments]

    def test_legend_matches_dominant_types(self, real_page, real_segments):
        expected = sorted({row["dominant_type"] for seg in real_segments
                           for row in seg["residues"]})
        assert _legend(_chart(real_page["residues"])) == expected
        assert expected == ["hydrogen_bond", "hydrophobic",
                            "metal_coordination", "water_bridge"]

    def test_no_water_molecule_in_segments(self, real_page):
        labels = [label for _, label in _segments(_chart(real_page["residues"]))]
        assert labels
        assert not any("SOL" in label for label in labels)

    def test_segment_chips_list_narrow_segments(self, real_page, real_segments):
        chips = CHIP_RE.findall(_chart(real_page["residues"]))
        assert chips == [f"{seg['molecule']} · {seg['n_residues']}"
                         for seg in real_segments]
        assert len(chips) == 6
        assert "seg_1_D927 · 1" in chips

    def test_type_page_water_bridge_single_colour(self, real_page,
                                                  real_water_segments):
        chart = _chart(real_page["water_bridge"])
        assert len(_segments(chart)) == len(real_water_segments) == 4
        assert len(_bars(chart)) == sum(seg["n_residues"]
                                        for seg in real_water_segments)
        assert _legend(chart) == ["water_bridge"]
        assert "SOL" not in chart

    def test_real_chart_before_table(self, real_page):
        for html in real_page.values():
            assert CHART_MARKER in html
            assert html.index("residue-chart") < html.index("hotspots-table")
