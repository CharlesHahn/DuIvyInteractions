# -*- coding: utf-8 -*-
"""残基热点 Web 页测试：/residues、/download/residues.csv、主页与类型页入口。

用 tmp_path 构造最小 h5；残基编号/名字全部来自 h5 的 residue_id/residue_name。
"""

import csv
import io
import re
from pathlib import Path
from urllib.parse import quote

import h5py
import numpy as np
import pytest
from fastapi.testclient import TestClient

from DuIvyInteractions.core.datas import AtomData, Group, Interaction
from DuIvyInteractions.io.h5 import save_interactions
from DuIvyInteractions.web import app as web_app

FRAMES = 4
MANY_PAIRS = 210            # 大数据集：用于验证全部行渲染（无截断）
MANY_RESIDUES = 2 * MANY_PAIRS
STATIC_CSS = (Path(__file__).resolve().parents[2]
              / "DuIvyInteractions" / "web" / "static" / "app.css")
MERGED_HB = "hydrogen_bond@hydrogen_bond.h5"   # 两批次同名类型的显示名
NAN = float("nan")
HOTSPOT_TABLE = '<table class="hotspots-table">'
PANEL_TABLE = '<table class="table hotspots-table">'
SUMMARY_RE = re.compile(r'<p class="hotspots-summary">(.*?)</p>', re.S)
TR_RE = re.compile(r"<tr>(.*?)</tr>", re.S)
TD_RE = re.compile(r"<td[^>]*>(.*?)</td>", re.S)
TAG_RE = re.compile(r"<[^>]+>")

# 六类相互作用：LYS10 参与其中五类（types 列折叠为 3 + "+2"）
HB_SPECS = [
    [("H_donor", "LYS", 10), ("H_acceptor", "ASP", 100)],
    [("H_donor", "LYS", 10), ("H_acceptor", "GLU", 200)],
    [("H_donor", "ARG", 40), ("H_acceptor", "ASP", 140)],
]
HB_EXIST = [[True, True, False, False],
            [False, True, True, False],
            [True, False, False, False]]
HB_DIST = [[2.0, 2.5, NAN, NAN],
           [NAN, 3.0, 3.4, NAN],
           [2.8, NAN, NAN, NAN]]
SB_SPECS = [
    [("charged_positive", "LYS", 10), ("charged_negative", "ASP", 100)],
    [("charged_positive", "LYS", 10), ("charged_negative", "GLU", 200)],
    [("charged_positive", "LYS", 30), ("charged_negative", "GLU", 300)],
]
SB_EXIST = [[True, False, True, False],
            [False, True, False, False],
            [False, False, False, True]]
SB_DIST = [[3.0, NAN, 3.2, NAN],
           [NAN, 4.0, NAN, NAN],
           [NAN, NAN, NAN, 5.0]]
WATER_SPECS = [[("H_donor", "LYS", 10), ("water", "SOL", 500),
                ("H_acceptor", "ASP", 100)]]
WATER_EXIST = [[True, False, False, False]]
WATER_DIST = [[3.0, NAN, NAN, NAN]]
HYDRO_SPECS = [[("hydrophobic", "LYS", 10), ("hydrophobic", "LEU", 50)]]
HYDRO_EXIST = [[True, True, True, True]]
HYDRO_DIST = [[3.5, 3.6, 3.4, 3.5]]
METAL_SPECS = [[("metal_binding", "LYS", 10), ("metal", "MG", 900)]]
METAL_EXIST = [[True, True, False, False]]
METAL_DIST = [[2.1, 2.2, NAN, NAN]]
PICATION_SPECS = [[("aromatic_ring", "PHE", 60), ("charged_positive", "LYS", 10)]]
PICATION_EXIST = [[False, True, True, False]]
PICATION_DIST = [[NAN, 4.0, 4.1, NAN]]


# ============================================================
# fixture 构造
# ============================================================

def _atom(gidx: int) -> AtomData:
    """占位原子（热点聚合不读原子）。"""
    return AtomData(atom_global_idx=gidx, atom_idx_in_residue=0, atom_name="C",
                    atom_type="", atom_element="C", atom_charge=0.0,
                    atom_mass=0.0)


def _group(gid: int, gtype: str, res: str, rid: int) -> Group:
    """按 (group_type, residue_name, residue_id) 造基团（氢键供体补 H 原子）。"""
    atoms = [_atom(gid)]
    if gtype == "H_donor":
        atoms.append(_atom(gid + 1))
    return Group(group_id=gid, group_type=gtype, molecule="test",
                 residue_name=res, residue_id=rid, atoms=atoms)


def _interaction(itype: str, specs, existence, distances) -> Interaction:
    """specs 每对为 [(group_type, residue_name, residue_id), ...]。"""
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


def _many_specs():
    """MANY_PAIRS 对唯一的 (LYS i) – (ASP 1000+i)。"""
    specs, existence = [], []
    for i in range(MANY_PAIRS):
        specs.append([("H_donor", "LYS", i), ("H_acceptor", "ASP", 1000 + i)])
        row = [False] * FRAMES
        row[0] = True
        existence.append(row)
    distances = [[3.0 if frame == 0 else NAN for frame in range(FRAMES)]
                 for _ in range(MANY_PAIRS)]
    return specs, existence, distances


@pytest.fixture(autouse=True)
def isolate_recent(monkeypatch, tmp_path):
    """隔离最近数据源文件，避免写入用户缓存。"""
    monkeypatch.setattr("DuIvyInteractions.visualizers.data.RECENT_PATH",
                        tmp_path / "recent.json")


@pytest.fixture
def dataset(tmp_path):
    """六类型小数据集（LYS10 参与五类；SOL500 只出现在水桥）。"""
    directory = tmp_path / "small"
    directory.mkdir()
    _save(directory, "hydrogen_bond", HB_SPECS, HB_EXIST, HB_DIST)
    _save(directory, "salt_bridge", SB_SPECS, SB_EXIST, SB_DIST)
    _save(directory, "water_bridge", WATER_SPECS, WATER_EXIST, WATER_DIST)
    _save(directory, "hydrophobic", HYDRO_SPECS, HYDRO_EXIST, HYDRO_DIST)
    _save(directory, "metal_coordination", METAL_SPECS, METAL_EXIST, METAL_DIST)
    _save(directory, "pi_cation", PICATION_SPECS, PICATION_EXIST, PICATION_DIST)
    return directory


@pytest.fixture
def many_dir(tmp_path):
    """单类型大数据集：2 * MANY_PAIRS 个残基。"""
    directory = tmp_path / "many"
    directory.mkdir()
    specs, existence, distances = _many_specs()
    _save(directory, "hydrogen_bond", specs, existence, distances)
    return directory


def _select_client(*paths) -> TestClient:
    """创建应用并逐个载入数据源。"""
    client = TestClient(web_app.create_app())
    for path in paths:
        response = client.post("/select", data={"path": str(path)},
                               follow_redirects=False)
        assert response.status_code == 303
        assert response.headers["location"] == "/"
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


def _rows(html: str, marker: str = HOTSPOT_TABLE):
    """marker 之后首个 tbody 的数据行 → 单元格文本列表。"""
    assert marker in html, f"table missing: {marker}"
    body = html.split(marker, 1)[1].split("<tbody>", 1)[1]
    body = body.split("</tbody>", 1)[0]
    return [[TAG_RE.sub("", cell).strip() for cell in TD_RE.findall(row)]
            for row in TR_RE.findall(body)]


def _labels(html: str, marker: str = HOTSPOT_TABLE):
    """表格首列（residue 标签）顺序。"""
    return [row[0] for row in _rows(html, marker)]


def _summary(html: str) -> str:
    """页面摘要行文本。"""
    match = SUMMARY_RE.search(html)
    assert match, "summary line missing"
    return match.group(1).strip()


def _csv(client: TestClient, ds: str, **params) -> str:
    """下载热点 CSV 并断言响应头。"""
    query = "&".join(f"{key}={quote(str(value), safe='')}"
                     for key, value in params.items())
    response = client.get(f"/download/residues.csv?ds="
                          f"{quote(ds, safe='/')}&{query}")
    assert response.status_code == 200, response.text
    assert response.headers["content-type"].startswith("text/csv")
    return response.text


def _rules(css: str, selector: str):
    """CSS 中某个选择器的声明块文本列表（.a 与 .a > b 不互相匹配）。"""
    pattern = re.compile(re.escape(selector) + r"\s*\{([^}]*)\}")
    return pattern.findall(css)


def _h5_labels(directory: Path):
    """直读 h5 groups/residue_id + residue_name → {NAME+id}。"""
    labels = set()
    for path in sorted(Path(directory).glob("*.h5")):
        with h5py.File(path, "r") as handle:
            group = handle["interaction_0/groups"]
            ids = group["residue_id"][:]
            names = [name.decode() if isinstance(name, bytes) else str(name)
                     for name in group["residue_name"][:]]
            labels.update(f"{name}{int(rid)}" for name, rid in zip(names, ids))
    return labels


# ============================================================
# 数据集热点页
# ============================================================

class TestDatasetPage:

    def test_status_and_residue_labels(self, dataset):
        rows = _rows(_get_page(_select_client(dataset), _ds(dataset)))
        labels = [row[0] for row in rows]
        assert rows
        assert "LYS10" in labels and "ARG40" in labels

    def test_summary_counts_and_water_excluded(self, dataset):
        html = _get_page(_select_client(dataset), _ds(dataset))
        summary = _summary(html)
        assert "10 residues" in summary
        assert "all interactions" in summary
        assert "water: excluded" in summary

    def test_occupancy_and_counts(self, dataset):
        rows = _rows(_get_page(_select_client(dataset), _ds(dataset)))
        lys = next(row for row in rows if row[0] == "LYS10")
        assert lys[1] == "100.0%"
        assert lys[2] == "4"          # frames_active（跨类型并集）
        assert lys[3] == "8"          # n_pairs（水桥对仍为蛋白残基计一次）
        assert lys[4] == "5"          # n_partners
        assert len(lys) == 6          # residue/occupancy/frames/pairs/partners/types

    def test_types_chips_rendered_in_full(self, dataset):
        html = _get_page(_select_client(dataset), _ds(dataset))
        types_cell = next(row for row in _rows(html) if row[0] == "LYS10")[5]
        for name in ("hydrogen_bond", "salt_bridge", "water_bridge",
                     "hydrophobic", "metal_coordination", "pi_cation"):
            assert name in types_cell
        assert "hotspots-types-more" not in html
        assert re.search(r">\+\d+<", types_cell) is None

    def test_no_distance_column_or_hint(self, dataset):
        html = _get_page(_select_client(dataset), _ds(dataset))
        assert "avg dist" not in html
        assert "D–OW" not in html


# ============================================================
# 类型过滤
# ============================================================

class TestTypeFilter:

    def test_only_selected_type(self, dataset):
        html = _get_page(_select_client(dataset), _ds(dataset),
                         type="salt_bridge")
        assert set(_labels(html)) == {"LYS10", "ASP100", "GLU200",
                                      "LYS30", "GLU300"}
        assert "ARG40" not in html
        assert "5 residues" in _summary(html)
        assert "salt_bridge" in _summary(html)

    def test_unknown_type_404(self, dataset):
        response = _select_client(dataset).get(
            _url(_ds(dataset), type="pi_stacking"))
        assert response.status_code == 404
        assert "interaction type not loaded" in response.text

    def test_no_distance_column_even_with_type_filter(self, dataset):
        html = _get_page(_select_client(dataset), _ds(dataset),
                         type="salt_bridge")
        rows = _rows(html)
        assert rows and all(len(row) == 6 for row in rows)
        assert "avg dist" not in html

    def test_water_bridge_has_no_distance_column(self, dataset):
        html = _get_page(_select_client(dataset), _ds(dataset),
                         type="water_bridge")
        rows = _rows(html)
        assert {row[0] for row in rows} == {"LYS10", "ASP100"}
        assert all(len(row) == 6 for row in rows)
        assert "avg dist" not in html


# ============================================================
# 数据集选择 / 404
# ============================================================

class TestDatasetSelection:

    def test_single_dataset_autoselected(self, dataset):
        html = _get_page(_select_client(dataset), "")
        assert "LYS10" in _labels(html)
        assert "pick a dataset" not in html

    def test_multiple_datasets_render_links(self, dataset, many_dir):
        client = _select_client(dataset, many_dir)
        response = client.get("/residues")
        assert response.status_code == 200
        html = response.text
        assert "pick a dataset" in html
        for name in ("small", "many"):
            assert f">{name}</a>" in html
        for path in (_ds(dataset), _ds(many_dir)):
            assert f'href="/residues?ds={quote(path, safe="/")}"' in html
        assert HOTSPOT_TABLE not in html
        assert "residue ids/names come from the h5" not in html

    def test_unknown_dataset_404(self, dataset):
        response = _select_client(dataset).get("/residues?ds=/nope/not-loaded")
        assert response.status_code == 404
        assert "dataset not loaded" in response.text

    def test_type_scoped_to_selected_dataset(self, dataset, many_dir):
        client = _select_client(dataset, many_dir)
        response = client.get(_url(_ds(many_dir), type="salt_bridge"))
        assert response.status_code == 404
        assert "interaction type not loaded" in response.text

    def test_merged_display_name_still_filters(self, dataset, many_dir):
        client = _select_client(dataset, many_dir)
        rows = _rows(_get_page(client, _ds(many_dir), type=MERGED_HB))
        assert len(rows) == MANY_RESIDUES
        assert rows[0][0] == "LYS0"

    def test_no_dataset_choice_page(self):
        html = TestClient(web_app.create_app()).get("/residues").text
        assert "pick a dataset" in html
        assert "no dataset loaded" in html


# ============================================================
# 水残基
# ============================================================

class TestWater:

    def test_default_excludes_water(self, dataset):
        html = _get_page(_select_client(dataset), _ds(dataset))
        assert "SOL500" not in html

    def test_include_water_adds_residue(self, dataset):
        html = _get_page(_select_client(dataset), _ds(dataset),
                         include_water="1")
        rows = _rows(html)
        sol = next(row for row in rows if row[0] == "SOL500")
        assert sol[3] == "1"
        assert sol[4] == "2"
        assert "water: included" in _summary(html)


# ============================================================
# 搜索与排序
# ============================================================

class TestSearchSort:

    def test_search_filters_rows(self, dataset):
        html = _get_page(_select_client(dataset), _ds(dataset),
                         type="salt_bridge", search="LYS")
        assert _labels(html) == ["LYS10", "LYS30"]

    def test_default_sort_is_occupancy_desc(self, dataset):
        labels = _labels(_get_page(_select_client(dataset), _ds(dataset),
                                   type="salt_bridge"))
        assert labels[0] == "LYS10"

    def test_occupancy_ascending(self, dataset):
        labels = _labels(_get_page(_select_client(dataset), _ds(dataset),
                                   type="salt_bridge", sort_by="occupancy",
                                   sort_dir="asc"))
        assert labels == ["LYS30", "GLU200", "GLU300", "ASP100", "LYS10"]

    def test_pairs_descending(self, dataset):
        labels = _labels(_get_page(_select_client(dataset), _ds(dataset),
                                   type="salt_bridge", sort_by="pairs",
                                   sort_dir="desc"))
        assert labels == ["LYS10", "ASP100", "LYS30", "GLU200", "GLU300"]

    def test_residue_ascending(self, dataset):
        labels = _labels(_get_page(_select_client(dataset), _ds(dataset),
                                   type="salt_bridge", sort_by="residue",
                                   sort_dir="asc"))
        assert labels == ["ASP100", "GLU200", "GLU300", "LYS10", "LYS30"]


# ============================================================
# 全量渲染（无行数上限）
# ============================================================

class TestFullRender:

    def test_all_rows_rendered(self, many_dir):
        html = _get_page(_select_client(many_dir), _ds(many_dir))
        rows = _rows(html)
        assert len(rows) == MANY_RESIDUES
        assert f"{MANY_RESIDUES} residues" in _summary(html)
        assert "showing top" not in html
        assert "refine the filter" not in html
        assert "residue ids/names come from the h5" not in html

    def test_no_render_cap_constants(self):
        assert not hasattr(web_app, "HOTSPOT_ROW_LIMIT")
        assert not hasattr(web_app, "HOTSPOT_TYPE_TOP")
        assert not hasattr(web_app, "HOTSPOT_TYPES_SHOWN")

    def test_dataset_table_uses_scroll_wrapper(self, many_dir):
        html = _get_page(_select_client(many_dir), _ds(many_dir))
        assert '<div class="table-wrap hotspots-scroll">' in html
        assert '<div class="table-wrap short">' not in html

    def test_scroll_wrapper_css(self):
        css = STATIC_CSS.read_text(encoding="utf-8")
        rules = _rules(css, ".hotspots-scroll")
        assert any("max-height" in rule for rule in rules)
        assert any("overflow-y: auto" in rule for rule in rules)
        for dead in (".hotspots-limit", ".hotspots-types-more",
                     ".hotspots-distance"):
            assert _rules(css, dead) == []


# ============================================================
# CSV 下载
# ============================================================

class TestCsvDownload:

    def test_header_and_row_count_match_total(self, dataset):
        client = _select_client(dataset)
        text = _csv(client, _ds(dataset))
        rows = list(csv.reader(io.StringIO(text)))
        assert rows[0] == ["residue", "occupancy", "frames_active", "n_pairs",
                           "n_partners", "n_types", "types"]
        assert len(rows) - 1 == 10

    def test_type_and_water_filters(self, dataset):
        client = _select_client(dataset)
        rows = list(csv.reader(io.StringIO(
            _csv(client, _ds(dataset), type="salt_bridge"))))
        assert len(rows) - 1 == 5
        rows = list(csv.reader(io.StringIO(
            _csv(client, _ds(dataset), type="water_bridge",
                 include_water="1"))))
        assert {row[0] for row in rows[1:]} == {"LYS10", "ASP100", "SOL500"}

    def test_attachment_filename(self, dataset):
        response = _select_client(dataset).get(
            f"/download/residues.csv?ds={quote(_ds(dataset), safe='/')}")
        assert 'filename="hotspots-small.csv"' in \
            response.headers["content-disposition"]

    def test_unknown_dataset_404(self, dataset):
        response = _select_client(dataset).get(
            "/download/residues.csv?ds=/nope/not-loaded")
        assert response.status_code == 404
        assert "dataset not loaded" in response.text


# ============================================================
# 入口：主页 + 类型页
# ============================================================

class TestEntries:

    def test_home_dataset_button_per_dataset(self, dataset, many_dir):
        html = _select_client(dataset, many_dir).get("/").text
        assert html.count('href="/residues?ds=') == 2
        for path in (_ds(dataset), _ds(many_dir)):
            assert f'href="/residues?ds={quote(path, safe="/")}"' in html
        assert ">hot spots</a>" in html

    def test_type_page_has_plain_panel(self, many_dir):
        html = _select_client(many_dir).get("/type/hydrogen_bond").text
        assert "top residues · this type" in html
        assert "<details" not in html
        assert "all residues for this dataset" not in html
        assert "residue ids/names come from the h5" not in html
        assert "<th>avg dist</th>" in html

    def test_type_page_panel_renders_all_residues(self, many_dir):
        html = _select_client(many_dir).get("/type/hydrogen_bond").text
        assert len(_rows(html, PANEL_TABLE)) == MANY_RESIDUES
        assert '<div class="table-wrap short">' in html

    def test_type_page_has_no_dataset_link(self, dataset):
        html = _select_client(dataset).get("/type/salt_bridge").text
        assert "all residues for this dataset" not in html
        link = (f'/residues?ds={quote(_ds(dataset), safe="/")}'
                f'&amp;type=salt_bridge')
        assert link not in html

    def test_type_page_panel_only_current_type(self, dataset):
        html = _select_client(dataset).get("/type/salt_bridge").text
        panel = html.split("top residues · this type", 1)[1].split(
            "active pairs per frame", 1)[0]
        labels = _labels(panel, PANEL_TABLE)
        assert set(labels) == {"LYS10", "ASP100", "GLU200", "LYS30", "GLU300"}
        assert "ARG40" not in panel

    def test_type_page_merged_name_panel(self, dataset, many_dir):
        client = _select_client(dataset, many_dir)
        html = client.get(f"/type/{quote(MERGED_HB, safe='')}").text
        assert "top residues · this type" in html
        assert "<details" not in html
        assert len(_rows(html, PANEL_TABLE)) == MANY_RESIDUES


# ============================================================
# 残基身份来源（h5 直读对照）
# ============================================================

class TestIdentity:

    def test_page_labels_equal_h5_residue_combo(self, dataset):
        html = _get_page(_select_client(dataset), _ds(dataset))
        expected = {label for label in _h5_labels(dataset)
                    if not label.startswith("SOL")}
        assert set(_labels(html)) == expected
        assert len(expected) == 10

    def test_no_pdb_wording(self, dataset):
        client = _select_client(dataset)
        assert "PDB" not in _get_page(client, _ds(dataset))
        assert "PDB" not in client.get("/type/salt_bridge").text
