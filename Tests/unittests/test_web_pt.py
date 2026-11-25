# -*- coding: utf-8 -*-
"""Web 可视化 π-堆积 P/T 类型分色与 meta 悬停的单元测试（合成数据，不依赖轨迹）。"""

import json
import re
from pathlib import Path
from typing import List, Optional

import numpy as np
import pytest
from fastapi.testclient import TestClient

from DuIvyInteractions.core.datas import AtomData, Group, Interaction
from DuIvyInteractions.io.h5 import save_interactions
from DuIvyInteractions.visualizers import (
    fig_pair, matrix_view, pair_drilldown, pair_types)
from DuIvyInteractions.visualizers import figures
from DuIvyInteractions.visualizers.data import MAX_SERIES_POINTS, PT_METRIC, PT_TYPES
from DuIvyInteractions.visualizers.figures import PT_COLORS, TYPE_RANGE
from DuIvyInteractions.web.app import create_app

FRAMES = 6
SMALL_BUCKET_POINTS = 50  # 缩小的降采样上限（桶宽 3 帧，便于构造众数回归）
PT_EXISTENCE = [[1, 1, 1, 1, 1, 0], [1, 0, 1, 1, 0, 1]]
PT_VALUES = [["P", "P", "T", "T", "N", "N"], ["N", "N", "P", "P", "T", "T"]]
TEMPLATE = Path(figures.__file__).parent.parent / "web" / "templates" / "type.html"
MATRIX_JS = Path(figures.__file__).parent.parent / "web" / "static" / "matrix.js"
APP_CSS = Path(figures.__file__).parent.parent / "web" / "static" / "app.css"
NPG_PALETTE = frozenset({
    "#e64b35", "#4dbbd5", "#00a087", "#3c5488", "#f39b7f",
    "#8491b4", "#91d1c2", "#dc0000", "#7e6148"})  # ggsci scale_color_npg 9 色


def _atom(gidx: int, name: str) -> AtomData:
    """最小原子元信息。"""
    return AtomData(atom_global_idx=gidx, atom_idx_in_residue=0, atom_name=name,
                    atom_type="ca", atom_element="C", atom_charge=0.0,
                    atom_mass=12.0)


def _group(gid: int, gidx: int, res: str, rid: int) -> Group:
    """芳香环基团（两原子即可生成 pair label）。"""
    return Group(group_id=gid, group_type="aromatic_ring", molecule="LIG",
                 residue_name=res, residue_id=rid,
                 atoms=[_atom(gidx, "C1"), _atom(gidx + 1, "C2")])


def _groups(n_pairs: int) -> List[tuple]:
    """每对的芳香环基团对。"""
    return [(_group(2 * i, 10 * i, "PHE", i),
             _group(2 * i + 1, 10 * i + 2, "TYR", 100 + i))
            for i in range(n_pairs)]


def _pi_interaction(existence: List[List[int]],
                    values: Optional[List[List[str]]]) -> Interaction:
    """合成 π-堆积 Interaction（values 为 None 时不带类型指标）。"""
    mask = np.asarray(existence, dtype=bool)
    metrics = {"distance": np.where(mask, 3.5, np.nan),
               "angle": np.where(mask, 45.0, np.nan)}
    if values is not None:
        metrics[PT_METRIC] = np.asarray(values, dtype="U1")
    return Interaction(interaction_type="pi_stacking", groups=_groups(mask.shape[0]),
                       existence=mask, metrics=metrics,
                       times=np.arange(mask.shape[1], dtype=float) * 10.0)


def _long_interaction(frames: int, values: List[str],
                      existence: Optional[np.ndarray] = None) -> Interaction:
    """超长 π-堆积 Interaction（可指定逐帧存在性，触发桶降采样）。"""
    mask = np.ones((1, frames), dtype=bool) if existence is None else existence
    return Interaction(
        interaction_type="pi_stacking", groups=_groups(1),
        existence=mask,
        metrics={"distance": np.full((1, frames), 3.0),
                 PT_METRIC: np.asarray([values], dtype="U1")},
        times=np.arange(frames, dtype=float))


def _hb_interaction() -> Interaction:
    """合成无字符串指标的氢键 Interaction。"""
    groups = [(_group(2 * i, 10 * i, "LYS", i),
               _group(2 * i + 1, 10 * i + 5, "ASP", i)) for i in range(2)]
    return Interaction(interaction_type="hydrogen_bond", groups=groups,
                       existence=np.ones((2, FRAMES), dtype=bool),
                       metrics={"distance": np.full((2, FRAMES), 2.8)},
                       times=np.arange(FRAMES, dtype=float) * 10.0)


def _row(interaction: Interaction, idx: int) -> dict:
    """矩阵视图中指定 pair 的行。"""
    return next(row for row in matrix_view(interaction)["rows"]
                if row["idx"] == idx)


def _typed_fig(interaction: Interaction) -> dict:
    """第一对的下钻图表 JSON。"""
    return fig_pair(pair_drilldown(interaction, 0))


def _run_length(runs: List[List[int]]) -> int:
    """区段列表覆盖的帧数。"""
    return sum(end - start for start, end in runs)


@pytest.fixture(autouse=True)
def isolate_recent(monkeypatch, tmp_path):
    """隔离最近数据源缓存。"""
    monkeypatch.setattr("DuIvyInteractions.visualizers.data.RECENT_PATH",
                        tmp_path / "recent.json")


@pytest.fixture
def pi_interaction() -> Interaction:
    """带 P/T/N 交替类型的 π-堆积数据。"""
    return _pi_interaction(PT_EXISTENCE, PT_VALUES)


@pytest.fixture
def plain_interaction() -> Interaction:
    """无类型指标的 π-堆积数据。"""
    return _pi_interaction(PT_EXISTENCE, None)


@pytest.fixture
def h5_dir(tmp_path) -> Path:
    """含 π-堆积与氢键两个 h5 的临时数据源目录。"""
    directory = tmp_path / "out"
    directory.mkdir()
    save_interactions([_pi_interaction(PT_EXISTENCE, PT_VALUES)],
                      str(directory / "pi_stacking.h5"))
    save_interactions([_hb_interaction()], str(directory / "hydrogen_bond.h5"))
    return directory


@pytest.fixture
def client(h5_dir: Path) -> TestClient:
    """已载入临时数据源的 TestClient。"""
    test_client = TestClient(create_app())
    resp = test_client.post("/select", data={"path": str(h5_dir)},
                            follow_redirects=False)
    assert resp.status_code == 303
    return test_client


@pytest.fixture
def pt_only_client(tmp_path) -> TestClient:
    """只含 P/T 活跃类型的 π-堆积类型页 TestClient。"""
    directory = tmp_path / "pt_only"
    directory.mkdir()
    save_interactions(
        [_pi_interaction([[1] * FRAMES, [1] * FRAMES],
                         [["P", "P", "T", "T", "P", "P"],
                          ["P", "T", "T", "T", "P", "P"]])],
        str(directory / "pi_stacking.h5"))
    test_client = TestClient(create_app())
    resp = test_client.post("/select", data={"path": str(directory)},
                            follow_redirects=False)
    assert resp.status_code == 303
    return test_client


class TestStructurePtBar:
    """3D 结构页的 π-堆积类型向量条。"""

    def test_pair_types_returns_raw_frames(self, pi_interaction):
        assert pair_types(pi_interaction, 0) == PT_VALUES[0]

    def test_pair_types_empty_without_metric(self, plain_interaction):
        assert pair_types(plain_interaction, 0) == []

    def test_pair_types_out_of_range(self, pi_interaction):
        assert pair_types(pi_interaction, 99) == []

    def test_structure_page_colors_pt_frames(self, client):
        html = client.get("/structure/pi_stacking/0").text
        assert 'class="on pt-p"' in html
        assert 'class="on pt-t"' in html
        assert "pt-legend" in html
        match = re.search(
            r'<script id="structureConfig"[^>]*>(.*?)</script>', html, re.S)
        assert match and json.loads(match.group(1))["pt_types"] == PT_VALUES[0]

    def test_structure_page_without_pt_has_no_legend(self, client):
        html = client.get("/structure/hydrogen_bond/0").text
        assert "pt-legend" not in html and "pt-p" not in html

    def test_structure_bar_colors_match_figures(self):
        css = APP_CSS.read_text(encoding="utf-8")
        for name, color in PT_COLORS.items():
            rule = ".structure-existence span.on.pt-%s { background: %s; }" % (
                name.lower(), color)
            assert rule in css


class TestMatrixTyped:
    """矩阵行的 typed 字段。"""

    def test_rows_carry_typed_keys(self, pi_interaction):
        row = _row(pi_interaction, 0)
        assert set(row["typed"]) == set(PT_TYPES)
        assert row["typed"]["T"] == [[2, 4]]

    def test_typed_splits_at_type_change(self, pi_interaction):
        assert _row(pi_interaction, 0)["typed"] == {
            "P": [[0, 2]], "T": [[2, 4]], "N": [[4, 5]]}
        assert _row(pi_interaction, 1)["typed"] == {
            "P": [[2, 4]], "T": [[5, 6]], "N": [[0, 1]]}

    def test_absent_type_is_empty_list(self):
        inter = _pi_interaction([[1, 1, 1, 1]], [["P", "P", "N", "N"]])
        assert _row(inter, 0)["typed"]["T"] == []

    def test_segment_lengths_sum_to_runs(self, pi_interaction):
        for idx in (0, 1):
            row = _row(pi_interaction, idx)
            typed_frames = sum(_run_length(runs) for runs in row["typed"].values())
            assert typed_frames == _run_length(row["runs"])

    def test_without_metric_no_typed_key(self, plain_interaction):
        assert "typed" not in _row(plain_interaction, 0)


class TestDrilldownTypes:
    """下钻数据的 types 序列。"""

    def test_values_align_with_times(self, pi_interaction):
        data = pair_drilldown(pi_interaction, 0)
        assert data["types"]["name"] == PT_METRIC
        assert len(data["types"]["values"]) == len(data["times"]) == FRAMES

    def test_per_frame_values(self, pi_interaction):
        assert pair_drilldown(pi_interaction, 0)["types"]["values"] == \
            ["P", "P", "T", "T", "N", "N"]

    def test_without_metric_no_types_key(self, plain_interaction):
        assert "types" not in pair_drilldown(plain_interaction, 0)

    def test_bucket_mode_and_tie_priority(self):
        values = ["T"] * 3000
        values[1], values[2] = "P", "T"
        values[4], values[5] = "N", "T"
        data = pair_drilldown(_long_interaction(3000, values), 0)
        assert data["downsampled"] is True
        assert len(data["types"]["values"]) == len(data["times"]) == MAX_SERIES_POINTS
        assert data["types"]["values"][0] == "T"
        assert data["types"]["values"][1] == "P"
        assert data["types"]["values"][3] == "T"

    def test_downsample_mode_uses_active_frames_only(self, monkeypatch):
        """降采样众数只统计活跃帧；桶内全非活跃时为空串。"""
        monkeypatch.setattr(
            "DuIvyInteractions.visualizers.data.MAX_SERIES_POINTS",
            SMALL_BUCKET_POINTS)
        frames = SMALL_BUCKET_POINTS + 100  # 桶宽 3 帧
        values = ["N"] * frames
        values[0] = "P"
        existence = np.zeros((1, frames), dtype=bool)
        existence[0, 0] = True
        data = pair_drilldown(_long_interaction(frames, values, existence), 0)
        assert data["downsampled"] is True
        assert len(data["types"]["values"]) == SMALL_BUCKET_POINTS
        assert data["types"]["values"][0] == "P"
        assert data["types"]["values"][1] == ""


class TestFiguresTyped:
    """fig_pair 的类型带子图。"""

    def test_type_traces_added(self, pi_interaction):
        assert [trace["name"] for trace in _typed_fig(pi_interaction)["data"]] == \
            ["existence", "type P", "type T", "type N", "distance", "angle"]

    def test_type_markers_use_pt_colors(self, pi_interaction):
        for trace in _typed_fig(pi_interaction)["data"]:
            if not trace["name"].startswith("type "):
                continue
            label = trace["name"].split(" ")[1]
            assert trace["marker"]["color"] == PT_COLORS[label]
            assert trace["marker"]["symbol"] == "square"
            assert trace["mode"] == "markers" and trace["showlegend"] is True

    def test_type_band_below_existence(self, pi_interaction):
        layout = _typed_fig(pi_interaction)["layout"]
        assert layout["yaxis2"]["range"] == TYPE_RANGE
        assert layout["yaxis2"]["domain"][1] <= layout["yaxis"]["domain"][0]
        assert layout["showlegend"] is True

    def test_metric_axes_shift_after_type_band(self, pi_interaction):
        fig = _typed_fig(pi_interaction)
        assert [trace["yaxis"] for trace in fig["data"][-2:]] == ["y3", "y4"]

    def test_type_trace_skips_inactive_frames(self, pi_interaction):
        n_trace = next(trace for trace in _typed_fig(pi_interaction)["data"]
                       if trace["name"] == "type N")
        assert n_trace["x"] == [40.0]

    def test_without_types_traces_unchanged(self, plain_interaction):
        fig = _typed_fig(plain_interaction)
        assert [trace["name"] for trace in fig["data"]] == \
            ["existence", "distance", "angle"]
        assert fig["layout"]["showlegend"] is False

    def test_without_types_axis_layout_unchanged(self, plain_interaction):
        layout = _typed_fig(plain_interaction)["layout"]
        assert "yaxis2" in layout and "yaxis3" in layout
        assert "yaxis4" not in layout and "xaxis4" not in layout


class TestTypePage:
    """类型页的 meta 悬停、色点、图例与矩阵配置。"""

    def test_meta_values_have_full_titles(self, client, h5_dir):
        html = client.get("/type/pi_stacking").text
        assert 'title="pi_stacking"' in html
        assert 'title="2"' in html and 'title="6"' in html
        assert 'title="0.0 → 50.0"' in html
        assert 'title="angle, distance, pistacking_type"' in html
        assert f'title="{h5_dir / "pi_stacking.h5"}"' in html

    def test_string_chips_have_pt_dots(self, client):
        html = client.get("/type/pi_stacking").text
        assert html.count('class="pt-dot"') >= 3
        assert "background: #3c5488" in html

    def test_matrix_legend(self, client):
        html = client.get("/type/pi_stacking").text
        assert 'class="pt-legend"' in html
        assert "pistacking type" in html
        assert ">P</span>" in html

    def test_matrix_legend_lists_present_types_only(self, pt_only_client):
        html = pt_only_client.get("/type/pi_stacking").text
        legend = html.split('class="pt-legend"', 1)[1].split("</div>", 1)[0]
        assert legend.count("pt-legend-item") == 2
        assert ">P</span>" in legend and ">T</span>" in legend
        assert ">N</span>" not in legend

    def test_matrix_config_pt_colors(self, client):
        html = client.get("/type/pi_stacking").text
        assert '"pt_colors"' in html
        assert '"P": "#3c5488"' in html
        assert '"T": "#e64b35"' in html
        assert '"N": "#8491b4"' in html

    def test_non_pt_type_has_no_pt_fields(self, client):
        html = client.get("/type/hydrogen_bond").text
        assert "pt-dot" not in html and "pt_colors" not in html
        assert "pt-legend" not in html
        assert 'title="hydrogen_bond"' in html


class TestColorContract:
    """P/T/N 颜色三元组的一致性约定。"""

    def test_pt_colors_triplet(self):
        assert PT_COLORS == {"P": "#3c5488", "T": "#e64b35", "N": "#8491b4"}

    def test_palette_is_npg(self):
        colors = list(figures.PT_COLORS.values()) + list(figures.CHART_COLORS)
        assert {color.lower() for color in colors} <= NPG_PALETTE
        chart = [color.lower() for color in figures.CHART_COLORS]
        assert len(set(chart)) == len(chart)

    def test_template_colors_match_figures(self):
        source = TEMPLATE.read_text(encoding="utf-8")
        for name, color in PT_COLORS.items():
            assert f'"{name}": "{color}"' in source
        pt_line = next(line for line in source.splitlines()
                       if line.lstrip().startswith("{% set pt_colors"))
        assert "#2563eb" not in source and "#d97706" not in source
        assert "#94a3b8" not in pt_line

    def test_matrix_js_colors_match_figures(self):
        source = MATRIX_JS.read_text(encoding="utf-8")
        for name, color in PT_COLORS.items():
            assert f'{name}: "{color}"' in source
        assert 'ACTIVE = "60,84,136"' in source
        assert 'COLOR_BAR_FG = "#8491b4"' in source

    def test_no_tailwind_colors_in_figures(self):
        source = Path(figures.__file__).read_text(encoding="utf-8")
        assert "#4f46e5" not in source
        assert "rgba(79,70,229" not in source

    def test_existence_trace_uses_npg(self, pi_interaction):
        trace = _typed_fig(pi_interaction)["data"][0]
        assert trace["name"] == "existence"
        assert trace["line"]["color"] == "#3c5488"
        assert trace["fillcolor"] == "rgba(60,84,136,0.16)"

    def test_matrix_tip_label_is_npg(self):
        css = (MATRIX_JS.parent / "app.css").read_text(encoding="utf-8")
        block = css.split(".matrix-tip .t-label", 1)[1].split("}", 1)[0]
        assert "#91d1c2" in block
