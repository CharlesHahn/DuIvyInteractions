# -*- coding: utf-8 -*-
"""Web 可视化单元测试：纯函数层（visualizers/）≥ 20 用例（spec SC-006）+ 路由冒烟。

不依赖 MDAnalysis 轨迹：用 tmp 目录构造 h5 fixture。
"""

from pathlib import Path

import numpy as np
import pytest

from DuIvyInteractions.core.datas import AtomData, Group, Interaction
from DuIvyInteractions.io.h5 import save_interactions
from DuIvyInteractions.visualizers import (
    add_recent, aggregate_table, cell_values, export_table, fig_pair,
    frame_activity, frame_view, get_exporter, interaction_meta, load_recent,
    matrix_view, metric_stats, overview, pair_drilldown, pair_frames,
    pair_summary, parse_datasource, runs_from_mask, string_stats,
)

FRAMES = 5
H5_DATA_DIR = Path(__file__).parent.parent / "interaction_h5data"


@pytest.fixture(autouse=True)
def isolate_recent(monkeypatch, tmp_path):
    monkeypatch.setattr("DuIvyInteractions.visualizers.data.RECENT_PATH",
                        tmp_path / "recent.json")


def _atom(gidx: int, name: str, element: str) -> AtomData:
    return AtomData(atom_global_idx=gidx, atom_idx_in_residue=0, atom_name=name,
                    atom_type="", atom_element=element, atom_charge=0.0,
                    atom_mass=0.0)


def _group(gid: int, gtype: str, gidx: int, name: str,
           res: str, rid: int, with_h: bool = False) -> Group:
    atoms = [_atom(gidx, name, "C")]
    if with_h:
        atoms.append(_atom(gidx, f"H{name}", "H"))
    return Group(group_id=gid, group_type=gtype, molecule="test",
                 residue_name=res, residue_id=rid, atoms=atoms)


def _interaction(itype: str, n_pairs: int, flags: bool) -> Interaction:
    existence = np.zeros((n_pairs, FRAMES), dtype=bool)
    existence[:, 0] = True
    if n_pairs > 1:
        existence[1, 1] = True
    rng = np.random.RandomState(42)
    metrics = {"distance": rng.uniform(2.0, 4.0, (n_pairs, FRAMES))}
    for i in range(n_pairs):
        metrics["distance"][i][~existence[i]] = np.nan
    groups = [(_group(2 * i, "H_donor", 10 * i, "C", "LYS", i, with_h=True),
               _group(2 * i + 1, "H_acceptor", 10 * i + 1, "O", "ASP", 100 + i))
              for i in range(n_pairs)]
    inter = Interaction(
        interaction_type=itype, groups=groups, existence=existence,
        metrics=metrics, times=np.arange(FRAMES, dtype=float) * 10.0)
    if flags:
        farr = np.zeros((n_pairs, FRAMES), dtype=bool)
        farr[0, 0] = True
        inter.flags["redundant_pi_stacking"] = farr
    return inter


@pytest.fixture
def h5_dir(tmp_path):
    """输出目录：2 个 h5（hydrogen_bond 无 flags / hydrophobic 有 flags）。"""
    d = tmp_path / "out"
    d.mkdir()
    save_interactions([_interaction("hydrogen_bond", 3, False)],
                      str(d / "hydrogen_bond.h5"))
    save_interactions([_interaction("hydrophobic", 2, True)],
                      str(d / "hydrophobic.h5"))
    (d / "corrupt.h5").write_text("not an h5")
    return d


@pytest.fixture
def hb(h5_dir):
    from DuIvyInteractions.io.h5 import load_interactions
    return load_interactions(str(h5_dir / "hydrogen_bond.h5"))[0]


@pytest.fixture
def hyd(h5_dir):
    from DuIvyInteractions.io.h5 import load_interactions
    return load_interactions(str(h5_dir / "hydrophobic.h5"))[0]


@pytest.fixture
def h5_data_dir():
    if not (H5_DATA_DIR / "pi_stacking.h5").exists():
        pytest.skip("real h5 test data not available")
    return H5_DATA_DIR


# ============================================================
# 数据源解析（T003）
# ============================================================

class TestParseDatasource:
    def test_directory_scans_all_valid_h5(self, h5_dir):
        items = parse_datasource(str(h5_dir))
        assert [it.interaction_type for _, _, it in items] == \
            ["hydrogen_bond", "hydrophobic"]

    def test_directory_skips_corrupt_with_warning(self, h5_dir, capsys):
        parse_datasource(str(h5_dir))
        assert "skip unreadable h5 corrupt.h5" in capsys.readouterr().out

    def test_single_file(self, h5_dir):
        items = parse_datasource(str(h5_dir / "hydrophobic.h5"))
        assert len(items) == 1 and items[0][0] == "hydrophobic"

    def test_nonexistent_path_raises(self, tmp_path):
        with pytest.raises(ValueError, match="path not found"):
            parse_datasource(str(tmp_path / "nope"))

    def test_structurally_invalid_hdf5_skipped(self, tmp_path):
        import h5py
        bad = tmp_path / "bad"
        bad.mkdir()
        with h5py.File(bad / "broken.h5", "w") as handle:
            handle.create_dataset("foo", data=[1])
        assert parse_datasource(str(bad)) == []

    def test_hdf5_missing_interaction_count_skipped(self, tmp_path):
        import h5py
        bad = tmp_path / "bad"
        bad.mkdir()
        with h5py.File(bad / "broken.h5", "w") as handle:
            handle.attrs["format_version"] = 1.0
        assert parse_datasource(str(bad)) == []

    def test_items_are_resolved_paths(self, h5_dir):
        items = parse_datasource(str(h5_dir))
        assert all(Path(p).is_absolute() for _, p, _ in items)


# ============================================================
# 总览（T004）
# ============================================================

class TestOverview:

    def test_card_counts(self, h5_dir):
        cards = overview(parse_datasource(str(h5_dir)))
        by_type = {c["type"]: c for c in cards}
        assert by_type["hydrogen_bond"]["n_pairs"] == 3
        assert by_type["hydrogen_bond"]["n_frames"] == FRAMES

    def test_card_flags_presence(self, h5_dir):
        cards = overview(parse_datasource(str(h5_dir)))
        by_type = {c["type"]: c for c in cards}
        assert by_type["hydrogen_bond"]["flags"] == []
        assert by_type["hydrophobic"]["flags"] == ["redundant_pi_stacking"]

    def test_card_occupancy_distribution(self, h5_dir):
        cards = overview(parse_datasource(str(h5_dir)))
        by_type = {c["type"]: c for c in cards}
        hist = by_type["hydrophobic"]["occupancy"]
        assert hist["count"] == 2 and sum(hist["hist"]) == 2


# ============================================================
# 类型表格（T005）
# ============================================================

class TestExportTable:

    def test_rows_and_columns(self, hb):
        table = export_table(hb)
        assert table["total"] == 3
        assert table["columns"] == ["pair_label", "occupancy",
                                    "avg_distance", "std_distance"]
        assert "avg_redundant_pi_stacking" not in table["columns"]

    def test_columns_include_flags_when_present(self, hyd):
        table = export_table(hyd)
        assert "avg_redundant_pi_stacking" in table["columns"]

    def test_flag_average_formatted(self, hyd):
        table = export_table(hyd)
        expected = get_exporter(hyd).get_pair_label(hyd, 0)
        row0 = next(r for r in table["rows"] if r["pair_label"] == expected)
        assert row0["avg_redundant_pi_stacking"] == "1.0000"

    def test_all_rows_returned(self, hb):
        assert len(export_table(hb)["rows"]) == 3

    def test_zero_occupancy_rows_skipped(self, hb):
        hb.existence[:, :] = False
        assert export_table(hb)["rows"] == []

    def test_search_filters_by_label(self, hb):
        table = export_table(hb, search="LYS0:")
        assert table["total"] == 1
        assert table["rows"][0]["pair_label"] == get_exporter(hb).get_pair_label(hb, 0)

    def test_search_case_insensitive(self, hb):
        assert export_table(hb, search="lys0:")["total"] == 1

    def test_sort_occupancy_desc(self, hb):
        rows = export_table(hb, sort_by="occupancy", sort_dir="desc")["rows"]
        assert rows[0]["occupancy"] == "0.4000"

    def test_sort_occupancy_asc(self, hb):
        rows = export_table(hb, sort_by="occupancy", sort_dir="asc")["rows"]
        assert rows[0]["occupancy"] == "0.2000"

    def test_unknown_sort_by_falls_back_to_index(self, hb):
        labels = [row["pair_label"] for row in export_table(hb, sort_by="nope")["rows"]]
        expected = [get_exporter(hb).get_pair_label(hb, index)
                    for index in (0, 1, 2)]
        assert labels == expected

    def test_label_format_matches_exporter(self, hb):
        rows = export_table(hb)["rows"]
        expected = get_exporter(hb).get_pair_label(hb, 0)
        assert rows[0]["pair_label"] == expected


# ============================================================
# 下钻与降采样（T006）
# ============================================================

class TestPairDrilldown:

    def test_series_lengths_match_frames(self, hb):
        data = pair_drilldown(hb, 0)
        assert data["existence"] == [1, 0, 0, 0, 0]
        for s in data["series"].values():
            assert len(s["times"]) == FRAMES and len(s["values"]) == FRAMES

    def test_metadata_and_groups(self, hb):
        data = pair_drilldown(hb, 0)
        assert data["pair_label"] == get_exporter(hb).get_pair_label(hb, 0)
        assert data["groups"][0]["residue_name"] == "LYS"

    def test_string_metrics_excluded_from_series(self, hb):
        it = hb
        it.metrics["pistacking_type"] = np.full(
            (it.n_pairs, FRAMES), "N", dtype="U1")
        data = pair_drilldown(it, 0)
        assert "pistacking_type" not in data["series"]

    def test_bucket_aggregation_keeps_single_frame_event(self):
        from DuIvyInteractions.visualizers.data import MAX_SERIES_POINTS
        frames = 3000
        existence = np.zeros((1, frames), dtype=bool)
        existence[0, 2] = True
        long = Interaction(
            interaction_type="hydrogen_bond",
            groups=[(_group(0, "H_donor", 0, "C", "LYS", 0, with_h=True),
                     _group(1, "H_acceptor", 1, "O", "ASP", 1))],
            existence=existence,
            metrics={"distance": np.full((1, frames), 3.0)},
            times=np.arange(frames, dtype=float))
        data = pair_drilldown(long, 0)
        assert data["downsampled"] is True
        assert len(data["times"]) == MAX_SERIES_POINTS
        assert sum(data["existence"]) == 1

    def test_bucket_metric_mean(self):
        frames = 3000
        values = np.zeros((1, frames))
        values[0, 1] = 4.0
        values[0, 2] = 6.0
        short = Interaction(
            interaction_type="hydrogen_bond",
            groups=[(_group(0, "H_donor", 0, "C", "LYS", 0, with_h=True),
                     _group(1, "H_acceptor", 1, "O", "ASP", 1))],
            existence=np.ones((1, frames), dtype=bool),
            metrics={"distance": values}, times=np.arange(frames, dtype=float))
        data = pair_drilldown(short, 0)
        assert data["series"]["distance"]["values"][1] == 5.0

    def test_long_series_downsampled_consistently(self):
        from DuIvyInteractions.visualizers.data import MAX_SERIES_POINTS
        frames = 3000
        long = Interaction(
            interaction_type="hydrogen_bond",
            groups=[(_group(0, "H_donor", 0, "C", "LYS", 0, with_h=True),
                     _group(1, "H_acceptor", 1, "O", "ASP", 1))],
            existence=np.ones((1, frames), dtype=bool),
            metrics={"distance": np.linspace(2.0, 4.0, frames).reshape(1, frames)},
            times=np.arange(frames, dtype=float))
        data = pair_drilldown(long, 0)
        assert len(data["times"]) == MAX_SERIES_POINTS
        assert len(data["existence"]) == MAX_SERIES_POINTS
        for series in data["series"].values():
            assert len(series["times"]) == MAX_SERIES_POINTS
            assert len(series["values"]) == MAX_SERIES_POINTS


# ============================================================
# 图表 JSON（T008）
# ============================================================

class TestFigures:

    def test_pair_figure_structure(self, hb):
        fig = fig_pair(pair_drilldown(hb, 0))
        assert set(fig) == {"data", "layout"}
        assert fig["data"][0]["name"] == "existence"
        assert fig["layout"]["hovermode"] == "x unified"

    def test_pair_figure_metric_subplot(self, hb):
        fig = fig_pair(pair_drilldown(hb, 0))
        assert fig["layout"]["yaxis2"]["anchor"] == "x2"
        assert fig["data"][1]["yaxis"] == "y2"

    def test_pair_figure_nan_becomes_null(self, hb):
        fig = fig_pair(pair_drilldown(hb, 0))
        assert fig["data"][1]["y"][2] is None

    def test_pair_figure_without_metrics(self, hb):
        hb.metrics.clear()
        fig = fig_pair(pair_drilldown(hb, 0))
        assert len(fig["data"]) == 1


# ============================================================
# 最近数据源（T007）
# ============================================================

class TestRecent:

    def test_add_and_load_roundtrip(self, monkeypatch, tmp_path):
        monkeypatch.setattr("DuIvyInteractions.visualizers.data.RECENT_PATH",
                            tmp_path / "recent.json")
        add_recent("/a")
        add_recent("/b")
        add_recent("/a")
        assert load_recent()[0] == "/a"
        assert len(load_recent()) == 2

    def test_corrupt_file_silently_rebuilt(self, monkeypatch, tmp_path):
        bad = tmp_path / "recent.json"
        bad.write_text("not json")
        monkeypatch.setattr("DuIvyInteractions.visualizers.data.RECENT_PATH", bad)
        assert load_recent() == []
        add_recent("/x")
        assert load_recent() == ["/x"]


# ============================================================
# 矩阵 / 单元格 / 帧 / 元数据（可视化数据层）
# ============================================================

class TestMatrix:

    def test_runs_from_mask(self):
        assert runs_from_mask([0, 1, 1, 0, 1]) == [[1, 3], [4, 5]]

    def test_runs_all_false(self):
        assert runs_from_mask([0, 0, 0]) == []

    def test_matrix_sort_by_occupancy(self, hb):
        view = matrix_view(hb)
        assert view["n_frames"] == FRAMES
        assert [row["idx"] for row in view["rows"]] == [1, 0, 2]

    def test_matrix_row_runs(self, hb):
        row = next(r for r in matrix_view(hb)["rows"] if r["idx"] == 1)
        assert row["runs"] == [[0, 2]]

    def test_matrix_filter_search(self, hb):
        assert len(matrix_view(hb, search="ASP100")["rows"]) == 1

    def test_matrix_filter_occupancy(self, hb):
        assert [r["idx"] for r in matrix_view(hb, min_occupancy=0.3)["rows"]] == [1]

    def test_matrix_sort_label(self, hb):
        labels = [r["label"] for r in
                  matrix_view(hb, sort_by="label", sort_dir="asc")["rows"]]
        assert labels == sorted(labels)

    def test_matrix_sort_direction(self, hb):
        asc = [r["idx"] for r in
               matrix_view(hb, sort_by="occupancy", sort_dir="asc")["rows"]]
        desc = [r["idx"] for r in
                matrix_view(hb, sort_by="occupancy", sort_dir="desc")["rows"]]
        assert asc == [0, 2, 1] and desc == [1, 0, 2]

    def test_cell_values(self, hb):
        cell = cell_values(hb, 0, 0)
        assert cell["exists"] is True and cell["time"] == 0.0
        assert set(cell["values"]) == {"distance"}

    def test_cell_inactive_nan_as_none(self, hb):
        assert cell_values(hb, 0, 1)["values"]["distance"] is None

    def test_cell_out_of_range(self, hb):
        with pytest.raises(ValueError):
            cell_values(hb, 99, 0)

    def test_frame_view(self, hb):
        view = frame_view(hb, 0)
        assert view["n_active"] == 3 and len(view["pairs"]) == 3

    def test_frame_view_zero_pairs(self):
        empty = Interaction(
            interaction_type="test", groups=[],
            existence=np.zeros((0, FRAMES), dtype=bool),
            metrics={"distance": np.zeros((0, FRAMES))},
            times=np.arange(FRAMES, dtype=float))
        view = frame_view(empty, 0)
        assert view["n_active"] == 0 and view["pairs"] == []


class TestMetadata:

    def test_meta_fields(self, hb):
        meta = interaction_meta(hb, "x.h5")
        assert meta["h5_path"] == "x.h5"
        assert meta["n_pairs"] == 3 and meta["n_frames"] == FRAMES
        assert meta["time_start"] == 0.0 and meta["time_end"] == 40.0
        assert meta["metrics"] == ["distance"]

    def test_occupancy_summary(self, hb):
        occ = interaction_meta(hb)["occupancy"]
        assert occ["count"] == 3 and sum(occ["hist"]) == 3

    def test_metric_stats(self, hb):
        stats = metric_stats(hb)
        assert stats[0]["name"] == "distance"
        assert sum(stats[0]["hist"]) == stats[0]["count"]
        assert stats[0]["min"] <= stats[0]["mean"] <= stats[0]["max"]

    def test_string_metric_skipped(self, hb):
        hb.metrics["pistacking_type"] = np.full((3, FRAMES), "N", dtype="U1")
        assert [s["name"] for s in metric_stats(hb)] == ["distance"]

    def test_metric_unit(self, hb):
        assert metric_stats(hb)[0]["unit"] == "Å"

    def test_file_info(self, hb, h5_dir):
        meta = interaction_meta(hb, str(h5_dir / "hydrogen_bond.h5"))
        assert meta["file"]["size_kb"] > 0
        assert meta["file"]["format_version"] == "1.0"
        assert meta["file"]["n_interactions"] == 1

    def test_frame_activity(self, hb):
        assert frame_activity(hb) == [3, 1, 0, 0, 0]


class TestPairFrames:

    def test_pagination(self, hb):
        page = pair_frames(hb, 0, start=0, limit=2)
        assert [r["frame"] for r in page["rows"]] == [0, 1]
        assert page["total"] == FRAMES
        assert [r["frame"] for r in pair_frames(hb, 0, start=4, limit=2)["rows"]] == [4]

    def test_flags_included(self, hyd):
        assert pair_frames(hyd, 0)["rows"][0]["flags"]["redundant_pi_stacking"] is True

    def test_bad_index(self, hb):
        with pytest.raises(ValueError):
            pair_frames(hb, 99)


# ============================================================
# 路由冒烟
# ============================================================

@pytest.fixture
def client(h5_dir):
    """已载入测试数据源的 TestClient。"""
    from fastapi.testclient import TestClient
    from DuIvyInteractions.web.app import create_app
    client = TestClient(create_app())
    resp = client.post("/select", data={"path": str(h5_dir)},
                       follow_redirects=False)
    assert resp.status_code == 303
    return client


class TestRoutes:

    def test_select_redirects(self, h5_dir):
        from fastapi.testclient import TestClient
        from DuIvyInteractions.web.app import create_app
        client = TestClient(create_app())
        resp = client.post("/select", data={"path": str(h5_dir)},
                           follow_redirects=False)
        assert resp.status_code == 303
        assert resp.headers["location"].startswith("/?warn=")
        assert "skipped" in client.get(resp.headers["location"]).text

    def test_select_invalid_path_reports_error(self):
        from fastapi.testclient import TestClient
        from DuIvyInteractions.web.app import create_app
        client = TestClient(create_app())
        resp = client.post("/select", data={"path": "/nope"}, follow_redirects=False)
        assert resp.status_code == 303 and "error=" in resp.headers["location"]
        assert "banner" in client.get(resp.headers["location"]).text

    def test_index_has_single_open_control(self, client):
        html = client.get("/").text
        assert html.count('class="source-form') == 1
        assert 'name="path"' in html and "browseBtn" in html

    def test_index_wires_browse_control(self, client):
        html = client.get("/").text
        assert "fetch('/api/pick'" in html
        assert "browseBtn" in html

    def test_type_page_links_back_to_source(self, client):
        html = client.get("/type/hydrogen_bond").text
        assert "change data source" in html
        assert "data-kind" not in html

    def test_type_page_panels(self, client):
        html = client.get("/type/hydrogen_bond").text
        assert "metadata" in html and "matrixCanvas" in html
        assert "frameSlider" in html and "occupancy distribution" in html
        assert "active pairs per frame" in html and "format version" in html
        assert "aggregates" in html and "apply filters" in html

    def test_type_page_unknown_type_404(self, client):
        assert client.get("/type/nope").status_code == 404

    def test_pair_page_renders_chart_and_frames(self, client):
        html = client.get("/pair/hydrogen_bond/0").text
        assert "pairChart" in html and "plotly-basic.min.js" in html
        assert "framesTable" in html and "pairFigure" in html

    def test_pair_page_out_of_range_404(self, client):
        assert client.get("/pair/hydrogen_bond/99").status_code == 404

    def test_openapi_disabled(self, client):
        assert client.get("/openapi.json").status_code == 404


class TestApi:

    def test_api_matrix(self, client):
        data = client.get("/api/matrix/hydrogen_bond").json()
        assert data["n_frames"] == FRAMES and len(data["rows"]) == 3

    def test_api_matrix_clamps_occupancy(self, client):
        resp = client.get("/api/matrix/hydrogen_bond?min_occupancy=nan")
        assert len(resp.json()["rows"]) == 3
        resp = client.get("/api/matrix/hydrogen_bond?min_occupancy=9")
        assert resp.json()["rows"] == []

    def test_api_cell(self, client):
        data = client.get("/api/cell/hydrogen_bond/0/0").json()
        assert data["exists"] is True and "distance" in data["values"]
        assert data["units"]["distance"] == "Å"

    def test_api_frame(self, client):
        data = client.get("/api/frame/hydrogen_bond/0").json()
        assert data["n_active"] == 3 and data["units"]["distance"] == "Å"

    def test_api_frame_limit(self, client):
        data = client.get("/api/frame/hydrogen_bond/0?limit=2").json()
        assert data["n_active"] == 3 and data["shown"] == 2
        assert data["limit_max"] == 2000

    def test_download_unknown_type_text_404(self, client):
        resp = client.get("/download/csv/nope")
        assert resp.status_code == 404
        assert resp.headers["content-type"].startswith("text/plain")

    def test_api_frames(self, client):
        data = client.get("/api/frames/hydrogen_bond/0?start=0").json()
        assert len(data["rows"]) == FRAMES
        assert data["units"]["distance"] == "Å"

    def test_api_cell_out_of_range_404(self, client):
        resp = client.get("/api/cell/hydrogen_bond/99/0")
        assert resp.status_code == 404
        assert resp.headers["content-type"].startswith("application/json")

    def test_invalid_page_parameter_renders_html(self, client):
        resp = client.get("/pair/hydrogen_bond/notint")
        assert resp.status_code == 422
        assert resp.headers["content-type"].startswith("text/html")

    def test_invalid_api_parameter_returns_json(self, client):
        resp = client.get("/api/cell/hydrogen_bond/notint/0")
        assert resp.status_code == 422
        assert resp.headers["content-type"].startswith("application/json")

    def test_pair_page_has_summary(self, client):
        html = client.get("/pair/hydrogen_bond/0").text
        assert "summary" in html and "occupancy" in html

    def test_csv_real_pi_stacking(self, h5_data_dir):
        from fastapi.testclient import TestClient
        from DuIvyInteractions.web.app import create_app
        client = TestClient(create_app())
        client.post("/select", data={"path": str(h5_data_dir)})
        resp = client.get("/download/csv/pi_stacking")
        assert resp.status_code == 200
        lines = resp.text.strip().splitlines()
        assert "avg_distance" in lines[0]
        assert "mode_" not in lines[0]
        assert len(lines) == 11  # header + 10 pairs（全量，非分页）

    def test_type_page_clamps_nan_occupancy(self, client):
        html = client.get("/type/hydrogen_bond?min_occupancy=nan").text
        assert '"min_occupancy": 0.0' in html

    def test_type_page_empty_occupancy_filter(self, client):
        assert client.get("/type/hydrogen_bond?min_occupancy=").status_code == 200

    def test_type_page_aggregates_respect_min_occupancy(self, client):
        html = client.get("/type/hydrogen_bond?min_occupancy=0.3").text
        assert "aggregates · 1 pairs" in html

    def test_invalid_sort_falls_back(self, client):
        html = client.get("/type/hydrogen_bond?sort_by=evil&sort_dir=evil").text
        assert '"sort_by": "occupancy"' in html
        assert '"sort_dir": "desc"' in html

    def test_type_page_pi_stacking_string_panel(self, h5_data_dir):
        from fastapi.testclient import TestClient
        from DuIvyInteractions.web.app import create_app
        client = TestClient(create_app())
        client.post("/select", data={"path": str(h5_data_dir)})
        html = client.get("/type/pi_stacking").text
        assert "string metrics" in html
        assert "pistacking_type" in html

    def test_csv_download_fields_and_rows(self, client):
        resp = client.get("/download/csv/hydrogen_bond")
        assert resp.headers["content-type"].startswith("text/csv")
        lines = resp.text.strip().split("\n")
        assert lines[0].startswith("pair_label,occupancy,avg_distance")
        assert len(lines) == 4  # header + 3 rows


# ============================================================
# 原生目录对话框（唯一的数据源入口）
# ============================================================

class TestPicker:

    @pytest.fixture
    def client(self):
        from fastapi.testclient import TestClient
        from DuIvyInteractions.web.app import create_app
        return TestClient(create_app())

    def test_pick_returns_selected_path(self, client, h5_dir, monkeypatch):
        monkeypatch.setattr("DuIvyInteractions.web.app.pick_directory",
                            lambda: str(h5_dir))
        resp = client.post("/api/pick")
        assert resp.status_code == 200 and resp.json()["path"] == str(h5_dir)

    def test_pick_cancel_returns_empty_path(self, client, monkeypatch):
        monkeypatch.setattr("DuIvyInteractions.web.app.pick_directory", lambda: "")
        assert client.post("/api/pick").json() == {"path": ""}

    def test_pick_failure_reports_error(self, client, monkeypatch):
        def boom():
            raise RuntimeError("no display name and no $DISPLAY")
        monkeypatch.setattr("DuIvyInteractions.web.app.pick_directory", boom)
        resp = client.post("/api/pick")
        assert resp.status_code == 400 and "DISPLAY" in resp.json()["error"]


# ============================================================
# 字符串指标 / 空数据 / 单对摘要
# ============================================================

class TestStringMetrics:

    def test_export_table_skips_string_metrics(self, hb):
        hb.metrics["pistacking_type"] = np.full((3, FRAMES), "N", dtype="U1")
        columns = export_table(hb)["columns"]
        assert "mode_pistacking_type" not in columns
        assert "avg_distance" in columns

    def test_string_stats(self, hb):
        values = np.array([["P", "N", "N", "N", "N"]] * 3, dtype="U1")
        hb.metrics["pistacking_type"] = values
        stats = string_stats(hb)
        assert stats[0]["name"] == "pistacking_type"
        assert stats[0]["top"][0]["count"] >= 1

    def test_frame_view_lists_string_metrics(self, hb):
        hb.metrics["pistacking_type"] = np.full((3, FRAMES), "N", dtype="U1")
        view = frame_view(hb, 0)
        assert "pistacking_type" in view["string_metrics"]
        assert view["pairs"][0]["values"]["pistacking_type"] == "N"

    def test_pair_summary_string_metrics(self, hb):
        values = np.array([["P", "N", "N", "N", "N"]] * 3, dtype="U1")
        hb.metrics["pistacking_type"] = values
        summary = pair_summary(hb, 0)
        assert summary["string_metrics"]["pistacking_type"][0]["value"] == "P"


class TestEmptyData:

    @staticmethod
    def _zero_frame() -> Interaction:
        groups = [(_group(2 * i, "H_donor", 10 * i, "C", "LYS", i, with_h=True),
                   _group(2 * i + 1, "H_acceptor", 10 * i + 1, "O", "ASP", i))
                  for i in range(2)]
        return Interaction(
            interaction_type="hydrogen_bond", groups=groups,
            existence=np.zeros((2, 0), dtype=bool),
            metrics={"distance": np.zeros((2, 0))}, times=np.zeros((0,)))

    def test_overview_zero_frames(self):
        card = overview([("hydrogen_bond", "x.h5", self._zero_frame())])[0]
        assert card["n_pairs"] == 2
        assert card["occupancy"]["count"] == 0

    def test_meta_zero_frames(self):
        meta = interaction_meta(self._zero_frame())
        assert meta["time_start"] is None
        assert meta["occupancy"]["count"] == 0

    def test_matrix_zero_frames(self):
        rows = matrix_view(self._zero_frame())["rows"]
        assert [row["runs"] for row in rows] == [[], []]

    def test_pair_frames_zero_frames(self):
        assert pair_frames(self._zero_frame(), 0)["rows"] == []


class TestPairSummary:

    def test_summary_counts_and_occupancy(self, hb):
        summary = pair_summary(hb, 0)
        assert summary["occupancy"] == 1 / FRAMES
        assert summary["metrics"]["distance"]["count"] == 1

    def test_summary_bad_index(self, hb):
        with pytest.raises(ValueError):
            pair_summary(hb, 99)

    def test_summary_flags(self, hyd):
        flags = pair_summary(hyd, 0)["flags"]
        assert flags["redundant_pi_stacking"]["active"] == 1
        assert flags["redundant_pi_stacking"]["frames"] == FRAMES


class TestAggregateTable:

    def test_rows_and_pagination(self, hb):
        table = aggregate_table(hb)
        assert table["total"] == 3 and len(table["rows"]) == 3
        assert table["pages"] == 1 and table["page"] == 1

    def test_sort_by_label(self, hb):
        rows = aggregate_table(hb, sort_by="label", sort_dir="asc")["rows"]
        labels = [row["pair_label"] for row in rows]
        assert labels == sorted(labels)

    def test_search_filter(self, hb):
        table = aggregate_table(hb, search="ASP100")
        assert table["total"] == 1

    def test_min_occupancy_filter(self, hb):
        table = aggregate_table(hb, min_occupancy=0.3)
        assert table["total"] == 1


class TestMultipleInteractions:

    def test_all_interactions_in_file_loaded(self, tmp_path):
        source = tmp_path / "multi"
        source.mkdir()
        save_interactions(
            [_interaction("hydrogen_bond", 1, False),
             _interaction("hydrophobic", 1, False)],
            str(source / "multi.h5"))
        assert {name for name, _, _ in parse_datasource(str(source))} == \
            {"hydrogen_bond", "hydrophobic"}

    def test_duplicate_type_gets_filename_suffix(self, tmp_path):
        source = tmp_path / "dup"
        source.mkdir()
        save_interactions([_interaction("hydrogen_bond", 1, False)],
                          str(source / "a.h5"))
        save_interactions([_interaction("hydrogen_bond", 1, False)],
                          str(source / "b.h5"))
        names = [name for name, _, _ in parse_datasource(str(source))]
        assert names == ["hydrogen_bond", "hydrogen_bond@b.h5"]

    def test_third_duplicate_gets_tilde_suffix(self, tmp_path):
        source = tmp_path / "triple"
        source.mkdir()
        save_interactions([_interaction("hydrogen_bond", 1, False)
                           for _ in range(3)], str(source / "multi.h5"))
        names = [name for name, _, _ in parse_datasource(str(source))]
        assert names == ["hydrogen_bond", "hydrogen_bond@multi.h5",
                         "hydrogen_bond@multi.h5~2"]
