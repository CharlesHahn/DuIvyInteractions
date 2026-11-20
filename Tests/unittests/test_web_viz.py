# -*- coding: utf-8 -*-
"""Web 可视化单元测试：纯函数层（visualizers/）≥ 20 用例（spec SC-006）+ 路由冒烟。

不依赖 MDAnalysis 轨迹：用 tmp 目录构造 h5 fixture。
"""

import shutil
from pathlib import Path

import numpy as np
import pytest

from DuIvyInteractions.core.datas import AtomData, Group, Interaction
from DuIvyInteractions.io.h5 import save_interactions
from DuIvyInteractions.visualizers import (
    add_recent, fig_existence_timeline, fig_metric_timeseries,
    fig_occupancy_hist, get_exporter, load_recent, overview,
    pair_drilldown, parse_datasource, type_table,
)

TEMP_DIR = Path(__file__).parent.parent.parent / "test_temp_viz"
FRAMES = 5


@pytest.fixture(autouse=True)
def setup_temp_dir():
    TEMP_DIR.mkdir(exist_ok=True)
    yield
    shutil.rmtree(TEMP_DIR, ignore_errors=True)


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
def h5_dir():
    """输出目录：2 个 h5（hydrogen_bond 无 flags / hydrophobic 有 flags）。"""
    d = TEMP_DIR / "out"
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
        assert "skip corrupt h5 corrupt.h5" in capsys.readouterr().out

    def test_single_file(self, h5_dir):
        items = parse_datasource(str(h5_dir / "hydrophobic.h5"))
        assert len(items) == 1 and items[0][0] == "hydrophobic"

    def test_nonexistent_path_raises(self):
        with pytest.raises(ValueError, match="path not found"):
            parse_datasource(str(TEMP_DIR / "nope"))

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

    def test_card_occupancies_length_and_values(self, h5_dir):
        cards = overview(parse_datasource(str(h5_dir)))
        by_type = {c["type"]: c for c in cards}
        occ = by_type["hydrophobic"]["occupancies"]
        assert len(occ) == 2 and occ[0] == 0.2 and occ[1] == 0.4


# ============================================================
# 类型表格（T005）
# ============================================================

class TestTypeTable:

    def test_rows_and_columns(self, hb):
        table = type_table(hb)
        assert table["total"] == 3
        assert table["columns"][0:3] == ["index", "pair_label", "occupancy"]
        assert "avg_distance" in table["columns"]
        assert "flag_redundant_pi_stacking" not in table["columns"]

    def test_columns_include_flags_when_present(self, hyd):
        table = type_table(hyd)
        assert "flag_redundant_pi_stacking" in table["columns"]

    def test_flag_fraction(self, hyd):
        table = type_table(hyd)
        row0 = next(r for r in table["rows"] if r["index"] == 0)
        assert row0["flag_redundant_pi_stacking"] == 1.0  # 1/1 活跃帧
        row1 = next(r for r in table["rows"] if r["index"] == 1)
        assert row1["flag_redundant_pi_stacking"] == 0.0  # 1/1 活跃帧未标记

    def test_pagination(self, hb):
        table = type_table(hb, page=1)
        assert table["pages"] == 1 and table["page"] == 1
        assert len(table["rows"]) == 3

    def test_page_clamp(self, hb):
        assert type_table(hb, page=99)["page"] == 1

    def test_search_filters_by_label(self, hb):
        table = type_table(hb, search="LYS0:")
        assert table["total"] == 1
        assert table["rows"][0]["pair_label"] == get_exporter(hb).get_pair_label(hb, 0)

    def test_search_case_insensitive(self, hb):
        assert type_table(hb, search="lys0:")["total"] == 1

    def test_sort_occupancy_desc(self, hb):
        rows = type_table(hb, sort_by="occupancy", sort_dir="desc")["rows"]
        assert rows[0]["occupancy"] == 2 / FRAMES

    def test_sort_occupancy_asc(self, hb):
        rows = type_table(hb, sort_by="occupancy", sort_dir="asc")["rows"]
        assert rows[0]["occupancy"] == 1 / FRAMES

    def test_unknown_sort_by_falls_back_to_index(self, hb):
        rows = type_table(hb, sort_by="nope")["rows"]
        assert [r["index"] for r in rows] == [0, 1, 2]

    def test_label_format_matches_exporter(self, hb):
        rows = type_table(hb)["rows"]
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

    def test_flags_included_when_present(self, hyd):
        assert pair_drilldown(hyd, 0)["flags"]["redundant_pi_stacking"] == [1, 0, 0, 0, 0]

    def test_string_metrics_excluded_from_series(self, hb):
        it = hb
        it.metrics["pistacking_type"] = np.full(
            (it.n_pairs, FRAMES), "N", dtype="U1")
        data = pair_drilldown(it, 0)
        assert "pistacking_type" not in data["series"]

    def test_downsample_caps_points(self):
        from DuIvyInteractions.visualizers.data import _downsample, MAX_SERIES_POINTS
        xs = list(range(5000))
        out_x, out_y = _downsample(xs, [v * 2.0 for v in xs], 1000)
        assert len(out_x) == 1000
        assert out_x[0] == 0 and out_x[-1] == 4999
        assert out_y[0] == 0.0 and out_y[-1] == 9998.0

    def test_downsample_small_input_unchanged(self):
        from DuIvyInteractions.visualizers.data import _downsample
        assert _downsample([0, 1], [1.0, 2.0], 2000) == ([0, 1], [1.0, 2.0])


# ============================================================
# 图表 JSON（T008）
# ============================================================

class TestFigures:

    def test_timeline_structure(self, hb):
        fig = fig_existence_timeline(list(hb.times),
                                     hb.existence[0].astype(int).tolist(), "L")
        assert list(fig.keys()) == ["data", "layout"]
        assert fig["data"][0]["type"] == "scatter"
        assert len(fig["data"][0]["y"]) == FRAMES

    def test_metric_masked_nan_to_none(self, hyd):
        values = hyd.metrics["distance"][0].tolist()
        values[2] = float("nan")
        fig = fig_metric_timeseries(list(hyd.times), values, "distance", "L")
        assert fig["data"][0]["y"][2] is None

    def test_occupancy_hist_bins(self, hb):
        occ = hb.occupancy().tolist()
        fig = fig_occupancy_hist(occ, "hydrogen_bond")
        assert fig["data"][0]["type"] == "histogram"
        assert fig["data"][0]["xbins"]["size"] == 0.05


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
# 路由冒烟（T012/T014）
# ============================================================

class TestRoutes:

    @pytest.fixture
    def client(self, h5_dir):
        from fastapi.testclient import TestClient
        from DuIvyInteractions.web.app import create_app
        app = create_app()
        client = TestClient(app)
        resp = client.post("/select", data={"path": str(h5_dir)},
                           follow_redirects=False)
        assert resp.status_code == 303
        return client

    def test_select_redirects(self, h5_dir):
        from fastapi.testclient import TestClient
        from DuIvyInteractions.web.app import create_app
        resp = TestClient(create_app()).post(
            "/select", data={"path": str(h5_dir)}, follow_redirects=False)
        assert resp.status_code == 303 and resp.headers["location"] == "/"

    def test_select_invalid_path_warns_and_redirects(self, capsys):
        from fastapi.testclient import TestClient
        from DuIvyInteractions.web.app import create_app
        client = TestClient(create_app(), raise_server_exceptions=False)
        resp = client.post("/select", data={"path": "/nope"}, follow_redirects=False)
        assert resp.status_code == 303  # redirect 回总览（FR-014：告警不崩溃）
        # 无半成品数据源
        client = TestClient(create_app())
        client.post("/select", data={"path": "/nope"}, follow_redirects=True)
        assert "尚无数据源" in client.get("/").text

    def test_index_shows_cards_after_select(self, client):
        html = client.get("/").text
        assert "hydrogen_bond" in html and "hydrophobic" in html

    def test_type_page_200_with_table(self, client):
        html = client.get("/type/hydrogen_bond").text
        assert "pair_label" in html and "ASP" in html

    def test_type_page_unknown_type_raises(self, client):
        from fastapi import HTTPException
        with pytest.raises(ValueError, match="not loaded"):
            client.get("/type/water_bridge")

    def test_pair_page_200(self, client):
        html = client.get("/pair/hydrogen_bond/0").text
        assert "figExistence" in html and "plotly-basic.min.js" in html

    def test_figures_json_structure(self, client):
        figs = client.get("/api/figures/hydrogen_bond/0.json").json()
        assert all(set(fig) >= {"data", "layout", "div"} for fig in figs)
        assert figs[0]["div"] == "figExistence"

    def test_csv_download_fields_and_rows(self, client):
        resp = client.get("/download/csv/hydrogen_bond")
        assert resp.headers["content-type"].startswith("text/csv")
        lines = resp.text.strip().split("\n")
        assert lines[0].startswith("index,pair_label,occupancy")
        assert len(lines) == 4  # header + 3 rows
