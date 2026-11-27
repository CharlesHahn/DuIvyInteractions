# -*- coding: utf-8 -*-
"""Web 多数据源测试：merge_batches 纯函数 + 批次追加/替换/移除/清空 + 跨批次重名。

不依赖真实轨迹：用 tmp_path 构造最小 h5。
"""

import re
from pathlib import Path
from urllib.parse import unquote

import numpy as np
import pytest
from fastapi.testclient import TestClient

from DuIvyInteractions.core.datas import AtomData, Group, Interaction
from DuIvyInteractions.io.h5 import save_interactions
from DuIvyInteractions.visualizers.data import merge_batches
from DuIvyInteractions.web import app as web_app

FRAMES = 4
SUFFIXED = "hydrogen_bond@hydrogen_bond.h5"
STRUCTURE_NOT_SELECTED = "structure not selected"


# ============================================================
# fixture 构造
# ============================================================

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


def _interaction(itype: str, n_pairs: int = 2) -> Interaction:
    """最小 Interaction：首帧全部活跃，其余帧不活跃。"""
    existence = np.zeros((n_pairs, FRAMES), dtype=bool)
    existence[:, 0] = True
    distance = np.full((n_pairs, FRAMES), 3.0)
    distance[~existence] = np.nan
    groups = [(_group(2 * i, "H_donor", 10 * i, "C", "LYS", i, with_h=True),
               _group(2 * i + 1, "H_acceptor", 10 * i + 1, "O", "ASP", 100 + i))
              for i in range(n_pairs)]
    return Interaction(
        interaction_type=itype, groups=groups, existence=existence,
        metrics={"distance": distance},
        times=np.arange(FRAMES, dtype=float) * 10.0)


def _write_dir(directory: Path, types) -> Path:
    """建数据源目录：每个类型写一个 <type>.h5。"""
    directory.mkdir()
    for name in types:
        save_interactions([_interaction(name)], str(directory / f"{name}.h5"))
    return directory


def _h5(directory: str, itype: str) -> str:
    """纯函数测试用的 h5 路径（不落盘）。"""
    return f"/data/{directory}/{itype}.h5"


def _items(directory: str, itype: str) -> list:
    """一个批次的最小 items（名称等同 interaction_type）。"""
    return [(itype, _h5(directory, itype), _interaction(itype))]


@pytest.fixture(autouse=True)
def isolate_recent(monkeypatch, tmp_path):
    """隔离最近数据源文件，避免写入用户缓存。"""
    monkeypatch.setattr("DuIvyInteractions.visualizers.data.RECENT_PATH",
                        tmp_path / "recent.json")


@pytest.fixture
def first_dir(tmp_path):
    """批次一：hydrogen_bond + hydrophobic。"""
    return _write_dir(tmp_path / "first", ("hydrogen_bond", "hydrophobic"))


@pytest.fixture
def second_dir(tmp_path):
    """批次二：hydrogen_bond（与批次一同名）+ pi_stacking。"""
    return _write_dir(tmp_path / "second", ("hydrogen_bond", "pi_stacking"))


@pytest.fixture
def make_client(monkeypatch):
    """创建 TestClient 并捕获 create_app 内部的会话状态。"""
    def build(**kwargs):
        captured: dict = {}
        original = web_app._register_pages

        def capture(app, state, browse_mode, root):
            captured["state"] = state
            original(app, state, browse_mode, root)

        monkeypatch.setattr(web_app, "_register_pages", capture)
        client = TestClient(web_app.create_app(**kwargs))
        return client, captured["state"]

    return build


@pytest.fixture
def loaded(make_client, first_dir, second_dir):
    """已连续载入两个批次的 (client, state)。"""
    client, state = make_client()
    assert _select(client, first_dir).headers["location"] == "/"
    assert _select(client, second_dir).headers["location"] == "/"
    return client, state


def _select(client, path):
    """POST /select（不跟随跳转）。"""
    return client.post("/select", data={"path": str(path)},
                       follow_redirects=False)


# ============================================================
# merge_batches 纯函数
# ============================================================

class TestMergeBatches:

    def test_single_batch_passthrough(self):
        batch = [("hydrogen_bond", _h5("a", "hydrogen_bond"),
                  _interaction("hydrogen_bond")),
                 ("hydrophobic", _h5("a", "hydrophobic"),
                  _interaction("hydrophobic"))]
        assert merge_batches([batch]) == batch

    def test_single_empty_batch(self):
        assert merge_batches([[]]) == []

    def test_no_batches(self):
        assert merge_batches([]) == []

    def test_cross_batch_duplicate_gets_filename_suffix(self):
        merged = merge_batches([_items("a", "hydrogen_bond"),
                                _items("b", "hydrogen_bond")])
        assert [name for name, _, _ in merged] == ["hydrogen_bond", SUFFIXED]

    def test_third_duplicate_gets_tilde_suffix(self):
        merged = merge_batches([_items(d, "hydrogen_bond")
                                for d in ("a", "b", "c")])
        assert [name for name, _, _ in merged] == \
            ["hydrogen_bond", SUFFIXED, f"{SUFFIXED}~2"]

    def test_batch_and_inner_order_preserved(self):
        merged = merge_batches([
            [("hydrogen_bond", _h5("a", "hydrogen_bond"),
              _interaction("hydrogen_bond")),
             ("hydrophobic", _h5("a", "hydrophobic"),
              _interaction("hydrophobic"))],
            [("hydrogen_bond", _h5("b", "hydrogen_bond"),
              _interaction("hydrogen_bond")),
             ("pi_stacking", _h5("b", "pi_stacking"),
              _interaction("pi_stacking"))],
        ])
        assert [name for name, _, _ in merged] == \
            ["hydrogen_bond", "hydrophobic", SUFFIXED, "pi_stacking"]
        assert [h5 for _, h5, _ in merged] == \
            [_h5("a", "hydrogen_bond"), _h5("a", "hydrophobic"),
             _h5("b", "hydrogen_bond"), _h5("b", "pi_stacking")]

    def test_distinct_types_keep_plain_names(self):
        merged = merge_batches([_items("a", "hydrogen_bond"),
                                _items("b", "salt_bridge")])
        assert [name for name, _, _ in merged] == \
            ["hydrogen_bond", "salt_bridge"]

    def test_names_globally_unique(self):
        merged = merge_batches([_items(d, "hydrogen_bond") for d in "abc"] +
                               [_items(d, "hydrophobic") for d in "abc"])
        names = [name for name, _, _ in merged]
        assert len(names) == len(set(names)) == 6


# ============================================================
# /select 追加、同路径替换
# ============================================================

class TestSelectAppend:

    def test_two_batches_appended(self, loaded, first_dir, second_dir):
        _, state = loaded
        assert [source["name"] for source in state["sources"]] == \
            ["first", "second"]
        assert [source["path"] for source in state["sources"]] == \
            [str(first_dir.resolve()), str(second_dir.resolve())]
        assert len(state["items"]) == 4

    def test_merged_items_ordered(self, loaded):
        _, state = loaded
        assert [name for name, _, _ in state["items"]] == \
            ["hydrogen_bond", "hydrophobic", SUFFIXED, "pi_stacking"]

    def test_overview_shows_two_card_groups(self, loaded):
        client, _ = loaded
        html = client.get("/").text
        assert html.count('class="card"') == 4
        assert html.count('class="card-label">hydrogen_bond') == 2

    def test_sidebar_lists_every_type(self, loaded):
        client, _ = loaded
        html = client.get("/").text
        aside = html.split('<aside class="sidebar">', 1)[1].split("</aside>", 1)[0]
        assert aside.count('class="side-item"') == 4
        base_label = '<span class="side-label">hydrogen_bond</span>'
        assert aside.count(base_label) == 2
        assert SUFFIXED not in aside
        hrefs = re.findall(r'data-t="hydrogen_bond" href="([^"]+)"', aside)
        assert len(hrefs) == 2 and len(set(hrefs)) == 2
        assert sum("%40" in href for href in hrefs) == 1
        assert aside.count('class="side-batch-count">2</span>') == 2

    def test_reselect_same_path_replaces_in_place(self, make_client,
                                                  first_dir, second_dir):
        client, state = make_client()
        _select(client, first_dir)
        _select(client, second_dir)
        save_interactions([_interaction("pi_cation")],
                          str(first_dir / "pi_cation.h5"))
        assert _select(client, first_dir).headers["location"] == "/"
        assert [source["name"] for source in state["sources"]] == \
            ["first", "second"]
        assert len(state["batch_items"][str(first_dir.resolve())]) == 3
        assert len(state["items"]) == 5
        names = [name for name, _, _ in state["items"]]
        assert "pi_cation" in names and len(names) == len(set(names))

    def test_trailing_slash_is_same_batch(self, make_client, first_dir):
        client, state = make_client()
        _select(client, first_dir)
        _select(client, f"{first_dir}/")
        assert len(state["sources"]) == 1
        assert len(state["items"]) == 2

    def test_invalid_select_keeps_batches(self, make_client, first_dir):
        client, state = make_client()
        _select(client, first_dir)
        resp = _select(client, "/definitely/not/here")
        assert resp.status_code == 303 and "error=" in resp.headers["location"]
        assert len(state["sources"]) == 1 and len(state["items"]) == 2

    def test_select_still_records_recent(self, loaded, first_dir, second_dir):
        from DuIvyInteractions.visualizers import load_recent
        recent = load_recent()
        assert recent[0] == str(second_dir)
        assert recent[1] == str(first_dir)


# ============================================================
# /select/remove 与 /select/clear
# ============================================================

class TestRemoveClear:

    def test_remove_one_batch(self, loaded, first_dir):
        client, state = loaded
        resp = client.post("/select/remove", data={"path": str(first_dir)},
                           follow_redirects=False)
        assert resp.status_code == 303 and resp.headers["location"] == "/"
        assert [source["name"] for source in state["sources"]] == ["second"]
        assert [name for name, _, _ in state["items"]] == \
            ["hydrogen_bond", "pi_stacking"]

    def test_remove_unknown_path_reports_error(self, loaded, tmp_path):
        client, state = loaded
        resp = client.post("/select/remove",
                           data={"path": str(tmp_path / "missing")},
                           follow_redirects=False)
        assert resp.status_code == 303
        location = unquote(resp.headers["location"])
        assert location.startswith("/?error=batch not loaded:")
        assert len(state["sources"]) == 2

    def test_remove_last_batch_returns_to_empty(self, make_client, first_dir):
        client, state = make_client()
        _select(client, first_dir)
        client.post("/select/remove", data={"path": str(first_dir)},
                    follow_redirects=False)
        assert state["sources"] == [] and state["items"] == []
        html = client.get("/").text
        assert 'class="card"' not in html
        main = html.split('<main class="content">', 1)[1].split("</main>", 1)[0]
        assert 'class="empty-brand"' in main
        aside = html.split('<aside class="sidebar">', 1)[1].split("</aside>", 1)[0]
        assert "+ load" in aside and "+ run" in aside
        assert 'id="loadModal"' in html

    def test_clear_all_batches(self, loaded):
        client, state = loaded
        resp = client.post("/select/clear", follow_redirects=False)
        assert resp.status_code == 303 and resp.headers["location"] == "/"
        assert state["sources"] == [] and state["items"] == []
        assert 'class="card"' not in client.get("/").text

    def test_clear_without_batches(self, make_client):
        client, state = make_client()
        resp = client.post("/select/clear", follow_redirects=False)
        assert resp.status_code == 303 and resp.headers["location"] == "/"
        assert state["sources"] == [] and state["items"] == []


# ============================================================
# 多批次下的既有路由
# ============================================================

class TestMergedRoutes:

    def test_type_page_plain_and_suffixed(self, loaded):
        client, _ = loaded
        assert client.get("/type/hydrogen_bond").status_code == 200
        resp = client.get(f"/type/{SUFFIXED}")
        assert resp.status_code == 200 and SUFFIXED in resp.text

    def test_pair_page_merged(self, loaded):
        assert loaded[0].get(f"/pair/{SUFFIXED}/0").status_code == 200

    def test_api_matrix_both_names(self, loaded):
        client, _ = loaded
        for name in ("hydrogen_bond", SUFFIXED):
            data = client.get(f"/api/matrix/{name}").json()
            assert data["n_frames"] == FRAMES and len(data["rows"]) == 2

    def test_api_cell_merged(self, loaded):
        data = loaded[0].get(f"/api/cell/{SUFFIXED}/0/0").json()
        assert data["exists"] is True
        assert data["values"]["distance"] == 3.0

    def test_api_frame_and_frames_merged(self, loaded):
        client, _ = loaded
        assert client.get(f"/api/frame/{SUFFIXED}/0").json()["n_active"] == 2
        assert len(client.get(f"/api/frames/{SUFFIXED}/0").json()["rows"]) == \
            FRAMES

    def test_download_csv_merged(self, loaded):
        resp = loaded[0].get(f"/download/csv/{SUFFIXED}")
        assert resp.status_code == 200
        assert resp.headers["content-type"].startswith("text/csv")
        lines = resp.text.strip().splitlines()
        assert lines[0].startswith("pair_label,occupancy,avg_distance")
        assert len(lines) == 3

    def test_structure_page_merged(self, loaded):
        resp = loaded[0].get(f"/structure/{SUFFIXED}/0")
        assert resp.status_code == 200 and SUFFIXED in resp.text

    def test_api_structure_interactions_without_structure(self, loaded):
        resp = loaded[0].get("/api/structure/interactions",
                             params={"type": SUFFIXED, "idx": 0})
        assert resp.status_code == 400
        assert resp.json()["detail"] == STRUCTURE_NOT_SELECTED

    def test_unknown_type_still_404(self, loaded):
        assert loaded[0].get("/api/matrix/nope").status_code == 404


# ============================================================
# 导航上下文
# ============================================================

class TestNavContext:

    def test_batches_context_for_multiple(self, loaded, first_dir,
                                          second_dir):
        _, state = loaded
        nav = web_app._nav(state, "hydrogen_bond")
        assert set(nav) == {"types", "datasource", "ds_name", "batches",
                            "active_type"}
        assert [batch["name"] for batch in nav["batches"]] == \
            ["first", "second"]
        assert [batch["path"] for batch in nav["batches"]] == \
            [str(first_dir.resolve()), str(second_dir.resolve())]
        assert [batch["n_types"] for batch in nav["batches"]] == [2, 2]
        assert nav["ds_name"] == "first +1"
        assert nav["active_type"] == "hydrogen_bond"
        assert len(nav["types"]) == 4
        assert " + " in nav["datasource"]

    def test_single_batch_ds_name(self, make_client, first_dir):
        client, state = make_client()
        _select(client, first_dir)
        nav = web_app._nav(state)
        assert nav["ds_name"] == "first"
        assert len(nav["batches"]) == 1
        assert nav["batches"][0]["n_types"] == 2
        assert nav["datasource"] == str(first_dir.resolve())

    def test_empty_nav(self, make_client):
        _, state = make_client()
        assert web_app._nav(state) == {
            "types": [], "datasource": None, "ds_name": "",
            "batches": [], "active_type": ""}

    def test_nav_source_without_batches_keeps_old_keys(self):
        nav = web_app._nav_source("/tmp/out", [], "hydrogen_bond")
        assert nav["ds_name"] == "out" and nav["batches"] == []
        assert nav["datasource"] == "/tmp/out"
        assert nav["active_type"] == "hydrogen_bond"

    def test_source_empty_state(self, make_client):
        _, state = make_client()
        assert web_app._source(state) == ("", [])
