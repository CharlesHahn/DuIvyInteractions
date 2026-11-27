# -*- coding: utf-8 -*-
"""P0 极简外壳测试：起始屏两入口、数据集分块、侧栏数据集树与导航上下文。

不依赖真实轨迹，也不 import 其它测试文件：用 tmp_path 构造最小 h5，
TestClient 走真实路由；样式断言读取 app.css 文本。
"""

import re
from pathlib import Path
from urllib.parse import quote

import numpy as np
import pytest
from fastapi.testclient import TestClient

from DuIvyInteractions.core.datas import AtomData, Group, Interaction
from DuIvyInteractions.io.h5 import save_interactions
from DuIvyInteractions.web import app as web_app

REPO_ROOT = Path(__file__).resolve().parents[2]
APP_CSS = REPO_ROOT / "DuIvyInteractions" / "web" / "static" / "app.css"
FRAMES = 4
SUFFIXED = "hydrogen_bond@hydrogen_bond.h5"
SUFFIXED_HREF = "/type/" + quote(SUFFIXED, safe="/")
FIRST_TYPES = ("hydrogen_bond", "hydrophobic")
SECOND_TYPES = ("hydrogen_bond", "pi_stacking")


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


@pytest.fixture(autouse=True)
def isolate_recent(monkeypatch, tmp_path):
    """隔离最近数据源文件，避免写入用户缓存。"""
    monkeypatch.setattr("DuIvyInteractions.visualizers.data.RECENT_PATH",
                        tmp_path / "recent.json")


@pytest.fixture
def first_dir(tmp_path):
    return _write_dir(tmp_path / "first", FIRST_TYPES)


@pytest.fixture
def second_dir(tmp_path):
    return _write_dir(tmp_path / "second", SECOND_TYPES)


@pytest.fixture
def make_client(monkeypatch):
    """创建 TestClient 并捕获 create_app 内部的会话状态。"""
    def build():
        captured = {}
        original = web_app._register_pages

        def capture(app, state, browse_mode, root):
            captured["state"] = state
            original(app, state, browse_mode, root)

        monkeypatch.setattr(web_app, "_register_pages", capture)
        return TestClient(web_app.create_app()), captured["state"]

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


def _aside(html: str) -> str:
    """侧栏片段（无侧栏时为空串）。"""
    if '<aside class="sidebar">' not in html:
        return ""
    return html.split('<aside class="sidebar">', 1)[1].split("</aside>", 1)[0]


def _blocks(html: str) -> list:
    """按数据集区块切分总览页（每块不含尾随的 add dataset 行）。"""
    return [part.split("</section>")[0]
            for part in html.split('<section class="ds-block">')[1:]]


def _block_named(html: str, name: str) -> str:
    """名称为 name 的数据集区块。"""
    matches = [block for block in _blocks(html)
               if f"<h2>{name}</h2>" in block]
    assert len(matches) == 1, "missing dataset block: %s" % name
    return matches[0]


# ============================================================
# 无数据：起始屏只有两个入口
# ============================================================

class TestEmptyLanding:

    def test_main_shows_brand_empty_state(self, make_client):
        """无数据时主区只渲染品牌空状态：纯展示，不含任何入口。"""
        client, _ = make_client()
        html = client.get("/").text
        for marker in ('class="start-choice"', "start-screen", "add-panel",
                       "+ add dataset", "ds-add", 'class="card"', "ds-block"):
            assert marker not in html, marker
        main = html.split('<main class="content">', 1)[1].split("</main>", 1)[0]
        assert main.count("<section") == 1
        brand = main.split('<section class="empty-brand">', 1)[1] \
                    .split("</section>", 1)[0]
        assert 'class="empty-brand-mark">DII</span>' in brand
        assert "<h1>DuIvyInteractions</h1>" in brand
        assert "non-covalent interaction analysis for MD trajectories" in brand
        assert brand.count('class="kbd"') == 2
        assert ">+ load</span>" in brand and ">+ run</span>" in brand
        assert not re.search(r"<a\b|<form\b|<button\b", brand)

    def test_no_empty_brand_with_data(self, loaded):
        assert "empty-brand" not in loaded[0].get("/").text

    def test_load_modal_is_rendered_hidden(self, make_client):
        client, _ = make_client()
        html = client.get("/").text
        assert 'class="load-modal hidden" id="loadModal"' in html
        assert 'class="source-form" id="loadForm" method="post" action="/select"' \
            in html

    def test_sidebar_is_always_present(self, make_client):
        client, _ = make_client()
        aside = _aside(client.get("/").text)
        assert "+ load" in aside and "+ run" in aside
        assert "data-load-open" in aside
        assert 'class="side-batch"' not in aside
        assert 'action="/select/clear"' not in aside

    def test_add_query_no_longer_renders_inline_form(self, make_client):
        client, _ = make_client()
        resp = client.get("/?add=1")
        assert resp.status_code == 200
        assert "add-panel" not in resp.text
        assert "<details" not in resp.text
        assert 'class="load-modal hidden" id="loadModal"' in resp.text


# ============================================================
# 顶栏：只留品牌
# ============================================================

class TestTopbarMinimal:

    def test_no_chips_and_no_add_link(self, loaded):
        html = loaded[0].get("/").text
        assert "ds-chip" not in html
        assert "add data source" not in html
        assert "topbar-spacer" not in html

    def test_brand_only(self, loaded):
        html = loaded[0].get("/").text
        assert html.count('class="topbar"') == 1
        assert '<span class="brand-mark">DII</span>' in html
        assert "<span>DuIvyInteractions</span>" in html


# ============================================================
# 有数据：按数据集分块
# ============================================================

class TestDatasetBlocks:

    def test_one_block_per_batch(self, loaded):
        html = loaded[0].get("/").text
        assert html.count('<section class="ds-block">') == 2
        assert html.count("ds-block-path") == 2
        assert "<h2>first</h2>" in html and "<h2>second</h2>" in html

    def test_cards_stay_inside_their_batch(self, loaded):
        html = loaded[0].get("/").text
        assert html.count('class="card"') == 4
        first = _block_named(html, "first")
        second = _block_named(html, "second")
        assert 'class="card-label">hydrogen_bond<' in first
        assert 'class="card-label">hydrophobic<' in first
        assert 'class="card-label">pi_stacking<' in second
        assert 'class="card-label">hydrophobic<' not in second

    def test_card_links_are_unique_per_batch(self, loaded):
        html = loaded[0].get("/").text
        assert 'href="/type/hydrogen_bond"' in _block_named(html, "first")
        assert SUFFIXED_HREF in _block_named(html, "second")

    def test_no_add_dataset_entry(self, loaded):
        html = loaded[0].get("/").text
        assert "+ add dataset" not in html
        assert "ds-add" not in html

    def test_add_query_keeps_single_load_form(self, loaded):
        html = loaded[0].get("/?add=1").text
        assert "add-panel" not in html
        assert '<details' not in html
        assert html.count('class="source-form"') == 1

    def test_remove_drops_block(self, loaded, first_dir):
        client, _ = loaded
        resp = client.post("/select/remove", data={"path": str(first_dir)})
        html = resp.text
        assert html.count('<section class="ds-block">') == 1
        assert "<h2>first</h2>" not in html
        assert "<h2>second</h2>" in html


# ============================================================
# 侧栏：数据集 → 类型 两级树
# ============================================================

class TestSidebarTree:

    def test_side_actions(self, loaded):
        aside = _aside(loaded[0].get("/").text)
        assert "+ load" in aside and "+ run" in aside
        assert 'href="/?add=1"' in aside and 'href="/run"' in aside

    def test_batch_rows_offer_remove(self, loaded):
        aside = _aside(loaded[0].get("/").text)
        assert aside.count('class="side-batch"') == 2
        assert aside.count('class="side-batch-remove"') == 2
        assert aside.count('action="/select/remove"') == 2
        assert 'action="/select/clear"' in aside

    def test_labels_show_base_links_are_unique(self, loaded):
        aside = _aside(loaded[0].get("/").text)
        assert aside.count('<span class="side-label">hydrogen_bond</span>') == 2
        assert SUFFIXED not in aside
        assert 'href="/type/hydrogen_bond"' in aside
        assert f'href="{SUFFIXED_HREF}"' in aside

    def test_sidebar_tree_on_type_page(self, loaded):
        html = loaded[0].get("/type/hydrogen_bond").text
        aside = _aside(html)
        assert aside.count('class="side-batch"') == 2
        assert aside.count('class="side-item active"') == 1
        assert "+ load" in aside and "+ run" in aside


# ============================================================
# 导航 / 总览上下文
# ============================================================

class TestNavAndGroups:

    def test_batch_summaries_expose_types(self, loaded):
        _, state = loaded
        batches = web_app._batch_summaries(state)
        assert [b["n_types"] for b in batches] == [2, 2]
        first = batches[0]["types"]
        assert [t["base"] for t in first] == list(FIRST_TYPES)
        assert [t["name"] for t in first] == list(FIRST_TYPES)
        assert [t["n_pairs"] for t in first] == [2, 2]
        second = batches[1]["types"]
        assert [t["base"] for t in second] == [SUFFIXED.split("@")[0],
                                              "pi_stacking"]
        assert [t["name"] for t in second] == [SUFFIXED, "pi_stacking"]

    def test_groups_split_by_batch(self, loaded):
        _, state = loaded
        groups = web_app._batch_groups(state)
        assert [g["name"] for g in groups] == ["first", "second"]
        assert [g["n_types"] for g in groups] == [2, 2]
        assert [c["data_type"] for c in groups[0]["cards"]] == list(FIRST_TYPES)
        assert [c["type"] for c in groups[1]["cards"]] == [SUFFIXED,
                                                          "pi_stacking"]

    def test_index_context_keeps_merged_cards(self, monkeypatch, loaded):
        client, _ = loaded
        captured = {}
        original = web_app.TEMPLATES.TemplateResponse

        def capture(request, name, context, **kwargs):
            captured.update(context)
            return original(request, name, context, **kwargs)

        monkeypatch.setattr(web_app.TEMPLATES, "TemplateResponse", capture)
        assert client.get("/").status_code == 200
        assert len(captured["cards"]) == 4
        assert len(captured["groups"]) == 2

    def test_index_context_has_shell_extras(self, monkeypatch, loaded):
        """主页上下文保留外壳所需的浏览方式与最近数据源（供 base 弹窗）。"""
        client, _ = loaded
        captured = {}
        original = web_app.TEMPLATES.TemplateResponse

        def capture(request, name, context, **kwargs):
            captured.update(context)
            return original(request, name, context, **kwargs)

        monkeypatch.setattr(web_app.TEMPLATES, "TemplateResponse", capture)
        assert client.get("/").status_code == 200
        assert captured["browse_mode"] == "picker"
        assert isinstance(captured["recent"], list)


# ============================================================
# 样式（仅末尾追加）
# ============================================================

@pytest.fixture(scope="module")
def css():
    return APP_CSS.read_text(encoding="utf-8")


class TestShellCss:

    @staticmethod
    def _rules(css, selector):
        rules = re.findall(re.escape(selector) + r"\s*\{([^}]*)\}", css)
        assert rules, "missing rule for %s" % selector
        return rules

    def test_load_modal_rules(self, css):
        modal = self._rules(css, ".load-modal")[0]
        assert "position: fixed" in modal
        assert "z-index: 90" in modal
        assert any("max-width" in rule
                   for rule in self._rules(css, ".load-panel"))
        assert "display: flex" in self._rules(css, ".load-row")[0]
        assert self._rules(css, ".load-recent")
        assert ".ds-add" not in css

    def test_dataset_and_sidebar_rules(self, css):
        assert any("text-overflow: ellipsis" in rule
                   for rule in self._rules(css, ".ds-block-path"))
        assert "margin-bottom" in self._rules(css, ".ds-block")[-1]
        assert "padding-left" in self._rules(css, ".side-batch .side-item")[-1]
        for selector in (".side-actions", ".side-action", ".side-batch",
                         ".side-batch-head", ".side-batch-remove",
                         ".side-clear"):
            assert self._rules(css, selector), selector
