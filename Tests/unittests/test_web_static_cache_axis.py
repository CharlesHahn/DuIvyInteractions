# -*- coding: utf-8 -*-
"""Web 静态资源缓存头 + 类型页柱状图 Y 轴。

缓存：/static/* 必须每次使用前向服务器校验（no-cache），避免浏览器启发式缓存
旧 JS 导致类型页矩阵空白。
Y 轴：active pairs per frame 与 occupancy distribution 两图左侧标注柱高计数
（顶 = 最大值，中 = 一半向下取整，底 = 0）；轴标题（active pairs / pairs）
竖排写在 Y 轴旁，不写进 panel 标题；柱子按比例填满绘图区宽度。
"""

import re
from pathlib import Path

import numpy as np
import pytest

from DuIvyInteractions.core.datas import AtomData, Group, Interaction
from DuIvyInteractions.io.h5 import save_interactions
from DuIvyInteractions.visualizers import frame_activity, interaction_meta

FRAMES = 12
TYPE = "hydrogen_bond"
N_PAIRS = 6
REPO_ROOT = Path(__file__).resolve().parents[2]
APP_CSS = REPO_ROOT / "DuIvyInteractions" / "web" / "static" / "app.css"
STATIC_ASSETS = ("matrix.js", "app.css", "structure.js")
CACHE_CONTROL = "no-cache"
YAXIS_RE = re.compile(r'<div class="bars-yaxis"[^>]*>(.*?)</div>', re.S)
YAXIS_TITLE_RE = re.compile(
    r'<div class="bars-yaxis-title"[^>]*>(.*?)</div>', re.S)
SPAN_RE = re.compile(r"<span>([^<]*)</span>")
INTEGER_RE = re.compile(r"^\d+$")
PLOT_MARKUP = 'class="bars-plot"'
YAXIS_MARKUP = 'class="bars-yaxis"'
YAXIS_TITLE_MARKUP = 'class="bars-yaxis-title"'


@pytest.fixture(autouse=True)
def isolate_recent(monkeypatch, tmp_path):
    monkeypatch.setattr("DuIvyInteractions.visualizers.data.RECENT_PATH",
                        tmp_path / "recent.json")


def _atom(gidx: int, name: str, element: str) -> AtomData:
    return AtomData(atom_global_idx=gidx, atom_idx_in_residue=0, atom_name=name,
                    atom_type="", atom_element=element, atom_charge=0.0,
                    atom_mass=0.0)


def _group(gid: int, gtype: str, gidx: int, res: str, rid: int,
           with_h: bool = False) -> Group:
    atoms = [_atom(gidx, "C", "C")]
    if with_h:
        atoms.append(_atom(gidx + 1, "H", "H"))
    return Group(group_id=gid, group_type=gtype, molecule="test",
                 residue_name=res, residue_id=rid, atoms=atoms)


def _existence() -> np.ndarray:
    """每帧恒 5 个活跃对；占用率 4+1+1 分布 → 最大 bin 4 对。"""
    existence = np.ones((N_PAIRS, FRAMES), dtype=bool)
    existence[0, 3] = False       # 第 0 对缺一帧
    existence[N_PAIRS - 1, :] = False
    existence[N_PAIRS - 1, 3] = True  # 末对只在该帧出现，补满 5 对
    return existence


def _interaction() -> Interaction:
    existence = _existence()
    metrics = {"distance": np.where(existence, 3.0, np.nan)}
    groups = [(_group(2 * i, "H_donor", 10 * i, "LYS", i, with_h=True),
               _group(2 * i + 1, "H_acceptor", 10 * i + 1, "ASP", 100 + i))
              for i in range(N_PAIRS)]
    return Interaction(interaction_type=TYPE, groups=groups,
                       existence=existence, metrics=metrics,
                       times=np.arange(FRAMES, dtype=float) * 10.0)


@pytest.fixture
def interaction() -> Interaction:
    return _interaction()


@pytest.fixture
def h5_dir(tmp_path, interaction):
    directory = tmp_path / "out"
    directory.mkdir()
    save_interactions([interaction], str(directory / (TYPE + ".h5")))
    return directory


@pytest.fixture
def client(h5_dir):
    from fastapi.testclient import TestClient
    from DuIvyInteractions.web.app import create_app
    client = TestClient(create_app())
    client.post("/select", data={"path": str(h5_dir)})
    return client


def _type_page(client):
    return client.get("/type/" + TYPE)


def _yaxis_blocks(html: str):
    """两处 bars-yaxis 的标签列表（顶/中/底）。"""
    return [SPAN_RE.findall(block) for block in YAXIS_RE.findall(html)]


def _yaxis_titles(html: str):
    """两处 bars-yaxis-title 的文本（active pairs / pairs）。"""
    return [text.strip() for text in YAXIS_TITLE_RE.findall(html)]


# ============================================================
# 缓存头：/static/* 每次校验；其它响应不因此改变
# ============================================================

class TestStaticCacheControl:

    @pytest.mark.parametrize("asset", STATIC_ASSETS)
    def test_static_asset_requires_revalidation(self, client, asset):
        response = client.get("/static/" + asset)
        assert response.status_code == 200
        assert response.headers["cache-control"] == CACHE_CONTROL

    def test_static_asset_keeps_etag(self, client):
        response = client.get("/static/app.css")
        assert response.headers["etag"]

    def test_conditional_request_returns_304(self, client):
        etag = client.get("/static/app.css").headers["etag"]
        response = client.get("/static/app.css",
                              headers={"If-None-Match": etag})
        assert response.status_code == 304
        assert response.headers["cache-control"] == CACHE_CONTROL

    def test_home_page_not_marked_no_cache(self, client):
        response = client.get("/")
        assert response.status_code == 200
        assert response.headers.get("cache-control") is None

    def test_type_page_still_reachable(self, client):
        response = _type_page(client)
        assert response.status_code == 200
        assert response.headers.get("cache-control") is None

    def test_api_matrix_not_marked(self, client):
        response = client.get("/api/matrix/" + TYPE)
        assert response.status_code == 200
        assert response.headers.get("cache-control") is None

    def test_unknown_static_path_still_404(self, client):
        assert client.get("/static/missing-asset.js").status_code == 404


# ============================================================
# 模板：两图各一处 Y 轴，数值为整数，X 轴标签保留
# ============================================================

class TestBarsYAxis:

    def test_two_plots_and_two_yaxes(self, client):
        html = _type_page(client).text
        assert html.count(PLOT_MARKUP) == 2
        assert html.count(YAXIS_MARKUP) == 2
        assert html.count(YAXIS_TITLE_MARKUP) == 2

    def test_yaxis_titles_sit_beside_axes(self, client):
        assert _yaxis_titles(_type_page(client).text) == ["active pairs",
                                                          "pairs"]

    def test_activity_yaxis_values(self, client, interaction):
        blocks = _yaxis_blocks(_type_page(client).text)
        top = max(frame_activity(interaction))
        assert top == 5
        assert blocks[0] == [str(top), str(top // 2), "0"]

    def test_occupancy_yaxis_values(self, client, interaction, h5_dir):
        blocks = _yaxis_blocks(_type_page(client).text)
        meta = interaction_meta(interaction, str(h5_dir / (TYPE + ".h5")))
        top = max(meta["occupancy"]["hist"])
        assert top == 4
        assert blocks[1] == [str(top), str(top // 2), "0"]

    def test_yaxis_labels_are_integers(self, client):
        for labels in _yaxis_blocks(_type_page(client).text):
            assert len(labels) == 3
            for label in labels:
                assert INTEGER_RE.match(label), label

    def test_plot_holds_original_bars(self, client):
        html = _type_page(client).text
        assert '<div class="bars big">' in html
        assert '<div class="bars big bars-occupancy">' in html

    def test_activity_x_labels_preserved(self, client):
        html = _type_page(client).text
        assert "<span>frame 0</span>" in html
        assert "<span>frame %d</span>" % (FRAMES - 1) in html

    def test_occupancy_x_labels_preserved(self, client):
        html = _type_page(client).text
        assert ("<span>0%</span><span>50%</span><span>100%</span>"
                in html)

    def test_titles_do_not_inline_axis_unit(self, client):
        html = _type_page(client).text
        assert "active pairs per frame" in html
        assert "occupancy distribution · %d pairs" % N_PAIRS in html
        assert "(y:" not in html


# ============================================================
# CSS：末尾追加 Y 轴样式，且不波及无 Y 轴的迷你柱状图
# ============================================================

class TestAxisCss:

    @pytest.fixture
    def css(self) -> str:
        return APP_CSS.read_text(encoding="utf-8")

    def test_bars_plot_rule(self, css):
        assert ".bars-plot {" in css
        assert "grid-template-columns: 16px 42px 1fr;" in css

    def test_bars_yaxis_title_rule(self, css):
        assert ".bars-yaxis-title {" in css
        assert "writing-mode: vertical-rl" in css

    def test_bars_yaxis_rule(self, css):
        assert ".bars-yaxis {" in css

    def test_bars_axis_padding_rule(self, css):
        assert ".bars-axis { padding-left: 70px; }" in css

    def test_big_bars_fill_plot_width(self, css):
        rules = re.findall(r"\.bars\.big span\s*\{([^}]*)\}", css)
        assert rules, "missing .bars.big span rule"
        final_rule = rules[-1]
        assert "flex: 1 1 0" in final_rule
        assert "max-width: none" in final_rule

    def test_mini_axis_alignment_restored(self, css):
        assert ".bars.mini + .bars-axis," in css
        assert ".hist-cell .bars-axis { padding-left: 0; }" in css

    def test_existing_bars_sizing_untouched(self, css):
        assert ".bars.big { height: 64px; gap: 3px; }" in css
        assert ".bars-axis {\n  display: flex;" in css
