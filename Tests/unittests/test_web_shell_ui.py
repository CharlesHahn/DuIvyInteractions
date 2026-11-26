# -*- coding: utf-8 -*-
"""Web 外壳 UI 测试：顶栏品牌、多批次 chips、宽屏自适应与 3D viewer 响应式高度。

不依赖 app.py：用独立 jinja2 环境直接渲染 base.html（该模板的 block 渲染为空，
可独立渲染）；样式断言直接读取 app.css 文本，取同名选择器的最后一条规则
（CSS 层叠下后写的同名规则生效）。
"""

import re
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
WEB_DIR = REPO_ROOT / "DuIvyInteractions" / "web"
TEMPLATE_DIR = WEB_DIR / "templates"
APP_CSS = WEB_DIR / "static" / "app.css"

BRAND_TEXT = "DuIvyInteractions"
CLASSIC_BLUE = "#38A7D0"
CLASSIC_PINK = "#F67088"


@pytest.fixture(scope="module")
def env():
    from jinja2 import Environment, FileSystemLoader, select_autoescape

    return Environment(loader=FileSystemLoader(str(TEMPLATE_DIR)),
                       autoescape=select_autoescape(["html", "xml"]))


@pytest.fixture(scope="module")
def css():
    return APP_CSS.read_text(encoding="utf-8")


def _render(env, title="overview", **context):
    return env.get_template("base.html").render(title=title, **context)


def _rules(css: str, selector: str):
    return re.findall(re.escape(selector) + r"\s*\{([^}]*)\}", css)


def _final_rule(css: str, selector: str) -> str:
    rules = _rules(css, selector)
    assert rules, "missing rule for %s" % selector
    return rules[-1]


def _batch(name: str, n_types: int) -> dict:
    path = "/data/%s" % name
    return {"path": path, "name": name, "n_types": n_types}


# ============================================================
# 品牌
# ============================================================


class TestBrand:
    def test_brand_text_is_single_word(self, env):
        html = _render(env)
        assert "<span>%s</span>" % BRAND_TEXT in html
        assert "DuIvy Interactions" not in html

    def test_brand_mark_text_is_dii(self, env):
        html = _render(env)
        assert '<span class="brand-mark">DII</span>' in html

    def test_page_title_keeps_dii_prefix(self, env):
        assert "<title>DII · overview</title>" in _render(env)
        assert "<title>DII · structure</title>" in _render(env, title="structure")

    def test_brand_mark_uses_duivy_classic_palette(self, css):
        rule = _final_rule(css, ".brand-mark")
        assert "linear-gradient(135deg, %s, %s)" % (CLASSIC_BLUE, CLASSIC_PINK) in rule
        assert "color: #fff" in rule
        assert "text-shadow" in rule


# ============================================================
# 宽屏自适应
# ============================================================


class TestContentWidth:
    def test_content_final_rule_removes_width_cap(self, css):
        """判定方式：末尾追加的 .content 规则含 max-width: none（CSS 后者生效）。

        历史规则里的 max-width: 1280px 保留不改，故不用“全文不再出现 1280px”判定。
        """
        rule = _final_rule(css, ".content")
        assert "max-width: none" in rule

    def test_wide_table_keeps_readability_cap(self, css):
        rule = _final_rule(css, ".content .table-wrap")
        assert "max-width: 1600px" in rule

    def test_matrix_and_viewer_are_not_width_capped(self, css):
        """核心要求：矩阵画布与 3D viewer 用满宽度，不套可读性上限。"""
        for selector in (".matrix-wrap", ".matrix-panel", ".structure-viewer"):
            for rule in _rules(css, selector):
                assert "max-width" not in rule, selector


# ============================================================
# 3D viewer 响应式高度
# ============================================================


class TestViewerHeight:
    def test_viewer_final_rule_is_responsive(self, css):
        """620px 是 min() 的旧浏览器回退值，随后的 min(72vh, 820px) 生效。"""
        rule = _final_rule(css, ".structure-viewer")
        assert "width: 100%" in rule
        assert "height: 620px" in rule
        assert "height: min(72vh, 820px)" in rule
        assert "position: relative" in rule

    def test_viewer_inner_div_fills_container(self, css):
        rule = _final_rule(css, ".structure-viewer > div")
        assert "width: 100%" in rule
        assert "height: 100%" in rule


# ============================================================
# 顶栏批次 chips
# ============================================================


class TestBatchChips:
    def test_three_batches_render_three_chips(self, env):
        batches = [_batch("web_run_test", 8), _batch("md_amber", 3),
                   _batch("md_charmm", 5)]
        html = _render(env, batches=batches)
        assert html.count('class="ds-chip"') == 3
        for b in batches:
            assert 'name="path" value="%s"' % b["path"] in html
            assert 'title="%s">%s</span>' % (b["path"], b["name"]) in html
            assert '<span class="ds-count">%d</span>' % b["n_types"] in html
        assert html.count('action="/select/remove"') == 3
        assert html.count('class="ds-remove"') == 3

    def test_single_batch_has_no_clear_all(self, env):
        html = _render(env, batches=[_batch("only_one", 2)])
        assert html.count('class="ds-chip"') == 1
        assert 'action="/select/clear"' not in html

    def test_multiple_batches_render_clear_all_once(self, env):
        html = _render(env, batches=[_batch("first", 1), _batch("second", 2)])
        assert html.count('action="/select/clear"') == 1
        assert '<form class="ds-clear" method="post" action="/select/clear">' in html
        assert "clear all" in html

    def test_empty_batches_falls_back_to_single_chip(self, env):
        html = _render(env, batches=[], ds_name="legacy_run",
                       datasource="/data/legacy_run")
        assert html.count('class="ds-chip"') == 1
        assert "legacy_run" in html
        assert 'action="/select/remove"' not in html
        assert 'action="/select/clear"' not in html
        assert 'class="ds-chips"' not in html

    def test_missing_batches_falls_back_to_single_chip(self, env):
        """其它页面/旧后端不传 batches，仍按 ds_name/datasource 渲染单 chip。"""
        html = _render(env, ds_name="old_ds", datasource="/data/old_ds")
        assert html.count('class="ds-chip"') == 1
        assert 'title="/data/old_ds"><span class="dot"></span>' in html

    def test_no_chip_without_any_datasource(self, env):
        html = _render(env)
        assert "ds-chip" not in html
        assert "ds-chips" not in html

    def test_chips_scroller_degrades_to_horizontal_scroll(self, css):
        """多批次布局方案：单行横向滚动（不换行），右侧按钮不被挤压。"""
        base = [rule for rule in _rules(css, ".ds-chips")
                if "display: flex" in rule]
        assert base, "missing flex .ds-chips rule"
        assert "overflow-x: auto" in base[0]
        assert "max-width: 46vw" in base[0]
