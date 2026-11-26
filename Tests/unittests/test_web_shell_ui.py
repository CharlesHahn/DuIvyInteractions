# -*- coding: utf-8 -*-
"""Web 外壳 UI 测试：顶栏品牌、侧栏数据集树、宽屏自适应与 3D viewer 响应式高度。

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


def _type(name: str, n_pairs: int = 1) -> dict:
    """侧栏类型条目：name 为全局唯一链接名，base 为显示的基础类型名。"""
    return {"name": name, "base": name.split("@")[0], "n_pairs": n_pairs}


def _batch(name: str, type_names) -> dict:
    """批次摘要：路径、名称、类型条数与类型明细（形状同 app._batch_summaries）。"""
    path = "/data/%s" % name
    types = [_type(type_name) for type_name in type_names]
    return {"path": path, "name": name, "n_types": len(types), "types": types}


def _render_sidebar(env, batches) -> str:
    """渲染带侧栏的外壳页：types 非空才渲染侧栏（形状同 type_summary）。"""
    types = [{"type": t["name"], "data_type": t["base"], "n_pairs": t["n_pairs"]}
             for batch in batches for t in batch["types"]]
    return _render(env, batches=batches, types=types)


def _aside(html: str) -> str:
    """侧栏片段（无侧栏时为空串）。"""
    if '<aside class="sidebar">' not in html:
        return ""
    return html.split('<aside class="sidebar">', 1)[1].split("</aside>", 1)[0]


def _side_batch(aside: str, name: str) -> str:
    """侧栏中名为 name 的数据集块（split 边界恰为下一块起点）。"""
    for chunk in aside.split('<div class="side-batch">')[1:]:
        head = '<span class="side-batch-name" title="/data/%s">%s</span>' % (name, name)
        if head in chunk:
            return chunk
    raise AssertionError("missing sidebar batch: %s" % name)


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
# 侧栏数据集树
# ============================================================


class TestSidebarTree:
    def test_three_batches_render_three_datasets(self, env):
        batches = [_batch("web_run_test", ["hydrogen_bond", "hydrophobic"]),
                   _batch("md_amber", ["salt_bridge"]),
                   _batch("md_charmm", ["pi_stacking", "halogen_bond"])]
        html = _render_sidebar(env, batches)
        assert html.count('class="side-batch"') == 3
        for batch in batches:
            assert 'name="path" value="%s"' % batch["path"] in html
            assert '<span class="side-batch-name" title="%s">%s</span>' % (
                batch["path"], batch["name"]) in html
            assert '<span class="side-batch-count">%d</span>' % \
                batch["n_types"] in html
        assert html.count('action="/select/remove"') == 3
        assert html.count('class="side-batch-remove"') == 3

    def test_sidebar_lists_batch_types(self, env):
        """side-label 显示基础名（无 @），href 用全局唯一名。"""
        batches = [_batch("first", ["hydrogen_bond",
                                    "hydrogen_bond@hydrogen_bond.h5"])]
        first = _side_batch(_aside(_render_sidebar(env, batches)), "first")
        label = '<span class="side-label">hydrogen_bond</span>'
        assert first.count('class="side-item"') == 2
        assert first.count(label) == 2
        assert "hydrogen_bond@" not in first
        assert 'href="/type/hydrogen_bond"' in first
        assert 'href="/type/hydrogen_bond%40hydrogen_bond.h5"' in first

    def test_multiple_batches_render_clear_all_once(self, env):
        html = _render_sidebar(env, [_batch("first", ["hydrogen_bond"]),
                                     _batch("second", ["salt_bridge"])])
        assert html.count('action="/select/clear"') == 1
        assert '<form method="post" action="/select/clear" class="side-clear">' \
            in html
        assert "clear all" in html

    def test_single_batch_renders_clear_all(self, env):
        """新规则：只要有数据集就提供清空（旧规则单批次隐藏 clear all）。"""
        html = _render_sidebar(env, [_batch("only_one", ["hydrogen_bond"])])
        assert html.count('class="side-batch"') == 1
        assert 'action="/select/clear"' in html

    def test_empty_or_missing_batches_hides_sidebar(self, env):
        """无批次（空列表 / 不传）时不渲染侧栏，也不再有单 chip 回退。"""
        empty = _render(env, batches=[])
        missing = _render(env)
        for html in (empty, missing):
            assert 'class="sidebar"' not in html
            assert "side-batch" not in html
            assert "ds-chip" not in html
            assert "ds-chips" not in html


# ============================================================
# 历史 chips 样式（样式表只追加不删，保留回归）
# ============================================================


class TestLegacyChipCss:
    def test_chips_scroller_degrades_to_horizontal_scroll(self, css):
        """历史 .ds-chips 规则仍在：单行横向滚动（不换行），按钮不被挤压。"""
        base = [rule for rule in _rules(css, ".ds-chips")
                if "display: flex" in rule]
        assert base, "missing flex .ds-chips rule"
        assert "overflow-x: auto" in base[0]
        assert "max-width: 46vw" in base[0]
