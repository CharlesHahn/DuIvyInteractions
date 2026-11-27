# -*- coding: utf-8 -*-
"""全局 load 弹窗测试：base 常驻弹窗、picker/server 分支、recent 与前端资源。

不依赖真实轨迹：用 tmp_path 构造最小 h5；recent 缓存用 fixture 隔离。
"""

import json
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest
from fastapi.testclient import TestClient

from DuIvyInteractions.core.datas import AtomData, Group, Interaction
from DuIvyInteractions.io.h5 import save_interactions
from DuIvyInteractions.web.app import create_app

REPO_ROOT = Path(__file__).resolve().parents[2]
WEB_DIR = REPO_ROOT / "DuIvyInteractions" / "web"
LOAD_JS = WEB_DIR / "static" / "load.js"
APP_CSS = WEB_DIR / "static" / "app.css"
BASE_HTML = WEB_DIR / "templates" / "base.html"
FRAMES = 4
RECENT_PATHS = ["/data/alpha", "/data/beta"]
HINT = "Browse selects a folder"
SERVER_HINT = "Server browsing is limited to the start-up directory"


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


@pytest.fixture(autouse=True)
def recent_path(monkeypatch, tmp_path):
    """隔离最近数据源文件并返回其路径（默认不存在 = 无 recent）。"""
    path = tmp_path / "recent.json"
    monkeypatch.setattr("DuIvyInteractions.visualizers.data.RECENT_PATH", path)
    return path


@pytest.fixture
def data_dir(tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    save_interactions([_interaction("hydrogen_bond")],
                      str(out / "hydrogen_bond.h5"))
    return out


@pytest.fixture
def picker():
    return TestClient(create_app())


@pytest.fixture
def server():
    return TestClient(create_app(browse_mode="server"))


@pytest.fixture
def loaded(picker, data_dir):
    assert picker.post("/select", data={"path": str(data_dir)},
                       follow_redirects=False).status_code == 303
    return picker


def _aside(html: str) -> str:
    if '<aside class="sidebar">' not in html:
        return ""
    return html.split('<aside class="sidebar">', 1)[1].split("</aside>", 1)[0]


def _main(html: str) -> str:
    return html.split('<main class="content">', 1)[1].split("</main>", 1)[0]


def _load_modal(html: str) -> str:
    """load 弹窗片段（缺失时为空串）。"""
    if '<div class="load-modal hidden" id="loadModal">' not in html:
        return ""
    return html.split('<div class="load-modal hidden" id="loadModal">', 1)[1] \
               .split('<div id="toast"', 1)[0]


def _node_ellipsis(cases):
    """用 Node 执行 load.js 里真实的 middleEllipsis（不复制实现）。"""
    node = shutil.which("node")
    if node is None:
        pytest.skip("node unavailable")
    script = """
const fs = require("fs");
const src = fs.readFileSync(process.argv[1], "utf8");
function block(name) {
  var at = src.indexOf("function " + name);
  if (at < 0) throw new Error("missing " + name);
  var depth = 0, i = src.indexOf("{", at);
  for (; i < src.length; i++) {
    if (src[i] === "{") depth++;
    else if (src[i] === "}") { depth--; if (depth === 0) return src.slice(at, i + 1); }
  }
  throw new Error("unterminated " + name);
}
eval(src.match(/var ELLIPSIS = [^;]+;/)[0]);
eval(src.match(/var HEAD_RATIO = [^;]+;/)[0]);
eval(block("middleEllipsis"));
var cases = JSON.parse(process.argv[2]);
console.log(JSON.stringify(cases.map(function (c) {
  return middleEllipsis(c[0], c[1]);
})));
"""
    result = subprocess.run([node, "-e", script, str(LOAD_JS),
                             json.dumps(cases)],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    return json.loads(result.stdout)


# ============================================================
# 所有页面都带全局弹窗
# ============================================================

class TestModalOnEveryPage:

    def test_home_has_global_load_modal(self, picker):
        html = picker.get("/").text
        assert 'id="loadModal"' in html
        assert "data-load-open" in _aside(html)
        assert 'class="source-form" id="loadForm" method="post" action="/select"' \
            in html

    def test_type_page_has_global_load_modal(self, loaded):
        html = loaded.get("/type/hydrogen_bond").text
        assert 'id="loadModal"' in html
        assert "data-load-open" in _aside(html)

    def test_pair_page_has_global_load_modal(self, loaded):
        html = loaded.get("/pair/hydrogen_bond/0").text
        assert 'id="loadModal"' in html
        assert "data-load-open" in _aside(html)

    def test_run_page_has_global_load_modal(self, picker):
        html = picker.get("/run").text
        assert 'id="loadModal"' in html
        assert "data-load-open" in _aside(html)

    def test_error_page_has_global_load_modal(self, picker):
        html = picker.get("/type/nope").text
        assert 'id="loadModal"' in html


# ============================================================
# 弹窗字段与文案
# ============================================================

class TestModalContent:

    def test_modal_starts_hidden(self, picker):
        assert 'class="load-modal hidden" id="loadModal"' in picker.get("/").text

    def test_fields_and_no_submit_button(self, picker):
        """弹窗只有输入框 + browse + ✕：提交手势是 Enter，不是按钮。"""
        modal = _load_modal(picker.get("/").text)
        assert 'id="loadPath"' in modal and 'name="path"' in modal
        assert 'id="loadForm"' in modal and 'action="/select"' in modal
        assert " required" in modal
        assert 'id="loadBrowse"' in modal and ">browse…</button>" in modal
        assert 'id="loadClose"' in modal
        assert ">load results</span>" in modal
        assert "type or paste a path, then press Enter" in modal
        assert ">Open<" not in modal
        assert "btn-primary" not in modal
        assert 'id="loadSubmit"' not in modal
        assert modal.count("<button") == 2

    def test_hint_matches_old_form(self, picker):
        html = picker.get("/").text
        assert HINT in html
        assert SERVER_HINT not in html

    def test_recent_links_rendered(self, recent_path, picker):
        recent_path.write_text(json.dumps(RECENT_PATHS))
        html = picker.get("/").text
        assert 'class="load-recent"' in html
        assert 'class="load-recent-title">recent</span>' in html
        assert html.count('<li><a href="#"') == len(RECENT_PATHS)
        for path in RECENT_PATHS:
            assert 'data-load-recent="%s"' % path in html
            assert 'title="%s"' % path in html
            assert ">%s</a>" % path in html

    def test_recent_limits_to_five(self, recent_path, picker):
        paths = ["/data/dataset_%d" % i for i in range(8)]
        recent_path.write_text(json.dumps(paths))
        html = picker.get("/").text
        block = html.split('<div class="load-recent">', 1)[1] \
                    .split("</div>", 1)[0]
        assert block.count("<li>") == 5
        assert "/data/dataset_0" in block
        assert "/data/dataset_5" not in block

    def test_no_recent_block_when_empty(self, picker):
        html = picker.get("/").text
        assert "load-recent" not in html


# ============================================================
# browse_mode 分支
# ============================================================

class TestBrowseModes:

    def test_picker_mode_has_no_browse_overlay(self, picker):
        html = picker.get("/").text
        assert 'data-browse-mode="picker"' in html
        assert 'id="browseModal"' not in html
        assert "/static/browse.js" not in html
        assert 'id="browseBtn"' not in html

    def test_server_mode_renders_browse_overlay(self, server):
        html = server.get("/").text
        assert 'data-browse-mode="server"' in html
        assert 'id="browseModal"' in html
        assert "/static/browse.js" in html
        assert 'id="browseBtn"' in html
        assert SERVER_HINT in html

    def test_load_js_is_served(self, picker):
        resp = picker.get("/static/load.js")
        assert resp.status_code == 200
        assert "loadModal" in resp.text


# ============================================================
# load.js 行为要点（静态断言）
# ============================================================

@pytest.fixture(scope="module")
def load_source():
    return LOAD_JS.read_text(encoding="utf-8")


class TestLoadJsSource:

    def test_is_es5_iife(self, load_source):
        assert "use strict" in load_source
        assert "=>" not in load_source
        assert "let " not in load_source and "const " not in load_source

    def test_open_and_close_targets(self, load_source):
        assert "loadModal" in load_source and "loadPath" in load_source
        assert "loadClose" in load_source and "focus()" in load_source
        assert "Escape" in load_source and "KEY_ESCAPE_CODE" in load_source

    def test_openers_and_recent(self, load_source):
        assert "data-load-open" in load_source
        assert "data-load-recent" in load_source
        assert "URLSearchParams" in load_source and "add=1" in load_source

    def test_browse_dispatch(self, load_source):
        assert "/api/pick" in load_source
        assert "mode=dir&kind=dataset" in load_source
        assert "#browseBtn" in load_source
        assert "server" in load_source

    def test_enter_submits_form(self, load_source):
        """Enter 走 requestSubmit（触发 required 校验），无按钮提交。"""
        assert "requestSubmit" in load_source
        assert 'FORM_ID = "loadForm"' in load_source
        assert "Enter" in load_source and "keydown" in load_source
        assert "submitForm" in load_source

    def test_recent_and_native_browse_submit(self, load_source):
        """recent 点击与原生浏览选完即提交，不再只填值。"""
        recent = load_source.split("function bindRecent", 1)[1] \
                            .split("function onPathKey", 1)[0]
        assert "submitForm()" in recent
        native = load_source.split("function lookupNative", 1)[1] \
                            .split("function browse", 1)[0]
        assert "submitForm()" in native

    def test_recent_middle_ellipsis(self, load_source):
        assert "middleEllipsis" in load_source
        assert "applyRecentLabels" in load_source
        assert 'getAttribute("data-load-recent")' in load_source
        assert 'addEventListener("resize"' in load_source
        assert "RESIZE_DEBOUNCE_MS" in load_source


# ============================================================
# 主页与 CSS
# ============================================================

class TestIndexAndCss:

    def test_main_shows_brand_empty_state(self, picker):
        """无数据时主区只渲染品牌空状态（纯展示，无入口）。"""
        html = picker.get("/").text
        main = _main(html)
        assert 'class="empty-brand"' in main
        assert "empty-brand-mark" in main and 'class="kbd"' in main
        assert "start-choice" not in html
        assert "add-panel" not in html
        assert "ds-add" not in html
        brand = main.split('<section class="empty-brand">', 1)[1] \
                    .split("</section>", 1)[0]
        assert "href=" not in brand and "<button" not in brand
        assert "<form" not in brand

    def test_empty_brand_css(self):
        css = APP_CSS.read_text(encoding="utf-8")
        brand = re.findall(r"\.empty-brand\s*\{([^}]*)\}", css)[-1]
        assert "min-height: calc(100vh - 230px)" in brand
        assert "flex-direction: column" in brand
        mark = re.findall(r"\.empty-brand-mark\s*\{([^}]*)\}", css)[-1]
        assert "linear-gradient(135deg, #38A7D0, #F67088)" in mark
        assert "width: 84px" in mark and "font-size: 26px" in mark

    def test_css_recent_is_vertical_and_truncated(self):
        css = APP_CSS.read_text(encoding="utf-8")
        ul = re.findall(r"\.load-recent ul\s*\{([^}]*)\}", css)[-1]
        assert "flex-direction: column" in ul
        links = re.findall(r"\.load-recent a\s*\{([^}]*)\}", css)
        assert any("text-overflow: ellipsis" in rule for rule in links)
        assert any("white-space: nowrap" in rule for rule in links)
        panel = re.findall(r"\.load-panel\s*\{([^}]*)\}", css)
        assert any("min-width: 0" in rule or "overflow: hidden" in rule
                   for rule in panel)

    def test_banner_close_css(self):
        css = APP_CSS.read_text(encoding="utf-8")
        assert re.search(r"\.banner-close\s*\{", css)
        banner = re.findall(r"\.banner\s*\{([^}]*)\}", css)[-1]
        assert "position: relative" in banner

    def test_index_has_no_add_dataset_entry(self, loaded):
        html = loaded.get("/").text
        assert "+ add dataset" not in html
        assert "<details" not in html
        assert '<section class="ds-block">' in html

    def test_css_has_load_modal_rules(self):
        css = APP_CSS.read_text(encoding="utf-8")
        for selector in (".load-modal", ".load-panel", ".load-head",
                         ".load-title", ".load-row", ".load-recent"):
            assert re.search(re.escape(selector) + r"\s*\{", css), selector
        assert "z-index: 90" in re.findall(r"\.load-modal\s*\{([^}]*)\}",
                                          css)[0]
        assert ".ds-add" not in css


# ============================================================
# recent 路径中间省略（Node 执行真实函数）
# ============================================================

class TestMiddleEllipsis:

    def test_short_path_unchanged(self):
        text = "/data/short.h5"
        assert _node_ellipsis([[text, 40]]) == [text]

    def test_long_path_keeps_both_ends(self):
        text = "/data/project_alpha/results/hydrogen_bond.h5"
        out = _node_ellipsis([[text, 20]])[0]
        assert out.count("…") == 1
        assert len(out) <= 20
        assert out[0] == text[0] and out[-1] == text[-1]
        assert 0 < out.index("…") < len(out) - 1

    def test_tiny_budget_has_head_and_tail(self):
        text = "/data/very/long/path/to/file.h5"
        out = _node_ellipsis([[text, 3]])[0]
        assert len(out) <= 3
        assert out[0] == text[0] and out[-1] == text[-1]
        assert "…" in out


# ============================================================
# 提示栏可关闭
# ============================================================

class TestBanners:

    def test_error_and_warning_banners_are_dismissible(self, picker):
        html = picker.get("/", params={"error": "boom",
                                       "warn": "2 file(s) skipped"}).text
        assert html.count('class="banner banner-dismissible"') == 2
        assert html.count('class="banner-close"') == 2
        assert 'aria-label="dismiss"' in html
        assert "2 file(s) skipped" in html

    def test_base_inline_script_hides_banner(self):
        source = BASE_HTML.read_text(encoding="utf-8")
        assert "querySelectorAll('.banner-close')" in source
        assert "classList.add('hidden')" in source
        assert "parentNode" in source
