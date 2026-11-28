# -*- coding: utf-8 -*-
"""应用内路由（pjax）与常驻 3D 坞：模板/脚本静态契约 + 服务端零改动回归。

无浏览器环境：不验证 WebGL canvas 搬移与脚本实际重放，只做可静态核对的契约断言。
"""

import re
from pathlib import Path

import pytest

H5_DATA_DIR = Path(__file__).resolve().parents[1] / "interaction_h5data"
WEB_DIR = Path(__file__).resolve().parents[2] / "DuIvyInteractions" / "web"
TEMPLATES = WEB_DIR / "templates"
ROUTER_JS = WEB_DIR / "static" / "router.js"
STRUCTURE_JS = WEB_DIR / "static" / "structure.js"
PAIR_JS = WEB_DIR / "static" / "pair.js"
BASE_HTML = TEMPLATES / "base.html"
PARTIAL = TEMPLATES / "_structure_view.html"
SCRIPT_SRC = re.compile(r'<script src="([^"]+)"')
PAGES = ("/", "/type/hydrogen_bond", "/pair/pi_stacking/0", "/residues",
         "/run")


@pytest.fixture(autouse=True)
def isolate_recent(monkeypatch, tmp_path):
    monkeypatch.setattr("DuIvyInteractions.visualizers.data.RECENT_PATH",
                        tmp_path / "recent.json")


@pytest.fixture
def client():
    from fastapi.testclient import TestClient
    from DuIvyInteractions.web.app import create_app
    client = TestClient(create_app())
    client.post("/select", data={"path": str(H5_DATA_DIR)})
    return client


def _router_js() -> str:
    return ROUTER_JS.read_text(encoding="utf-8")


def _structure_js() -> str:
    return STRUCTURE_JS.read_text(encoding="utf-8")


def _function_body(text: str, name: str) -> str:
    match = re.search(r"  function %s\(.*?\n  \}\n" % re.escape(name), text,
                      re.S)
    assert match, f"missing function {name}"
    return match.group(0)


# ============================================================
# 常驻坞：全站唯一，位于 main 之外
# ============================================================


class TestViewerDock:
    def test_base_declares_single_dock(self):
        html = BASE_HTML.read_text(encoding="utf-8")
        for element_id in ("viewerDock", "molstarViewer", "existenceBar",
                           "frameLabel", "groupLegend", "structureStatus"):
            assert 'id="%s"' % element_id in html, element_id
        assert 'class="viewer-dock"' in html and " hidden" in html
        assert 'id="pageScripts"' in html and "data-page-scripts" in html
        assert '<script src="/static/router.js"></script>' in html

    @pytest.mark.parametrize("url", PAGES)
    def test_every_page_has_exactly_one_viewer(self, client, url):
        html = client.get(url).text
        assert html.count('id="molstarViewer"') == 1
        assert html.count('id="viewerDock"') == 1
        assert html.count('id="existenceBar"') == 1

    def test_dock_outside_main_slot_inside_main(self, client):
        html = client.get("/pair/pi_stacking/0").text
        main_start = html.index('<main class="content">')
        main_end = html.index("</main>")
        slot = html.index('id="structureDockSlot"')
        dock = html.index('id="viewerDock"')
        assert main_start < slot < main_end
        assert dock > main_end

    def test_partial_has_slot_without_viewer(self):
        text = PARTIAL.read_text(encoding="utf-8")
        assert 'id="structureDockSlot"' in text
        assert 'id="molstarViewer"' not in text
        assert 'id="existenceBar"' not in text

    def test_structure_config_carries_structure_key(self, client):
        html = client.get("/pair/pi_stacking/0").text
        config = re.search(
            r'<script id="structureConfig"[^>]*>(.*?)</script>', html, re.S)
        assert config and '"topology"' in config.group(1)


# ============================================================
# 无 JS 兜底与脚本装配
# ============================================================


class TestScriptWiring:
    def test_navigation_links_stay_plain_anchors(self, client):
        html = client.get("/pair/pi_stacking/1").text
        assert re.search(r'<a[^>]*href="/pair/pi_stacking/0[^"]*"', html)
        assert re.search(r'<a[^>]*href="/pair/pi_stacking/2[^"]*"', html)
        type_page = client.get("/type/hydrogen_bond").text
        assert 'class="side-item' in type_page
        assert '<a class="side-item' in type_page

    def test_page_scripts_wrapper_outside_main(self, client):
        html = client.get("/pair/pi_stacking/0").text
        wrapper = html.index('id="pageScripts"')
        assert wrapper > html.index("</main>")

    def test_heavy_bundles_marked_data_once(self, client):
        html = client.get("/pair/pi_stacking/0").text
        scripts = SCRIPT_SRC.findall(html)
        page = html[html.index('id="pageScripts"'):]
        for src in ("/static/vendor/molstar.js", "/static/structure.js",
                    "/static/plotly-basic.min.js"):
            assert src in scripts
            assert re.search(r'<script src="%s" data-once>' % re.escape(src),
                             page), src
        assert 'src="/static/pair.js" data-once' not in page


# ============================================================
# router.js 契约
# ============================================================


class TestRouterScript:
    def test_exports_public_api(self):
        text = _router_js()
        for token in ("window.DiiRouter = {", "navigate: navigate",
                      "currentUrl: currentUrl", "onAfterSwap: onAfterSwap"):
            assert token in text, token

    def test_uses_history_and_dom_parsing(self):
        text = _router_js()
        for token in ("history.pushState", 'DOMParser',
                      "replaceChildren", "popstate", "X-DII-PJAX",
                      "onAfterSwap", "data-once"):
            assert token in text, token
        assert "innerHTML" not in text

    def test_skips_api_download_static_and_modified_clicks(self):
        text = _router_js()
        for token in ('"/api/"', '"/download/"', '"/static/"',
                      "defaultPrevented", "metaKey", "target",
                      "hasAttribute(\"download\")"):
            assert token in text, token

    def test_script_replay_dedupes_once_bundles(self):
        text = _router_js()
        assert "loadedOnce" in text
        assert "hasScript(src)" in text
        assert "script.async = false;" in text

    def test_nested_navigation_guard_and_fallback(self):
        text = _router_js()
        assert "var pending = false;" in text
        assert "if (pending) return" in text
        assert "location.assign(url)" in text


# ============================================================
# structure.js 坞生命周期
# ============================================================


class TestStructureDockLifecycle:
    def test_exports_attach_detach_and_is_mounted(self):
        text = _structure_js()
        for token in ("attach: attach", "detach: detach",
                      "isMounted: function () { return !!viewer; }",
                      "mount: mount", "updatePair: updatePair"):
            assert token in text, token

    def test_singleton_guard_prevents_second_viewer(self):
        text = _structure_js()
        assert ("if (window.DiiStructure && "
                "typeof window.DiiStructure.mount === \"function\")") in text
        assert "if (mounted) return;" in _function_body(text, "mount")

    def test_attach_moves_dock_and_compares_structure_key(self):
        body = _function_body(_structure_js(), "attach")
        assert "structureDockSlot" in body
        assert "sameStructure(config)" in body
        assert "reloadStructure()" in body
        helper = _function_body(_structure_js(), "sameStructure")
        assert "currentStructure.topology" in helper
        assert "currentStructure.trajectory" in helper

    def test_detach_keeps_viewer_alive(self):
        text = _structure_js()
        body = _function_body(text, "detach")
        assert "dock.hidden = true;" in body
        assert "appendChild(dock)" in body
        assert "viewer.dispose" not in text
        assert "dispose(" not in body

    def test_after_swap_hook_refreshes_pair_data(self):
        text = _structure_js()
        hook = _function_body(text, "registerSwapHook")
        assert "DiiRouter.onAfterSwap" in hook
        assert "attach(next)" in hook
        refresh = _function_body(text, "refreshPairData")
        assert "/api/pair/" in refresh
        assert "updatePair(payload)" in refresh




# ============================================================
# 交换前钩子：先摘坞再换 main
# ============================================================


class TestBeforeSwap:
    def test_router_exports_and_orders_before_swap_hook(self):
        text = _router_js()
        assert "onBeforeSwap: onBeforeSwap" in text
        body = _function_body(text, "navigate")
        assert "runBeforeSwap(url);" in body
        assert "swapMain(doc);" in body
        assert body.index("runBeforeSwap(url);") < body.index("swapMain(doc);")
        assert text.index("runBeforeSwap(url);") < text.index(
            "target.replaceChildren")

    def test_structure_registers_before_swap_detach(self):
        text = _structure_js()
        hook = _function_body(text, "registerSwapHook")
        assert "onBeforeSwap(detachFromMain)" in hook
        guard = _function_body(text, "detachFromMain")
        assert "dock.parentNode !== document.body" in guard
        assert "detach()" in guard
        detach = _function_body(text, "detach")
        assert "dock.hidden = true;" in detach
        assert "appendChild(dock)" in detach
        assert "viewer.dispose" not in detach

    def test_router_syncs_sidebar_active(self):
        text = _router_js()
        body = _function_body(text, "syncSidebar")
        assert '".side-item"' in body
        assert 'classList.add("active")' in body
        assert "syncSidebar(url);" in _function_body(text, "navigate")



# ============================================================
# 首屏脚本登记 + attach 诊断（R6）
# ============================================================


class TestPjaxDiagnostics:
    def test_router_seeds_initial_once_scripts(self):
        text = _router_js()
        assert "function seedLoadedScripts" in text
        body = _function_body(text, "seedLoadedScripts")
        assert "PAGE_SCRIPTS_ID" in body
        assert "script[src]" in body
        assert "ONCE_ATTR" in body
        assert "loadedOnce[src] = true;" in body
        assert "hasScript" not in body

    def test_router_calls_seed_at_iife_end(self):
        text = _router_js()
        definition = text.index("function seedLoadedScripts")
        call = text.rindex("seedLoadedScripts();")
        popstate = text.index('window.addEventListener("popstate", onPopState);')
        assert definition < call
        assert call > popstate
        assert text.rstrip().endswith("})();")
        assert "DOMContentLoaded" in text

    def test_attach_logs_structured_decision(self):
        body = _function_body(_structure_js(), "attach")
        assert 'console.warn("[dii] dock attach", {' in body
        for field in ("same: same", "ready: ready", "hasViewer: hasViewer",
                      "reloaded: reloaded", "from: currentStructure",
                      "to: structureKey(config)"):
            assert field in body, field

    def test_same_structure_empty_fingerprint_guards(self):
        body = _function_body(_structure_js(), "sameStructure")
        assert "if (!next.topology && !next.trajectory) return true;" in body
        assert "if (!currentStructure.topology && !currentStructure.trajectory)" in body
        assert "return !!viewer && ready;" in body

    def test_singleton_guard_still_first(self):
        text = _structure_js()
        head = text[:text.index("var cfgEl =")]
        assert "if (window.DiiStructure &&" in head
        assert "typeof window.DiiStructure.mount === \"function\"" in head

# ============================================================
# 服务端零改动回归
# ============================================================


class TestServerUnchanged:
    def test_pair_and_type_pages_still_server_rendered(self, client):
        pair = client.get("/pair/pi_stacking/0")
        assert pair.status_code == 200
        for token in ("occupancy", 'id="pairChart"', 'id="pairFigure"',
                      'id="framesTable"', 'id="pairSummaryBody"'):
            assert token in pair.text, token
        type_page = client.get("/type/hydrogen_bond")
        assert type_page.status_code == 200
        assert 'id="matrixCanvas"' in type_page.text

    def test_pair_script_binds_head_container(self):
        text = PAIR_JS.read_text(encoding="utf-8")
        body = _function_body(text, "bindPairNav")
        assert 'document.getElementById("pairHead")' in body
        assert "document.addEventListener" not in body
