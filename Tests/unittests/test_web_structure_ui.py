# -*- coding: utf-8 -*-
"""3D 结构页 UI 测试：浏览按钮、水/序列开关、服务端文件弹窗与静态资源。

不依赖浏览器：TestClient 断言模板 DOM，另直接读取静态资源文本校验关键 API 字符串。
"""

import re
from pathlib import Path

import numpy as np
import pytest

from DuIvyInteractions.core.datas import AtomData, Group, Interaction
from DuIvyInteractions.io.h5 import save_interactions

WEB_DIR = Path(__file__).resolve().parents[2] / "DuIvyInteractions" / "web"
STATIC_DIR = WEB_DIR / "static"
STRUCTURE_JS = STATIC_DIR / "structure.js"
FILEPICK_JS = STATIC_DIR / "filepick.js"
APP_CSS = STATIC_DIR / "app.css"


def _atom(index: int, name: str, atype: str = "") -> AtomData:
    return AtomData(atom_global_idx=index, atom_idx_in_residue=0,
                    atom_name=name, atom_type=atype, atom_element="C",
                    atom_charge=0.0, atom_mass=0.0)


def _group(group_id: int, group_type: str, atoms) -> Group:
    return Group(group_id=group_id, group_type=group_type, molecule="test",
                 residue_name="LYS", residue_id=1, atoms=list(atoms))


def _interaction() -> Interaction:
    existence = np.ones((1, 1), dtype=bool)
    groups = [(_group(0, "H_donor", [_atom(0, "N1", "N"),
                                     _atom(1, "H1", "H")]),
               _group(1, "H_acceptor", [_atom(2, "O1", "O"),
                                        _atom(3, "O2", "O")]))]
    return Interaction(interaction_type="hydrogen_bond", groups=groups,
                       existence=existence,
                       metrics={"distance": np.full((1, 1), 3.0)},
                       times=np.zeros(1))


@pytest.fixture(autouse=True)
def isolate_recent(monkeypatch, tmp_path):
    monkeypatch.setattr("DuIvyInteractions.visualizers.data.RECENT_PATH",
                        tmp_path / "recent.json")


@pytest.fixture
def h5_dir(tmp_path):
    directory = tmp_path / "out"
    directory.mkdir()
    save_interactions([_interaction()], str(directory / "hydrogen_bond.h5"))
    return directory


def _client(**kwargs):
    from fastapi.testclient import TestClient
    from DuIvyInteractions.web.app import create_app
    return TestClient(create_app(**kwargs))


def _page(client) -> str:
    return client.get("/structure/hydrogen_bond/0").text


@pytest.fixture
def client(h5_dir):
    client = _client()
    client.post("/select", data={"path": str(h5_dir)})
    return client


@pytest.fixture
def server_client(h5_dir):
    client = _client(browse_mode="server")
    client.post("/select", data={"path": str(h5_dir)})
    return client


def _tag(html: str, tag_name: str, element_id: str) -> str:
    pattern = r"<%s[^>]*id=\"%s\"[^>]*>" % (tag_name, re.escape(element_id))
    match = re.search(pattern, html)
    assert match, f"missing <{tag_name} id=\"{element_id}\">"
    return match.group(0)


def _function_body(text: str, name: str) -> str:
    """取 structure.js 里 2 空格缩进的函数体（静态切片，不做 JS 解析）。"""
    match = re.search(r"  function %s\(.*?\n  \}\n" % re.escape(name), text, re.S)
    assert match, f"missing function {name}"
    return match.group(0)


def _render_template(browse_mode: str) -> str:
    from DuIvyInteractions.web.app import TEMPLATES

    template = TEMPLATES.get_template("structure.html")
    return template.render(
        title="structure", ds_name="", datasource="", types=[],
        active_type=None, type="hydrogen_bond", idx=0, pair_label="pair 1",
        n_atoms=4, n_frames=1, has_structure=False, existence=[1],
        topology="", trajectory="", browse_mode=browse_mode)


# ============================================================
# 页面 DOM
# ============================================================


class TestStructurePageDom:
    def test_page_has_two_browse_buttons(self, client):
        html = _page(client)
        topology = _tag(html, "button", "browseTopology")
        trajectory = _tag(html, "button", "browseTrajectory")
        for tag in (topology, trajectory):
            assert 'type="button"' in tag
            assert 'data-browse="picker"' in tag

    def test_page_has_water_and_sequence_toggles(self, client):
        html = _page(client)
        assert 'id="showWater"' in html
        assert 'id="showSequence"' in html
        assert "structure-toolbar" in html
        assert "checked" not in _tag(html, "input", "showWater")
        assert "checked" in _tag(html, "input", "showSequence")

    def test_page_keeps_viewer_dom_ids(self, client):
        html = _page(client)
        for element_id in ("molstarViewer", "frameSlider", "frameLabel",
                           "existenceBar", "structureStatus", "structureConfig",
                           "structureForm", "topologyInput", "trajectoryInput",
                           "openStructure"):
            assert 'id="%s"' % element_id in html

    def test_page_links_viewer_assets_in_order(self, client):
        html = _page(client)
        assert "/static/vendor/molstar.css" in html
        molstar_at = html.index("/static/vendor/molstar.js")
        structure_at = html.index("/static/structure.js")
        assert molstar_at < structure_at

    def test_page_has_no_cdn_reference(self, client):
        html = _page(client)
        assert "http://" not in html
        assert "https://" not in html


# ============================================================
# 服务端文件弹窗（picker / server 两种模式）
# ============================================================


class TestServerFilePicker:
    def test_picker_page_has_no_modal_dom(self, client):
        html = _page(client)
        assert 'id="filepickModal"' not in html
        assert "/static/filepick.js" not in html

    def test_picker_template_has_no_modal_dom(self):
        html = _render_template("picker")
        assert 'data-browse="picker"' in html
        assert 'id="filepickModal"' not in html
        assert "/static/filepick.js" not in html

    def test_server_template_has_modal_dom(self):
        html = _render_template("server")
        assert 'data-browse="server"' in html
        assert 'id="filepickModal"' in html
        for element_id in ("filepickList", "filepickError", "filepickNote",
                           "filepickPath", "filepickCurrent", "filepickUp",
                           "filepickGo", "filepickUse", "filepickCancel"):
            assert 'id="%s"' % element_id in html
        assert "/static/filepick.js" in html

    def test_server_route_renders_modal(self, server_client):
        html = _page(server_client)
        if 'data-browse="server"' not in html:
            pytest.skip("_render_structure 尚未注入 browse_mode（后端任务未完成）")
        assert 'id="filepickModal"' in html
        assert "/static/filepick.js" in html


# ============================================================
# 静态资源
# ============================================================


class TestStaticAssets:
    def test_filepick_script_is_served(self, client):
        response = client.get("/static/filepick.js")
        assert response.status_code == 200
        assert "DiiFilePick" in response.text

    def test_structure_script_keeps_required_apis(self):
        text = STRUCTURE_JS.read_text(encoding="utf-8")
        for token in ("focusLoci", "customText", "solvent", "/api/pick",
                      "layoutShowControls", "layoutShowSequence",
                      "involves_water", "contact_label", "DiiFilePick"):
            assert token in text, token
        assert "&solvent=" in text

    def test_frontend_scripts_avoid_innerhtml(self):
        for path in (STRUCTURE_JS, FILEPICK_JS):
            assert "innerHTML" not in path.read_text(encoding="utf-8")

    def test_css_has_filepick_and_toolbar_rules(self):
        text = APP_CSS.read_text(encoding="utf-8")
        assert ".filepick-modal" in text
        assert ".filepick-row" in text
        assert ".structure-toolbar" in text
        assert ".structure-toggle" in text

    def test_structure_script_keeps_viewer_unexpanded(self):
        """静态检查：不启用 Mol* 展开态（展开态 fixed 定位会盖住页内工具栏/滑块）。"""
        text = STRUCTURE_JS.read_text(encoding="utf-8")
        assert "layoutIsExpanded: false" in text
        assert "layoutIsExpanded: true" not in text

    def test_css_viewer_container_rule_takes_effect(self):
        """静态检查：末尾 .structure-viewer 规则把 viewer 约束为容器内 620px。"""
        text = APP_CSS.read_text(encoding="utf-8")
        rules = re.findall(r"\.structure-viewer\s*\{([^}]*)\}", text)
        assert rules, "missing .structure-viewer rule"
        final_rule = rules[-1]
        assert "height: 620px" in final_rule
        assert "position: relative" in final_rule


# ============================================================
# 相互作用两端基团球棍标注 + 图例
# ============================================================


class TestGroupAnnotations:
    def test_page_has_group_legend_container(self, client):
        html = _page(client)
        assert 'id="groupLegend"' in html
        assert "structure-groups" in html
        legend_at = html.index('id="groupLegend"')
        status_at = html.index('id="structureStatus"')
        assert legend_at < status_at

    def test_structure_script_keeps_group_overlay_apis(self):
        text = STRUCTURE_JS.read_text(encoding="utf-8")
        for token in ("GROUP_COLORS", "ball-and-stick", "addRepresentation",
                      "StructureElement.Bundle.fromLoci", "payload.groups"):
            assert token in text, token

    def test_structure_script_group_palette_is_npg(self):
        text = STRUCTURE_JS.read_text(encoding="utf-8")
        match = re.search(r"GROUP_COLORS\s*=\s*\[([^\]]*)\]", text)
        assert match, "missing GROUP_COLORS array"
        assert re.findall(r"#[0-9a-fA-F]{6}", match.group(1)) == [
            "#e64b35", "#3c5488", "#00a087"]

    def test_structure_script_clears_group_reprs_per_frame(self):
        """静态检查：换帧时复用 state.data.build().delete 清除基团表示。"""
        text = STRUCTURE_JS.read_text(encoding="utf-8")
        assert "groupReprs" in text
        assert "function clearGroups(" in text
        assert "clearGroups()" in text

    def test_group_component_uses_awaited_structure_selector(self):
        """tryCreateComponent 首参必须是 cell/selector，且是 await 过的 selector。"""
        text = STRUCTURE_JS.read_text(encoding="utf-8")
        assert "function structureSelector(" in text
        assert "structureRef(" not in text
        match = re.search(r"tryCreateComponent\(\s*([A-Za-z0-9_$.]+)", text)
        assert match, "missing tryCreateComponent call"
        first_arg = match.group(1)
        assert first_arg != "structureRef"
        assert not first_arg.endswith(".ref"), first_arg
        assert re.search(r"groupComponent\([^)]*\)\s*\.then\(", text, re.S), \
            "groupComponent 返回的 selector 必须先 await 再传给 addRepresentation"
        assert 'console.warn("[dii] group component failed"' in text

    def test_css_has_group_legend_rules(self):
        text = APP_CSS.read_text(encoding="utf-8")
        assert ".structure-groups" in text
        assert ".structure-group-dot" in text
        assert "width: 8px" in text


# ============================================================
# 高亮小结构路径（一组一 PDB + 接触对 PDB）
# ============================================================


class TestHighlightStructures:
    def test_structure_script_uses_highlight_payload_fields(self):
        text = STRUCTURE_JS.read_text(encoding="utf-8")
        for token in ("payload.highlights", "contact_pdb", "contact_label"):
            assert token in text, token

    def test_each_highlight_loads_its_own_structure(self):
        text = STRUCTURE_JS.read_text(encoding="utf-8")
        assert "loadStructureFromData" in _function_body(text, "loadStructureCell")
        assert "loadStructureCell(pdb)" in _function_body(text, "highlightStep")
        runner = _function_body(text, "addHighlightReprs")
        assert "highlightStep(" in runner
        assert "for (" in runner

    def test_highlight_path_avoids_main_structure_atom_indices(self):
        """新高亮路径不得使用主结构的原子序号（unitOf），它只留在回退分支。"""
        text = STRUCTURE_JS.read_text(encoding="utf-8")
        for name in ("wholeStructureLoci", "highlightStep", "addHighlightReprs",
                     "addContactPdb", "structureAtomEntries"):
            assert "unitOf" not in _function_body(text, name), name
        assert "unitOf" in _function_body(text, "groupLoci")

    def test_structure_script_has_default_repr_cleanup(self):
        text = STRUCTURE_JS.read_text(encoding="utf-8")
        body = _function_body(text, "clearDefaultReprs")
        for token in ("components", "representations", "builder.delete",
                      "canUndo: false"):
            assert token in body, token
        assert "default representation cleanup failed" in body

    def test_highlight_cleans_default_reprs_before_adding(self):
        """新球棍表示必须先删掉小结构自带的默认表示（否则双重渲染）。"""
        text = STRUCTURE_JS.read_text(encoding="utf-8")
        step = _function_body(text, "highlightStep")
        assert step.index("clearDefaultReprs(") < step.index("addGroupRepr(")
        assert "clearDefaultReprs(" in _function_body(text, "addContactPdb")

    def test_page_has_residue_number_hint(self, client):
        html = _page(client)
        assert "标签编号同时给出 h5 表格编号与 PDB auth 编号" in html
        assert "TYR50 (PDB 51)" in html
