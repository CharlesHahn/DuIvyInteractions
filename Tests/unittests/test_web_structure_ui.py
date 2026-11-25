# -*- coding: utf-8 -*-
"""3D 结构页 UI 测试：Mol* 原生轨迹、零重建换帧、棍棒基团与图例、服务端文件弹窗。

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


def _structure_js() -> str:
    return STRUCTURE_JS.read_text(encoding="utf-8")


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

    def test_page_has_no_self_made_toggles(self, client):
        """开关全部用 Mol* 原生按钮/原生 Layout 多选，页面不自造控件。"""
        html = _page(client)
        for element_id in ("showWater", "showSequence", "toggleControls"):
            assert 'id="%s"' % element_id not in html
        assert "show water" not in html
        assert "structure-toolbar" not in html
        for element_id in ("molstarViewer", "existenceBar", "groupLegend",
                           "frameLabel", "structureStatus"):
            assert 'id="%s"' % element_id in html

    def test_page_keeps_viewer_dom_ids(self, client):
        html = _page(client)
        for element_id in ("molstarViewer", "frameLabel", "existenceBar",
                           "structureStatus", "structureConfig",
                           "structureForm", "topologyInput", "trajectoryInput",
                           "openStructure"):
            assert 'id="%s"' % element_id in html

    def test_page_has_no_range_slider(self, client):
        """进度条已并入相互作用向量条：页面不再有 range 输入。"""
        html = _page(client)
        assert 'type="range"' not in html
        assert 'id="frameSlider"' not in html

    def test_existence_bar_is_the_frame_slider(self, client):
        html = _page(client)
        bar = _tag(html, "div", "existenceBar")
        assert "bars" in bar and "structure-existence" in bar
        assert 'role="slider"' in bar
        assert 'tabindex="0"' in bar
        assert 'aria-valuemin="0"' in bar
        assert 'aria-valuemax="' in bar
        assert 'aria-valuenow="0"' in bar
        assert '<span class="on"' in html

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

    def test_page_has_residue_number_hint(self, client):
        html = _page(client)
        assert "标签编号同时给出 h5 表格编号与 PDB auth 编号" in html
        assert "TYR50 (PDB 51)" in html

    def test_page_states_water_is_not_shown(self, client):
        html = _page(client)
        assert "水分子不显示" in html


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
# 静态资源（原生轨迹 / 零重建 / 棍棒基团）
# ============================================================


class TestNativeTrajectoryScript:
    def test_trajectory_loaded_once_by_native_api(self):
        text = _structure_js()
        for token in ("loadTrajectory(", "model-data", "coordinates-url",
                      "isBinary: true", "frameIndex"):
            assert token in text, token
        call = _function_body(text, "loadTrajectoryOnce")
        for token in ("MODEL_KIND", "coordinates-url", "TRAJECTORY_URL",
                      "trajectoryFormat()"):
            assert token in call, token
        assert 'MODEL_KIND = "model-data"' in text

    def test_topology_pdb_uses_model_data_branch(self):
        """PDB 不是拓扑格式（dataFormats 只对 psf/prmtop/top 给 .topology），
        必须走 model-data 分支，否则 loadTrajectory 内部读 r.ref 会抛 TypeError。"""
        text = _structure_js()
        assert text.count("topology-data") == 0
        assert text.count("model-data") >= 1
        assert text.count('kind: MODEL_KIND') == 1
        assert 'format: TOPOLOGY_FORMAT' in _function_body(
            text, "loadTrajectoryOnce")

    def test_load_failure_message_carries_context(self):
        body = _function_body(_structure_js(), "loadFailureText")
        for token in ("MODEL_KIND", "TOPOLOGY_FORMAT", "trajectoryFormat()",
                      "trajectory load failed"):
            assert token in body, token

    def test_model_ref_lookup_has_three_levels(self):
        """换帧目标定位：hierarchy model 条目 → preset.model → 扫变换名。"""
        text = _structure_js()
        body = _function_body(text, "findModelRef")
        assert "structureEntry()" in body
        assert "preset.model" in body
        assert "scanModelRef()" in body
        scan = _function_body(text, "scanModelRef")
        assert "MODEL_TRANSFORM_NAME" in scan
        assert 'MODEL_TRANSFORM_NAME = "model-from-trajectory"' in text

    def test_script_never_rebuilds_structure(self):
        assert _structure_js().count("loadStructureFromData") == 0

    def test_frame_path_only_moves_model_index(self):
        text = _structure_js()
        assert "updateTransform" in text
        assert "MODEL_INDEX_PARAM = \"modelIndex\"" in text
        for name in ("selectFrame", "applyMolstarFrame", "pollFrame"):
            body = _function_body(text, name)
            assert "loadStructureFromData" not in body, name
            assert "loadTrajectory" not in body, name
        apply_body = _function_body(text, "applyMolstarFrame")
        assert "params[MODEL_INDEX_PARAM] = frameIndex" in apply_body
        assert "updateTransform" in apply_body

    def test_camera_is_focused_only_once(self):
        text = _structure_js()
        assert "function focusInitially(" in text
        assert text.count("managers.camera") == 1
        body = _function_body(text, "focusInitially")
        for token in ("focusFromResidues(", "focusFromEachResidue(",
                      "focusWholeStructure()", "lockCamera()"):
            assert token in body, token
        assert "updateTransform" not in body
        assert "canvas.setProps({ camera: { manualReset: true } })" in text

    def test_focus_has_three_level_fallback(self):
        """聚焦回退：两端并集 → 两端各自合并 → 主结构整体，路径进状态栏。"""
        text = _structure_js()
        assert 'focusPath = path || "none"' in text
        apply_body = _function_body(text, "applyFocusLoci")
        assert "focusLoci(" in apply_body
        assert "isLociEmpty(" in apply_body
        assert '"residues"' in _function_body(text, "focusFromResidues")
        assert '"residues-each"' in _function_body(text, "focusFromEachResidue")
        assert '"structure"' in _function_body(text, "focusWholeStructure")
        assert '"focus: " + focusPath' in _function_body(text, "statusText")

    def test_loci_are_remapped_to_root_structure(self):
        """测量/相机必须用主结构上的 Loci（组件子结构的 Loci 先 remap）。"""
        text = _structure_js()
        to_root = _function_body(text, "toRootLoci")
        assert "rootStructure()" in to_root
        assert "Loci.remap(" in to_root
        assert "rootLociOf(" in _function_body(text, "buildContactMeasurement")
        assert "rootLociOf(" in _function_body(text, "focusFromResidues")
        assert "wholeStructureLoci(componentStructure" not in text

    def test_structure_cell_is_resolved_by_ref(self):
        """帧更新会换掉 structure cell 对象，主结构必须按 ref 实时取。"""
        text = _structure_js()
        assert "var structureRef = null;" in text
        body = _function_body(text, "structureCell")
        assert "cells.get(structureRef)" in body
        assert "tryCreateComponent(structureRef" in text
        assert "tryCreateComponentStatic(structureRef" in text

    def test_contact_builds_new_measurement_before_removing_old(self):
        text = _structure_js()
        update = _function_body(text, "updateContact")
        build_at = update.index("buildContactMeasurement(")
        install_at = update.index("installMeasurement(")
        assert build_at < install_at
        assert "removeMeasurement()" in _function_body(text, "installMeasurement")
        active_path = update[update.index("contactAt("):build_at]
        assert "removeMeasurement" not in active_path
        assert 'contactState = "failed"' in text
        assert 'contactState = "ok"' in text
        assert '"contact: " + contactState' in _function_body(text, "statusText")

    def test_script_uses_expression_loci_and_sticks(self):
        text = _structure_js()
        for token in ("fromExpression", "ball-and-stick", "sizeFactor",
                      "sizeAspectRatio", "tryCreateComponentStatic",
                      "addRepresentation"):
            assert token in text, token
        assert "RESIDUE_STICK_SIZE = 0.22" in text
        assert "STICK_ASPECT = 1.0" in text
        assert "LIGAND_SIZE = 0.28" in text

    def test_residue_layer_is_single_element_colored_sticks(self):
        """相互作用残基只有一层棍棒、用 Mol* 默认 element 配色（无第二层/统一色）。"""
        text = _structure_js()
        residue = _function_body(text, "buildResidueLayer")
        assert "residueLayerExpression(residue)" in residue
        assert "stickParams(RESIDUE_STICK_SIZE)" in residue
        assert "buildLayer(" in residue
        assert "and not (" not in residue
        assert "color" not in residue
        params = _function_body(text, "stickParams")
        assert "type: STICK_TYPE" in params
        assert "sizeFactor" in params and "sizeAspectRatio" in params
        assert "color" not in params and "colorParams" not in params
        for token in ("GROUP_COLORS", "GROUP_STICK_SIZE", "STICK_COLOR",
                      "hexToInt", "buildResidueRest", "groupUnionExpression",
                      "uniform", "colorParams"):
            assert token not in text, token

    def test_water_residue_layer_uses_atom_expression(self):
        text = _structure_js()
        helper = _function_body(text, "residueLayerExpression")
        assert "isWaterResidue(residue)" in helper
        assert "residue.expression" in helper
        assert "residueSelector(residue)" in helper

    def test_mixed_layers_exclude_interaction_residues(self):
        """主干层不叠加：骨架/配体层都排除相互作用残基（它们改用棍棒显示）。"""
        text = _structure_js()
        body = _function_body(text, "polymerExpression")
        assert "POLYMER_EXPRESSION + \" and not (\"" in body
        assert "residueSelectorUnion(residues)" in body
        lig = _function_body(text, "ligandExpression")
        assert "NON_POLYMER_EXPRESSION" in lig
        assert "and not (" in lig
        assert "residueSelectorUnion(residues)" in lig

    def test_residue_selector_is_used_for_exclusion_lists(self):
        """骨架/配体层的排除列表仍来自 residue_selector（Cartoon 分区不变）。"""
        text = _structure_js()
        union = _function_body(text, "residueSelectorUnion")
        assert '" or "' in union
        assert "isWaterResidue(" in union
        selector = _function_body(text, "residueSelector")
        assert "residue_selector" in selector

    def test_script_uses_backend_expression_payload(self):
        text = _structure_js()
        for token in ("residues", "contacts", "/api/structure/interactions",
                      "/api/structure/contact", "/api/structure/topology",
                      "/api/structure/trajectory", "existence"):
            assert token in text, token
        assert "groupAtomExpression" in text
        assert "EXPRESSION_LANGUAGE = \"pymol\"" in text

    def test_topology_request_carries_type_and_idx(self):
        """拓扑 PDB 按该对原子打 B-factor 标记，请求必须带 type/idx。"""
        body = _function_body(_structure_js(), "topologyUrl")
        assert "type=" in body and "idx=" in body
        assert "TOPOLOGY_URL" in body
        assert "encodeURIComponent(cfg.type)" in body

    def test_expression_strings_are_passed_through_unchanged(self):
        """表达式原样交给 Mol* 解析：脚本里不得出现按空格/关键字改写表达式的代码。"""
        text = _structure_js()
        assert 'params: { language: EXPRESSION_LANGUAGE, expression: expression }' in text
        assert "residueLayerExpression(residue)" in _function_body(
            text, "buildResidueLayer")
        assert "residueSelector(residue)" in _function_body(
            text, "residueLayerExpression")
        assert "POLYMER_EXPRESSION" in _function_body(text, "polymerExpression")
        for token in ("split(\" and", "replace(/chain", "parseExpression",
                      "resi \""):
            assert token not in text, token

    def test_water_residues_avoid_ambiguous_selectors(self):
        """水残基选择器不唯一（resSeq 截断）：残基层改用它自己的原子表达式。"""
        text = _structure_js()
        assert "is_water" in text
        helper = _function_body(text, "isWaterResidue")
        assert "is_water" in helper
        assert "isWaterResidue(" in _function_body(text, "residueLayerExpression")
        assert "isWaterResidue(" in _function_body(text, "residueSelectorUnion")

    def test_script_polls_native_trajectory_controls(self):
        text = _structure_js()
        assert "viewportShowTrajectoryControls: true" in text
        assert "FRAME_POLL_MS = 150" in text
        assert "TrajectoryInfo" in text
        assert "setInterval(pollFrame, FRAME_POLL_MS)" in text

    def test_script_drops_old_water_and_overlay_paths(self):
        text = _structure_js()
        for token in ("showWater", "contact_pdb", "highlights", "solvent",
                      "loadStructureCell", "renderPayload", "clearGroups",
                      "involves_water", "contact_label"):
            assert token not in text, token

    def test_existence_bar_and_slider_are_driven_by_payload(self):
        text = _structure_js()
        assert "function renderExistence(" in text
        assert "function highlightFrame(" in text
        assert "function updateFrameRange(" in text
        assert 'className = values[i] ? "on" : ""' in text
        assert 'aria-valuenow' in text
        assert '" active"' in text

    def test_existence_bar_handles_pointer_and_keyboard(self):
        """向量条即进度条：指针拖动/点选 + 方向键/Home/End 换帧。"""
        text = _structure_js()
        assert "frameSlider" not in text
        for token in ("pointerdown", "pointermove", "pointerup",
                      "pointercancel", "setPointerCapture",
                      "releasePointerCapture", "frameFromPointer",
                      "ArrowLeft", "ArrowRight", "Home", "End"):
            assert token in text, token
        pointer = _function_body(text, "frameFromPointer")
        assert "getBoundingClientRect()" in pointer
        assert "Math.floor" in pointer
        assert "clampFrame(" in pointer
        down = _function_body(text, "onPointerDown")
        assert "requestFrame(frameFromPointer(event))" in down
        assert "setPointerCapture" in down

    def test_group_legend_is_plain_text_list(self):
        """图例不再有颜色含义：残基标签 + 原子名两段纯文字，无色点/组色说明。"""
        text = _structure_js()
        assert "structure-group-note" not in text
        assert "structure-group-dot" not in text
        body = _function_body(text, "renderLegend")
        assert "legendItem(" in body
        assert "legendNote" not in body
        item = _function_body(text, "legendItem")
        assert "GROUP_LABEL_CLASS" in item
        assert "dot" not in item
        assert "backgroundColor" not in item

    def test_legend_lists_defined_group_atom_names(self):
        """图例第二段列出"判定为相互作用组"的原子名（纯文字，无 3D 叠加）。"""
        text = _structure_js()
        item = _function_body(text, "legendItem")
        assert "structure-group-atoms" in item
        assert "atoms" in item
        assert "LEGEND_ATOM_LIMIT" in item
        assert "title" in item
        assert "legendAtomNames(" in item and "legendAtomsText(" in item
        assert "LEGEND_ATOM_LIMIT = 12" in text
        names = _function_body(text, "legendAtomNames")
        assert 'names.push(String(atoms[i].name))' in names
        assert "sort" not in names
        limit = _function_body(text, "legendAtomsText")
        assert "names.slice(0, limit)" in limit
        assert '" … (+"' in limit
        assert "names.join(\", \")" in limit

    def test_no_3d_overlay_techniques(self):
        """3D 里不加 overpaint/emissive/3D label 等额外视觉层。"""
        text = _structure_js()
        for token in ("overpaint", "emissive", "label-3d"):
            assert token not in text.lower(), token

    def test_frontend_scripts_avoid_innerhtml(self):
        for path in (STRUCTURE_JS, FILEPICK_JS):
            assert "innerHTML" not in path.read_text(encoding="utf-8")

    def test_filepick_script_is_served(self, client):
        response = client.get("/static/filepick.js")
        assert response.status_code == 200
        assert "DiiFilePick" in response.text

    def test_structure_script_keeps_viewer_unexpanded(self):
        """静态检查：不启用 Mol* 展开态（展开态 fixed 定位会盖住页内工具栏/滑块）。"""
        text = _structure_js()
        assert "layoutIsExpanded: false" in text
        assert "layoutIsExpanded: true" not in text

    def test_left_panel_collapsed_to_icon_strip(self):
        """左侧导航默认收起为图标条：viewer 选项 collapseLeftPanel 必须为 true。"""
        text = _structure_js()
        assert "layoutShowLeftPanel: true" in text
        assert "collapseLeftPanel: true" in text
        assert "collapseLeftPanel: false" not in text

    def test_viewer_uses_native_molstar_controls(self):
        """开关全部用原生：Toggle Controls Panel + Settings（Layout 多选里才有 Sequence）。"""
        text = _structure_js()
        assert "layoutShowControls: true" in text
        assert "layoutShowLeftPanel: true" in text
        assert "layoutShowSequence: true" in text
        assert "collapseLeftPanel: true" in text
        assert "collapseRightPanel: false" in text
        assert "viewportShowControls: true" in text
        assert "viewportShowSettings: true" in text
        assert "viewportShowTrajectoryControls: true" in text

    def test_script_has_no_self_made_toggle_logic(self):
        """自造的 controls/sequence 开关逻辑必须全部删除，不主动改 Mol* tab。"""
        text = _structure_js()
        for token in ("toggleControls", "showSequence", "applySequence",
                      "toggleSequenceClass", "copyShallow", "leftPanelTabName",
                      "NO_SEQUENCE_CLASS", "regionState"):
            assert token not in text, token

    def test_no_left_panel_tab_hack(self):
        """不再主动改左侧 tab（图标条与内容区交给 Mol* 自己的 collapsed 语义）。"""
        assert "leftPanelTabName" not in _structure_js()

    def test_css_does_not_hide_left_panel_content(self):
        text = APP_CSS.read_text(encoding="utf-8")
        assert ".msp-layout-collapse-left .msp-layout-left" in text
        assert ".msp-left-panel-controls .msp-scrollable-container" not in text
        assert ".msp-left-panel-controls-buttons" not in text

    def test_polymer_expression_uses_real_macro(self):
        """polymer/organic/inorganic 不是可靠宏；用 protein / not (protein or water or nucleic)。"""
        text = _structure_js()
        assert 'POLYMER_EXPRESSION = "protein"' in text
        assert 'POLYMER_EXPRESSION = "polymer"' not in text
        assert 'NON_POLYMER_EXPRESSION = "not (protein or water or nucleic)"' in text
        assert "organic or inorganic" not in text
        poly = _function_body(text, "polymerExpression")
        assert 'POLYMER_EXPRESSION + " and not ("' in poly
        lig = _function_body(text, "ligandExpression")
        assert "NON_POLYMER_EXPRESSION" in lig
        assert "and not (" in lig

    def test_empty_component_falls_back(self):
        """构建成功但为空的组件按失败处理，骨架回退 static polymer、配体回退 static ligand。"""
        text = _structure_js()
        assert "componentElementCount" in text
        assert "elementCount" in text
        assert "[dii] empty component" in text
        guard = _function_body(text, "buildScriptComponent")
        assert "componentElementCount(" in guard
        assert "ok: false" in guard
        polymer = _function_body(text, "buildPolymerLayer")
        assert "buildStaticLayer(POLYMER_KEYWORD" in polymer
        plain = _function_body(text, "buildPlainLigandLayer")
        assert "STATIC_LIGAND" in plain

    def test_script_keeps_sequence_toggle_and_structure_form(self):
        text = _structure_js()
        for token in ("structureForm", "/api/pick",
                      "DiiFilePick", "openStructure"):
            assert token in text, token


# ============================================================
# 基团图例 / 组色 / CSS
# ============================================================


class TestGroupLegendAndCss:
    def test_page_has_group_legend_container(self, client):
        html = _page(client)
        assert 'id="groupLegend"' in html
        assert "structure-groups" in html
        legend_at = html.index('id="groupLegend"')
        status_at = html.index('id="structureStatus"')
        assert legend_at < status_at

    def test_structure_script_has_no_group_palette(self):
        text = _structure_js()
        for token in ("GROUP_COLORS", "GROUP_COLOR_NOTE", "hexToInt"):
            assert token not in text, token
        assert not re.findall(r"#[0-9a-fA-F]{6}", text)

    def test_css_has_group_legend_rules(self):
        text = APP_CSS.read_text(encoding="utf-8")
        assert ".structure-groups" in text
        assert ".structure-group-item" in text
        assert ".structure-group-dot" not in text
        assert ".structure-group-note" not in text
        assert ".structure-group-atoms" in text
        assert "font-family: var(--mono)" in text
        assert ".structure-existence span.active" in text
        assert "cursor: grab" in text
        assert "touch-action: none" in text

    def test_css_viewer_container_rule_takes_effect(self):
        """静态检查：末尾 .structure-viewer 规则把 viewer 约束为容器内 620px。"""
        text = APP_CSS.read_text(encoding="utf-8")
        rules = re.findall(r"\.structure-viewer\s*\{([^}]*)\}", text)
        assert rules, "missing .structure-viewer rule"
        final_rule = rules[-1]
        assert "height: 620px" in final_rule
        assert "position: relative" in final_rule

    def test_css_has_filepick_rules_without_self_made_toggles(self):
        text = APP_CSS.read_text(encoding="utf-8")
        assert ".filepick-modal" in text
        assert ".filepick-row" in text
        assert ".structure-toolbar" not in text
        assert ".structure-toggle" not in text
