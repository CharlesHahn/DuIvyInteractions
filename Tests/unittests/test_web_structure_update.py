# -*- coding: utf-8 -*-
"""structure.js 原地切 pair 接口的静态断言（无浏览器依赖）。

检查：window.DiiStructure 契约、mount 幂等、updatePair 不重载结构/轨迹、
group_expression 优先与 b 区间回退、失败软处理。
"""

import re
from pathlib import Path

WEB_DIR = Path(__file__).resolve().parents[2] / "DuIvyInteractions" / "web"
STRUCTURE_JS = WEB_DIR / "static" / "structure.js"
HEX_COLOR = re.compile(r"#[0-9a-fA-F]{6}")


def _source() -> str:
    return STRUCTURE_JS.read_text(encoding="utf-8")


def _function_body(text: str, name: str) -> str:
    """取 structure.js 里 2 空格缩进的函数体（静态切片，不做 JS 解析）。"""
    match = re.search(r"  function %s\(.*?\n  \}\n" % re.escape(name), text,
                      re.S)
    assert match, f"missing function {name}"
    return match.group(0)


class TestUpdateApi:
    def test_exposes_window_contract(self):
        text = _source()
        assert "window.DiiStructure = {" in text
        for token in ("mount: mount", "updatePair: updatePair",
                      "isMounted:", "return !!viewer;"):
            assert token in text, token

    def test_mount_is_idempotent(self):
        body = _function_body(_source(), "mount")
        assert "if (mounted) return;" in body
        assert "mounted = true;" in body

    def test_script_self_mounts(self):
        """模板未改：脚本末尾自挂载，外层后续可幂等再调 mount()。"""
        text = _source()
        assert text.index("window.DiiStructure = {") < text.index(
            "window.DiiStructure.mount()")

    def test_mount_keeps_bootstrap_branches(self):
        body = _function_body(_source(), "mount")
        for token in ("highlightFrame(frame)", "shouldAutoLoad()",
                      "autoLoadStructure()", "cfg.has_structure",
                      "initViewer()", "select topology and trajectory"):
            assert token in body, token


class TestUpdatePair:
    def test_does_not_reload_topology_or_trajectory(self):
        body = _function_body(_source(), "updatePair")
        for token in ("TOPOLOGY_URL", "TRAJECTORY_URL", "loadTrajectoryOnce(",
                      "loadScene(", "fetch(", "clearDefaultComponents("):
            assert token not in body, token

    def test_rebuilds_pair_parts(self):
        body = _function_body(_source(), "updatePair")
        for token in ("clearResidueLayers()", "clearMeasurements()",
                      "buildResidueLayers()", "renderExistence(",
                      "updateFrameRange()", "renderLegend(", "renderTypeLegend("):
            assert token in body, token

    def test_follows_payload_contract(self):
        body = _function_body(_source(), "updatePair")
        for token in ("next.groups", "next.existence", "next.pt_types",
                      "next.n_frames", "next.pair_label", "next.frame"):
            assert token in body, token

    def test_refocuses_and_applies_frame(self):
        body = _function_body(_source(), "updatePair")
        assert "focusInitially()" in body
        assert "selectFrame(clampFrame(target))" in body

    def test_failure_is_soft(self):
        body = _function_body(_source(), "updatePair")
        assert 'console.warn("[dii] pair update failed"' in body
        assert "pair update failed: " in body
        assert "return Promise.resolve();" in body
        assert "if (!ready || !structureRef) {" in body

    def test_cleanup_only_targets_residue_layers(self):
        body = _function_body(_source(), "clearResidueLayers")
        assert "var stale = residueRefs;" in body
        assert "residueRefs = [];" in body
        assert "builder.delete(stale[i])" in body
        assert "deleteMeasurementCells(stale)" in _function_body(
            _source(), "clearMeasurements")


class TestGroupSelection:
    def test_group_expression_first_expression_fallback(self):
        body = _function_body(_source(), "groupSelection")
        assert "residue.group_expression" in body
        assert "residue.expression" in body

    def test_residue_layer_uses_group_selection_for_water(self):
        body = _function_body(_source(), "residueLayerExpression")
        assert "isWaterResidue(residue)" in body
        assert "groupSelection(residue)" in body

    def test_measurements_use_group_selection(self):
        body = _function_body(_source(), "groupLociStep")
        assert "groupSelection(residue)" in body
        assert "residue.expression" in body

    def test_no_innerhtml_and_no_hex_palette(self):
        text = _source()
        assert "innerHTML" not in text
        assert not HEX_COLOR.findall(text)

    def test_focus_uses_group_selection(self):
        """重聚焦也必须与 pair 无关（b 区间只对首屏 PDB 有效）。"""
        for name in ("residueUnionExpression", "focusFromEachResidue"):
            assert "groupSelection(" in _function_body(_source(), name), name
