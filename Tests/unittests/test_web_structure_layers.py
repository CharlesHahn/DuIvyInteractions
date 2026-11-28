# -*- coding: utf-8 -*-
"""结构背景层实现契约：内置 static 组件 + 诊断去重 + 残基表达式保留。

背景层（卡通/配体/离子）不再拼 protein/water/nucleic 选择语言表达式；
残基棍棒、测量、聚焦仍用后端 chain/resi/name 表达式。
"""

import re
from pathlib import Path

STATIC_ROOT = Path(__file__).resolve().parents[2] / "DuIvyInteractions" / "web" / "static"
STRUCTURE_JS = STATIC_ROOT / "structure.js"
BACKGROUND_NAMES = ("buildPolymerLayer", "buildLigandLayer",
                    "buildStaticLayer", "buildBackgroundLayers")


def _source() -> str:
    return STRUCTURE_JS.read_text(encoding="utf-8")


def _function_body(text: str, name: str) -> str:
    match = re.search(r"  function %s\(.*?\n  \}\n" % re.escape(name), text,
                      re.S)
    assert match, f"missing function {name}"
    return match.group(0)


def test_background_expression_paths_removed():
    text = _source()
    for token in ("POLYMER_EXPRESSION", "NON_POLYMER_EXPRESSION",
                  "polymerExpression(", "ligandExpression(",
                  "residueSelectorUnion(", "buildPlainLigandLayer"):
        assert token not in text, token
    assert "protein and not (" not in text
    assert "not (protein or water or nucleic)" not in text


def test_polymer_layer_uses_static_keyword():
    body = _function_body(_source(), "buildPolymerLayer")
    assert 'buildStaticLayer("polymer"' in body
    assert "residues" not in body
    assert "excluded" not in body
    static = _function_body(_source(), "buildStaticLayer")
    assert "tryCreateComponentStatic(structureRef, keyword" in static


def test_ligand_layer_falls_back_to_ion():
    body = _function_body(_source(), "buildLigandLayer")
    assert 'buildStaticLayer("ligand"' in body
    assert 'buildStaticLayer("ion"' in body
    assert body.index('buildStaticLayer("ligand"') < body.index(
        'buildStaticLayer("ion"')


def test_background_layers_built_only_on_mount_and_reload():
    text = _source()
    body = _function_body(text, "buildBackgroundLayers")
    assert "buildPolymerLayer()" in body
    assert "buildLigandLayer()" in body
    assert "residues" not in body
    assert "buildBackgroundLayers()" in _function_body(text, "buildLayers")
    assert 'loadScene("reload")' in _function_body(text, "reloadStructure")
    assert "buildBackgroundLayers()" not in _function_body(text, "updatePair")


def test_clear_functions_have_disjoint_responsibilities():
    text = _source()
    residue = _function_body(text, "clearResidueLayers")
    background = _function_body(text, "clearBackgroundLayers")
    assert "backgroundRefs" not in residue
    assert "residueRefs" in residue
    assert "backgroundRefs" in background
    assert "residueRefs" not in background
    reload_body = _function_body(text, "reloadStructure")
    assert "clearResidueLayers()" in reload_body
    assert "clearBackgroundLayers()" in reload_body


def test_update_pair_reports_timing_and_keeps_background():
    body = _function_body(_source(), "updatePair")
    assert 'console.warn("[dii] pair updated in"' in body
    assert "backgroundKept: true" in body
    assert "clearResidueLayers()" in body
    assert "clearMeasurements()" in body
    assert "buildResidueLayers()" in body
    assert "buildBackgroundLayers()" not in body


def test_dock_attach_uses_console_warn():
    text = _source()
    assert 'console.warn("[dii] dock attach", {' in text
    assert 'console.info("[dii] dock attach"' not in text


def test_empty_component_warns_once_with_set():
    text = _source()
    helper = _function_body(text, "warnOnce")
    assert "new Set()" in text
    assert "warnedOnce.has(key)" in helper
    assert "warnedOnce.add(key)" in helper
    body = _function_body(text, "buildScriptComponent")
    assert "warnOnce(" in body
    assert "[dii] empty component" in body


def test_static_failure_warns_once_without_throw():
    body = _function_body(_source(), "buildStaticLayer")
    assert "warnOnce(" in body
    assert "[dii] static component failed" in body
    assert "throw" not in body
    assert "noteFailure(" not in body


def test_static_components_are_tracked_for_rebuild():
    body = _function_body(_source(), "buildStaticLayer")
    assert "trackBackgroundRef(component)" in body
    clear = _function_body(_source(), "clearBackgroundLayers")
    assert "backgroundRefs" in clear


def test_residue_sticks_still_use_backend_selectors():
    text = _source()
    layer = _function_body(text, "residueLayerExpression")
    assert "isWaterResidue(residue)" in layer
    assert "residueSelector(residue)" in layer
    assert "group_expression" in _function_body(text, "groupSelection")
    assert "residue_selector" in _function_body(text, "residueSelector")


def test_measurements_and_focus_still_use_group_selection():
    text = _source()
    assert "groupSelection(residue)" in _function_body(text, "groupLociStep")
    assert "groupSelection(" in _function_body(text, "residueUnionExpression")
    assert "groupSelection(" in _function_body(text, "focusFromEachResidue")


def test_pymol_expression_language_only_for_residue_layers():
    text = _source()
    assert 'EXPRESSION_LANGUAGE = "pymol"' in text
    script = _function_body(text, "buildScriptComponent")
    assert "language: EXPRESSION_LANGUAGE" in script
    for name in BACKGROUND_NAMES:
        assert "buildScriptComponent" not in _function_body(text, name), name
