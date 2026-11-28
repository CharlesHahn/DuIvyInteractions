# -*- coding: utf-8 -*-
"""结构复用判据与加载诊断：structure_version 优先，loadScene(reason) 可追踪。

静态断言（无浏览器）。目标：同一已选结构下切 interaction 类型时
attach 判定 same → 不重载轨迹；若重载，Console 的 [dii] loadScene 行指出原因。
"""

import re
from pathlib import Path

STRUCTURE_JS = (Path(__file__).resolve().parents[2] / "DuIvyInteractions"
                / "web" / "static" / "structure.js")


def _source() -> str:
    return STRUCTURE_JS.read_text(encoding="utf-8")


def _function_body(name: str) -> str:
    text = _source()
    match = re.search(r"  function %s\(.*?\n  \}\n" % re.escape(name), text,
                      re.S)
    assert match, f"missing function {name}"
    return match.group(0)


def test_structure_key_includes_version():
    body = _function_body("structureKey")
    assert "version:" in body
    assert "config.structure_version" in body
    assert "topology:" in body and "trajectory:" in body


def test_same_structure_compares_version_before_paths():
    body = _function_body("sameStructure")
    assert "if (next.version && currentStructure.version)" in body
    assert "return next.version === currentStructure.version;" in body
    version_guard = body.index("next.version && currentStructure.version")
    path_compare = body.index("next.topology === currentStructure.topology")
    assert version_guard < path_compare


def test_empty_fingerprint_guards_retained():
    body = _function_body("sameStructure")
    assert "if (!next.topology && !next.trajectory) return true;" in body
    assert "if (!currentStructure.topology && !currentStructure.trajectory)" in body
    assert "return !!viewer && ready;" in body


def test_load_scene_logs_reason_and_structure_fields():
    body = _function_body("loadScene")
    assert 'console.warn("[dii] loadScene", reason, {' in body
    for field in ("topology: cfg.topology", "trajectory: cfg.trajectory",
                  "version: cfg.structure_version"):
        assert field in body, field
    assert "setStatus(\"loading topology and trajectory…\", false)" in body


def test_load_scene_records_current_structure_before_fetch():
    body = _function_body("loadScene")
    assert "currentStructure = structureKey(cfg);" in body
    record = body.index("currentStructure = structureKey(cfg)")
    fetch = body.index("fetchTopologyText()")
    assert record < fetch


def test_every_load_scene_call_has_reason():
    text = _source()
    calls = re.findall(r"loadScene\(([^)]*)\)", text)
    call_args = [value.strip() for value in calls[1:]]
    assert '"mount"' in call_args
    assert '"reload"' in call_args
    assert all(value for value in call_args), call_args


def test_mount_path_starts_load_scene():
    body = _function_body("initViewer")
    assert 'loadScene("mount")' in body


def test_reload_path_clears_all_layers_then_reloads():
    body = _function_body("reloadStructure")
    assert "clearMeasurements()" in body
    assert "clearResidueLayers()" in body
    assert "clearBackgroundLayers()" in body
    assert 'loadScene("reload")' in body
    order = [body.index("clearMeasurements()"),
             body.index("clearResidueLayers()"),
             body.index("clearBackgroundLayers()"),
             body.index('loadScene("reload")')]
    assert order == sorted(order)


def test_attach_reloads_only_when_structure_differs():
    body = _function_body("attach")
    assert "sameStructure(config)" in body
    assert "reloadStructure()" in body
    guard = body.index("if (same)")
    reload = body.index("reloadStructure()")
    assert guard < reload


def test_update_pair_still_keeps_background_and_skips_load_scene():
    body = _function_body("updatePair")
    assert "buildBackgroundLayers()" not in body
    assert "loadScene(" not in body
    assert "clearResidueLayers()" in body
    assert "buildResidueLayers()" in body
