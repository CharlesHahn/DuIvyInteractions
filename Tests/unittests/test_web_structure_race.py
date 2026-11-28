# -*- coding: utf-8 -*-
"""3D 结构页竞态回归：loadTrajectory resolve 时 Mol* 结构树可能尚未更新。

本文件只做源码静态切片断言 + 可选的 Node 逻辑测试；无浏览器环境无法复现
真实竞态，故不启动 Mol*，也不使用 TestClient。
"""

import json
import re
import shutil
import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
STRUCTURE_JS = (REPO_ROOT / "DuIvyInteractions" / "web" / "static"
                / "structure.js")
WAIT_INTERVAL_MS = 100
WAIT_TIMEOUT_MS = 5000
DIAGNOSTIC_HEAD_CHARS = 40
# 假拓扑：换行落在前 40 字符内，用来验证诊断 head 的转义
TOPOLOGY = ("data_test\n#\nloop_\n_atom_site.id\n""ATOM 1 N N ASN A 1 1 1 ASN A N ? 0.000 0.000 0.000 1.00 0.00\n")
CFG = {"has_structure": True, "n_frames": 101}

NODE_SCRIPT = r"""
const fs = require("fs");
const source = fs.readFileSync(process.argv[1], "utf8");
const payload = JSON.parse(process.argv[2]);

function extract(name) {
  const match = source.match(
    new RegExp("  function " + name + "\\([\\s\\S]*?\\n  \\}"));
  if (!match) throw new Error("missing function " + name);
  return match[0];
}

function constant(name) {
  const match = source.match(new RegExp("var " + name + " = (\\d+);"));
  if (!match) throw new Error("missing constant " + name);
  return Number(match[1]);
}

function buildWait(structureEntry, setTimeoutImpl, timeoutText) {
  const factory = new Function(
    "structureEntry", "setTimeout", "STRUCTURE_WAIT_INTERVAL_MS",
    "STRUCTURE_WAIT_TIMEOUT_MS", "structureTimeoutText",
    extract("waitForStructure") + "\nreturn waitForStructure;");
  return factory(structureEntry, setTimeoutImpl,
    constant("STRUCTURE_WAIT_INTERVAL_MS"),
    constant("STRUCTURE_WAIT_TIMEOUT_MS"), timeoutText);
}

function buildTimeoutText(cfg, topologySnapshot) {
  const factory = new Function(
    "cfg", "topologySnapshot", "STRUCTURE_WAIT_TIMEOUT_MS",
    "DIAGNOSTIC_HEAD_CHARS",
    extract("diagnosticHead") + "\n" + extract("structureDiagnostics") + "\n" +
      extract("structureTimeoutText") + "\nreturn structureTimeoutText;");
  return factory(cfg, topologySnapshot,
    constant("STRUCTURE_WAIT_TIMEOUT_MS"),
    constant("DIAGNOSTIC_HEAD_CHARS"));
}

(async () => {
  const result = {
    interval: constant("STRUCTURE_WAIT_INTERVAL_MS"),
    timeoutMs: constant("STRUCTURE_WAIT_TIMEOUT_MS"),
    headChars: constant("DIAGNOSTIC_HEAD_CHARS")
  };

  // A. 前两次结构树为空、第三次出现：轮询应 resolve 同一个 entry
  const queue = [];
  let calls = 0;
  const entry = {cell: {}};
  const waitSuccess = buildWait(
    () => { calls += 1; return calls >= 3 ? entry : null; },
    (fn) => { queue.push(fn); return queue.length; },
    () => "unexpected timeout");
  const pending = waitSuccess();
  while (queue.length) queue.shift()();
  const value = await pending;
  result.success = {
    sameEntry: value === entry,
    polls: calls,
    pendingTimers: queue.length
  };

  // B. 永不出现：每 100ms 一次轮询，5000ms 上限后应 reject
  const timeoutQueue = [];
  let timeoutPolls = 0;
  let timeoutTimers = 0;
  const waitNever = buildWait(
    () => { timeoutPolls += 1; return null; },
    (fn) => {
      timeoutTimers += 1;
      timeoutQueue.push(fn);
      return timeoutTimers;
    },
    () => "timeout-diagnostic");
  const settled = waitNever().then(
    () => ({rejected: false}),
    (error) => ({rejected: true, message: error.message}));
  while (timeoutQueue.length) timeoutQueue.shift()();
  result.timeout = await settled;
  result.timeout.polls = timeoutPolls;
  result.timeout.timers = timeoutTimers;

  // C. 用假拓扑 / 假配置构造真实超时文案
  const timeoutText = buildTimeoutText(payload.cfg, payload.topology);
  result.message = timeoutText();
  result.messageTopologyChars = payload.topology.length;

  process.stdout.write(JSON.stringify(result));
})().catch((error) => {
  process.stderr.write(String((error && error.stack) || error));
  process.exit(1);
});
"""

NODE_BIN = shutil.which("node")
needs_node = pytest.mark.skipif(NODE_BIN is None, reason="node is unavailable")


def _source() -> str:
    return STRUCTURE_JS.read_text(encoding="utf-8")


def _function(name: str) -> str:
    """取 2 空格缩进的函数体（静态切片，不做 JS 解析）。"""
    match = re.search(r"  function %s\([\s\S]*?\n  \}" % re.escape(name),
                      _source())
    assert match, f"missing function {name}"
    return match.group(0)


def _constant(name: str) -> int:
    match = re.search(r"var %s = (\d+);" % re.escape(name), _source())
    assert match, f"missing constant {name}"
    return int(match.group(1))


def _run_node(payload: dict) -> dict:
    completed = subprocess.run(
        [NODE_BIN, "-e", NODE_SCRIPT, str(STRUCTURE_JS),
         json.dumps(payload)],
        capture_output=True, text=True, timeout=60)
    assert completed.returncode == 0, completed.stderr
    return json.loads(completed.stdout)


def _expected_timeout_message(topology: str, cfg: dict) -> str:
    head = (topology[:DIAGNOSTIC_HEAD_CHARS]
            .replace("\\", "\\\\")
            .replace("\r", "\\r")
            .replace("\n", "\\n")
            .replace('"', '\\"'))
    return ('no structure after %dms (topology %d chars, head "%s", '
            'has_structure=%s, n_frames=%s)' % (
                WAIT_TIMEOUT_MS, len(topology), head,
                str(cfg["has_structure"]).lower(), cfg["n_frames"]))


# ============================================================
# 常量与轮询等待
# ============================================================


class TestWaitConstants:
    def test_wait_constants_are_centralized(self):
        text = _source()
        assert _constant("STRUCTURE_WAIT_INTERVAL_MS") == WAIT_INTERVAL_MS
        assert _constant("STRUCTURE_WAIT_TIMEOUT_MS") == WAIT_TIMEOUT_MS
        assert _constant("DIAGNOSTIC_HEAD_CHARS") == DIAGNOSTIC_HEAD_CHARS
        wait_at = text.index("function waitForStructure")
        assert text.index("var STRUCTURE_WAIT_INTERVAL_MS") < wait_at
        assert text.index("var STRUCTURE_WAIT_TIMEOUT_MS") < wait_at


class TestWaitForStructure:
    def test_polls_structure_tree_with_interval(self):
        body = _function("waitForStructure")
        assert "structureEntry()" in body
        assert "resolve(entry)" in body
        assert "waited += STRUCTURE_WAIT_INTERVAL_MS" in body
        assert "setTimeout(poll, STRUCTURE_WAIT_INTERVAL_MS)" in body
        # 禁止魔法数字
        assert "100" not in body
        assert "5000" not in body

    def test_timeout_branch_rejects_with_diagnostics(self):
        body = _function("waitForStructure")
        assert re.search(r"if \(waited >= STRUCTURE_WAIT_TIMEOUT_MS\)", body)
        assert "reject(new Error(structureTimeoutText()))" in body


# ============================================================
# afterLoad 走等待，不再一次性读取后立刻抛错
# ============================================================


class TestAfterLoad:
    def test_after_load_waits_for_structure(self):
        body = _function("afterLoad")
        assert "waitForStructure()" in body
        assert "refOf(entry.cell)" in body
        assert "modelRef = findModelRef();" in body

    def test_after_load_keeps_post_load_order(self):
        body = _function("afterLoad")
        order = ["clearDefaultComponents()", "buildLayers()", "focusInitially()",
                 "selectFrame(frame)", "startPolling()"]
        positions = [body.index(token) for token in order]
        assert positions == sorted(positions)

    def test_legacy_one_shot_throw_is_removed(self):
        text = _source()
        assert "trajectory produced no structure" not in text
        body = _function("afterLoad")
        assert not re.search(r"if\s*\(\s*!entry\s*\)\s*throw", body)
        assert "var entry = structureEntry();" not in body


# ============================================================
# 拓扑响应内容校验与超时诊断
# ============================================================


class TestTopologyValidation:
    def test_topology_snapshot_is_module_state(self):
        text = _source()
        assert re.search(r'var topologySnapshot = "";', text)
        assert text.index("var topologySnapshot") < text.index(
            "function fetchTopologyText")

    def test_fetch_topology_validates_mmcif(self):
        body = _function("fetchTopologyText")
        assert "isCifText(text)" in body
        assert "topology response is not mmCIF" in body
        assert "text.length" in body
        assert "topologySnapshot = text;" in body
        assert "return response.text();" in body
        assert "return text;" in body
        check = _function("isCifText")
        assert 'data_' in check and '_atom_site.' in check

    def test_diagnostic_head_is_escaped_single_line(self):
        body = _function("diagnosticHead")
        assert "slice(0, DIAGNOSTIC_HEAD_CHARS)" in body
        assert body.count(".replace(") >= 3
        assert "\\n" in body


class TestTimeoutDiagnostics:
    def test_timeout_message_carries_diagnostics(self):
        timeout_body = _function("structureTimeoutText")
        diagnostics = _function("structureDiagnostics")
        assert "STRUCTURE_WAIT_TIMEOUT_MS" in timeout_body
        assert "no structure after " in timeout_body
        for token in ("topologySnapshot.length", " chars", "head",
                      "has_structure", "n_frames"):
            assert token in diagnostics, token
        assert "diagnosticHead(topologySnapshot)" in diagnostics

    def test_load_scene_surfaces_failure_text(self):
        body = _function("loadScene")
        assert "setStatus(loadFailureText(error), true)" in body
        assert "error.message" in _function("loadFailureText")


# ============================================================
# Node 逻辑测试：直接跑源码里的 waitForStructure（IIFE 私有函数靠提取）
# ============================================================


@needs_node
class TestWaitForStructureLogic:
    def test_polling_resolves_and_times_out(self):
        result = _run_node({"topology": TOPOLOGY, "cfg": CFG})
        assert result["interval"] == WAIT_INTERVAL_MS
        assert result["timeoutMs"] == WAIT_TIMEOUT_MS
        assert result["success"] == {"sameEntry": True, "polls": 3,
                                     "pendingTimers": 0}
        assert result["timeout"]["rejected"] is True
        assert result["timeout"]["message"] == "timeout-diagnostic"
        assert result["timeout"]["timers"] == WAIT_TIMEOUT_MS // WAIT_INTERVAL_MS
        assert result["timeout"]["polls"] == (
            WAIT_TIMEOUT_MS // WAIT_INTERVAL_MS + 1)

    def test_timeout_message_uses_fake_topology(self):
        result = _run_node({"topology": TOPOLOGY, "cfg": CFG})
        assert result["messageTopologyChars"] == len(TOPOLOGY)
        assert result["message"] == _expected_timeout_message(TOPOLOGY, CFG)
        assert "\n" not in result["message"]
        assert "\\n" in result["message"]
        assert result["headChars"] == DIAGNOSTIC_HEAD_CHARS
