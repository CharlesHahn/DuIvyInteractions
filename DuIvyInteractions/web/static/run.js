/* W3：Web 发起计算（run.html 表单页 + job.html 任务页共用；按页面元素存在与否分支） */
(function () {
  "use strict";

  var PICK_ENDPOINT = "/api/pick";
  var RUN_ENDPOINT = "/api/run/";
  var SELECT_URL = "/select";
  var POLL_INTERVAL_MS = 800;
  var PICK_MODE_FILE = "file";
  var PICK_MODE_DIR = "dir";
  var PICK_KIND = "structure";
  var BROWSE_PICKER = "picker";
  var BROWSE_SERVER = "server";
  var BROWSE_TARGET_ATTR = "data-target";
  var BROWSE_MODE_ATTR = "data-browse";
  var BROWSE_KIND_ATTR = "data-kind";
  var JOB_ACTIVE = {pending: true, running: true};
  var JOB_TERMINAL = {done: true, failed: true, cancelled: true};
  var JOB_STATES = ["pending", "running", "done", "failed", "cancelled"];
  var PLACEHOLDER = {
    topology: "<topology>",
    trajectory: "<trajectory>",
    output: "<output>",
    ff: "<ff>",
    strategy: "<strategy>",
    interactions: "<interactions>"
  };
  var PAIRING_ANY = "any";
  var STRATEGY_PER_TUPLE = "per_tuple";
  var HYDROPHOBIC = "hydrophobic";
  var INTERACTION_SELECTOR = 'input[name="interactions"]';
  var INTERACTION_CHECKED_SELECTOR = 'input[name="interactions"]:checked';
  var BROWSE_BUTTON_SELECTOR = "button[data-browse]";

  function byId(id) { return document.getElementById(id); }

  function setStatus(node, text, isError) {
    if (!node) return;
    node.textContent = text || "";
    node.classList.toggle("hidden", !text);
    node.classList.toggle("error", !!text && !!isError);
  }

  function show(node) { if (node) node.classList.remove("hidden"); }
  function hide(node) { if (node) node.classList.add("hidden"); }

  function fieldValue(node) {
    return node && node.value ? node.value.trim() : "";
  }

  function setInteractions(form, checked) {
    var boxes = form.querySelectorAll(INTERACTION_SELECTOR);
    for (var i = 0; i < boxes.length; i++) boxes[i].checked = checked;
  }

  function checkedInteractions(form) {
    var boxes = form.querySelectorAll(INTERACTION_CHECKED_SELECTOR);
    var names = [];
    for (var i = 0; i < boxes.length; i++) names.push(boxes[i].value);
    return names;
  }

  function selectLines(text) {
    var raw = (text || "").split("\n");
    var lines = [];
    for (var i = 0; i < raw.length; i++) {
      var line = raw[i].trim();
      if (line) lines.push(line);
    }
    return lines;
  }

  function quoteArg(text) {
    return '"' + text.replace(/"/g, '\\"') + '"';
  }

  function buildCommand(refs) {
    var selections = selectLines(refs.selectInput.value);
    var interactions = checkedInteractions(refs.form);
    var parts = [
      "dii run",
      "-t", fieldValue(refs.topology) || PLACEHOLDER.topology,
      "-f", fieldValue(refs.trajectory) || PLACEHOLDER.trajectory,
      "-o", fieldValue(refs.output) || PLACEHOLDER.output,
      "--ff", fieldValue(refs.ff) || PLACEHOLDER.ff,
      "--strategy", fieldValue(refs.strategy) || PLACEHOLDER.strategy,
      "--interactions", interactions.length ? interactions.join(",") : PLACEHOLDER.interactions
    ];
    if (selections.length) {
      parts.push("--select");
      for (var i = 0; i < selections.length; i++) parts.push(quoteArg(selections[i]));
    }
    if (refs.pairing.value && refs.pairing.value !== PAIRING_ANY) {
      parts.push("--pairing", refs.pairing.value);
    }
    if (refs.markRedundant.checked) parts.push("--mark-redundant");
    return parts.join(" ");
  }

  function updateNotes(refs) {
    if (refs.perTupleNote) {
      refs.perTupleNote.classList.toggle(
        "hidden", refs.strategy.value !== STRATEGY_PER_TUPLE);
    }
    if (refs.markNote) {
      var boxes = refs.form.querySelectorAll(INTERACTION_SELECTOR);
      var hydrophobic = false;
      for (var i = 0; i < boxes.length; i++) {
        if (boxes[i].value === HYDROPHOBIC && boxes[i].checked) hydrophobic = true;
      }
      refs.markNote.classList.toggle(
        "hidden", !(refs.markRedundant.checked && !hydrophobic));
    }
  }

  function updatePreview(refs) {
    refs.preview.textContent = buildCommand(refs);
    updateNotes(refs);
  }

  function applyPickedPath(refs, input, path) {
    input.value = path;
    updatePreview(refs);
    setStatus(refs.status, "selected " + path, false);
  }

  function readResponse(response) {
    return response.json().then(function (data) {
      return {ok: response.ok, status: response.status, data: data};
    }, function () {
      return {ok: false, status: response.status, data: null};
    });
  }

  function pickNative(refs, input, kind) {
    var mode = kind === PICK_MODE_DIR ? PICK_MODE_DIR : PICK_MODE_FILE;
    setStatus(refs.status, "waiting for the system dialog…", false);
    fetch(PICK_ENDPOINT, {
      method: "POST",
      headers: {"Content-Type": "application/x-www-form-urlencoded"},
      body: "mode=" + encodeURIComponent(mode) + "&kind=" + encodeURIComponent(PICK_KIND)
    }).then(readResponse).then(function (result) {
      var data = result.data || {};
      if (data.path) {
        applyPickedPath(refs, input, String(data.path));
        return;
      }
      setStatus(refs.status, data.error
        ? "browse unavailable: " + data.error : "no path selected", true);
    }, function () {
      setStatus(refs.status, "browse failed: network error", true);
    });
  }

  function openServerBrowser(refs, input, kind) {
    if (!window.DiiFilePick || typeof window.DiiFilePick.open !== "function") {
      setStatus(refs.status, "server file browser is not available", true);
      return;
    }
    window.DiiFilePick.open({
      kind: kind,
      onPick: function (path) { applyPickedPath(refs, input, String(path)); }
    });
  }

  function browseKind(button) {
    return button.getAttribute(BROWSE_KIND_ATTR) === PICK_MODE_DIR
      ? PICK_MODE_DIR : PICK_MODE_FILE;
  }

  function onBrowseClick(refs, event) {
    var button = event.currentTarget;
    var input = byId(button.getAttribute(BROWSE_TARGET_ATTR));
    if (!input) {
      setStatus(refs.status, "browse target input is missing", true);
      return;
    }
    var mode = button.getAttribute(BROWSE_MODE_ATTR);
    var kind = browseKind(button);
    if (mode === BROWSE_PICKER) pickNative(refs, input, kind);
    else if (mode === BROWSE_SERVER) openServerBrowser(refs, input, kind);
    else setStatus(refs.status, "browsing is not available", true);
  }

  function initFormPage() {
    var form = byId("runForm");
    var preview = byId("commandPreview");
    if (!form || !preview) return;
    var refs = {
      form: form,
      preview: preview,
      topology: byId("topologyInput"),
      trajectory: byId("trajectoryInput"),
      output: byId("outputInput"),
      ff: byId("ffSelect"),
      strategy: byId("strategySelect"),
      selectInput: byId("selectInput"),
      pairing: byId("pairingSelect"),
      markRedundant: byId("markRedundant"),
      perTupleNote: byId("perTupleNote"),
      markNote: byId("markNote"),
      status: byId("runStatus")
    };
    if (!refs.topology || !refs.trajectory || !refs.output || !refs.ff ||
        !refs.strategy || !refs.selectInput || !refs.pairing ||
        !refs.markRedundant) return;
    var buttons = form.querySelectorAll(BROWSE_BUTTON_SELECTOR);
    for (var i = 0; i < buttons.length; i++) {
      buttons[i].addEventListener("click", function (event) {
        onBrowseClick(refs, event);
      });
    }
    var allBtn = byId("interactionsAll");
    var noneBtn = byId("interactionsNone");
    if (allBtn) allBtn.addEventListener("click", function () {
      setInteractions(form, true);
      updatePreview(refs);
    });
    if (noneBtn) noneBtn.addEventListener("click", function () {
      setInteractions(form, false);
      updatePreview(refs);
    });
    form.addEventListener("input", function () { updatePreview(refs); });
    form.addEventListener("change", function () { updatePreview(refs); });
    updatePreview(refs);
  }

  function isTerminal(status) { return !!JOB_TERMINAL[status]; }

  function isActive(status) { return !!JOB_ACTIVE[status]; }

  function elapsedText(value) {
    var number = Number(value);
    return isFinite(number) ? number.toFixed(1) + "s" : String(value);
  }

  function readSince(logEl) {
    var value = parseInt(logEl.getAttribute("data-since"), 10);
    return isFinite(value) && value > 0 ? value : 0;
  }

  function updateJobStatus(state, payload) {
    var status = payload.status || state.statusEl.getAttribute("data-status");
    if (status) state.statusEl.setAttribute("data-status", status);
    if (state.stateEl) {
      for (var i = 0; i < JOB_STATES.length; i++) {
        state.stateEl.classList.remove("job-state-" + JOB_STATES[i]);
      }
      state.stateEl.classList.add("job-state-" + status);
      state.stateEl.textContent = status;
    }
    if (state.elapsedEl && payload.elapsed !== undefined &&
        payload.elapsed !== null) {
      state.elapsedEl.textContent = elapsedText(payload.elapsed);
    }
    if (state.cancelBtn) state.cancelBtn.disabled = !isActive(status);
  }

  function splitLog(logEl) {
    var text = logEl.textContent;
    if (!text) return [];
    var lines = text.split("\n");
    if (lines.length && lines[lines.length - 1] === "") lines.pop();
    return lines;
  }

  function sameBlock(left, right) {
    if (left.length !== right.length) return false;
    for (var i = 0; i < left.length; i++) {
      if (left[i] !== right[i]) return false;
    }
    return true;
  }

  function overlapLength(existing, incoming) {
    var limit = Math.min(existing.length, incoming.length);
    for (var k = limit; k > 0; k--) {
      if (sameBlock(existing.slice(existing.length - k),
                    incoming.slice(0, k))) return k;
    }
    return 0;
  }

  function containedAt(block, lines) {
    if (!block.length || block.length > lines.length) return -1;
    for (var start = 0; start + block.length <= lines.length; start++) {
      if (sameBlock(block, lines.slice(start, start + block.length))) {
        return start;
      }
    }
    return -1;
  }

  function appendLines(logEl, lines, reconcile) {
    if (!lines || !lines.length) return;
    var fresh = lines;
    if (reconcile) {
      var existing = splitLog(logEl);
      var overlap = overlapLength(existing, fresh);
      if (overlap) {
        fresh = fresh.slice(overlap);
      } else if (existing.length) {
        var at = containedAt(existing, fresh);
        if (at >= 0) fresh = fresh.slice(at + existing.length);
      }
      if (!fresh.length) return;
    }
    var text = fresh.join("\n");
    if (logEl.textContent &&
        logEl.textContent.charAt(logEl.textContent.length - 1) !== "\n") {
      text = "\n" + text;
    }
    logEl.appendChild(document.createTextNode(text + "\n"));
    logEl.scrollTop = logEl.scrollHeight;
  }

  function renderH5(listEl, files) {
    if (!listEl) return;
    while (listEl.firstChild) listEl.removeChild(listEl.firstChild);
    var names = files || [];
    if (!names.length) {
      var empty = document.createElement("li");
      empty.className = "mono";
      empty.textContent = "no h5 files";
      listEl.appendChild(empty);
      return;
    }
    for (var i = 0; i < names.length; i++) {
      var item = document.createElement("li");
      item.className = "mono";
      item.textContent = String(names[i]);
      listEl.appendChild(item);
    }
  }

  function showResult(state, payload) {
    if (state.timer) {
      clearTimeout(state.timer);
      state.timer = null;
    }
    var status = payload.status || state.statusEl.getAttribute("data-status");
    if (state.resultEl) {
      state.resultEl.classList.remove("hidden");
      if (payload.output_dir) {
        state.resultEl.setAttribute("data-output-dir",
          String(payload.output_dir));
      }
    }
    hide(state.doneEl);
    hide(state.failedEl);
    hide(state.cancelledEl);
    if (status === "done") {
      if (payload.h5_files && payload.h5_files.length !== undefined) {
        renderH5(state.h5El, payload.h5_files);
      }
      show(state.doneEl);
    } else if (status === "failed") {
      if (payload.error && state.failedErrorEl) {
        state.failedErrorEl.textContent = String(payload.error);
        state.failedErrorEl.classList.remove("hidden");
      }
      show(state.failedEl);
    } else if (status === "cancelled") {
      show(state.cancelledEl);
    }
    updateJobStatus(state, payload);
  }

  function schedulePoll(state) {
    if (state.timer) clearTimeout(state.timer);
    state.timer = setTimeout(function () {
      state.timer = null;
      pollJob(state);
    }, POLL_INTERVAL_MS);
  }

  function pollJob(state) {
    fetch(RUN_ENDPOINT + encodeURIComponent(state.jobId) +
          "?since=" + state.since).then(readResponse).then(function (result) {
      if (!result.ok) {
        var detail = result.data && result.data.detail
          ? String(result.data.detail) : "HTTP " + result.status;
        setStatus(state.errorEl, "poll failed: " + detail, true);
        if (result.status === 404 || result.status === 400) return;
        schedulePoll(state);
        return;
      }
      var data = result.data || {};
      setStatus(state.errorEl, "", false);
      appendLines(state.logEl, data.lines, state.firstPoll);
      state.firstPoll = false;
      if (typeof data.next === "number") state.since = data.next;
      updateJobStatus(state, data);
      if (isTerminal(data.status)) {
        showResult(state, data);
        return;
      }
      schedulePoll(state);
    }, function () {
      setStatus(state.errorEl, "poll failed: network error", true);
      schedulePoll(state);
    });
  }

  function cancelJob(state) {
    if (state.cancelBtn) state.cancelBtn.disabled = true;
    fetch(RUN_ENDPOINT + encodeURIComponent(state.jobId) + "/cancel", {
      method: "POST"
    }).then(readResponse).then(function (result) {
      if (result.ok && result.data && result.data.status) {
        var status = String(result.data.status);
        if (isTerminal(status)) showResult(state, {status: status});
        else updateJobStatus(state, {status: status});
        return;
      }
      if (state.cancelBtn) state.cancelBtn.disabled = false;
      var detail = result.data && result.data.detail
        ? String(result.data.detail) : "HTTP " + result.status;
      setStatus(state.errorEl, "cancel failed: " + detail, true);
    }, function () {
      if (state.cancelBtn) state.cancelBtn.disabled = false;
      setStatus(state.errorEl, "cancel failed: network error", true);
    });
  }

  function openResults(state) {
    var outputDir = state.resultEl
      ? state.resultEl.getAttribute("data-output-dir") : "";
    if (!outputDir) {
      setStatus(state.errorEl, "output directory is not available", true);
      return;
    }
    var form = document.createElement("form");
    form.method = "post";
    form.action = SELECT_URL;
    var input = document.createElement("input");
    input.type = "hidden";
    input.name = "path";
    input.value = outputDir;
    form.appendChild(input);
    document.body.appendChild(form);
    form.submit();
  }

  function initJobPage() {
    var statusEl = byId("jobStatus");
    var logEl = byId("jobLog");
    if (!statusEl || !logEl) return;
    var jobId = statusEl.getAttribute("data-job-id");
    if (!jobId) return;
    var state = {
      jobId: jobId,
      statusEl: statusEl,
      logEl: logEl,
      stateEl: byId("jobState"),
      elapsedEl: byId("jobElapsed"),
      errorEl: byId("jobStatusError"),
      cancelBtn: byId("cancelRun"),
      resultEl: byId("jobResult"),
      doneEl: byId("jobDone"),
      failedEl: byId("jobFailed"),
      cancelledEl: byId("jobCancelled"),
      h5El: byId("jobH5List"),
      failedErrorEl: byId("jobFailedError"),
      openBtn: byId("openResults"),
      since: readSince(logEl),
      firstPoll: true,
      timer: null
    };
    if (state.cancelBtn) {
      state.cancelBtn.addEventListener("click", function () { cancelJob(state); });
    }
    if (state.openBtn) {
      state.openBtn.addEventListener("click", function () { openResults(state); });
    }
    var initial = statusEl.getAttribute("data-status");
    if (isTerminal(initial)) {
      showResult(state, {status: initial});
      return;
    }
    schedulePoll(state);
  }

  initFormPage();
  initJobPage();
})();
