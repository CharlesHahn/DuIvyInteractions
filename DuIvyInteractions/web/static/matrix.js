/* pairs × frames 矩阵：canvas 绘制 + 悬停逐帧取值 + 滚动虚拟化 */
(function () {
  "use strict";
  var cfgEl = document.getElementById("matrixConfig");
  if (!cfgEl) return;
  var cfg = JSON.parse(cfgEl.textContent);
  var canvas = document.getElementById("matrixCanvas");
  var wrap = document.querySelector(".matrix-wrap");
  var scroller = document.getElementById("matrixScroll");
  var spacer = document.getElementById("matrixSpacer");
  var tip = document.getElementById("matrixTip");
  var info = document.getElementById("matrixInfo");
  var slider = document.getElementById("frameSlider");
  var frameInfo = document.getElementById("frameInfo");
  var frameBody = document.getElementById("frameBody");
  var frameInput = document.getElementById("frameInput");
  var framePrev = document.getElementById("framePrev");
  var frameNext = document.getElementById("frameNext");
  var viewStartInput = document.getElementById("viewStart");
  var viewEndInput = document.getElementById("viewEnd");
  var applyWindowBtn = document.getElementById("applyWindow");
  var windowInfo = document.getElementById("windowInfo");
  var ctx = canvas.getContext("2d");

  var ROW_H = 16;        // 行高（px）
  var LABEL_W_MAX = 320;   // 标签列最大宽度
  var LABEL_W_MIN = 150;   // 标签列最小宽度
  var LABEL_FRACTION = 0.4;
  var AXIS_H = 22;       // 底部时间轴高
  var OCC_W = 46;        // 标签列内占有率条宽
  var PAD = 8;
  var ACTIVE = "60,84,136";      // #3c5488 的 RGB 分量（矩阵活跃块，非 π-堆积类型）
  var CACHE_LIMIT = 2000;
  var MIN_COLUMN_WIDTH = 1e-6;  // 防止 0 帧时除零
  var ALPHA_MIN = 0.3, ALPHA_SPAN = 0.7;
  var MIN_RUN_WIDTH = 0.6;
  var TIP_OFFSET = 14, TIP_MARGIN = 8;
  var VALUE_DIGITS = 3;
  var COLOR_BG = "#ffffff";
  var COLOR_ZEBRA = "#fbfbfd";
  var COLOR_LABEL = "#475569";
  var COLOR_MUTED = "#7b8794";
  var COLOR_BAR_BG = "#eceef3";
  var COLOR_BAR_FG = "#8491b4";  // 行内 occupancy 小条
  var COLOR_AXIS_BG = "#f8fafc";
  var COLOR_AXIS_LINE = "#e5e7eb";
  var COLOR_HOVER_ROW = "rgba(15,23,42,0.06)";
  var COLOR_HOVER_COL = "rgba(15,23,42,0.07)";
  var COLOR_FRAME_LINE = "rgba(15,23,42,0.55)";
  var COLOR_FRAME_COL = "rgba(60,84,136,0.10)";  // 选中帧列底色（与 hover 同层，不改变活跃块配色）
  // π-堆积类型配色：与 figures.PT_COLORS 保持一致
  var PT_COLORS = {P: "#3c5488", T: "#e64b35", N: "#8491b4"};
  var ptColors = cfg.pt_colors || PT_COLORS;
  var LABEL_FONT = "11.5px ui-monospace, Consolas, monospace";
  var AXIS_FONT = "10.5px ui-monospace, Consolas, monospace";
  var EMPTY_FONT = "13px system-ui, sans-serif";
  var AXIS_LABEL_MARGIN = 56;
  var MAX_VISIBLE_FRAMES = 600;
  var FRAME_PAGE = 200;

  function labelColWidth() {
    return Math.min(LABEL_W_MAX,
                    Math.max(LABEL_W_MIN, canvas.clientWidth * LABEL_FRACTION));
  }

  var rows = [], times = [], nFrames = cfg.n_frames;
  var totalPairs = 0, loaded = false, matrixError = "";
  var hoverRow = -1, hoverCol = -1, selectedFrame = 0;
  var viewStart = 0, viewEnd = 0;
  var tipSeq = 0, frameSeq = 0;
  var cellCache = {};
  var sliderTimer = null;
  var dpr = window.devicePixelRatio || 1;

  function resize() {
    var w = canvas.clientWidth, h = canvas.clientHeight;
    canvas.width = Math.round(w * dpr);
    canvas.height = Math.round(h * dpr);
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    if (rows.length) {
      spacer.style.height = Math.max(0, rows.length * ROW_H - plotHeight()) + "px";
    }
    draw();
  }

  function plotHeight() {
    return canvas.clientHeight - AXIS_H;
  }

  function columnWidth() {
    return Math.max(MIN_COLUMN_WIDTH,
                    (canvas.clientWidth - labelColWidth()) /
                    Math.max(viewEnd - viewStart, 1));
  }

  function columnX(frame) {
    return labelColWidth() + (frame - viewStart) * columnWidth();
  }

  function visibleRange() {
    var top = scroller.scrollTop;
    var first = Math.max(0, Math.floor(top / ROW_H));
    var last = Math.min(rows.length, Math.ceil((top + plotHeight()) / ROW_H));
    return [first, last, top];
  }

  function draw() {
    var w = canvas.clientWidth, h = canvas.clientHeight;
    ctx.clearRect(0, 0, w, h);
    ctx.fillStyle = COLOR_BG;
    ctx.fillRect(0, 0, w, h);
    if (!rows.length || nFrames === 0) {
      drawEmptyState();
      return;
    }
    var range = visibleRange(), colW = columnWidth();
    drawRows(range, colW);
    drawOverlays(range, colW);
    drawSelectedFrame();
    drawAxis();
  }

  function drawEmptyState() {
    ctx.fillStyle = COLOR_MUTED;
    ctx.font = EMPTY_FONT;
    var text = "loading…";
    if (matrixError) {
      text = "matrix failed to load";
    } else if (loaded) {
      if (nFrames === 0) {
        text = "this h5 contains no frames";
      } else {
        text = totalPairs === 0
          ? "this h5 contains no pairs"
          : "no pairs match the current filter";
      }
    }
    ctx.fillText(text, labelColWidth() + 24, 40);
  }

  function drawRows(range, colW) {
    for (var r = range[0]; r < range[1]; r++) {
      var row = rows[r];
      var y = r * ROW_H - range[2];
      if (r % 2 === 0) {
        ctx.fillStyle = COLOR_ZEBRA;
        ctx.fillRect(0, y, canvas.clientWidth, ROW_H);
      }
      drawRuns(row, y, colW);
      drawLabel(row, y);
    }
  }

  function drawRuns(row, y, colW) {
    var alpha = ALPHA_MIN + ALPHA_SPAN * row.occupancy;
    if (!row.typed) {
      drawSegments(row.runs, y, colW, "rgba(" + ACTIVE + "," + alpha.toFixed(3) + ")");
      return;
    }
    Object.keys(row.typed).forEach(function (type) {
      drawSegments(row.typed[type], y, colW, typeFill(type, alpha));
    });
  }

  function drawSegments(runs, y, colW, fill) {
    ctx.fillStyle = fill;
    runs.forEach(function (run) {
      var start = Math.max(run[0], viewStart);
      var end = Math.min(run[1], viewEnd);
      if (end <= start) return;
      var x = labelColWidth() + (start - viewStart) * colW;
      var width = Math.max(MIN_RUN_WIDTH, (end - start) * colW);
      ctx.fillRect(x, y + 1, width, ROW_H - 2);
    });
  }

  function typeFill(type, alpha) {
    var hex = ptColors[type] || PT_COLORS.N;
    var value = parseInt(hex.slice(1), 16);
    var rgb = [value >> 16 & 255, value >> 8 & 255, value & 255];
    return "rgba(" + rgb.join(",") + "," + alpha.toFixed(3) + ")";
  }

  function drawOverlays(range, colW) {
    if (hoverRow >= 0) {
      var hy = hoverRow * ROW_H - range[2];
      ctx.fillStyle = COLOR_HOVER_ROW;
      ctx.fillRect(0, hy, canvas.clientWidth, ROW_H);
    }
    if (hoverCol >= 0) {
      ctx.fillStyle = COLOR_HOVER_COL;
      ctx.fillRect(columnX(hoverCol), 0, Math.max(1, colW), plotHeight());
    }
  }

  function drawLabel(row, y) {
    ctx.save();
    ctx.beginPath();
    ctx.rect(0, y, labelColWidth() - OCC_W - 2 * PAD, ROW_H);
    ctx.clip();
    ctx.fillStyle = COLOR_LABEL;
    ctx.font = LABEL_FONT;
    ctx.textBaseline = "middle";
    ctx.fillText(row.label, PAD, y + ROW_H / 2);
    ctx.restore();
    var barY = y + ROW_H / 2 - 2;
    ctx.fillStyle = COLOR_BAR_BG;
    ctx.fillRect(labelColWidth() - OCC_W - PAD, barY, OCC_W, 4);
    ctx.fillStyle = COLOR_BAR_FG;
    ctx.fillRect(labelColWidth() - OCC_W - PAD, barY, OCC_W * row.occupancy, 4);
  }

  function drawSelectedFrame() {
    if (selectedFrame < viewStart || selectedFrame >= viewEnd) return;
    var x = columnX(selectedFrame);
    ctx.fillStyle = COLOR_FRAME_COL;
    ctx.fillRect(x, 0, Math.max(1, columnWidth()), plotHeight());
    ctx.fillStyle = COLOR_FRAME_LINE;
    ctx.fillRect(x, 0, 1, plotHeight());
  }

  function drawAxis() {
    var y = plotHeight();
    ctx.fillStyle = COLOR_AXIS_BG;
    ctx.fillRect(0, y, canvas.clientWidth, AXIS_H);
    ctx.strokeStyle = COLOR_AXIS_LINE;
    ctx.beginPath();
    ctx.moveTo(0, y + 0.5);
    ctx.lineTo(canvas.clientWidth, y + 0.5);
    ctx.stroke();
    var step = Math.max(1, Math.ceil((viewEnd - viewStart) / 8));
    ctx.fillStyle = COLOR_MUTED;
    ctx.font = AXIS_FONT;
    ctx.textBaseline = "middle";
    for (var f = viewStart; f < viewEnd; f += step) {
      var x = columnX(f);
      if (x > canvas.clientWidth - AXIS_LABEL_MARGIN) break;
      ctx.fillText(formatTime(times[f]), x + 4, y + AXIS_H / 2);
    }
  }

  function formatTime(value) {
    if (typeof value !== "number") return "";
    return Math.abs(value) >= 1000 ? value.toFixed(0) : value.toFixed(1);
  }

  function frameAt(x) {
    var frame = viewStart + Math.floor((x - labelColWidth()) / columnWidth());
    return Math.min(Math.max(frame, viewStart), viewEnd - 1);
  }

  function hitTest(event) {
    var rect = canvas.getBoundingClientRect();
    var x = event.clientX - rect.left;
    var y = event.clientY - rect.top;
    var inPlot = y >= 0 && y < plotHeight();
    // 底部时间轴条（AXIS_H）作为列头热区：顶部是数据行，放热区会抢首行点击
    var inAxis = rows.length > 0 && y >= plotHeight() &&
      y <= canvas.clientHeight;
    var row = Math.floor((y + scroller.scrollTop) / ROW_H);
    var inColumns = viewEnd > viewStart && x > labelColWidth();
    return {
      row: inPlot && row >= 0 && row < rows.length ? row : -1,
      col: inPlot && inColumns ? frameAt(x) : -1,
      axisCol: inAxis && inColumns ? frameAt(x) : -1,
      x: x, y: y
    };
  }

  function showTip(hit) {
    var row = rows[hit.row];
    var seq = ++tipSeq;
    tip.classList.remove("hidden");
    tip.innerHTML = "";
    appendTipTitle(row);
    addTipRow("occupancy", (row.occupancy * 100).toFixed(1) + "%");
    if (hit.col < 0) {
      placeTip(hit);
      return;
    }
    addTipRow("frame", hit.col + " · " + formatTime(times[hit.col]) + " ps");
    placeTip(hit);
    fetchTipValues(row.idx, hit.col, seq, hit);
  }

  function fetchTipValues(pairIdx, frame, seq, hit) {
    fetchCell(pairIdx, frame).then(function (cell) {
      if (seq !== tipSeq) return;
      addTipValues(cell);
      placeTip(hit);
    }).catch(function () {
      if (seq === tipSeq) addTipRow("values", "unavailable");
    });
  }

  function appendTipTitle(row) {
    var title = document.createElement("div");
    title.className = "t-label";
    title.textContent = row.label;
    tip.appendChild(title);
  }

  function addTipRow(key, value) {
    var line = document.createElement("div");
    line.className = "t-row";
    var k = document.createElement("span");
    k.className = "t-key";
    k.textContent = key;
    var v = document.createElement("span");
    v.className = "t-val";
    v.textContent = value;
    line.appendChild(k);
    line.appendChild(v);
    tip.appendChild(line);
  }

  function addTipValues(cell) {
    addTipRow("exists", cell.exists ? "true" : "false");
    Object.keys(cell.values).forEach(function (name) {
      var unit = (cell.units || {})[name];
      addTipRow(unit ? name + " (" + unit + ")" : name,
                cell.values[name] === null ? "—" : String(cell.values[name]));
    });
    Object.keys(cell.flags || {}).forEach(function (name) {
      addTipRow(name, cell.flags[name] ? "true" : "false");
    });
  }

  function placeTip(hit) {
    var maxX = wrap.clientWidth - tip.offsetWidth - TIP_MARGIN;
    var left = Math.min(hit.x + TIP_OFFSET, Math.max(TIP_MARGIN, maxX));
    var top = Math.min(hit.y + TIP_OFFSET,
                       wrap.clientHeight - tip.offsetHeight - TIP_MARGIN);
    tip.style.left = left + "px";
    tip.style.top = Math.max(TIP_MARGIN, top) + "px";
  }

  function readJson(resp) {
    if (!resp.ok) throw new Error("HTTP " + resp.status);
    return resp.json();
  }

  function fetchCell(pair, frame) {
    var key = pair + ":" + frame;
    if (cellCache[key]) return Promise.resolve(cellCache[key]);
    if (Object.keys(cellCache).length > CACHE_LIMIT) cellCache = {};
    return fetch("/api/cell/" + encodeURIComponent(cfg.type) + "/" + pair + "/" + frame)
      .then(readJson)
      .then(function (data) { cellCache[key] = data; return data; });
  }

  function onMove(event) {
    var hit = hitTest(event);
    var col = hit.axisCol >= 0 ? hit.axisCol : hit.col;
    var changed = hit.row !== hoverRow || col !== hoverCol;
    hoverRow = hit.row;
    hoverCol = col;
    if (!changed) return;
    canvas.style.cursor = hit.row >= 0 || hit.axisCol >= 0 ? "pointer" : "default";
    draw();
    if (hit.row >= 0) {
      showTip(hit);
    } else {
      tipSeq += 1;
      tip.classList.add("hidden");
    }
  }

  function onClick(event) {
    var hit = hitTest(event);
    if (hit.axisCol >= 0) {
      loadFrame(hit.axisCol);
      return;
    }
    if (hit.row < 0) return;
    openPair(rows[hit.row].idx, selectedFrame);
  }

  function loadMatrix() {
    fetch(matrixUrl())
      .then(readJson)
      .then(applyMatrix)
      .catch(handleMatrixError);
  }

  function matrixUrl() {
    return cfg.url + "?search=" + encodeURIComponent(cfg.search || "") +
      "&min_occupancy=" + (cfg.min_occupancy || 0) +
      "&sort_by=" + encodeURIComponent(cfg.sort_by || "occupancy") +
      "&sort_dir=" + encodeURIComponent(cfg.sort_dir || "desc");
  }

  function applyMatrix(data) {
    rows = data.rows;
    times = data.times;
    nFrames = data.n_frames;
    totalPairs = data.total_pairs;
    loaded = true;
    matrixError = "";
    resetHover();
    resetWindow();
    updateScrollHeight();
    info.textContent = rows.length + " / " + data.total_pairs + " pairs";
    resize();
    if (nFrames > 0) {
      loadFrame(initialFrame());
    } else {
      showNoFrames();
    }
  }

  function initialFrame() {
    var frame = parseInt(cfg.frame, 10);
    if (!isFinite(frame) || frame < 0) return 0;
    return Math.min(frame, nFrames - 1);
  }

  function resetHover() {
    hoverRow = -1;
    hoverCol = -1;
    tip.classList.add("hidden");
  }

  function updateScrollHeight() {
    var maxScroll = nFrames
      ? Math.max(0, rows.length * ROW_H - plotHeight())
      : 0;
    spacer.style.height = maxScroll + "px";
    if (scroller.scrollTop > maxScroll) scroller.scrollTop = maxScroll;
  }

  function resetWindow() {
    viewStart = 0;
    viewEnd = Math.min(nFrames, MAX_VISIBLE_FRAMES);
    syncWindow();
  }

  function syncWindow() {
    if (viewStartInput) viewStartInput.value = viewStart;
    if (viewEndInput) viewEndInput.value = Math.max(viewStart, viewEnd - 1);
    if (windowInfo) {
      windowInfo.textContent = nFrames
        ? "frames " + viewStart + "–" + (viewEnd - 1) + " of " + nFrames
        : "no frames";
    }
  }

  function applyWindow(start, last) {
    if (!isFinite(start) || !isFinite(last) || nFrames === 0) {
      showWindowError();
      return;
    }
    var range = clampWindow(start, last);
    viewStart = range[0];
    viewEnd = range[1];
    hoverRow = -1;
    hoverCol = -1;
    tip.classList.add("hidden");
    syncWindow();
    draw();
    if (selectedFrame < viewStart || selectedFrame >= viewEnd) {
      loadFrame(viewStart);
    }
  }

  function clampWindow(start, last) {
    start = Math.max(0, Math.min(Math.floor(start), nFrames - 1));
    last = Math.max(start, Math.min(Math.floor(last), nFrames - 1));
    if (last - start + 1 > MAX_VISIBLE_FRAMES) {
      last = Math.min(nFrames - 1, start + MAX_VISIBLE_FRAMES - 1);
    }
    return [start, last + 1];
  }

  function showWindowError() {
    if (windowInfo) windowInfo.textContent = "invalid window input";
  }

  function ensureWindow(frame) {
    if (frame >= viewStart && frame < viewEnd) return;
    var half = Math.floor(MAX_VISIBLE_FRAMES / 2);
    applyWindow(frame - half, frame + half);
  }

  function showNoFrames() {
    syncFrameNav(0);
    frameInfo.textContent = "no frames in this h5";
    frameBody.textContent = "no frames in this h5";
  }

  function handleMatrixError(error) {
    loaded = true;
    matrixError = String(error);
    info.textContent = "failed to load matrix: " + error;
    frameInfo.textContent = "matrix unavailable";
    frameBody.textContent = "matrix unavailable";
    draw();
  }

  function clampFrame(value) {
    if (nFrames <= 0) return 0;
    if (!isFinite(value)) return selectedFrame;
    return Math.min(Math.max(Math.round(value), 0), nFrames - 1);
  }

  // 帧导航控件与矩阵列头选帧共用同一 selectedFrame，任何入口后同步四态
  function syncFrameNav(frame) {
    if (frameInput) frameInput.value = frame;
    if (framePrev) framePrev.disabled = frame <= 0;
    if (frameNext) frameNext.disabled = !nFrames || frame >= nFrames - 1;
  }

  function stepFrame(delta) {
    var frame = selectedFrame + delta;
    if (frame < 0 || frame >= nFrames) return;
    loadFrame(frame);
  }

  function commitFrameInput() {
    if (!frameInput) return;
    var frame = clampFrame(parseInt(frameInput.value, 10));
    if (frame === selectedFrame) {
      frameInput.value = selectedFrame;
      return;
    }
    loadFrame(frame);
  }

  function loadFrame(frame, limit) {
    var seq = ++frameSeq;
    frame = clampFrame(frame);
    selectedFrame = frame;
    if (slider) slider.value = frame;
    syncFrameNav(frame);
    ensureWindow(frame);
    draw();
    frameInfo.textContent = "frame " + frame + " · loading…";
    fetch("/api/frame/" + encodeURIComponent(cfg.type) + "/" + frame +
          "?limit=" + (limit || FRAME_PAGE))
      .then(readJson)
      .then(function (data) {
        if (seq !== frameSeq) return;
        syncUrl(frame);
        renderFrame(data);
      })
      .catch(function (error) {
        if (seq === frameSeq) {
          frameInfo.textContent = "failed to load frame " + frame + ": " + error;
          frameBody.textContent = "failed to load frame " + frame;
        }
      });
  }

  // 把当前帧号写回 URL（保留 search/min_occupancy/sort_by/sort_dir/page 等既有参数）
  function syncUrl(frame) {
    if (!window.history || !window.history.replaceState) return;
    var url = new URL(window.location.href);
    url.searchParams.set("frame", frame);
    window.history.replaceState(null, "", url.pathname + url.search + url.hash);
  }

  function renderFrame(data) {
    frameInfo.textContent = "frame " + data.frame + " · " + formatTime(data.time) +
      " ps · " + data.n_active + " active";
    frameBody.classList.remove("muted");
    frameBody.innerHTML = "";
    if (!data.pairs.length) {
      frameBody.textContent = "no interaction active at this frame";
      return;
    }
    frameBody.appendChild(buildFrameTable(data));
    appendLoadMore(data);
  }

  function buildFrameTable(data) {
    var table = document.createElement("table");
    table.appendChild(frameHead(data));
    var tbody = document.createElement("tbody");
    data.pairs.forEach(function (pair) {
      var tr = frameRow(pair, data);
      tr.className = "clickable";
      tr.addEventListener("click", function () { openPair(pair.idx); });
      tbody.appendChild(tr);
    });
    table.appendChild(tbody);
    return table;
  }

  function frameHead(data) {
    var thead = document.createElement("thead");
    var head = document.createElement("tr");
    var units = data.units || {};
    var metrics = data.metrics.map(function (name) {
      return units[name] ? name + " (" + units[name] + ")" : name;
    });
    var columns = ["pair", "occupancy"].concat(metrics)
      .concat(data.string_metrics || []).concat(data.flags);
    columns.forEach(function (name) {
      var th = document.createElement("th");
      th.textContent = name;
      head.appendChild(th);
    });
    thead.appendChild(head);
    return thead;
  }

  function frameRow(pair, data) {
    var tr = document.createElement("tr");
    addCell(tr, pair.label, "pair-cell");
    addCell(tr, (pair.occupancy * 100).toFixed(1) + "%");
    data.metrics.forEach(function (name) { addCell(tr, formatValue(pair.values[name])); });
    (data.string_metrics || []).forEach(function (name) {
      addCell(tr, formatValue(pair.values[name]));
    });
    data.flags.forEach(function (name) { addCell(tr, pair.flags[name] ? "yes" : "no"); });
    return tr;
  }

  function appendLoadMore(data) {
    if (data.shown >= data.n_active) return;
    if (data.limit_max && data.shown >= data.limit_max) {
      var note = document.createElement("p");
      note.className = "muted";
      note.textContent = "showing first " + data.shown + " of " +
        data.n_active + " active pairs (max per request)";
      frameBody.appendChild(note);
      return;
    }
    var button = document.createElement("button");
    button.className = "btn btn-sm";
    button.textContent = "load more (" + (data.n_active - data.shown) + " remaining)";
    button.addEventListener("click", function () {
      loadFrame(data.frame, data.shown + FRAME_PAGE);
    });
    frameBody.appendChild(button);
  }

  /* 站内跳转优先走 pjax 路由；无 router 时保留整页兜底 */
  function navigateTo(url) {
    if (window.DiiRouter && typeof window.DiiRouter.navigate === "function") {
      window.DiiRouter.navigate(url);
      return;
    }
    window.location.href = url;
  }

  function openPair(idx, frame) {
    var url = "/pair/" + encodeURIComponent(cfg.type) + "/" + idx;
    if (typeof frame === "number" && frame >= 0) url += "?frame=" + frame;
    navigateTo(url);
  }

  function addCell(tr, text, cls) {
    var td = document.createElement("td");
    td.textContent = text;
    if (cls) td.className = cls;
    tr.appendChild(td);
  }

  function formatValue(value) {
    if (value === null || value === undefined) return "—";
    return typeof value === "number" ? value.toFixed(VALUE_DIGITS) : String(value);
  }

  canvas.addEventListener("mousemove", onMove);
  canvas.addEventListener("mouseleave", function () {
    hoverRow = -1;
    hoverCol = -1;
    tipSeq += 1;
    tip.classList.add("hidden");
    canvas.style.cursor = "default";
    draw();
  });
  canvas.addEventListener("click", onClick);
  scroller.addEventListener("scroll", draw);
  if (applyWindowBtn) {
    applyWindowBtn.addEventListener("click", function () {
      applyWindow(parseInt(viewStartInput.value, 10),
                  parseInt(viewEndInput.value, 10));
    });
  }
  if (slider) {
    slider.addEventListener("input", function () {
      if (sliderTimer) clearTimeout(sliderTimer);
      sliderTimer = setTimeout(function () {
        loadFrame(parseInt(slider.value, 10));
      }, 80);
    });
  }
  if (framePrev) {
    framePrev.addEventListener("click", function () { stepFrame(-1); });
  }
  if (frameNext) {
    frameNext.addEventListener("click", function () { stepFrame(1); });
  }
  if (frameInput) {
    frameInput.addEventListener("keydown", function (event) {
      if (event.key === "Enter") {
        event.preventDefault();
        commitFrameInput();
      } else if (event.key === "ArrowLeft") {
        event.preventDefault();
        stepFrame(-1);
      } else if (event.key === "ArrowRight") {
        event.preventDefault();
        stepFrame(1);
      }
    });
    frameInput.addEventListener("blur", commitFrameInput);
  }
  window.addEventListener("resize", resize);
  resize();
  loadMatrix();
})();
