/* 单对页：Plotly 对齐曲线 + 逐帧明细表分页 + 页内切 pair（不整页跳转） */
(function () {
  "use strict";
  var figEl = document.getElementById("pairFigure");
  var cfgEl = document.getElementById("pairConfig");
  if (!figEl || !cfgEl) return;
  var figure = JSON.parse(figEl.textContent);
  var cfg = JSON.parse(cfgEl.textContent);
  var chart = document.getElementById("pairChart");
  var VALUE_DIGITS = 3;
  var TIME_DIGITS = 2;
  var PLOT_CONFIG = {
    responsive: true,
    displaylogo: false,
    scrollZoom: false,
    modeBarButtonsToRemove: ["select2d", "lasso2d"]
  };
  var frames = null;      /* {el, state}：帧表分页上下文（切 pair 后重置） */
  var switching = false;  /* 防重复点击 */

  renderChart();
  renderFrames();
  bindPairNav();

  function renderChart() {
    if (typeof Plotly === "undefined") {
      chart.textContent = "plotly.js failed to load";
      return;
    }
    try {
      createChart();
    } catch (error) {
      chart.textContent = "failed to render chart: " + error;
    }
  }

  function createChart() {
    var promise = Plotly.newPlot(chart, figure.data, figure.layout,
      PLOT_CONFIG);
    if (promise && promise.catch) {
      promise.catch(function (error) {
        chart.textContent = "failed to render chart: " + error;
      });
    }
  }

  function renderFrames() {
    var el = frameElements();
    var state = {start: 0, next: null, pageSize: 1};
    frames = {el: el, state: state};
    el.prev.addEventListener("click", function () { goPrev(el, state); });
    el.next.addEventListener("click", function () { goNext(el, state); });
    loadFramePage(el, state);
  }

  function frameElements() {
    return {
      table: document.getElementById("framesTable"),
      info: document.getElementById("framesInfo"),
      prev: document.getElementById("framesPrev"),
      next: document.getElementById("framesNext")
    };
  }

  function readJson(resp) {
    if (!resp.ok) throw new Error("HTTP " + resp.status);
    return resp.json();
  }

  function loadFramePage(el, state) {
    el.prev.disabled = true;
    el.next.disabled = true;
    fetch("/api/frames/" + encodeURIComponent(cfg.type) + "/" + cfg.idx +
          "?start=" + state.start)
      .then(readJson)
      .then(function (data) { fillFrameTable(el, state, data); })
      .catch(function (error) {
        el.table.textContent = "failed to load frames: " + error;
      });
  }

  function goPrev(el, state) {
    if (state.start <= 0) return;
    state.start = Math.max(0, state.start - state.pageSize);
    loadFramePage(el, state);
  }

  function goNext(el, state) {
    if (state.next === null) return;
    state.start = state.next;
    loadFramePage(el, state);
  }

  function fillFrameTable(el, state, data) {
    state.pageSize = data.page_size || data.rows.length || 1;
    state.next = data.next_start === undefined ? null : data.next_start;
    el.info.textContent = frameInfoText(data);
    el.prev.disabled = data.start <= 0;
    el.next.disabled = state.next === null;
    el.table.replaceChildren();
    if (!data.rows.length) {
      el.table.textContent = "no frames";
      return;
    }
    el.table.appendChild(buildTable(data));
  }

  function frameInfoText(data) {
    if (!data.total) return "no frames";
    return "frames " + (data.start + 1) + "–" +
      (data.start + data.rows.length) + " of " + data.total;
  }

  function buildTable(data) {
    var table = document.createElement("table");
    var thead = document.createElement("thead");
    thead.appendChild(headerRow(data));
    table.appendChild(thead);
    var tbody = document.createElement("tbody");
    data.rows.forEach(function (row) { tbody.appendChild(bodyRow(row, data)); });
    table.appendChild(tbody);
    return table;
  }

  function headerRow(data) {
    var tr = document.createElement("tr");
    var units = data.units || {};
    var metrics = data.metrics.map(function (name) {
      return units[name] ? name + " (" + units[name] + ")" : name;
    });
    ["frame", "time (ps)", "exists"].concat(metrics)
      .concat(data.string_metrics || []).concat(data.flags)
      .forEach(function (name) {
        var th = document.createElement("th");
        th.textContent = name;
        tr.appendChild(th);
      });
    return tr;
  }

  function bodyRow(row, data) {
    var tr = document.createElement("tr");
    if (!row.exists) tr.className = "inactive";
    addCell(tr, row.frame);
    addCell(tr, formatNumber(row.time, TIME_DIGITS));
    addCell(tr, row.exists ? "yes" : "no");
    data.metrics.forEach(function (name) { addCell(tr, formatValue(row.values[name])); });
    (data.string_metrics || []).forEach(function (name) {
      addCell(tr, formatValue(row.values[name]));
    });
    data.flags.forEach(function (name) { addCell(tr, row.flags[name] ? "yes" : "no"); });
    return tr;
  }

  function addCell(tr, value) {
    var td = document.createElement("td");
    td.textContent = value;
    tr.appendChild(td);
  }

  function formatValue(value) {
    if (value === null || value === undefined) return "—";
    return typeof value === "number" ? formatNumber(value, VALUE_DIGITS) : String(value);
  }

  function formatNumber(value, digits) {
    return typeof value === "number" ? value.toFixed(digits) : "—";
  }

  /* ---------------- 页内切 pair（结构/轨迹常驻，不整页跳转） ---------------- */

  function bindPairNav() {
    var head = document.getElementById("pairHead");
    if (head) head.addEventListener("click", onPairNavClick);
  }

  function onPairNavClick(event) {
    var trigger = event.target && event.target.closest
      ? event.target.closest("[data-pair-nav]") : null;
    if (!trigger) return;
    event.preventDefault();
    switchPair(cfg.idx + Number(trigger.getAttribute("data-pair-nav")));
  }

  function panelUrl(idx) {
    return "/api/pair/" + encodeURIComponent(cfg.type) + "/" + idx +
      "/panels";
  }

  function dataUrl(idx) {
    return "/api/pair/" + encodeURIComponent(cfg.type) + "/" + idx;
  }

  function switchPair(target) {
    target = parseInt(target, 10);
    if (switching || isNaN(target)) return;
    if (target < 0 || target >= cfg.n_pairs) return;
    switching = true;
    Promise.all([
      fetch(panelUrl(target)).then(readJson),
      fetch(dataUrl(target)).then(readJson)
    ]).then(function (results) {
      applyPanels(results[0]);
      return applyPairData(results[1], target);
    }).catch(function (error) {
      console.warn("[dii] pair switch failed", error);
      toast("failed to switch pair: " + (error && error.message
        ? error.message : error));
    }).then(function () {
      switching = false;
    });
  }

  function applyPanels(panels) {
    replaceSection("pairHead", panels.head);
    replaceSection("pairSummary", panels.summary);
    replaceSection("pairGroups", panels.groups);
  }

  /* 用 DOMParser 解析出的节点替换容器内容（只操作 DOM 节点） */
  function replaceSection(elementId, html) {
    var host = document.getElementById(elementId);
    if (!host || !html) return;
    var body = new DOMParser().parseFromString(html, "text/html").body;
    var nodes = Array.prototype.slice.call(body.childNodes);
    host.replaceChildren.apply(host, nodes);
  }

  function applyPairData(payload, target) {
    cfg.idx = target;
    if (payload.n_pairs) cfg.n_pairs = payload.n_pairs;
    if (payload.pair_label) cfg.pair_label = payload.pair_label;
    if (frames) {
      frames.state.start = 0;
      frames.state.next = null;
      loadFramePage(frames.el, frames.state);
    }
    reactChart(payload.figure);
    syncHistory(target);
    syncNav();
    return updateStructure(payload);
  }

  function reactChart(nextFigure) {
    if (!nextFigure || typeof Plotly === "undefined") return;
    try {
      Plotly.react(chart, nextFigure.data, nextFigure.layout, PLOT_CONFIG);
    } catch (error) {
      console.warn("[dii] chart update failed", error);
    }
  }

  /* 3D 原地更新：未挂载时跳过（页面其它部分仍已切换） */
  function updateStructure(payload) {
    var api = window.DiiStructure;
    if (!api || typeof api.updatePair !== "function") {
      console.warn("[dii] DiiStructure is unavailable; 3D keeps the previous pair");
      return Promise.resolve();
    }
    if (typeof api.isMounted === "function" && !api.isMounted()) {
      console.warn("[dii] DiiStructure is not mounted; 3D keeps the previous pair");
      return Promise.resolve();
    }
    return Promise.resolve(api.updatePair(payload)).catch(function (error) {
      console.warn("[dii] 3D pair update failed", error);
    });
  }

  function syncHistory(target) {
    if (!window.history || typeof history.replaceState !== "function") return;
    var url = "/pair/" + encodeURIComponent(cfg.type) + "/" + target;
    try {
      history.replaceState(null, "", url);
    } catch (error) {
      /* 浏览器限制 replaceState 频率时忽略；页面内数据已更新 */
    }
  }

  /* 首/末 pair 的链接服务端不渲染；这里再兜底同步禁用态 */
  function syncNav() {
    var links = document.querySelectorAll("[data-pair-nav]");
    for (var i = 0; i < links.length; i++) {
      var step = Number(links[i].getAttribute("data-pair-nav"));
      var target = cfg.idx + step;
      if (target < 0 || target >= cfg.n_pairs) {
        links[i].setAttribute("aria-disabled", "true");
      } else {
        links[i].removeAttribute("aria-disabled");
      }
    }
  }

  function toast(text) {
    if (typeof window.showToast === "function") window.showToast(text, true);
  }
})();
