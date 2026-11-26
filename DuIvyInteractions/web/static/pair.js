/* 单对页：Plotly 对齐曲线 + 逐帧明细表分页 */
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

  renderChart();
  renderFrames();

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
    var promise = Plotly.newPlot(chart, figure.data, figure.layout, {
      responsive: true,
      displaylogo: false,
      scrollZoom: false,
      modeBarButtonsToRemove: ["select2d", "lasso2d"]
    });
    if (promise && promise.catch) {
      promise.catch(function (error) {
        chart.textContent = "failed to render chart: " + error;
      });
    }
  }

  function renderFrames() {
    var el = frameElements();
    var state = {start: 0, next: null, pageSize: 1};
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
    el.table.innerHTML = "";
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
})();
