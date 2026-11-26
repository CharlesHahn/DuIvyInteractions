/* aggregates 表格按需加载：fetch 下一批行并追加，不整页刷新 */
(function () {
  "use strict";
  var cfgEl = document.getElementById("rowsConfig");
  var body = document.getElementById("rowsBody");
  if (!cfgEl || !body) return;
  var cfg = JSON.parse(cfgEl.textContent);
  var status = document.getElementById("rowsStatus");
  var button = document.getElementById("rowsMoreBtn");
  var info = document.getElementById("rowsInfo");
  var limit = parseInt(cfg.limit, 10);
  if (!isFinite(limit) || limit < 1) limit = 50;
  var total = parseInt(cfg.total, 10);
  if (!isFinite(total) || total < 0) total = 0;
  // 服务端可能经 ?page=N 首屏渲染了非第一批行，offset 以渲染起点为基准
  var baseOffset = parseInt(cfg.offset, 10);
  if (!isFinite(baseOffset) || baseOffset < 0) baseOffset = 0;
  var baseUrl = cfg.url ||
    "/api/rows/" + encodeURIComponent(cfg.type || "");
  var rendered = body.querySelectorAll("tr").length;
  var loading = false;

  function rowsUrl() {
    return baseUrl + "?search=" + encodeURIComponent(cfg.search || "") +
      "&min_occupancy=" + (cfg.min_occupancy || 0) +
      "&sort_by=" + encodeURIComponent(cfg.sort_by || "occupancy") +
      "&sort_dir=" + encodeURIComponent(cfg.sort_dir || "desc") +
      "&offset=" + (baseOffset + rendered) + "&limit=" + limit;
  }

  function readJson(resp) {
    if (!resp.ok) throw new Error("HTTP " + resp.status);
    return resp.json();
  }

  function setStatus(text, failed) {
    if (!status) return;
    status.textContent = text;
    status.classList.toggle("error", Boolean(failed));
  }

  function hideButton() {
    if (button) button.classList.add("hidden");
  }

  function updateProgress() {
    if (info) {
      info.textContent = "showing " + (baseOffset + rendered) +
        " of " + total;
    }
    if (rendered >= total) {
      hideButton();
      setStatus("all " + total + " rows loaded", false);
      return;
    }
    if (button) {
      button.classList.remove("hidden");
      button.textContent = "load more (" + (total - rendered) + " remaining)";
    }
    setStatus("", false);
  }

  function textCell(text) {
    var td = document.createElement("td");
    td.textContent = text;
    return td;
  }

  function pairCell(entry) {
    var td = document.createElement("td");
    var link = document.createElement("a");
    link.className = "pair-link";
    link.href = "/pair/" + encodeURIComponent(cfg.type) + "/" +
      encodeURIComponent(entry.idx);
    link.textContent = entry.label;
    td.appendChild(link);
    return td;
  }

  function buildRow(entry, columns) {
    var tr = document.createElement("tr");
    entry.values.forEach(function (value, position) {
      if (columns[position] === "pair_label") {
        tr.appendChild(pairCell(entry));
      } else {
        tr.appendChild(textCell(value));
      }
    });
    return tr;
  }

  function appendRows(data) {
    var columns = data.columns || [];
    (data.rows || []).forEach(function (entry) {
      body.appendChild(buildRow(entry, columns));
    });
    rendered = body.querySelectorAll("tr").length;
    total = parseInt(data.total, 10);
    if (!isFinite(total)) total = rendered;
    updateProgress();
  }

  function loadMore() {
    if (loading) return;
    loading = true;
    setStatus("loading…", false);
    if (button) button.disabled = true;
    fetch(rowsUrl())
      .then(readJson)
      .then(function (data) {
        loading = false;
        if (button) button.disabled = false;
        if (!data.rows || !data.rows.length) {
          hideButton();
          setStatus("no more rows", false);
          return;
        }
        appendRows(data);
      })
      .catch(function (error) {
        loading = false;
        if (button) button.disabled = false;
        setStatus("failed to load rows: " + error, true);
      });
  }

  if (button) button.addEventListener("click", loadMore);
  updateProgress();
})();
