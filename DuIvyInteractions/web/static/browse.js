/* 服务端目录浏览弹窗：/api/list 列目录 → 进入目录 / 选中 .h5 → POST /select 跳转；无弹窗节点时静默退出 */
(function () {
  "use strict";

  var MODAL_ID = "browseModal";
  var OPEN_BUTTON_ID = "browseBtn";
  var LIST_ENDPOINT = "/api/list";
  var SELECT_ENDPOINT = "/select";
  var PATH_FIELD_SELECTOR = '.source-form input[name="path"]';
  var H5_SUFFIX = ".h5";
  var TRUNCATED_TEXT = "list truncated (limit reached)";
  var LOADING_TEXT = "loading…";
  var EMPTY_TEXT = "no sub-folders or .h5 files here";
  var FAILED_TEXT = "could not load this folder";
  var UP_ICON = "↑";
  var DIR_ICON = "▸";
  var FILE_ICON = "·";
  var UP_LABEL = "..";
  var SIZE_DIGITS = 1;
  var SIZE_SUFFIX = " KB";
  var ENTER_KEY = "Enter";
  var ESCAPE_KEY = "Escape";
  var KEY_ENTER_CODE = 13;
  var KEY_ESCAPE_CODE = 27;

  var modal = document.getElementById(MODAL_ID);
  if (!modal) return;
  var listEl = document.getElementById("browseList");
  var errorEl = document.getElementById("browseError");
  var noteEl = document.getElementById("browseNote");
  var pathInput = document.getElementById("browsePathInput");
  var currentEl = document.getElementById("browseCurrent");
  var upBtn = document.getElementById("browseUp");
  var goBtn = document.getElementById("browseGo");
  var useBtn = document.getElementById("browseUse");
  var closeBtn = document.getElementById("browseClose");
  if (!listEl || !errorEl || !noteEl || !pathInput || !currentEl ||
      !upBtn || !goBtn || !useBtn || !closeBtn) {
    return;
  }

  var state = {path: "", rel: "", parent: null, loading: false, listed: false};

  function toggleText(node, text) {
    node.textContent = text || "";
    if (text) {
      node.classList.remove("hidden");
    } else {
      node.classList.add("hidden");
    }
  }

  function showError(text) {
    toggleText(errorEl, text);
  }

  function showNote(text) {
    toggleText(noteEl, text);
  }

  function clearNode(node) {
    while (node.firstChild) node.removeChild(node.firstChild);
  }

  function messageRow(text) {
    var row = document.createElement("p");
    row.className = "browse-message";
    row.textContent = text;
    return row;
  }

  function buildRow(icon, label, sizeText) {
    var row = document.createElement("button");
    row.type = "button";
    row.className = "browse-row";
    var iconEl = document.createElement("span");
    iconEl.className = "browse-icon";
    iconEl.textContent = icon;
    var nameEl = document.createElement("span");
    nameEl.className = "browse-name";
    nameEl.textContent = label;
    row.appendChild(iconEl);
    row.appendChild(nameEl);
    if (sizeText) {
      var sizeEl = document.createElement("span");
      sizeEl.className = "browse-size";
      sizeEl.textContent = sizeText;
      row.appendChild(sizeEl);
    }
    return row;
  }

  function markSelected(row) {
    var rows = listEl.querySelectorAll(".browse-row");
    for (var i = 0; i < rows.length; i++) {
      rows[i].classList.remove("selected");
    }
    row.classList.add("selected");
  }

  function pathSeparator(base) {
    if (base.indexOf("/") >= 0) return "/";
    if (base.indexOf("\\") >= 0) return "\\";
    return "/";
  }

  function joinPath(base, name) {
    if (!base) return name;
    var lastChar = base.charAt(base.length - 1);
    if (lastChar === "/" || lastChar === "\\") return base + name;
    return base + pathSeparator(base) + name;
  }

  function isH5File(name) {
    var suffixAt = name.length - H5_SUFFIX.length;
    return suffixAt >= 0 && name.slice(suffixAt).toLowerCase() === H5_SUFFIX;
  }

  function formatSize(sizeKb) {
    if (typeof sizeKb !== "number" || !isFinite(sizeKb)) return "";
    return sizeKb.toFixed(SIZE_DIGITS) + SIZE_SUFFIX;
  }

  function pathField() {
    return document.querySelector(PATH_FIELD_SELECTOR);
  }

  function submitPath(path) {
    var field = pathField();
    if (field) field.value = path;
    var form = document.createElement("form");
    form.method = "post";
    form.action = SELECT_ENDPOINT;
    form.acceptCharset = "UTF-8";
    var input = document.createElement("input");
    input.type = "hidden";
    input.name = "path";
    input.value = path;
    form.appendChild(input);
    document.body.appendChild(form);
    form.submit();
  }

  function updateControls() {
    upBtn.disabled = state.loading || !state.parent;
    goBtn.disabled = state.loading;
    useBtn.disabled = state.loading || !state.listed;
  }

  function renderListing(data) {
    clearNode(listEl);
    var dirs = data.dirs || [];
    var files = data.files || [];
    var rows = 0;
    if (state.parent) {
      listEl.appendChild(buildUpRow());
      rows += 1;
    }
    for (var i = 0; i < dirs.length; i++) {
      listEl.appendChild(buildDirRow(dirs[i]));
      rows += 1;
    }
    for (var j = 0; j < files.length; j++) {
      if (!isH5File(files[j].name)) continue;
      listEl.appendChild(buildFileRow(files[j]));
      rows += 1;
    }
    if (!rows) listEl.appendChild(messageRow(EMPTY_TEXT));
  }

  function buildUpRow() {
    var row = buildRow(UP_ICON, UP_LABEL, "");
    row.addEventListener("click", function () {
      markSelected(row);
      loadDir(state.parent);
    });
    return row;
  }

  function buildDirRow(name) {
    var row = buildRow(DIR_ICON, name, "");
    row.addEventListener("click", function () {
      markSelected(row);
      loadDir(joinPath(state.path, name));
    });
    return row;
  }

  function buildFileRow(file) {
    var row = buildRow(FILE_ICON, file.name, formatSize(file.size_kb));
    row.addEventListener("click", function () {
      markSelected(row);
      submitPath(joinPath(state.path, file.name));
    });
    return row;
  }

  function readJson(resp) {
    return resp.json().then(function (data) {
      return {ok: resp.ok, status: resp.status, data: data, error: null};
    }, function () {
      return {ok: false, status: resp.status, data: null, error: null};
    });
  }

  function requestList(path) {
    var url = LIST_ENDPOINT + "?path=" + encodeURIComponent(path);
    return fetch(url).then(readJson, function (error) {
      return {ok: false, status: 0, data: null, error: error};
    });
  }

  function failureMessage(result) {
    if (result.error) return "browse failed: " + result.error;
    if (result.data && result.data.detail) return String(result.data.detail);
    return "HTTP " + result.status;
  }

  function applyListing(data) {
    state.path = data.path || "";
    state.rel = data.rel || "";
    state.parent = data.parent || null;
    state.listed = true;
    pathInput.value = state.path;
    currentEl.textContent = state.rel || state.path;
    renderListing(data);
    showNote(data.truncated ? TRUNCATED_TEXT : "");
    updateControls();
  }

  function showFailure(result) {
    state.listed = false;
    pathInput.value = state.path;
    currentEl.textContent = state.rel || state.path;
    clearNode(listEl);
    listEl.appendChild(messageRow(FAILED_TEXT));
    showError(failureMessage(result));
    updateControls();
  }

  function loadDir(path) {
    if (state.loading) return;
    state.loading = true;
    state.listed = false;
    updateControls();
    showError("");
    showNote("");
    clearNode(listEl);
    listEl.appendChild(messageRow(LOADING_TEXT));
    requestList(path).then(function (result) {
      state.loading = false;
      if (result.ok && result.data) {
        applyListing(result.data);
      } else {
        showFailure(result);
      }
    });
  }

  function isOpen() {
    return !modal.classList.contains("hidden");
  }

  function openModal() {
    var field = pathField();
    var initialPath = field && field.value ? field.value : "";
    modal.classList.remove("hidden");
    closeBtn.focus();
    loadDir(initialPath);
  }

  function closeModal() {
    modal.classList.add("hidden");
    showError("");
    showNote("");
  }

  function onDocumentKey(event) {
    var isEscape = event.key === ESCAPE_KEY || event.keyCode === KEY_ESCAPE_CODE;
    if (!isEscape || !isOpen()) return;
    event.preventDefault();
    closeModal();
  }

  function onPathKey(event) {
    var isEnter = event.key === ENTER_KEY || event.keyCode === KEY_ENTER_CODE;
    if (!isEnter) return;
    event.preventDefault();
    loadDir(pathInput.value);
  }

  function onModalClick(event) {
    if (event.target === modal) closeModal();
  }

  function onUpClick() {
    if (state.parent) loadDir(state.parent);
  }

  function onGoClick() {
    loadDir(pathInput.value);
  }

  function onUseClick() {
    if (state.listed && state.path) submitPath(state.path);
  }

  var openButtons = document.querySelectorAll("#" + OPEN_BUTTON_ID);
  for (var i = 0; i < openButtons.length; i++) {
    openButtons[i].addEventListener("click", openModal);
  }
  closeBtn.addEventListener("click", closeModal);
  upBtn.addEventListener("click", onUpClick);
  goBtn.addEventListener("click", onGoClick);
  useBtn.addEventListener("click", onUseClick);
  pathInput.addEventListener("keydown", onPathKey);
  modal.addEventListener("click", onModalClick);
  document.addEventListener("keydown", onDocumentKey);
  updateControls();
})();
