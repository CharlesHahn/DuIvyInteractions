/* 服务端目录浏览弹窗（结构文件）：/api/list 列目录 → 点击文件回调；无弹窗 DOM 时静默退出 */
(function () {
  "use strict";
  var ENDPOINT = "/api/list?kind=structure&path=";
  var SUFFIXES = [".tpr", ".xtc", ".trr", ".gro", ".pdb", ".dcd", ".nc",
    ".itp", ".psf", ".prmtop", ".parm7", ".mol2"];
  var TRUNCATED_TEXT = "list truncated (limit reached)";
  var LOADING_TEXT = "loading…", EMPTY_TEXT = "no sub-folders or files here";
  var FAILED_TEXT = "could not load this folder";
  var NO_FILE_TEXT = "select a structure file (click one in the list)";
  var ICONS = {up: "↑", dir: "▸", file: "·"}, UP_LABEL = "..";
  var SIZE_DIGITS = 1, SIZE_SUFFIX = " KB";
  var ENTER_KEY = "Enter", ESCAPE_KEY = "Escape";
  var ENTER_CODE = 13, ESCAPE_CODE = 27;
  var modal = document.getElementById("filepickModal");
  if (!modal) return;
  var listEl = byId("filepickList"), errorEl = byId("filepickError");
  var noteEl = byId("filepickNote"), pathInput = byId("filepickPath");
  var currentEl = byId("filepickCurrent"), upBtn = byId("filepickUp");
  var goBtn = byId("filepickGo"), useBtn = byId("filepickUse");
  var closeBtn = byId("filepickClose"), cancelBtn = byId("filepickCancel");
  if (!listEl || !errorEl || !noteEl || !pathInput || !currentEl ||
      !upBtn || !goBtn || !useBtn || !closeBtn || !cancelBtn) return;
  var state = {path: "", parent: null, loading: false, listed: false}, onPick = null;
  function byId(id) { return document.getElementById(id); }
  function el(tag, className, text) {
    var node = document.createElement(tag);
    node.className = className;
    node.textContent = text;
    return node;
  }
  function clear(node) { while (node.firstChild) node.removeChild(node.firstChild); }
  function setText(node, text) {
    node.textContent = text || "";
    node.classList.toggle("hidden", !text);
  }
  function note(text) { return el("p", "filepick-message", text); }
  function sizeText(sizeKb) {
    return typeof sizeKb === "number" && isFinite(sizeKb)
      ? sizeKb.toFixed(SIZE_DIGITS) + SIZE_SUFFIX : "";
  }
  function row(icon, label, size, onClick) {
    var button = el("button", "filepick-row", "");
    button.type = "button";
    button.appendChild(el("span", "filepick-icon", icon));
    button.appendChild(el("span", "filepick-name", label));
    if (size) button.appendChild(el("span", "filepick-size", size));
    button.addEventListener("click", onClick);
    return button;
  }
  function joinPath(base, name) {
    if (!base) return name;
    var last = base.charAt(base.length - 1);
    if (last === "/" || last === "\\") return base + name;
    return base + (base.indexOf("/") >= 0 ? "/" : "\\") + name;
  }
  function isStructureFile(name) {
    var lower = (name || "").toLowerCase();
    for (var i = 0; i < SUFFIXES.length; i++) {
      var at = lower.length - SUFFIXES[i].length;
      if (at > 0 && lower.slice(at) === SUFFIXES[i]) return true;
    }
    return false;
  }
  function updateControls() {
    upBtn.disabled = state.loading || !state.parent;
    goBtn.disabled = state.loading;
    useBtn.disabled = state.loading || !isStructureFile(pathInput.value.trim());
  }
  function pick(path) {
    if (!path) {
      setText(errorEl, NO_FILE_TEXT);
      return;
    }
    var callback = onPick;
    closeModal();
    if (callback) callback(path);
  }
  function renderListing(data) {
    clear(listEl);
    if (state.parent) listEl.appendChild(row(ICONS.up, UP_LABEL, "",
      function () { loadDir(state.parent); }));
    (data.dirs || []).forEach(function (name) {
      listEl.appendChild(row(ICONS.dir, name, "", function () {
        loadDir(joinPath(state.path, name));
      }));
    });
    (data.files || []).forEach(function (file) {
      if (!isStructureFile(file.name)) return;
      listEl.appendChild(row(ICONS.file, file.name, sizeText(file.size_kb),
        function () { pick(joinPath(state.path, file.name)); }));
    });
    if (!listEl.children.length) listEl.appendChild(note(EMPTY_TEXT));
  }
  function requestList(path) {
    return fetch(ENDPOINT + encodeURIComponent(path)).then(function (resp) {
      return resp.json().then(function (data) {
        return {ok: resp.ok, status: resp.status, data: data};
      }, function () { return {ok: false, status: resp.status, data: null}; });
    }, function () { return {ok: false, status: 0, data: null}; });
  }
  function loadDir(path) {
    if (state.loading) return;
    state.loading = true;
    state.listed = false;
    updateControls();
    setText(errorEl, "");
    setText(noteEl, "");
    clear(listEl);
    listEl.appendChild(note(LOADING_TEXT));
    requestList(path).then(function (result) {
      state.loading = false;
      if (!result.ok || !result.data) {
        pathInput.value = state.path;
        clear(listEl);
        listEl.appendChild(note(FAILED_TEXT));
        var data = result.data;
        setText(errorEl, data && data.detail ? String(data.detail)
          : result.status ? "HTTP " + result.status : "browse failed");
        updateControls();
        return;
      }
      var listing = result.data;
      state.path = listing.path || "";
      state.parent = listing.parent || null;
      state.listed = true;
      pathInput.value = state.path;
      currentEl.textContent = listing.rel || state.path;
      renderListing(listing);
      setText(noteEl, listing.truncated ? TRUNCATED_TEXT : "");
      updateControls();
    });
  }
  function closeModal() {
    modal.classList.add("hidden");
    setText(errorEl, "");
    setText(noteEl, "");
  }
  function open(options) {
    onPick = options && options.onPick ? options.onPick : null;
    modal.classList.remove("hidden");
    closeBtn.focus();
    if (!state.listed && !state.loading) loadDir("");
  }
  function loadPath() { loadDir(pathInput.value); }
  closeBtn.addEventListener("click", closeModal);
  cancelBtn.addEventListener("click", closeModal);
  upBtn.addEventListener("click", function () {
    if (state.parent) loadDir(state.parent);
  });
  goBtn.addEventListener("click", loadPath);
  useBtn.addEventListener("click", function () { pick(pathInput.value.trim()); });
  pathInput.addEventListener("input", updateControls);
  pathInput.addEventListener("keydown", function (event) {
    if (event.key !== ENTER_KEY && event.keyCode !== ENTER_CODE) return;
    event.preventDefault();
    loadPath();
  });
  modal.addEventListener("click", function (event) {
    if (event.target === modal) closeModal();
  });
  document.addEventListener("keydown", function (event) {
    var escape = event.key === ESCAPE_KEY || event.keyCode === ESCAPE_CODE;
    if (!escape || modal.classList.contains("hidden")) return;
    event.preventDefault();
    closeModal();
  });
  updateControls();
  window.DiiFilePick = {open: open};
})();
