/* 全局 load 弹窗：侧栏 + load 打开、recent 填充、Browse 分派（原生 picker / 服务端浏览） */
(function () {
  "use strict";

  var MODAL_ID = "loadModal";
  var PATH_ID = "loadPath";
  var FORM_ID = "loadForm";
  var CLOSE_ID = "loadClose";
  var BROWSE_ID = "loadBrowse";
  var SERVER_TRIGGER_SELECTOR = "#browseBtn";
  var OPEN_SELECTOR = "[data-load-open]";
  var RECENT_SELECTOR = "[data-load-recent]";
  var PICK_BODY = "mode=dir&kind=dataset";
  var SERVER_MODE = "server";
  var ENTER_KEY = "Enter";
  var ESCAPE_KEY = "Escape";
  var KEY_ENTER_CODE = 13;
  var KEY_ESCAPE_CODE = 27;
  var ELLIPSIS = "…";
  var CHAR_PX = 6.6;          /* 11.5px 等宽字体的近似单字符宽度 */
  var MIN_BUDGET = 24;
  var HEAD_RATIO = 0.6;       /* 保留头部比例，其余给尾部 */
  var PANEL_FALLBACK_WIDTH = 60;
  var PANEL_HORIZONTAL_PADDING = 48;
  var RESIZE_DEBOUNCE_MS = 100;

  var modal = document.getElementById(MODAL_ID);
  if (!modal) return;
  var pathInput = document.getElementById(PATH_ID);
  var form = document.getElementById(FORM_ID);
  var closeBtn = document.getElementById(CLOSE_ID);
  var browseBtn = document.getElementById(BROWSE_ID);

  function showToast(text) {
    if (typeof window.showToast === "function") {
      window.showToast(text);
    }
  }

  function submitForm() {
    if (!form) return;
    if (form.requestSubmit) {
      form.requestSubmit();
      return;
    }
    form.submit();
  }

  function middleEllipsis(text, budget) {
    if (!text || text.length <= budget) return text;
    var keep = budget - 1;
    var head = Math.ceil(keep * HEAD_RATIO);
    var tail = keep - head;
    if (tail < 1) {
      tail = 1;
      head = Math.max(1, keep - tail);
    }
    if (head < 1) {
      head = 1;
      tail = Math.max(1, keep - head);
    }
    return text.slice(0, head) + ELLIPSIS + text.slice(-tail);
  }

  function recentBudget() {
    var panel = document.querySelector(".load-panel");
    var width = panel && panel.clientWidth ? panel.clientWidth
                                           : PANEL_FALLBACK_WIDTH;
    return Math.max(MIN_BUDGET,
                    Math.floor((width - PANEL_HORIZONTAL_PADDING) / CHAR_PX));
  }

  function applyRecentLabels() {
    var nodes = document.querySelectorAll(RECENT_SELECTOR);
    if (!nodes.length) return;
    var budget = recentBudget();
    for (var i = 0; i < nodes.length; i++) {
      var full = nodes[i].getAttribute("data-load-recent") || "";
      nodes[i].textContent = middleEllipsis(full, budget);
    }
  }

  function onResize() {
    if (onResize.timer) clearTimeout(onResize.timer);
    onResize.timer = setTimeout(applyRecentLabels, RESIZE_DEBOUNCE_MS);
  }

  function openModal() {
    modal.classList.remove("hidden");
    if (pathInput) pathInput.focus();
    applyRecentLabels();
  }

  function closeModal() {
    modal.classList.add("hidden");
  }

  function isServerMode() {
    return document.body.dataset.browseMode === SERVER_MODE;
  }

  function lookupNative() {
    fetch('/api/pick', {
      method: "POST",
      headers: {"Content-Type": "application/x-www-form-urlencoded"},
      body: PICK_BODY
    })
      .then(function (resp) { return resp.json(); })
      .then(function (data) {
        if (data.path) {
          if (pathInput) pathInput.value = data.path;
          submitForm();
          return;
        }
        if (data.error) showToast("browse unavailable: " + data.error);
      })
      .catch(function (error) { showToast("browse failed: " + error); });
  }

  function browse() {
    if (!isServerMode()) {
      lookupNative();
      return;
    }
    var trigger = document.querySelector(SERVER_TRIGGER_SELECTOR);
    if (trigger) trigger.click();
  }

  function bindOpeners() {
    var nodes = document.querySelectorAll(OPEN_SELECTOR);
    for (var i = 0; i < nodes.length; i++) {
      nodes[i].addEventListener("click", function (event) {
        event.preventDefault();
        openModal();
      });
    }
  }

  function bindRecent() {
    var nodes = document.querySelectorAll(RECENT_SELECTOR);
    for (var i = 0; i < nodes.length; i++) {
      nodes[i].addEventListener("click", function (event) {
        event.preventDefault();
        if (pathInput) pathInput.value = this.getAttribute("data-load-recent") || "";
        submitForm();
      });
    }
  }

  function onPathKey(event) {
    var isEnter = event.key === ENTER_KEY || event.keyCode === KEY_ENTER_CODE;
    if (!isEnter) return;
    event.preventDefault();
    submitForm();
  }

  function onKey(event) {
    var isEscape = event.key === ESCAPE_KEY || event.keyCode === KEY_ESCAPE_CODE;
    if (!isEscape || modal.classList.contains("hidden")) return;
    event.preventDefault();
    closeModal();
  }

  /* 入口链接为 /?add=1：打开页面即展开弹窗（服务端 ?add=1 已由页面链接给出） */
  var params = new URLSearchParams(window.location.search);
  var wantsOpen = params.get("add") === "1";

  bindOpeners();
  bindRecent();
  applyRecentLabels();
  if (closeBtn) closeBtn.addEventListener("click", closeModal);
  if (browseBtn) browseBtn.addEventListener("click", browse);
  if (pathInput) pathInput.addEventListener("keydown", onPathKey);
  modal.addEventListener("click", function (event) {
    if (event.target === modal) closeModal();
  });
  document.addEventListener("keydown", onKey);
  window.addEventListener("resize", onResize);
  if (wantsOpen) openModal();
})();
