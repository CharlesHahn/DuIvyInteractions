/* 应用内路由（pjax）：只替换 <main>，页面脚本按需重放，viewer 坞常驻不销毁。
   无依赖 ES5 IIFE；未经浏览器端到端验证（见任务报告风险说明）。 */
(function () {
  "use strict";
  var PAGE_SCRIPTS_ID = "pageScripts";
  var MAIN_SELECTOR = "main.content";
  var PJAX_HEADER = { "X-DII-PJAX": "1" };
  var ONCE_ATTR = "data-once";
  var ROUTING_CLASS = "routing";
  var EXCLUDED_PREFIXES = ["/api/", "/download/", "/static/"];
  var afterSwap = [];
  var beforeSwap = [];
  var loadedOnce = {};
  var pending = false;

  function navigate(url, options) {
    var opts = options || {};
    if (pending) return Promise.resolve(false);
    pending = true;
    document.body.classList.add(ROUTING_CLASS);
    return fetch(url, { headers: PJAX_HEADER })
      .then(function (response) {
        if (!response.ok) throw new Error("HTTP " + response.status);
        return response.text();
      })
      .then(function (html) {
        var doc = new DOMParser().parseFromString(html, "text/html");
        runBeforeSwap(url);
        swapMain(doc);
        swapTitle(doc);
        replayPageScripts(doc);
        if (opts.push !== false) history.pushState({}, "", url);
        window.scrollTo(0, 0);
        syncSidebar(url);
        runAfterSwap(doc, url);
        return true;
      })
      .catch(function (error) {
        console.warn("[dii] pjax navigate failed; falling back", error);
        window.location.assign(url);
        return false;
      })
      .then(function (result) {
        pending = false;
        document.body.classList.remove(ROUTING_CLASS);
        return result;
      });
  }

  function swapMain(doc) {
    var target = document.querySelector(MAIN_SELECTOR);
    var incoming = doc.querySelector(MAIN_SELECTOR);
    if (!target || !incoming) throw new Error("missing main.content");
    target.replaceChildren.apply(
      target, Array.prototype.slice.call(incoming.childNodes));
  }

  function swapTitle(doc) {
    var title = doc.querySelector("title");
    if (title) document.title = title.textContent;
  }

  /* 重放 #pageScripts：先清空容器（容器本身保留），再按文档顺序重建。
     src 命中 data-once 且已加载过的脚本跳过（molstar/structure/plotly）。 */
  /* 首屏整页加载的脚本没经过重放：启动时把已有 data-once 脚本登记为已加载，
     否则第二次导航会因容器被清空而重复执行 molstar/structure/plotly。 */
  function seedLoadedScripts() {
    var host = document.getElementById(PAGE_SCRIPTS_ID);
    if (!host) return;
    var scripts = host.querySelectorAll("script[src][" + ONCE_ATTR + "]");
    for (var i = 0; i < scripts.length; i++) {
      var src = scripts[i].getAttribute("src");
      if (src) loadedOnce[src] = true;
    }
  }

  function replayPageScripts(doc) {
    var host = document.getElementById(PAGE_SCRIPTS_ID);
    var incoming = doc.getElementById(PAGE_SCRIPTS_ID);
    if (!host || !incoming) return;
    host.replaceChildren();
    copyStylesheets(incoming, host);
    var scripts = incoming.querySelectorAll("script");
    for (var i = 0; i < scripts.length; i++) {
      replayScript(scripts[i], host);
    }
  }

  function copyStylesheets(incoming, host) {
    var links = incoming.querySelectorAll('link[rel="stylesheet"]');
    for (var i = 0; i < links.length; i++) {
      var href = links[i].getAttribute("href");
      if (!href || hasStylesheet(href)) continue;
      var link = document.createElement("link");
      link.rel = "stylesheet";
      link.href = href;
      host.appendChild(link);
    }
  }

  function hasStylesheet(href) {
    var links = document.querySelectorAll("link");
    for (var i = 0; i < links.length; i++) {
      if (links[i].getAttribute("href") === href) return true;
    }
    return false;
  }

  function replayScript(source, host) {
    var src = source.getAttribute("src");
    var once = source.hasAttribute(ONCE_ATTR);
    if (src && once && (loadedOnce[src] || hasScript(src))) return;
    var script = document.createElement("script");
    var attrs = source.attributes;
    for (var i = 0; i < attrs.length; i++) {
      script.setAttribute(attrs[i].name, attrs[i].value);
    }
    if (src) {
      script.async = false;
      loadedOnce[src] = true;
    } else {
      script.textContent = source.textContent;
    }
    host.appendChild(script);
  }

  function hasScript(src) {
    var scripts = document.querySelectorAll("script[src]");
    for (var i = 0; i < scripts.length; i++) {
      if (scripts[i].getAttribute("src") === src) return true;
    }
    return false;
  }

  /* 交换前钩子：常驻坞等 main 内资源需先摘回 body，避免随 main 一起离场 */
  function runBeforeSwap(url) {
    for (var i = 0; i < beforeSwap.length; i++) {
      try {
        beforeSwap[i](url);
      } catch (error) {
        console.warn("[dii] beforeSwap callback failed", error);
      }
    }
  }

  /* 侧栏在 main 外不随 pjax 替换：按当前路径同步 active 高亮 */
  function syncSidebar(url) {
    var path;
    try {
      path = new URL(url, window.location.href).pathname;
    } catch (error) {
      return;
    }
    var items = document.querySelectorAll(".side-item");
    for (var i = 0; i < items.length; i++) {
      var href = items[i].getAttribute("href");
      if (href === path) items[i].classList.add("active");
      else items[i].classList.remove("active");
    }
  }

  function runAfterSwap(doc, url) {
    for (var i = 0; i < afterSwap.length; i++) {
      try {
        afterSwap[i](doc, url);
      } catch (error) {
        console.warn("[dii] afterSwap callback failed", error);
      }
    }
  }

  function shouldIntercept(link, event) {
    if (event.defaultPrevented) return false;
    if (event.button !== 0) return false;
    if (event.metaKey || event.ctrlKey || event.shiftKey || event.altKey) {
      return false;
    }
    if (!link || !link.getAttribute) return false;
    var href = link.getAttribute("href");
    if (!href || href.charAt(0) === "#") return false;
    if (link.target || link.hasAttribute("download")) return false;
    var url;
    try {
      url = new URL(link.href, window.location.href);
    } catch (error) {
      return false;
    }
    if (url.origin !== window.location.origin) return false;
    for (var i = 0; i < EXCLUDED_PREFIXES.length; i++) {
      if (url.pathname.indexOf(EXCLUDED_PREFIXES[i]) === 0) return false;
    }
    return true;
  }

  function onDocumentClick(event) {
    var link = event.target && event.target.closest
      ? event.target.closest("a[href]") : null;
    if (!shouldIntercept(link, event)) return;
    event.preventDefault();
    navigate(link.href);
  }

  function onPopState() {
    navigate(window.location.pathname + window.location.search, { push: false });
  }

  function currentUrl() {
    return window.location.pathname + window.location.search;
  }

  function onAfterSwap(fn) {
    if (typeof fn === "function") afterSwap.push(fn);
  }

  function onBeforeSwap(fn) {
    if (typeof fn === "function") beforeSwap.push(fn);
  }

  window.DiiRouter = {
    navigate: navigate,
    currentUrl: currentUrl,
    onAfterSwap: onAfterSwap,
    onBeforeSwap: onBeforeSwap
  };

  document.addEventListener("click", onDocumentClick);
  window.addEventListener("popstate", onPopState);

  /* 脚本位于 #pageScripts 之前，立即执行时可能还查不到容器；DOMContentLoaded 再补一次 */
  seedLoadedScripts();
  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", seedLoadedScripts);
  }
})();
