/* 3D 结构页：Mol* 原生轨迹（拓扑 PDB + 原始轨迹文件）只载入一次；
   帧切换只更新 modelIndex（不重建结构、不动相机），相互作用残基以棍棒呈现。 */
(function () {
  "use strict";

  /* pjax 重放/重复引入时复用既有实例：绝不创建第二个 viewer */
  if (window.DiiStructure && typeof window.DiiStructure.mount === "function") {
    return;
  }

  var cfgEl = document.getElementById("structureConfig");
  if (!cfgEl) return;
  var cfg = JSON.parse(cfgEl.textContent);

  var host = document.getElementById("molstarViewer");
  /* 单对页曲线视图不引入本条例；退化页面里没有 viewer 时安静退出 */
  if (!host) return;
  var form = document.getElementById("structureForm");
  var topologyInput = document.getElementById("topologyInput");
  var trajectoryInput = document.getElementById("trajectoryInput");
  var openButton = document.getElementById("openStructure");
  var frameLabel = document.getElementById("frameLabel");
  var existenceBar = document.getElementById("existenceBar");
  var statusEl = document.getElementById("structureStatus");
  var groupLegendEl = document.getElementById("groupLegend");
  var browseButtons = document.querySelectorAll("button[data-browse]");

  var TOPOLOGY_URL = "/api/structure/topology";
  var TRAJECTORY_URL = "/api/structure/trajectory";
  var INTERACTIONS_URL = "/api/structure/interactions";
  var SELECT_URL = "/api/structure/select";
  /* h5 的 run/ 组自动载入：auto_attempt 为真时先试后端，失败才回退手选 */
  var AUTO_URL_KEY = "auto_url";
  var AUTO_TRIED_FLAG = "auto";
  var AUTO_SKIP_VALUE = "0";
  var AUTO_CONFLICT_STATUS = 409;
  var AUTO_LOADING_TEXT = "loading h5 recorded structure…";
  var AUTO_UNAVAILABLE_TEXT = "h5 recorded trajectory unavailable";
  var AUTO_CONFIRM_HINT = " — confirm the paths";
  var AUTO_FAILED_PREFIX = "auto structure failed: HTTP ";
  /* 会话版本号：结构每次重选都自增，URL 变则缓存键变（绕开浏览器/代理缓存） */
  var VERSION_PARAM = "v";
  var PICK_ENDPOINT = "/api/pick";
  var PICK_BODY = "mode=file&kind=structure";

  var TOPOLOGY_FORMAT = "mmcif";
  var MODEL_KIND = "model-data";
  var DEFAULT_TRAJECTORY_FORMAT = "xtc";
  var BINARY_FORMATS = { xtc: true, trr: true, dcd: true };
  /* 拓扑与轨迹原子数不一致：多半命中了上一次会话的轨迹缓存，提示重选轨迹 */
  var MISMATCH_MARKER = "element count mismatch";
  var MISMATCH_HINT = " (topology and trajectory atom counts differ — " +
    "re-select the trajectory)";

  /* 后端给的表达式是 PyMOL 风格选择语言；Mol* 只有 script 类型的组件会按
     language 解析字符串（表达式类型需要编译好的 AST，字符串会被当常量）。 */
  var EXPRESSION_LANGUAGE = "pymol";
  var COMPONENT_SCRIPT = "script";
  /* 背景层用 Mol* 内置 static 组件关键字（polymer/ligand/ion），
     不拼 protein/water/nucleic 选择语言表达式 */
  var STICK_TYPE = "ball-and-stick";
  var CARTOON_TYPE = "cartoon";

  /* 相互作用残基只有一层棍棒：sizeFactor 与 Mol* 原生 click-focus 目标残基
     的 ball-and-stick（sizeFactor:.22 / sizeAspectRatio:.73）取同量级，
     颜色不指定 → 用 ball-and-stick 的默认主题 element-symbol（元素配色） */
  var RESIDUE_STICK_SIZE = 0.22;
  var STICK_ASPECT = 1.0;           /* 1.0 = 棍棒（无球） */
  var LIGAND_SIZE = 0.28;           /* 配体/离子球棍 */
  var DASH_LENGTH = 0.12;
  var LINE_SIZE = 0.06;
  var DEBOUNCE_MS = 60;             /* 滑块拖动到换帧的防抖 */
  var FRAME_POLL_MS = 150;          /* 轮询 Mol* 原生轨迹控件的当前帧 */
  /* loadTrajectory resolve 与 hierarchy 更新之间存在竞态：轮询等结构树出现 */
  var STRUCTURE_WAIT_INTERVAL_MS = 100;
  var STRUCTURE_WAIT_TIMEOUT_MS = 5000;
  var DIAGNOSTIC_HEAD_CHARS = 40;    /* 诊断信息里展示的拓扑开头字符数 */
  var FRAME_PARAM = "frame";        /* 地址栏帧号参数（视图之间贯通） */

  var TAG_POLYMER = "dii-polymer";
  var TAG_LIGAND = "dii-ligand";
  var TAG_RESIDUE = "dii-residue-";
  var TAG_FOCUS = "dii-focus";
  var TAG_FOCUS_ITEM = "dii-focus-";
  var TAG_MEASURE = "dii-measure-";

  var LEGEND_ATOM_LIMIT = 12;      /* 图例里最多列出的原子名个数，其余用 (+N) */
  var GROUP_ITEM_CLASS = "structure-group-item";
  var GROUP_LABEL_CLASS = "structure-group-label";
  var LEGEND_CLASS = "pt-legend";
  var LEGEND_ITEM_CLASS = "pt-legend-item";
  var LEGEND_DOT_CLASS = "pt-dot";
  /* 与后端 figures.PT_COLORS / 服务端模板中的 P/T 色一致（写 rgb 以避免自造色板） */
  var PT_LEGEND_TYPES = ["P", "T"];
  var PT_DOT_COLORS = { P: "rgb(60, 84, 136)", T: "rgb(230, 75, 53)" };
  var STATUS_CLASS = "structure-status";
  var STATUS_ERROR_CLASS = "structure-status error";
  var BROWSE_TARGET_ATTR = "data-target";
  var BROWSE_MODE_ATTR = "data-browse";
  var BROWSE_PICKER = "picker";
  var BROWSE_SERVER = "server";
  var MODEL_INDEX_PARAM = "modelIndex";
  var MODEL_INDEX_LABEL = "Model Index";
  var MODEL_TRANSFORM_NAME = "model-from-trajectory";

  /* 尽量贴近 Mol* 默认：开关全部用原生按钮/原生设置项，本页不自造控件 */
  var VIEWER_OPTIONS = {
    layoutIsExpanded: false,
    layoutShowControls: true,
    layoutShowLeftPanel: true,
    /* 左侧原生 32px 图标条（图标按钮原生保留） */
    collapseLeftPanel: true,
    /* 右侧保持可见：原生 Controls→Settings→Layout 多选（Sequence/Log/Left/Right）
       是 Sequence 的官方开关，面板不可达用户就没法用它 */
    collapseRightPanel: false,
    layoutShowSequence: true,
    layoutShowLog: false,
    layoutShowRemoteState: false,
    /* 原生视口按钮：Toggle Controls Panel / Settings & Controls Info */
    viewportShowControls: true,
    viewportShowSettings: true,
    viewportShowSelectionMode: false,
    /* 原生动画/播放按钮（轨迹逐帧播放）：Mol* 默认开启，显式打开 */
    viewportShowAnimation: true,
    viewportShowTrajectoryControls: true
  };

  var viewer = null;
  var ready = false;
  var trajectory = null;        /* loadTrajectory 的返回值 {model, coords, preset} */
  var structureRef = null;      /* 主结构 cell 的 ref：帧更新后 cell 对象会换，始终按 ref 取 */
  var modelRef = null;          /* model-from-trajectory 的 ref（改帧用） */
  var data = null;              /* /api/structure/interactions 的 JSON */
  var topologySnapshot = "";    /* 最近一次通过校验的拓扑文本（超时诊断用） */
  var frame = clampFrame(cfg.frame || 0);  /* 后端注入的初始帧（缺省 0） */
  var measurements = [];        /* 当前帧的虚线测量 [{selection, representation}] */
  var residueRefs = [];         /* 当前 pair 的残基棍棒组件 ref（updatePair 整批替换） */
  var backgroundRefs = [];      /* 卡通/配体/离子背景组件 ref（updatePair 一并重建） */
  var mounted = false;          /* mount 幂等标记 */
  var currentStructure = { topology: "", trajectory: "" };  /* 已载入结构指纹 */
  var focusPath = "";           /* 状态栏：focus: residues|residues-each|structure */
  var measureSeq = 0;
  var frameSeq = 0;             /* 已请求的换帧序号（滑块） */
  var appliedSeq = 0;           /* 已完成的换帧序号（轮询据此避让） */
  var failures = 0;
  var pendingTimer = null;
  var pollTimer = null;
  var dragging = false;         /* existence 条按住拖动中 */

  function setStatus(text, isError) {
    if (!statusEl) return;
    statusEl.textContent = text || "";
    statusEl.className = isError ? STATUS_ERROR_CLASS : STATUS_CLASS;
  }

  function noteFailure(label, error) {
    failures += 1;
    console.warn("[dii] " + label, error);
  }

  /* 同一表达式/组件失败只提示一次：重建、换帧不刷屏 */
  var warnedOnce = new Set();

  function warnOnce(key) {
    if (warnedOnce.has(key)) return;
    warnedOnce.add(key);
    console.warn.apply(console, Array.prototype.slice.call(arguments, 1));
  }

  function readJson(response) {
    return response.json().then(function (payload) {
      if (!response.ok) {
        var detail = payload && (payload.detail || payload.error)
          ? (payload.detail || payload.error) : "HTTP " + response.status;
        throw new Error(String(detail));
      }
      return payload;
    });
  }

  /* ---------------- 帧范围与状态栏 ---------------- */

  function existenceValues() {
    if (data && data.existence) return data.existence;
    return cfg.existence || [];
  }

  function typeValues() {
    var types = data && data.pt_types ? data.pt_types : (cfg.pt_types || []);
    var values = existenceValues();
    if (values.length && types.length !== values.length) return [];
    return types;
  }

  function frameLimit() {
    var total = data && typeof data.n_frames === "number"
      ? data.n_frames : (cfg.n_frames || 0);
    var values = data && data.existence ? data.existence : null;
    if (values && values.length && (!total || values.length < total)) {
      return values.length;
    }
    return total;
  }

  function clampFrame(value) {
    var limit = frameLimit();
    var number = parseInt(value, 10);
    if (isNaN(number) || number < 0) return 0;
    if (limit && number > limit - 1) return limit - 1;
    return number;
  }

  function updateFrameLabel(value) {
    if (!frameLabel) return;
    var limit = frameLimit();
    frameLabel.textContent = limit
      ? "frame " + value + " / " + limit
      : "no structure";
  }

  function highlightFrame(value) {
    if (!existenceBar) return;
    var values = existenceValues();
    var types = typeValues();
    var bars = existenceBar.children;
    for (var i = 0; i < bars.length; i++) {
      var active = !!values[i];
      var name = active ? types[i] : "";
      bars[i].className = (active ? "on" : "") +
        (name ? " pt-" + String(name).toLowerCase() : "") +
        (i === value ? " active" : "");
    }
    existenceBar.setAttribute("aria-valuenow", String(value));
  }

  function renderExistence(values, types) {
    if (!existenceBar || !values || !values.length) return;
    if (types && types.length !== values.length) types = [];
    while (existenceBar.firstChild) {
      existenceBar.removeChild(existenceBar.firstChild);
    }
    for (var i = 0; i < values.length; i++) {
      var bar = document.createElement("span");
      var active = !!values[i];
      var name = active && types ? types[i] : "";
      bar.className = (active ? "on" : "") +
        (name ? " pt-" + String(name).toLowerCase() : "");
      bar.title = "frame " + i + ": " +
        (active ? "active" : "inactive") +
        (name ? " (" + name + ")" : "");
      existenceBar.appendChild(bar);
    }
  }

  function updateFrameRange() {
    var limit = frameLimit();
    if (existenceBar && limit > 0) {
      existenceBar.setAttribute("aria-valuemax", String(limit - 1));
    }
    frame = clampFrame(frame);
    updateFrameLabel(frame);
    highlightFrame(frame);
  }

  /* 换帧后把帧号写回地址栏（保留 view=3d 等其它参数），刷新/切换视图不丢帧 */
  function syncFrameUrl(value) {
    if (!window.history || typeof window.history.replaceState !== "function") {
      return;
    }
    var url = new URL(window.location.href);
    if (url.searchParams.get(FRAME_PARAM) === String(value)) return;
    url.searchParams.set(FRAME_PARAM, String(value));
    try {
      window.history.replaceState(null, "",
        url.pathname + url.search + url.hash);
    } catch (error) {
      /* 浏览器对 replaceState 有频率限制：写不进去不影响换帧 */
    }
  }

  function statusText() {
    var parts = ["frame " + frame + " / " + frameLimit(),
                 "trajectory loaded"];
    var source = runSourceLabel();
    if (source) parts.push(source);
    var residues = data && data.residues ? data.residues : [];
    if (residues.length) parts.push(residues.length + " groups");
    if (focusPath) parts.push("focus: " + focusPath);
    if (failures) parts.push(failures + " warnings");
    return parts.join(" · ");
  }

  /* 已载入结构且 h5 带 run/ 记录时，状态栏标出来源（只提示，不重复载入） */
  function runSourceLabel() {
    var meta = cfg.run_meta;
    if (!cfg.has_structure || !meta) return "";
    if (!meta.topology || !meta.trajectory) return "";
    return "source: run (" + basename(meta.topology) + " + " +
      basename(meta.trajectory) + ")";
  }

  function basename(path) {
    var parts = String(path || "").split("/");
    return parts[parts.length - 1];
  }

  /* ---------------- 数据请求 ---------------- */

  function interactionUrl() {
    return INTERACTIONS_URL + "?type=" + encodeURIComponent(cfg.type) +
      "&idx=" + cfg.idx;
  }

  /* 拓扑 mmCIF 由后端按该对原子打 B-factor 标记（表达式据此圈选），必须带 type/idx */
  function topologyUrl() {
    return TOPOLOGY_URL + "?type=" + encodeURIComponent(cfg.type) +
      "&idx=" + cfg.idx + "&" + VERSION_PARAM + "=" +
      encodeURIComponent(cfg.structure_version || 0);
  }

  /* 原始轨迹文件同样按会话版本号取：URL 相同会复用上一次会话缓存的轨迹 */
  function trajectoryUrl() {
    return TRAJECTORY_URL + "?" + VERSION_PARAM + "=" +
      encodeURIComponent(cfg.structure_version || 0);
  }

  /* 后端 400 的 detail 才是真实原因（如原子名/原子数不一致），状态栏直接展示 */
  function fetchTopologyText() {
    return fetch(topologyUrl()).then(function (response) {
      if (!response.ok) {
        return response.json().then(function (payload) {
          var detail = payload && (payload.detail || payload.error);
          throw new Error(detail ? String(detail)
            : "topology HTTP " + response.status);
        }, function () {
          throw new Error("topology HTTP " + response.status);
        });
      }
      return response.text();
    }).then(function (text) {
      if (!isCifText(text)) {
        throw new Error("topology response is not mmCIF (" + text.length +
          " chars, head \"" + diagnosticHead(text) + "\")");
      }
      topologySnapshot = text;
      return text;
    });
  }

  /* 错误页 / JSON 错误不能当拓扑喂给 Mol*：必须是含 atom_site 的 mmCIF */
  function isCifText(text) {
    return !!text && text.indexOf("data_") !== -1 &&
      text.indexOf("_atom_site.") !== -1;
  }

  /* 单行化诊断文本：换行/引号转义，避免多行内容挤爆状态栏 */
  function diagnosticHead(text) {
    return String(text || "").slice(0, DIAGNOSTIC_HEAD_CHARS)
      .replace(/\\/g, "\\\\")
      .replace(/\r/g, "\\r")
      .replace(/\n/g, "\\n")
      .replace(/"/g, "\\\"");
  }

  function fetchInteractions() {
    return fetch(interactionUrl()).then(readJson).then(function (payload) {
      data = payload || null;
      renderExistence(data && data.existence, data && data.pt_types);
      updateFrameRange();
      renderLegend(data && data.residues ? data.residues : []);
      return data;
    }).catch(function (error) {
      noteFailure("interaction data failed", error);
      data = null;
      return null;
    });
  }

  function trajectoryFormat() {
    var path = String(cfg.trajectory || "").trim();
    if (!path && trajectoryInput) path = trajectoryInput.value;
    var parts = String(path).toLowerCase().split(".");
    var ext = parts.length > 1 ? parts[parts.length - 1] : "";
    return BINARY_FORMATS[ext] ? ext : DEFAULT_TRAJECTORY_FORMAT;
  }

  /* ---------------- 一次载入轨迹 ---------------- */

  /* PDB 不是拓扑格式：dataFormats 只有 psf/prmtop/top 会给出 .topology，
     所以拓扑文字必须以 model-data 走 parseTrajectory + createModel。 */
  function loadTrajectoryOnce(topologyText) {
    return viewer.loadTrajectory({
      model: {
        kind: MODEL_KIND,
        data: topologyText,
        format: TOPOLOGY_FORMAT
      },
      coordinates: {
        kind: "coordinates-url",
        url: trajectoryUrl(),
        isBinary: true,
        format: trajectoryFormat()
      }
    });
  }

  function loadScene() {
    currentStructure = structureKey(cfg);
    setStatus("loading topology and trajectory…", false);
    return Promise.all([fetchTopologyText(), fetchInteractions()])
      .then(function (results) {
        return loadTrajectoryOnce(results[0]);
      })
      .then(function (result) {
        trajectory = result;
        return afterLoad();
      })
      .catch(function (error) {
        setStatus(loadFailureText(error), true);
      });
  }

  function loadFailureText(error) {
    var message = error && error.message ? error.message : String(error);
    var text = "trajectory load failed (model=" + MODEL_KIND + ", format=" +
      TOPOLOGY_FORMAT + ", coords=" + trajectoryFormat() + "): " + message;
    /* 原子数不一致多半是旧轨迹缓存：追加可操作提示 */
    if (String(message).toLowerCase().indexOf(MISMATCH_MARKER) >= 0) {
      text += MISMATCH_HINT;
    }
    return text;
  }

  function afterLoad() {
    return waitForStructure().then(function (entry) {
      structureRef = refOf(entry.cell);
      modelRef = findModelRef();
      return clearDefaultComponents()
        .then(function () { return buildLayers(); })
        .then(function () { return focusInitially(); })
        .then(function () { return selectFrame(frame); })
        .then(function () { startPolling(); });
    });
  }

  /* 轮询等待 Mol* 结构树出现（loadTrajectory resolve 与 hierarchy 更新之间存在竞态） */
  function waitForStructure() {
    var waited = 0;
    return new Promise(function (resolve, reject) {
      (function poll() {
        var entry = structureEntry();
        if (entry) return resolve(entry);
        if (waited >= STRUCTURE_WAIT_TIMEOUT_MS) {
          return reject(new Error(structureTimeoutText()));
        }
        waited += STRUCTURE_WAIT_INTERVAL_MS;
        setTimeout(poll, STRUCTURE_WAIT_INTERVAL_MS);
      })();
    });
  }

  /* 超时诊断：拓扑规模/开头 + 后端注入的结构状态，区分"没结构"与"结构树没更新" */
  function structureTimeoutText() {
    return "no structure after " + STRUCTURE_WAIT_TIMEOUT_MS + "ms (" +
      structureDiagnostics() + ")";
  }

  function structureDiagnostics() {
    return "topology " + topologySnapshot.length + " chars, head \"" +
      diagnosticHead(topologySnapshot) + "\", has_structure=" +
      String(!!cfg.has_structure) + ", n_frames=" + String(cfg.n_frames);
  }

  /* 改帧目标是 model-from-trajectory 变换：优先 hierarchy 的 model 条目，
     其次 loadTrajectory 返回值里的 preset（model-data 分支下是 To 构建器，
     只有 .ref 没有 .cell，故统一用 refOf 取 ref），最后在 state 里按变换名兜底扫。 */
  function findModelRef() {
    var entry = structureEntry();
    var ref = entry && entry.model ? refOf(entry.model) : null;
    if (ref) return ref;
    var preset = trajectory ? trajectory.preset : null;
    ref = preset && preset.model ? refOf(preset.model) : null;
    if (ref) return ref;
    return scanModelRef();
  }

  function scanModelRef() {
    var found = null;
    viewer.plugin.state.data.cells.forEach(function (cell) {
      var transformer = cell.transform ? cell.transform.transformer : null;
      if (found || !transformer || transformer.name !== MODEL_TRANSFORM_NAME) {
        return;
      }
      found = cell.transform.ref;
    });
    return found;
  }

  /* ---------------- 结构 / 组件定位 ---------------- */

  function structures() {
    if (!viewer || !viewer.plugin) return [];
    var hierarchy = viewer.plugin.managers.structure.hierarchy;
    return hierarchy.current.structures;
  }

  function structureEntry() {
    var entries = structures();
    return entries.length ? entries[0] : null;
  }

  /* 当前主结构 cell（按 ref 实时取，帧切换后不会拿到旧 cell 对象） */
  function structureCell() {
    if (!structureRef || !viewer || !viewer.plugin) return null;
    return viewer.plugin.state.data.cells.get(structureRef) || null;
  }

  function rootStructure() {
    var cell = structureCell();
    return cell && cell.obj ? cell.obj.data : null;
  }

  function refOf(target) {
    if (!target) return null;
    if (typeof target === "string") return target;
    if (target.transform && target.transform.ref) return target.transform.ref;
    if (target.ref) return target.ref;
    if (target.cell && target.cell.transform) return target.cell.transform.ref;
    return null;
  }

  /* 整条结构的 Loci：fromExpression 只认编译好的 AST 或查询函数（字符串会被
     ml() 当常量），因此用函数形式拿 MolScript builder 构造 "all"。 */
  function wholeStructureLoci(structure) {
    if (!structure) return null;
    return molstar.lib.structure.StructureElement.Loci.fromExpression(
      structure, structureExpression);
  }

  function structureExpression(builder) {
    return builder.struct.generator.all();
  }

  /* 组件子结构上的 Loci → 主结构上的 Loci。StructureComponent 的 unit 由
     Unit.getChild 生成、保留同一个 unit.id，Loci.remap 正是按 unit.id +
     elements 行号重映射；测量/相机都要求 Loci 落在主结构上。 */
  function toRootLoci(loci) {
    var root = rootStructure();
    if (!loci || !root) return null;
    if (loci.structure === root) return loci;
    var remapped;
    try {
      remapped = molstar.lib.structure.StructureElement.Loci.remap(loci, root);
    } catch (error) {
      noteFailure("loci remap failed", error);
      return null;
    }
    return isLociEmpty(remapped) ? null : remapped;
  }

  function rootLociOf(structure) {
    return toRootLoci(wholeStructureLoci(structure));
  }

  function isLociEmpty(loci) {
    if (!loci) return true;
    return molstar.lib.structure.StructureElement.Loci.isEmpty(loci);
  }

  function unionLoci(first, second) {
    if (!first) return second;
    if (!second) return first;
    try {
      return molstar.lib.structure.StructureElement.Loci.union(first, second);
    } catch (error) {
      noteFailure("loci union failed", error);
      return first;
    }
  }

  /* ---------------- 表示重建（只做一次） ---------------- */

  function clearDefaultComponents() {
    var builder = viewer.plugin.state.data.build();
    var entries = structures();
    var removed = 0;
    for (var i = 0; i < entries.length; i++) {
      var components = entries[i].components || [];
      for (var j = 0; j < components.length; j++) {
        var ref = refOf(components[j].cell);
        if (!ref) continue;
        builder.delete(ref);
        removed += 1;
      }
    }
    if (!removed) return Promise.resolve();
    return Promise.resolve(builder.commit({ canUndo: false }))
      .catch(function (error) {
        noteFailure("default component cleanup failed", error);
      });
  }

  function structureBuilders() {
    return viewer.plugin.builders.structure;
  }

  /* script 组件：params = {language, expression}，由 Mol* 按语言解析字符串 */
  function buildScriptComponent(expression, tag, label) {
    if (!expression || !structureRef) {
      return Promise.resolve({ ok: false, component: null });
    }
    var spec = {
      type: {
        name: COMPONENT_SCRIPT,
        params: { language: EXPRESSION_LANGUAGE, expression: expression }
      },
      nullIfEmpty: true,
      label: label || ""
    };
    return Promise.resolve(
      structureBuilders().tryCreateComponent(structureRef, spec, tag)
    ).then(function (component) {
      if (!componentElementCount(component)) {
        warnOnce("empty:" + tag + "|" + expression,
          "[dii] empty component:", tag, "expr=" + expression);
        return { ok: false, component: null };
      }
      return { ok: true, component: component };
    }).catch(function (error) {
      noteFailure("component failed: " + tag, error);
      return { ok: false, component: null };
    });
  }

  /* 组件构建"成功但为空"也要当失败：否则 cartoon/配体层会静默消失 */
  function componentElementCount(component) {
    var structure = component && component.data ? component.data : null;
    return structure && typeof structure.elementCount === "number"
      ? structure.elementCount : 0;
  }

  function addLayerRepresentation(component, params, tag) {
    if (!component) return Promise.resolve(null);
    var builders = structureBuilders().representation;
    return Promise.resolve(
      builders.addRepresentation(component, params, { tag: tag })
    ).catch(function (error) {
      noteFailure("representation failed: " + tag, error);
      return null;
    });
  }

  function buildLayer(expression, tag, label, params) {
    return buildScriptComponent(expression, tag, label)
      .then(function (result) {
        trackResidueRef(tag, result);
        return addLayerRepresentation(result.component, params,
          tag + "-repr");
      });
  }

  /* 记录残基棍棒组件 ref：updatePair 时与背景层一起清掉后重建 */
  function trackResidueRef(tag, result) {
    var component = result ? result.component : null;
    if (!component || String(tag).indexOf(TAG_RESIDUE) !== 0) return;
    var ref = refOf(component.cell);
    if (ref) residueRefs.push(ref);
  }

  /* 记录卡通 polymer / 配体 ligand / 离子组件 ref（updatePair 重建背景层用） */
  function trackBackgroundRef(component) {
    if (!component) return;
    var ref = refOf(component.cell);
    if (ref) backgroundRefs.push(ref);
  }

  /* tryCreateComponentStatic(parent, keyword, {label})：keyword 是 Mol* 内置
     static 关键字（polymer/ligand/ion/water），不经过选择语言解析。 */
  function buildStaticLayer(keyword, label, params) {
    var builders = structureBuilders();
    return Promise.resolve(
      builders.tryCreateComponentStatic(structureRef, keyword,
        { label: label })
    ).then(function (component) {
      if (!component) return null;
      trackBackgroundRef(component);
      return addLayerRepresentation(component, params,
        "static-" + keyword + "-repr").then(function () {
          return component;
        });
    }).catch(function (error) {
      warnOnce("static:" + keyword,
        "[dii] static component failed:", keyword, error);
      return null;
    });
  }

  /* 不传任何颜色参数：ball-and-stick 的 defaultColorTheme 是 element-symbol，
     天然与 Mol* 原生棍棒（点击聚焦后那层）颜色一致 */
  function stickParams(sizeFactor) {
    return {
      type: STICK_TYPE,
      typeParams: {
        sizeFactor: sizeFactor,
        sizeAspectRatio: STICK_ASPECT
      }
    };
  }

  function cartoonParams() {
    return { type: CARTOON_TYPE };
  }

  function ligandParams() {
    return {
      type: STICK_TYPE,
      typeParams: { sizeFactor: LIGAND_SIZE }
    };
  }

  /* 主结构按残基分区，不再叠加额外几何：骨架与配体层都排除相互作用残基。
     整页首屏与 updatePair 共用：背景层（卡通 polymer + 配体 ligand）+ 棍棒层。 */
  function buildLayers() {
    residueRefs = [];
    backgroundRefs = [];
    return buildBackgroundLayers()
      .then(function () { return buildResidueLayers(); });
  }

  /* 背景层：卡通 polymer + 配体/离子 ligand（排除列表取当前 data.residues） */
  function buildBackgroundLayers() {
    return buildPolymerLayer()
      .then(function () { return buildLigandLayer(); });
  }

  /* 只重建 pair 相关的残基棍棒层；骨架/配体由 buildBackgroundLayers 单独重建 */
  function buildResidueLayers() {
    var residues = data && data.residues ? data.residues : [];
    var chain = Promise.resolve();
    for (var i = 0; i < residues.length; i++) {
      chain = chain.then(residueStep(residues[i], i));
    }
    return chain;
  }

  /* 卡通背景：Mol* 内置 static 组件，不拼选择语言表达式 */
  function buildPolymerLayer() {
    return buildStaticLayer("polymer", "polymer", cartoonParams());
  }

  function residueStep(residue, index) {
    return function () {
      return buildResidueLayer(residue, index);
    };
  }

  /* 相互作用残基：整残基一层棍棒，元素配色（没有第二层、没有统一色） */
  function buildResidueLayer(residue, index) {
    var expression = residueLayerExpression(residue);
    if (!expression) return Promise.resolve(null);
    var label = residue.label || "group " + (index + 1);
    return buildLayer(expression, TAG_RESIDUE + index, label,
      stickParams(RESIDUE_STICK_SIZE));
  }

  /* 水残基的 residue_selector 可能因 PDB resSeq 4 列截断而命中其它水分子；
     group_expression（后端给出，整水分子即 residue_selector）缺失时回退 b 区间。 */
  function residueLayerExpression(residue) {
    if (isWaterResidue(residue)) {
      return groupSelection(residue) ||
        (residue && residue.expression ? String(residue.expression) : null);
    }
    return residueSelector(residue);
  }

  /* 组级选择器：优先后端 group_expression（chain/resi + 原子名，与 pair 无关）；
     缺失（无结构旧 payload）时回退 b 区间表达式。 */
  function groupSelection(residue) {
    if (!residue) return null;
    if (residue.group_expression) return String(residue.group_expression);
    return residue.expression ? String(residue.expression) : null;
  }

  /* 残基级选择器（后端 residue_selector）；缺失时该残基不建表示（不猜几何） */
  function residueSelector(residue) {
    return residue && residue.residue_selector
      ? String(residue.residue_selector) : null;
  }

  /* 水的 resSeq 被 PDB 4 列截断后不唯一（15688→5689），残基级选择器会命中
     其它水分子；水分子 3 个原子本来全在组里，故跳过"其余原子"层与排除列表。 */
  function isWaterResidue(residue) {
    return !!residue && residue.is_water === true;
  }

  function residueUnionExpression(residues) {
    var parts = [];
    for (var i = 0; i < residues.length; i++) {
      var expression = groupSelection(residues[i]);
      if (expression) parts.push("(" + expression + ")");
    }
    return parts.length ? parts.join(" or ") : null;
  }

  /* 配体背景：内置 static ligand；为空再试 ion；都为空则跳过（背景可选） */
  function buildLigandLayer() {
    return buildStaticLayer("ligand", "ligand", ligandParams())
      .then(function (component) {
        if (component) return component;
        return buildStaticLayer("ion", "ion", ligandParams());
      });
  }

  /* ---------------- 图例 ---------------- */

  function renderLegend(residues) {
    if (!groupLegendEl) return;
    while (groupLegendEl.firstChild) {
      groupLegendEl.removeChild(groupLegendEl.firstChild);
    }
    for (var i = 0; i < residues.length; i++) {
      groupLegendEl.appendChild(legendItem(residues[i], i));
    }
  }

  /* π-堆积类型图例：有类型数据时显示 P/T，否则移除（updatePair 重建） */
  function renderTypeLegend(values) {
    var legend = document.querySelector("." + LEGEND_CLASS);
    var types = values && values.length ? values : [];
    if (!types.length) {
      if (legend && legend.parentNode) legend.parentNode.removeChild(legend);
      return;
    }
    if (!legend) {
      if (!groupLegendEl || !groupLegendEl.parentNode) return;
      legend = document.createElement("div");
      legend.className = LEGEND_CLASS;
      groupLegendEl.parentNode.insertBefore(legend, groupLegendEl);
    }
    while (legend.firstChild) legend.removeChild(legend.firstChild);
    for (var i = 0; i < PT_LEGEND_TYPES.length; i++) {
      legend.appendChild(typeLegendItem(PT_LEGEND_TYPES[i]));
    }
  }

  function typeLegendItem(name) {
    var item = document.createElement("span");
    item.className = LEGEND_ITEM_CLASS;
    var dot = document.createElement("span");
    dot.className = LEGEND_DOT_CLASS;
    dot.style.backgroundColor = PT_DOT_COLORS[name];
    item.appendChild(dot);
    item.appendChild(document.createTextNode(name));
    return item;
  }

  /* 图例：残基标签 + 判定为相互作用组的原子名（纯文字，3D 里不做任何叠加） */
  function legendItem(residue, index) {
    var item = document.createElement("span");
    item.className = GROUP_ITEM_CLASS;
    var names = legendAtomNames(residue);
    var label = document.createElement("span");
    label.className = GROUP_LABEL_CLASS;
    label.textContent = residue && residue.label
      ? residue.label : "group " + (index + 1);
    var atoms = document.createElement("span");
    atoms.className = "structure-group-atoms";
    atoms.textContent = legendAtomsText(names, LEGEND_ATOM_LIMIT);
    item.title = names.join(", ");
    item.appendChild(label);
    item.appendChild(atoms);
    return item;
  }

  /* 原子名顺序沿用后端（= h5 group.atoms 顺序），不排序 */
  function legendAtomNames(residue) {
    var atoms = residue && residue.atoms ? residue.atoms : [];
    var names = [];
    for (var i = 0; i < atoms.length; i++) {
      if (atoms[i] && atoms[i].name) names.push(String(atoms[i].name));
    }
    return names;
  }

  function legendAtomsText(names, limit) {
    if (!names.length) return "";
    if (!limit || names.length <= limit) return names.join(", ");
    return names.slice(0, limit).join(", ") +
      " … (+" + (names.length - limit) + ")";
  }

  /* ---------------- 相机（只在首次载入后聚焦一次） ---------------- */

  /* 三级回退：两端并集 → 两端各自合并 → 主结构整体；路径写进状态栏便于反馈 */
  function focusInitially() {
    var residues = data && data.residues ? data.residues : [];
    return focusFromResidues(residues)
      .then(function (path) { return path || focusFromEachResidue(residues); })
      .then(function (path) { return path || focusWholeStructure(); })
      .then(function (path) {
        focusPath = path || "none";
        if (!path) noteFailure("camera focus failed", "no loci");
        lockCamera();
        console.warn("[dii] focus:", focusPath);
      });
  }

  function focusFromResidues(residues) {
    return buildScriptComponent(residueUnionExpression(residues), TAG_FOCUS,
      "focus").then(function (result) {
      return applyFocusLoci(rootLociOf(componentStructure(result)),
        "residues");
    });
  }

  function focusFromEachResidue(residues) {
    var merged = null;
    var chain = Promise.resolve();
    for (var i = 0; i < residues.length; i++) {
      chain = chain.then(residueFocusStep(residues[i], i));
    }
    function residueFocusStep(residue, index) {
      return function () {
        var expression = groupSelection(residue);
        if (!expression) return null;
        return buildScriptComponent(expression,
          TAG_FOCUS_ITEM + index, "focus " + (index + 1)
        ).then(function (result) {
          merged = unionLoci(merged,
            rootLociOf(componentStructure(result)));
          return null;
        });
      };
    }
    return chain.then(function () {
      return applyFocusLoci(merged, "residues-each");
    });
  }

  function focusWholeStructure() {
    return applyFocusLoci(wholeStructureLoci(rootStructure()), "structure");
  }

  function applyFocusLoci(loci, path) {
    if (!loci || isLociEmpty(loci)) return null;
    focusLoci(loci);
    return path;
  }

  function focusLoci(loci) {
    var camera = viewer && viewer.plugin ? viewer.plugin.managers.camera : null;
    if (!camera || !loci) return;
    try {
      camera.focusLoci(loci);
    } catch (error) {
      noteFailure("camera focus failed", error);
    }
  }

  /* Mol* 在 manualReset 为 false 时每次 canvas commit 都会整场景重置相机，
     会把这次聚焦覆盖掉；聚焦后锁住相机，后续帧切换就不再动它。 */
  function lockCamera() {
    var canvas = viewer && viewer.plugin ? viewer.plugin.canvas3d : null;
    if (!canvas || typeof canvas.setProps !== "function") return;
    try {
      canvas.setProps({ camera: { manualReset: true } });
    } catch (error) {
      noteFailure("camera lock failed", error);
    }
  }

  /* ---------------- 逐帧：只换坐标帧 + 虚线 ---------------- */

  function modelCell() {
    var entry = structureEntry();
    return entry && entry.model && entry.model.cell ? entry.model.cell : null;
  }

  function currentFrameIndex() {
    var cell = modelCell();
    if (!cell || !cell.obj) return null;
    var info = molstar.lib.structure.Model.TrajectoryInfo.get(cell.obj.data);
    if (info && typeof info.index === "number") return info.index;
    var values = cell.params && cell.params.values;
    if (values && typeof values.modelIndex === "number") {
      return values.modelIndex;
    }
    return null;
  }

  /* 官方改帧方式：update model-from-trajectory 的 modelIndex（不碰结构） */
  function applyMolstarFrame(frameIndex) {
    if (!modelRef) {
      noteFailure("model cell is missing", frameIndex);
      return Promise.resolve();
    }
    var params = {};
    params[MODEL_INDEX_PARAM] = frameIndex;
    return viewer.plugin.state.updateTransform(viewer.plugin.state.data,
      modelRef, params, MODEL_INDEX_LABEL).catch(function (error) {
      noteFailure("frame switch failed", error);
    });
  }

  function selectFrame(value) {
    var target = clampFrame(value);
    var seq = ++frameSeq;
    frame = target;
    highlightFrame(target);
    updateFrameLabel(target);
    if (!ready) return Promise.resolve();
    return applyMolstarFrame(target)
      .then(function () {
        if (seq !== frameSeq) return null;
        return updateMeasurements(target);
      })
      .then(function () {
        if (seq === frameSeq) appliedSeq = seq;
        setStatus(statusText(), failures > 0);
      });
  }

  function activeAt(value) {
    var values = data && data.existence ? data.existence : null;
    if (!values) return true;
    return !!values[value];
  }

  /* 每条虚线连两个基团的几何中心：把"该基团全部原子"的 Loci 交给 addDistance，
     Mol* 的 StructureSelectionsDistance3D 自己取包围球球心并测中心距 */
  function updateMeasurements(value) {
    var seq = ++measureSeq;
    if (!data || !activeAt(value)) return replaceMeasurements([]);
    return buildMeasurements().then(function (next) {
      if (!next) return null;
      if (seq !== measureSeq) return dropMeasurements(next);
      return replaceMeasurements(next);
    }).catch(function (error) {
      noteFailure("distance overlay failed", error);
    });
  }

  /* 先按组取 root Loci，再按相邻段画线（水桥 3 组 → 2 段） */
  function buildMeasurements() {
    var groups = data && data.residues ? data.residues : [];
    var pairs = segmentPairs(groups.length);
    if (!pairs.length) return Promise.resolve([]);
    var loci = [];
    var chain = Promise.resolve();
    for (var i = 0; i < groups.length; i++) {
      chain = chain.then(groupLociStep(groups, i, loci));
    }
    return chain.then(function () { return drawSegments(loci, pairs); });
  }

  function groupLociStep(groups, index, loci) {
    return function () {
      var residue = groups[index];
      var expression = groupSelection(residue) ||
        (residue && residue.expression ? String(residue.expression) : null);
      if (!expression) return null;
      return buildScriptComponent(expression, TAG_MEASURE + index,
        "measure " + (index + 1)).then(function (result) {
        loci[index] = rootLociOf(componentStructure(result));
        return null;
      });
    };
  }

  /* 相邻基团两两成段：2 组 → 1 段；水桥 3 组 → 2 段 */
  function segmentPairs(count) {
    var pairs = [];
    for (var i = 0; i + 1 < count; i++) pairs.push([i, i + 1]);
    return pairs;
  }

  function drawSegments(loci, pairs) {
    var drawn = [];
    var chain = Promise.resolve();
    for (var i = 0; i < pairs.length; i++) {
      chain = chain.then(segmentStep(loci, pairs[i], i, drawn));
    }
    return chain.then(function () {
      return drawn.length ? drawn : null;
    });
  }

  /* 只有第一段带 customText（= pair_label）；水桥第二段留给 Mol* 显示距离 */
  function segmentStep(loci, pair, index, drawn) {
    return function () {
      var first = loci[pair[0]];
      var second = loci[pair[1]];
      if (!first || !second) return null;
      var text = index === 0 ? pairLabel() : "";
      return addDistance(first, second, text).then(function (result) {
        if (result) drawn.push(result);
        return null;
      });
    };
  }

  /* 与表格/CSV 完全一致的相互作用名，逐字使用（不做任何拼接） */
  function pairLabel() {
    return data && data.pair_label ? String(data.pair_label) : "";
  }

  function componentStructure(result) {
    var component = result ? result.component : null;
    return component && component.data ? component.data : null;
  }

  function addDistance(first, second, text) {
    var options = {
      lineParams: { dashLength: DASH_LENGTH, linesSize: LINE_SIZE }
    };
    if (text) options.customText = text;
    var measurementManager = viewer.plugin.managers.structure.measurement;
    return Promise.resolve(measurementManager.addDistance(first, second,
      options)).then(function (result) {
      if (!result) console.warn("[dii] distance not created");
      return result || null;
    }).catch(function (error) {
      noteFailure("distance failed", error);
      return null;
    });
  }

  function replaceMeasurements(next) {
    var stale = measurements;
    measurements = next || [];
    return deleteMeasurementCells(stale);
  }

  function dropMeasurements(cells) {
    return deleteMeasurementCells(cells);
  }

  function deleteMeasurementCells(cells) {
    if (!cells || !cells.length || !ready) return Promise.resolve();
    var builder = viewer.plugin.state.data.build();
    for (var i = 0; i < cells.length; i++) {
      if (cells[i].selection) builder.delete(refOf(cells[i].selection));
      if (cells[i].representation) {
        builder.delete(refOf(cells[i].representation));
      }
    }
    return Promise.resolve(builder.commit({ canUndo: false }))
      .catch(function (error) {
        noteFailure("measurement cleanup failed", error);
      });
  }

  /* ---------------- 原地切 pair（不重载结构/轨迹） ---------------- */

  /* 只清残基棍棒组件；背景层由 clearBackgroundLayers 单独管理 */
  function clearResidueLayers() {
    var stale = residueRefs;
    residueRefs = [];
    if (!stale.length || !ready) return Promise.resolve();
    var builder = viewer.plugin.state.data.build();
    for (var i = 0; i < stale.length; i++) builder.delete(stale[i]);
    return Promise.resolve(builder.commit({ canUndo: false }))
      .catch(function (error) {
        noteFailure("residue cleanup failed", error);
      });
  }

  /* 只清卡通/配体/离子背景组件；切 pair 时保留，换结构时才调用 */
  function clearBackgroundLayers() {
    var stale = backgroundRefs;
    backgroundRefs = [];
    if (!stale.length || !ready) return Promise.resolve();
    var builder = viewer.plugin.state.data.build();
    for (var i = 0; i < stale.length; i++) builder.delete(stale[i]);
    return Promise.resolve(builder.commit({ canUndo: false }))
      .catch(function (error) {
        noteFailure("background cleanup failed", error);
      });
  }

  function clearMeasurements() {
    var stale = measurements;
    measurements = [];
    return deleteMeasurementCells(stale);
  }

  /* performance 在旧环境可能缺失；仅用于诊断日志 */
  function nowMs() {
    return typeof performance !== "undefined" && performance.now
      ? performance.now() : Date.now();
  }

  /* 切 pair：只换残基棍棒/测量/帧数据；卡通+配体背景与 structure/trajectory
     均常驻不重建。payload = GET /api/pair/{type}/{idx} 的 JSON。 */
  function updatePair(payload) {
    var next = payload || {};
    if (!ready || !structureRef) {
      setStatus("3D viewer is not ready; pair not updated", true);
      return Promise.resolve();
    }
    var t0 = nowMs();
    data = {
      residues: next.groups || [],
      existence: next.existence || [],
      pt_types: next.pt_types || [],
      n_frames: typeof next.n_frames === "number"
        ? next.n_frames : cfg.n_frames,
      pair_label: next.pair_label ? String(next.pair_label) : ""
    };
    var target = typeof next.frame === "number" ? next.frame : frame;
    return clearResidueLayers()
      .then(function () { return clearMeasurements(); })
      .then(function () { return buildResidueLayers(); })
      .then(function () {
        renderExistence(data.existence, data.pt_types);
        updateFrameRange();
        renderLegend(data.residues);
        renderTypeLegend(data.pt_types);
        return null;
      })
      .then(function () { return focusInitially(); })
      .then(function () { return selectFrame(clampFrame(target)); })
      .then(function () {
        console.warn("[dii] pair updated in", Math.round(nowMs() - t0),
          "ms", {
            stickCount: data.residues.length,
            backgroundKept: true
          });
      })
      .catch(function (error) {
        console.warn("[dii] pair update failed", error);
        setStatus("pair update failed: " + (error && error.message
          ? error.message : error), true);
      });
  }

  /* ---------------- 常驻坞：attach / detach（pjax 页面切换复用 viewer） ---------------- */

  function structureKey(config) {
    return {
      topology: String(config && config.topology ? config.topology : ""),
      trajectory: String(config && config.trajectory ? config.trajectory : "")
    };
  }

  function sameStructure(config) {
    var next = structureKey(config);
    /* 页面没给结构信息：视为相同，避免无谓重载 */
    if (!next.topology && !next.trajectory) return true;
    /* 指纹缺失但 viewer 已就绪：旧轨迹仍在，同样不重载 */
    if (!currentStructure.topology && !currentStructure.trajectory) {
      return !!viewer && ready;
    }
    return next.topology === currentStructure.topology &&
      next.trajectory === currentStructure.trajectory;
  }

  function attachInitialDock() {
    var dock = document.getElementById("viewerDock");
    var slot = document.getElementById("structureDockSlot");
    if (!dock || !slot) return;
    slot.appendChild(dock);
    dock.hidden = false;
  }

  /* 把坞搬进当前页槽位；结构相同只刷新帧条/帧号，不同才重载结构 */
  function attach(config) {
    var dock = document.getElementById("viewerDock");
    var slot = document.getElementById("structureDockSlot");
    var same = !!viewer && ready && sameStructure(config);
    var hasViewer = !!viewer;
    var reloaded = false;
    var result;
    if (!dock || !slot) {
      result = Promise.resolve({ attached: false, reloaded: false });
    } else {
      if (dock.parentNode !== slot) slot.appendChild(dock);
      dock.hidden = false;
      cfg = config;
      if (same) {
        frame = clampFrame(frame);
        renderExistence(cfg.existence || [], cfg.pt_types || []);
        updateFrameRange();
        renderTypeLegend(cfg.pt_types || []);
        resizeViewer();
        result = Promise.resolve({ attached: true, reloaded: false });
      } else if (!hasViewer) {
        reloaded = true;
        mount();
        resizeViewer();
        result = Promise.resolve({ attached: true, reloaded: true });
      } else if (!ready) {
        result = Promise.resolve({ attached: true, reloaded: false });
      } else {
        reloaded = true;
        result = reloadStructure().then(function () {
          resizeViewer();
          return { attached: true, reloaded: true };
        });
      }
    }
    console.warn("[dii] dock attach", {
      same: same,
      ready: ready,
      hasViewer: hasViewer,
      reloaded: reloaded,
      from: currentStructure,
      to: structureKey(config)
    });
    return result;
  }

  /* 拓扑/轨迹变了：残基、背景、测量全部清掉后重新载入（viewer 容器不销毁） */
  function reloadStructure() {
    frame = clampFrame(cfg.frame || 0);
    return clearMeasurements()
      .then(function () { return clearResidueLayers(); })
      .then(function () { return clearBackgroundLayers(); })
      .then(function () {
        data = null;
        return loadScene();
      });
  }

  /* 移回 body 末尾并保持 hidden；viewer/trajectory 不销毁、不重置 */
  function detach() {
    var dock = document.getElementById("viewerDock");
    if (!dock) return;
    dock.hidden = true;
    if (dock.parentNode !== document.body) document.body.appendChild(dock);
  }

  /* pjax 交换前：坞若还挂在 main 里，先摘回 body，避免随 main 一起离场 */
  function detachFromMain() {
    var dock = document.getElementById("viewerDock");
    if (!dock) return;
    if (dock.parentNode && dock.parentNode !== document.body) {
      detach();
    }
  }

  /* DOM 搬移后 WebGL 画布需要重新量尺寸 */
  function resizeViewer() {
    if (!viewer || typeof viewer.handleResize !== "function") return;
    setTimeout(function () {
      try {
        viewer.handleResize();
      } catch (error) {
        console.warn("[dii] viewer resize failed", error);
      }
    }, 0);
  }

  /* pjax 交换后：有槽位+配置 → attach 并按新 pair 刷新 3D 数据；否则隐藏坞 */
  function registerSwapHook() {
    if (!window.DiiRouter || typeof window.DiiRouter.onAfterSwap !== "function") {
      return;
    }
    if (typeof window.DiiRouter.onBeforeSwap === "function") {
      window.DiiRouter.onBeforeSwap(detachFromMain);
    }
    window.DiiRouter.onAfterSwap(function (doc) {
      var slot = document.getElementById("structureDockSlot") ||
        (doc && doc.getElementById("structureDockSlot"));
      var configEl = document.getElementById("structureConfig") ||
        (doc && doc.getElementById("structureConfig"));
      if (!slot || !configEl) {
        detach();
        return null;
      }
      var next;
      try {
        next = JSON.parse(configEl.textContent);
      } catch (error) {
        noteFailure("structure config parse failed", error);
        return null;
      }
      return attach(next).then(function (result) {
        if (!result || result.reloaded) return null;
        return refreshPairData();
      }).catch(function (error) {
        console.warn("[dii] structure attach failed", error);
      });
    });
  }

  /* 结构未变时用 /api/pair 数据原地更新 pair 相关表示 */
  function refreshPairData() {
    var pairCfgEl = document.getElementById("pairConfig");
    if (!pairCfgEl || !ready) return Promise.resolve();
    var pairCfg;
    try {
      pairCfg = JSON.parse(pairCfgEl.textContent);
    } catch (error) {
      console.warn("[dii] pair config parse failed", error);
      return Promise.resolve();
    }
    return fetch("/api/pair/" + encodeURIComponent(pairCfg.type) + "/" +
      pairCfg.idx)
      .then(readJson)
      .then(function (payload) { return updatePair(payload); })
      .catch(function (error) {
        console.warn("[dii] pair data update failed", error);
      });
  }

  /* 原生轨迹控件（首帧/上一帧/下一帧）改帧后把当前帧同步回滑块与 existence 条；
     滑块拖动（pendingTimer）或换帧未落地（frameSeq != appliedSeq）时避让。 */
  function pollFrame() {
    if (pendingTimer || frameSeq !== appliedSeq) return;
    var index = currentFrameIndex();
    if (index === null || index < 0 || index === frame) return;
    frame = index;
    syncFrameUrl(index);
    highlightFrame(index);
    updateFrameLabel(index);
    updateMeasurements(index).then(function () {
      setStatus(statusText(), failures > 0);
    });
  }

  function startPolling() {
    if (pollTimer) clearInterval(pollTimer);
    pollTimer = setInterval(pollFrame, FRAME_POLL_MS);
  }

  /* ---------------- existence 条即进度条：指针拖动/点选 + 键盘 ---------------- */

  function requestFrame(value) {
    frame = clampFrame(value);
    highlightFrame(frame);
    updateFrameLabel(frame);
    if (pendingTimer) clearTimeout(pendingTimer);
    pendingTimer = setTimeout(function () {
      pendingTimer = null;
      syncFrameUrl(frame);
      selectFrame(frame);
    }, DEBOUNCE_MS);
  }

  function frameFromPointer(event) {
    var box = existenceBar.getBoundingClientRect();
    var limit = frameLimit();
    if (!box.width || !limit) return 0;
    var ratio = (event.clientX - box.left) / box.width;
    return clampFrame(Math.floor(ratio * limit));
  }

  function onPointerDown(event) {
    if (!ready || !frameLimit()) return;
    dragging = true;
    if (existenceBar.setPointerCapture) {
      existenceBar.setPointerCapture(event.pointerId);
    }
    requestFrame(frameFromPointer(event));
    event.preventDefault();
  }

  function onPointerMove(event) {
    if (!dragging) return;
    requestFrame(frameFromPointer(event));
  }

  function onPointerUp(event) {
    if (!dragging) return;
    dragging = false;
    if (existenceBar.releasePointerCapture &&
        existenceBar.hasPointerCapture &&
        existenceBar.hasPointerCapture(event.pointerId)) {
      existenceBar.releasePointerCapture(event.pointerId);
    }
  }

  function onFrameKey(event) {
    var limit = frameLimit();
    if (!limit) return;
    if (event.key === "ArrowLeft") requestFrame(frame - 1);
    else if (event.key === "ArrowRight") requestFrame(frame + 1);
    else if (event.key === "Home") requestFrame(0);
    else if (event.key === "End") requestFrame(limit - 1);
    else return;
    event.preventDefault();
  }

  /* ---------------- 序列开关 / 浏览 / 选择结构 ---------------- */

  function pickNative(input) {
    setStatus("waiting for file selection…", false);
    fetch(PICK_ENDPOINT, {
      method: "POST",
      headers: { "Content-Type": "application/x-www-form-urlencoded" },
      body: PICK_BODY
    }).then(readJson).then(function (payload) {
      var path = payload && payload.path ? String(payload.path) : "";
      if (!path) {
        setStatus("no file selected", false);
        return;
      }
      input.value = path;
      setStatus("selected " + path, false);
    }).catch(function (error) {
      setStatus(error.message, true);
    });
  }

  function openServerBrowser(input) {
    if (!window.DiiFilePick || typeof window.DiiFilePick.open !== "function") {
      setStatus("server file browser is not available", true);
      return;
    }
    window.DiiFilePick.open({
      onPick: function (path) {
        input.value = path;
        setStatus("selected " + path, false);
      }
    });
  }

  function onBrowseClick(event) {
    var button = event.currentTarget;
    var input = document.getElementById(
      button.getAttribute(BROWSE_TARGET_ATTR));
    if (!input) return;
    var mode = button.getAttribute(BROWSE_MODE_ATTR);
    if (mode === BROWSE_PICKER) {
      pickNative(input);
    } else if (mode === BROWSE_SERVER) {
      openServerBrowser(input);
    } else {
      setStatus("structure browsing is not available", true);
    }
  }

  function selectStructure(event) {
    event.preventDefault();
    var topology = topologyInput.value.trim();
    var trajectoryPath = trajectoryInput.value.trim();
    if (!topology || !trajectoryPath) {
      setStatus("topology and trajectory are required", true);
      return;
    }
    openButton.disabled = true;
    setStatus("opening structure…", false);
    var body = "topology=" + encodeURIComponent(topology) +
      "&trajectory=" + encodeURIComponent(trajectoryPath);
    fetch(SELECT_URL, {
      method: "POST",
      headers: { "Content-Type": "application/x-www-form-urlencoded" },
      body: body
    }).then(readJson).then(function () {
      /* 手动选择优先：本次会话不再让自动尝试覆盖用户显式指定的结构 */
      markAutoSkipped();
      window.location.reload();
    }).catch(function (error) {
      openButton.disabled = false;
      setStatus(error.message, true);
    });
  }

  /* ---------------- 自动载入 h5 记录的拓扑/轨迹 ---------------- */

  /* auto=0 表示本次会话已尝试过（失败后写入地址栏，刷新不再反复请求） */
  function shouldAutoLoad() {
    return !!cfg.auto_attempt && !!cfg[AUTO_URL_KEY] && !autoSkipRequested();
  }

  function autoSkipRequested() {
    var url = new URL(window.location.href);
    return url.searchParams.get(AUTO_TRIED_FLAG) === AUTO_SKIP_VALUE;
  }

  function markAutoSkipped() {
    if (!window.history || typeof window.history.replaceState !== "function") {
      return;
    }
    var url = new URL(window.location.href);
    if (url.searchParams.get(AUTO_TRIED_FLAG) === AUTO_SKIP_VALUE) return;
    url.searchParams.set(AUTO_TRIED_FLAG, AUTO_SKIP_VALUE);
    try {
      window.history.replaceState(null, "", url.pathname + url.search +
        url.hash);
    } catch (error) {
      /* 浏览器限制 replaceState 频率时忽略：本次已回退手选 */
    }
  }

  function autoRequestBody() {
    return "type=" + encodeURIComponent(cfg.type) + "&idx=" + cfg.idx;
  }

  /* 409 时 h5 记录的路径仍值得展示：预填输入框让用户只需确认 */
  function prefillAutoPaths(payload) {
    if (payload && payload.topology && topologyInput) {
      topologyInput.value = String(payload.topology);
    }
    if (payload && payload.trajectory && trajectoryInput) {
      trajectoryInput.value = String(payload.trajectory);
    }
  }

  function autoFallback(status, payload) {
    var body = payload || {};
    if (status === AUTO_CONFLICT_STATUS) {
      prefillAutoPaths(body);
      setStatus(AUTO_UNAVAILABLE_TEXT + " (" + (body.reason || body.detail ||
        status) + ")" + AUTO_CONFIRM_HINT, true);
    } else {
      setStatus(body.detail || body.error || (AUTO_FAILED_PREFIX + status),
        true);
    }
    markAutoSkipped();
  }

  /* 200 后整页刷新：后端已写入结构状态，重载的页面直接进正常载入流程 */
  function autoLoadStructure() {
    setStatus(AUTO_LOADING_TEXT, false);
    fetch(cfg[AUTO_URL_KEY], {
      method: "POST",
      headers: { "Content-Type": "application/x-www-form-urlencoded" },
      body: autoRequestBody()
    }).then(function (response) {
      if (response.ok) {
        window.location.reload();
        return null;
      }
      return response.json().then(function (payload) {
        autoFallback(response.status, payload);
      }).catch(function () {
        autoFallback(response.status, null);
      });
    }).catch(function (error) {
      autoFallback(0, { detail: error.message });
    });
  }

  function initViewer() {
    if (!window.molstar || !molstar.Viewer) {
      setStatus("Mol* viewer is not available", true);
      return;
    }
    molstar.Viewer.create(host, VIEWER_OPTIONS).then(function (instance) {
      viewer = instance;
      ready = true;
      return loadScene();
    }).catch(function (error) {
      setStatus("viewer init failed: " + error.message, true);
    });
  }

  if (existenceBar) {
    existenceBar.addEventListener("pointerdown", onPointerDown);
    existenceBar.addEventListener("pointermove", onPointerMove);
    existenceBar.addEventListener("pointerup", onPointerUp);
    existenceBar.addEventListener("pointercancel", onPointerUp);
    existenceBar.addEventListener("keydown", onFrameKey);
  }

  for (var i = 0; i < browseButtons.length; i++) {
    browseButtons[i].addEventListener("click", onBrowseClick);
  }

  if (form) form.addEventListener("submit", selectStructure);

  /* ---------------- 挂载入口：幂等，外层可在切 pair 时复用同一 viewer ---------------- */

  function mount() {
    if (mounted) return;
    mounted = true;
    attachInitialDock();
    highlightFrame(frame);
    /* 自动尝试优先于 has_structure：切到另一数据集时旧结构与该 h5 记录不一致，
       必须先按记录重载，否则会用旧宇宙算新拓扑（400）。 */
    if (shouldAutoLoad()) {
      autoLoadStructure();
    } else if (cfg.has_structure) {
      updateFrameLabel(frame);
      initViewer();
    } else {
      setStatus("select topology and trajectory to load the structure", false);
    }
  }

  window.DiiStructure = {
    mount: mount,
    updatePair: updatePair,
    attach: attach,
    detach: detach,
    isMounted: function () { return !!viewer; }
  };

  registerSwapHook();

  /* 首屏自挂载（幂等）；pjax 后续走 DiiRouter.onAfterSwap 回调 */
  if (window.DiiStructure) window.DiiStructure.mount();
})();
