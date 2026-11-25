/* 3D 结构页：Mol* 原生轨迹（拓扑 PDB + 原始轨迹文件）只载入一次；
   帧切换只更新 modelIndex（不重建结构、不动相机），相互作用残基以棍棒呈现。 */
(function () {
  "use strict";

  var cfgEl = document.getElementById("structureConfig");
  if (!cfgEl) return;
  var cfg = JSON.parse(cfgEl.textContent);

  var host = document.getElementById("molstarViewer");
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
  var CONTACT_URL = "/api/structure/contact";
  var SELECT_URL = "/api/structure/select";
  var PICK_ENDPOINT = "/api/pick";
  var PICK_BODY = "mode=file&kind=structure";

  var TOPOLOGY_FORMAT = "pdb";
  var MODEL_KIND = "model-data";
  var DEFAULT_TRAJECTORY_FORMAT = "xtc";
  var BINARY_FORMATS = { xtc: true, trr: true, dcd: true };

  /* 后端给的表达式是 PyMOL 风格选择语言；Mol* 只有 script 类型的组件会按
     language 解析字符串（表达式类型需要编译好的 AST，字符串会被当常量）。 */
  var EXPRESSION_LANGUAGE = "pymol";
  var COMPONENT_SCRIPT = "script";
  var POLYMER_KEYWORD = "polymer";        /* static 组件关键字（兜底用） */
  /* 宏：protein / water / nucleic 是选择语言里真实存在的宏；
     polymer 只是数值属性、organic/inorganic 不是宏，裸用会解析成空选择 */
  var POLYMER_EXPRESSION = "protein";
  var NON_POLYMER_EXPRESSION = "not (protein or water or nucleic)";
  var STATIC_LIGAND = "ligand";
  var STATIC_ION = "ion";
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

  var TAG_POLYMER = "dii-polymer";
  var TAG_LIGAND = "dii-ligand";
  var TAG_RESIDUE = "dii-residue-";
  var TAG_FOCUS = "dii-focus";
  var TAG_FOCUS_ITEM = "dii-focus-";
  var TAG_CONTACT_A = "dii-contact-a";
  var TAG_CONTACT_B = "dii-contact-b";

  var LEGEND_ATOM_LIMIT = 12;      /* 图例里最多列出的原子名个数，其余用 (+N) */
  var GROUP_ITEM_CLASS = "structure-group-item";
  var GROUP_LABEL_CLASS = "structure-group-label";
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
    viewportShowAnimation: false,
    viewportShowTrajectoryControls: true
  };

  var viewer = null;
  var ready = false;
  var trajectory = null;        /* loadTrajectory 的返回值 {model, coords, preset} */
  var structureRef = null;      /* 主结构 cell 的 ref：帧更新后 cell 对象会换，始终按 ref 取 */
  var modelRef = null;          /* model-from-trajectory 的 ref（改帧用） */
  var data = null;              /* /api/structure/interactions 的 JSON */
  var frame = 0;
  var measurement = null;       /* 当前唯一一条虚线测量 {selection, representation} */
  var contactState = "";        /* 状态栏：contact: ok|failed */
  var focusPath = "";           /* 状态栏：focus: residues|residues-each|structure */
  var contactSeq = 0;
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
    var bars = existenceBar.children;
    for (var i = 0; i < bars.length; i++) {
      bars[i].className = (values[i] ? "on" : "") +
        (i === value ? " active" : "");
    }
    existenceBar.setAttribute("aria-valuenow", String(value));
  }

  function renderExistence(values) {
    if (!existenceBar || !values || !values.length) return;
    while (existenceBar.firstChild) {
      existenceBar.removeChild(existenceBar.firstChild);
    }
    for (var i = 0; i < values.length; i++) {
      var bar = document.createElement("span");
      bar.className = values[i] ? "on" : "";
      bar.title = "frame " + i + ": " +
        (values[i] ? "active" : "inactive");
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

  function statusText() {
    var parts = ["frame " + frame + " / " + frameLimit(),
                 "trajectory loaded"];
    var residues = data && data.residues ? data.residues : [];
    if (residues.length) parts.push(residues.length + " groups");
    if (contactState) parts.push("contact: " + contactState);
    if (focusPath) parts.push("focus: " + focusPath);
    if (failures) parts.push(failures + " warnings");
    return parts.join(" · ");
  }

  /* ---------------- 数据请求 ---------------- */

  function interactionUrl() {
    return INTERACTIONS_URL + "?type=" + encodeURIComponent(cfg.type) +
      "&idx=" + cfg.idx;
  }

  /* 拓扑 PDB 由后端按该对原子打 B-factor 标记（表达式据此圈选），必须带 type/idx */
  function topologyUrl() {
    return TOPOLOGY_URL + "?type=" + encodeURIComponent(cfg.type) +
      "&idx=" + cfg.idx;
  }

  function contactUrl(value) {
    return CONTACT_URL + "?type=" + encodeURIComponent(cfg.type) +
      "&idx=" + cfg.idx + "&frame=" + value;
  }

  function fetchTopologyText() {
    return fetch(topologyUrl()).then(function (response) {
      if (!response.ok) {
        throw new Error("topology HTTP " + response.status);
      }
      return response.text();
    });
  }

  function fetchInteractions() {
    return fetch(interactionUrl()).then(readJson).then(function (payload) {
      data = payload || null;
      renderExistence(data && data.existence);
      updateFrameRange();
      renderLegend(data && data.residues ? data.residues : []);
      return data;
    }).catch(function (error) {
      noteFailure("interaction data failed", error);
      data = null;
      return null;
    });
  }

  function fetchContact(value) {
    return fetch(contactUrl(value)).then(readJson).then(function (payload) {
      return payload && payload.contact ? payload.contact : null;
    }).catch(function (error) {
      noteFailure("contact request failed", error);
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
        url: TRAJECTORY_URL,
        isBinary: true,
        format: trajectoryFormat()
      }
    });
  }

  function loadScene() {
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
    return "trajectory load failed (model=" + MODEL_KIND + ", format=" +
      TOPOLOGY_FORMAT + ", coords=" + trajectoryFormat() + "): " + message;
  }

  function afterLoad() {
    var entry = structureEntry();
    if (!entry) throw new Error("trajectory produced no structure");
    structureRef = refOf(entry.cell);
    modelRef = findModelRef();
    return clearDefaultComponents()
      .then(function () { return buildLayers(); })
      .then(function () { return focusInitially(); })
      .then(function () { return selectFrame(0); })
      .then(function () { startPolling(); });
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
        console.warn("[dii] empty component:", tag, "expr=" + expression);
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
        return addLayerRepresentation(result.component, params,
          tag + "-repr");
      });
  }

  /* tryCreateComponentStatic(parent, keyword, {label}, tags)：第三个参数是
     label/tags，不是标签字符串（第四个参数必须是数组，否则展开会报错）。 */
  function buildStaticLayer(keyword, label, params) {
    var builders = structureBuilders();
    return Promise.resolve(
      builders.tryCreateComponentStatic(structureRef, keyword,
        { label: label })
    ).then(function (component) {
      return addLayerRepresentation(component, params,
        "static-" + keyword + "-repr");
    }).catch(function (error) {
      noteFailure("static component failed: " + keyword, error);
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

  /* 主结构按残基分区，不再叠加额外几何：骨架与配体层都排除相互作用残基，
     相互作用残基本身拆成"其余原子（灰）+ 定义基团原子（组色）"两块互补棍棒。 */
  function buildLayers() {
    var residues = data && data.residues ? data.residues : [];
    var chain = buildPolymerLayer(residues);
    chain = chain.then(function () { return buildLigandLayer(residues); });
    for (var i = 0; i < residues.length; i++) {
      chain = chain.then(residueStep(residues[i], i));
    }
    return chain;
  }

  function buildPolymerLayer(residues) {
    return buildScriptComponent(polymerExpression(residues), TAG_POLYMER,
      "polymer").then(function (result) {
      if (result.component) {
        return addLayerRepresentation(result.component, cartoonParams(),
          TAG_POLYMER + "-repr");
      }
      if (result.ok) return null;
      return buildStaticLayer(POLYMER_KEYWORD, "polymer", cartoonParams());
    });
  }

  function polymerExpression(residues) {
    var excluded = residueSelectorUnion(residues);
    if (!excluded) return POLYMER_EXPRESSION;
    return POLYMER_EXPRESSION + " and not (" + excluded + ")";
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

  /* 水残基的 residue_selector 不唯一（resSeq 4 列截断），改用精确的原子并集 */
  function residueLayerExpression(residue) {
    if (isWaterResidue(residue)) {
      return residue && residue.expression
        ? String(residue.expression) : null;
    }
    return residueSelector(residue);
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

  function residueSelectorUnion(residues) {
    var parts = [];
    for (var i = 0; i < residues.length; i++) {
      if (isWaterResidue(residues[i])) continue;
      var selector = residueSelector(residues[i]);
      if (selector) parts.push("(" + selector + ")");
    }
    return parts.length ? parts.join(" or ") : null;
  }

  /* 配体/离子层：非蛋白、非水、非核酸，并排除相互作用残基（它们已由棍棒显示） */
  function ligandExpression(residues) {
    var excluded = residueSelectorUnion(residues);
    if (!excluded) return NON_POLYMER_EXPRESSION;
    return "(" + NON_POLYMER_EXPRESSION + ") and not (" + excluded + ")";
  }

  function residueUnionExpression(residues) {
    var parts = [];
    for (var i = 0; i < residues.length; i++) {
      if (residues[i] && residues[i].expression) {
        parts.push("(" + residues[i].expression + ")");
      }
    }
    return parts.length ? parts.join(" or ") : null;
  }

  function buildLigandLayer(residues) {
    return buildScriptComponent(ligandExpression(residues), TAG_LIGAND,
      "ligand").then(function (result) {
      if (result.component) {
        return addLayerRepresentation(result.component, ligandParams(),
          TAG_LIGAND + "-repr");
      }
      if (result.ok) return null;
      return buildPlainLigandLayer();
    });
  }

  function buildPlainLigandLayer() {
    return buildScriptComponent(NON_POLYMER_EXPRESSION, TAG_LIGAND, "ligand")
      .then(function (result) {
        if (result.component) {
          return addLayerRepresentation(result.component, ligandParams(),
            TAG_LIGAND + "-repr");
        }
        if (result.ok) return null;
        return buildStaticLayer(STATIC_LIGAND, "ligand", ligandParams())
          .then(function () {
            return buildStaticLayer(STATIC_ION, "ion", ligandParams());
          });
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
        if (!residue || !residue.expression) return null;
        return buildScriptComponent(residue.expression,
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
        return updateContact(target);
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

  function contactAt(value) {
    if (data && data.contacts) {
      return Promise.resolve(data.contacts[value] || null);
    }
    return fetchContact(value);
  }

  /* 先建新测量、成功后才删旧的：建失败时旧虚线保留，状态栏 contact: failed */
  function updateContact(value) {
    var seq = ++contactSeq;
    if (!data || !activeAt(value)) return removeMeasurement();
    return contactAt(value).then(function (contact) {
      if (seq !== contactSeq || !contact) return null;
      return buildContactMeasurement(contact);
    }).then(function (next) {
      if (!next) return null;
      if (seq !== contactSeq) return dropMeasurement(next);
      return installMeasurement(next);
    }).catch(function (error) {
      contactState = "failed";
      noteFailure("contact overlay failed", error);
    });
  }

  function buildContactMeasurement(contact) {
    var residues = data && data.residues ? data.residues : [];
    var first = groupAtomExpression(residues[contact[0]], contact[1]);
    var second = groupAtomExpression(residues[contact[2]], contact[3]);
    if (!first || !second) return null;
    return Promise.all([
      buildScriptComponent(first, TAG_CONTACT_A, "contact a"),
      buildScriptComponent(second, TAG_CONTACT_B, "contact b")
    ]).then(function (results) {
      var left = rootLociOf(componentStructure(results[0]));
      var right = rootLociOf(componentStructure(results[1]));
      if (!left || !right) {
        contactState = "failed";
        console.warn("[dii] contact overlay failed: no root loci", contact);
        return null;
      }
      return addDistance(left, right, contactLabel(contact));
    });
  }

  function groupAtomExpression(residue, position) {
    var atoms = residue && residue.atoms ? residue.atoms : [];
    if (!(position >= 0 && position < atoms.length)) return null;
    return atoms[position].expression || null;
  }

  function contactLabel(contact) {
    var residues = data && data.residues ? data.residues : [];
    var first = residues[contact[0]];
    var second = residues[contact[2]];
    var left = first && first.label ? first.label : "group " + contact[0];
    var right = second && second.label ? second.label : "group " + contact[2];
    return left + " ┄ " + right;
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
      if (!result) {
        contactState = "failed";
        console.warn("[dii] contact overlay failed: distance not created");
        return null;
      }
      return result;
    }).catch(function (error) {
      contactState = "failed";
      noteFailure("distance failed", error);
      return null;
    });
  }

  function installMeasurement(next) {
    return removeMeasurement().then(function () {
      measurement = next;
      contactState = "ok";
      return next;
    });
  }

  function dropMeasurement(cells) {
    return deleteMeasurementCells(cells);
  }

  function removeMeasurement() {
    var current = measurement;
    measurement = null;
    return deleteMeasurementCells(current);
  }

  function deleteMeasurementCells(cells) {
    if (!cells || !ready) return Promise.resolve();
    var builder = viewer.plugin.state.data.build();
    if (cells.selection) builder.delete(refOf(cells.selection));
    if (cells.representation) builder.delete(refOf(cells.representation));
    return Promise.resolve(builder.commit({ canUndo: false }))
      .catch(function (error) {
        noteFailure("measurement cleanup failed", error);
      });
  }

  /* 原生轨迹控件（首帧/上一帧/下一帧）改帧后把当前帧同步回滑块与 existence 条；
     滑块拖动（pendingTimer）或换帧未落地（frameSeq != appliedSeq）时避让。 */
  function pollFrame() {
    if (pendingTimer || frameSeq !== appliedSeq) return;
    var index = currentFrameIndex();
    if (index === null || index < 0 || index === frame) return;
    frame = index;
    highlightFrame(index);
    updateFrameLabel(index);
    updateContact(index).then(function () {
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
      window.location.reload();
    }).catch(function (error) {
      openButton.disabled = false;
      setStatus(error.message, true);
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

  highlightFrame(0);
  if (cfg.has_structure) {
    updateFrameLabel(0);
    initViewer();
  } else {
    setStatus("select topology and trajectory to load the structure", false);
  }
})();
