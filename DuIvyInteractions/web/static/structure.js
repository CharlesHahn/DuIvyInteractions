/* 3D 结构页：Mol* 渲染单帧 PDB，默认视角聚焦相互作用，并提供水/序列开关与路径浏览 */
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
  var slider = document.getElementById("frameSlider");
  var frameLabel = document.getElementById("frameLabel");
  var existenceBar = document.getElementById("existenceBar");
  var statusEl = document.getElementById("structureStatus");
  var showWater = document.getElementById("showWater");
  var showSequence = document.getElementById("showSequence");
  var browseButtons = document.querySelectorAll("button[data-browse]");

  var DEBOUNCE_MS = 80;        // 滑块输入到请求的防抖间隔
  var CACHE_LIMIT = 5;         // 前端 PDB 文本 LRU 帧数（键含 solvent）
  var PREFETCH_AHEAD = 1;      // 预取当前帧之后的一帧
  var PDB_FORMAT = "pdb";
  var DASH_LENGTH = 0.12;      // Mol* 距离测量虚线间距
  var LINE_SIZE = 0.06;        // Mol* 距离测量线宽
  var SOLVENT_HIDE = "hide";
  var SOLVENT_FULL = "full";
  var REGION_FULL = "full";
  var REGION_HIDDEN = "hidden";
  var NO_SEQUENCE_CLASS = "structure-viewer--no-sequence";
  var BROWSE_TARGET_ATTR = "data-target";
  var BROWSE_MODE_ATTR = "data-browse";
  var BROWSE_PICKER = "picker";
  var BROWSE_SERVER = "server";
  var PICK_ENDPOINT = "/api/pick";
  var PICK_BODY = "mode=file&kind=structure";
  var VIEWER_OPTIONS = {
    layoutIsExpanded: false,
    layoutShowControls: true,
    layoutShowLeftPanel: true,
    layoutShowSequence: true,
    layoutShowLog: false,
    layoutShowRemoteState: false,
    viewportShowControls: true,
    viewportShowSelectionMode: false,
    viewportShowAnimation: false
  };
  var STATUS_CLASS = "structure-status";
  var STATUS_ERROR_CLASS = "structure-status error";
  var GROUP_COLORS = ["#e64b35", "#3c5488", "#00a087"];   // NPG 配色，按 group 顺序循环取色
  var GROUP_TYPE = "ball-and-stick";    // 两端基团用球棍，蛋白质骨架保持 cartoon
  var GROUP_SIZE_FACTOR = 0.3;          // ball-and-stick 球半径系数
  var GROUP_SIZE_ASPECT = 0.6;          // ball-and-stick 键粗细比
  var GROUP_COLOR_THEME = "uniform";    // Mol* 颜色主题：纯色
  var GROUP_BUNDLE_TYPE = "bundle";     // StructureComponent 只接受 static/script/expression/bundle
  var GROUP_TAG_PREFIX = "dii-group-";
  var GROUP_ITEM_CLASS = "structure-group-item";
  var GROUP_LABEL_CLASS = "structure-group-label";
  var GROUP_DOT_CLASS = "structure-group-dot";
  var GROUP_LEGEND_ID = "groupLegend";
  var GROUP_NOTE = " · group overlay unavailable";
  var CONTACT_NOTE = " · contact overlay unavailable";
  var HIGHLIGHT_FALLBACK_NOTE = " · highlight fallback";
  var HIGHLIGHT_FAILED_NOTE = " · highlight unavailable";
  var groupLegendEl = document.getElementById(GROUP_LEGEND_ID);

  var viewer = null;
  var ready = false;
  var latestFrame = 0;
  var requestSeq = 0;
  var pendingTimer = null;
  var renderChain = Promise.resolve();
  var cache = {};
  var cacheOrder = [];
  var contacts = [];
  var groupReprs = [];
  var mainStructure = null;    // 主结构 cell（cartoon），回退路径与聚焦用
  var highlightGroups = 0;     // 本帧成功加球棍的高亮组数（状态栏展示）
  var waterAutoChecked = false;

  function setStatus(text, isError) {
    if (!statusEl) return;
    statusEl.textContent = text || "";
    statusEl.className = isError ? STATUS_ERROR_CLASS : STATUS_CLASS;
  }

  function solventMode() {
    return showWater && showWater.checked ? SOLVENT_FULL : SOLVENT_HIDE;
  }

  function cacheKey(frame, solvent) {
    return frame + ":" + solvent;
  }

  function touchCache(key) {
    var index = cacheOrder.indexOf(key);
    if (index >= 0) cacheOrder.splice(index, 1);
    cacheOrder.push(key);
  }

  function cacheGet(key) {
    if (!Object.prototype.hasOwnProperty.call(cache, key)) return null;
    touchCache(key);
    return cache[key];
  }

  function cachePut(key, payload) {
    cache[key] = payload;
    touchCache(key);
    while (cacheOrder.length > CACHE_LIMIT) {
      delete cache[cacheOrder.shift()];
    }
  }

  function clearCache() {
    cache = {};
    cacheOrder = [];
  }

  function highlightFrame(frame) {
    if (!existenceBar) return;
    var bars = existenceBar.children;
    for (var i = 0; i < bars.length; i++) {
      var active = cfg.existence && cfg.existence[i];
      bars[i].className = (active ? "on" : "") + (i === frame ? " current" : "");
    }
  }

  function updateFrameLabel(frame) {
    if (!frameLabel) return;
    frameLabel.textContent = cfg.n_frames
      ? "frame " + frame + " / " + cfg.n_frames
      : "no structure";
  }

  function requestFrame(frame, solvent) {
    var url = "/api/structure/view?type=" + encodeURIComponent(cfg.type) +
      "&idx=" + cfg.idx + "&frame=" + frame + "&solvent=" + solvent;
    return fetch(url).then(function (response) {
      return response.json().then(function (data) {
        if (!response.ok) {
          var detail = data && data.detail ? data.detail : "HTTP " + response.status;
          throw new Error(String(detail));
        }
        return data;
      });
    });
  }

  function queueRender(payload) {
    renderChain = renderChain.then(function () {
      if (payload.frame !== latestFrame) return null;
      return renderPayload(payload);
    }).catch(function (error) {
      setStatus("render failed: " + error.message, true);
    });
    return renderChain;
  }

  function autoEnableWater(payload) {
    if (waterAutoChecked) return false;
    waterAutoChecked = true;
    if (payload.involves_water !== true || !showWater || showWater.checked) {
      return false;
    }
    showWater.checked = true;
    return true;
  }

  function loadFrame(frame) {
    if (!ready || frame < 0 || frame >= cfg.n_frames) return;
    latestFrame = frame;
    updateFrameLabel(frame);
    var solvent = solventMode();
    var key = cacheKey(frame, solvent);
    var cached = cacheGet(key);
    if (cached) {
      queueRender(cached);
      prefetch(frame + PREFETCH_AHEAD);
      return;
    }
    var sequence = ++requestSeq;
    setStatus("loading frame " + frame + "…", false);
    requestFrame(frame, solvent).then(function (payload) {
      cachePut(key, payload);
      if (sequence !== requestSeq) return;
      if (autoEnableWater(payload)) {
        loadFrame(frame);
        return;
      }
      queueRender(payload);
      prefetch(frame + PREFETCH_AHEAD);
    }).catch(function (error) {
      if (sequence === requestSeq) setStatus(error.message, true);
    });
  }

  function prefetch(frame) {
    var solvent = solventMode();
    var key = cacheKey(frame, solvent);
    if (frame < 0 || frame >= cfg.n_frames || cacheGet(key)) return;
    requestFrame(frame, solvent).then(function (payload) {
      cachePut(key, payload);
    }).catch(function () { /* 预取失败保持静默，滚动时再请求 */ });
  }

  function renderPayload(payload) {
    var overlay = { note: "", failed: false };
    setStatus("rendering frame " + payload.frame + "…", false);
    highlightFrame(payload.frame);
    return clearContacts().then(function () {
      return clearGroups();
    }).then(function () {
      return removeStructures();
    }).then(function () {
      return loadStructureCell(payload.pdb);
    }).then(function (cell) {
      mainStructure = cell;
      return drawOverlay(payload);
    }).then(function (result) {
      overlay = result || overlay;
    }).then(function () {
      setStatus(statusText(payload) + overlay.note, !!overlay.failed);
    });
  }

  /* 默认走高亮小结构路径；highlights 缺失时退回主结构（旧 groups/contact 路径） */
  function drawOverlay(payload) {
    renderGroupLegend(payload.highlights || payload.groups || []);
    if (!payload.highlights || !payload.highlights.length) {
      return drawMainOverlay(payload).then(function (result) {
        result.note = HIGHLIGHT_FALLBACK_NOTE + result.note;
        return result;
      });
    }
    return drawHighlightOverlay(payload).catch(function (error) {
      console.warn("[dii] highlight overlay failed", error);
      return { note: HIGHLIGHT_FAILED_NOTE, failed: true };
    });
  }

  /* 每条 highlight PDB 单独载入、整条结构加球棍；contact_pdb 两个原子画虚线 */
  function drawHighlightOverlay(payload) {
    var state = { note: "", failed: false };
    var focus = null;
    return addHighlightReprs(payload).then(function (groupFailed) {
      if (groupFailed) {
        state.note += GROUP_NOTE;
        state.failed = true;
      }
      return addContactPdb(payload);
    }).then(function (contact) {
      focus = contact.loci;
      if (contact.failed) {
        state.note += CONTACT_NOTE;
        state.failed = true;
      }
      if (!contact.loci && payload.contact) {
        return drawContact(structureDataOf(mainStructure), payload.contact,
          contactText(payload)).catch(function () {
            state.note += CONTACT_NOTE;
            state.failed = true;
          });
      }
      return null;
    }).then(function () {
      if (focus) {
        focusLoci(focus);
      } else {
        focusLoci(wholeStructureLoci(structureDataOf(mainStructure)));
      }
      return state;
    });
  }

  /* 回退：主结构 + groups.atoms + payload.contact（旧行为，保留） */
  function drawMainOverlay(payload) {
    var state = { note: "", failed: false };
    var structure = structureData();
    return addGroupReprs(payload).then(function (groupFailed) {
      if (groupFailed) {
        state.note += GROUP_NOTE;
        state.failed = true;
      }
      if (!(payload.active && payload.contact)) return null;
      return drawContact(structure, payload.contact, contactText(payload))
        .catch(function () {
          state.note += CONTACT_NOTE;
          state.failed = true;
        });
    }).then(function () {
      focusMainOverlay(payload, structure);
      return state;
    });
  }

  function contactText(payload) {
    return payload.contact_label ? String(payload.contact_label) : "";
  }

  function statusText(payload) {
    var parts = ["frame " + payload.frame, payload.active ? "active" : "inactive"];
    if (typeof payload.n_atoms_shown === "number" &&
        typeof payload.n_atoms_total === "number") {
      parts.push("atoms " + payload.n_atoms_shown + "/" + payload.n_atoms_total);
    }
    if (payload.involves_water) parts.push("water involved");
    if (payload.solvent === SOLVENT_FULL) parts.push("water shown");
    var groups = payload.groups || payload.highlights || [];
    if (groups.length) parts.push(groups.length + " groups");
    if (highlightGroups) parts.push("highlight: " + highlightGroups + " groups");
    return parts.join(" · ");
  }

  function clearContacts() {
    if (!contacts.length || !viewer) return Promise.resolve();
    var builder = viewer.plugin.state.data.build();
    for (var i = 0; i < contacts.length; i++) {
      if (contacts[i].selection) builder.delete(contacts[i].selection);
      if (contacts[i].representation) builder.delete(contacts[i].representation);
    }
    contacts = [];
    return Promise.resolve(builder.commit({ canUndo: false })).catch(function () {
      return null;
    });
  }

  function clearGroups() {
    if (!groupReprs.length || !viewer) return Promise.resolve();
    var builder = viewer.plugin.state.data.build();
    for (var i = 0; i < groupReprs.length; i++) {
      if (groupReprs[i].representation) {
        builder.delete(groupReprs[i].representation);
      }
      if (groupReprs[i].component) builder.delete(groupReprs[i].component);
    }
    groupReprs = [];
    return Promise.resolve(builder.commit({ canUndo: false })).catch(function () {
      return null;
    });
  }

  function removeStructures() {
    mainStructure = null;
    highlightGroups = 0;
    var hierarchy = viewer.plugin.managers.structure.hierarchy;
    var roots = hierarchy.current.trajectories;
    if (!roots || !roots.length) return Promise.resolve();
    return Promise.resolve(hierarchy.remove(roots, false));
  }

  function cellOfStructure(structure) {
    return structure.transform && structure.transform.cell
      ? structure.transform.cell : structure.cell;
  }

  function structureCell() {
    var structures = viewer.plugin.managers.structure.hierarchy.current.structures;
    if (!structures.length) return null;
    return cellOfStructure(structures[structures.length - 1]);
  }

  function structureData() {
    return structureDataOf(structureCell());
  }

  function structureDataOf(cell) {
    return cell && cell.obj ? cell.obj.data : null;
  }

  function structureCellRefs() {
    var structures = viewer.plugin.managers.structure.hierarchy.current.structures;
    var refs = [];
    for (var i = 0; i < structures.length; i++) {
      refs.push(structures[i].cell.transform.ref);
    }
    return refs;
  }

  function newestStructureCell(before) {
    var structures = viewer.plugin.managers.structure.hierarchy.current.structures;
    for (var i = structures.length - 1; i >= 0; i--) {
      var cell = cellOfStructure(structures[i]);
      if (cell && before.indexOf(cell.transform.ref) < 0) return cell;
    }
    return null;
  }

  /* 载入一条结构（主结构或高亮小结构）并返回它的 cell：
     实际新增的那条（等同 hierarchy 最后一项），不依赖条目排序 */
  function loadStructureCell(data) {
    var before = structureCellRefs();
    return viewer.loadStructureFromData(data, PDB_FORMAT).then(function () {
      var cell = newestStructureCell(before);
      if (!cell) throw new Error("structure was not added");
      return cell;
    });
  }

  /* StateBuilder.to() 只认 string / cell（.transform+.parent+.status）/ 带 .ref 的
     selector；hierarchy 条目是 {kind:"structure", cell, ...}，没有 .ref，
     故把 cell 本身作为新组件的父级（与 Mol* 内部 applyTheme 传 I.cell 一致）。 */
  function structureSelector() {
    return structureCell();
  }

  function structureEntryOf(cell) {
    var ref = cell && cell.transform ? cell.transform.ref : null;
    if (!ref) return null;
    var structures = viewer.plugin.managers.structure.hierarchy.current.structures;
    for (var i = 0; i < structures.length; i++) {
      if (structures[i].cell.transform.ref === ref) return structures[i];
    }
    return null;
  }

  /* 小结构经默认 preset 自带元素配色表示；只删表示、保留结构本身（Loci 仍要用） */
  function clearDefaultReprs(cell) {
    var entry = structureEntryOf(cell);
    if (!entry || !entry.components) return Promise.resolve();
    var builder = viewer.plugin.state.data.build();
    var removed = 0;
    for (var i = 0; i < entry.components.length; i++) {
      var reprs = entry.components[i].representations || [];
      for (var j = 0; j < reprs.length; j++) {
        if (!reprs[j]) continue;
        builder.delete(reprs[j].cell || reprs[j]);
        removed += 1;
      }
    }
    if (!removed) return Promise.resolve();
    return Promise.resolve(builder.commit({ canUndo: false })).catch(function (error) {
      console.warn("[dii] default representation cleanup failed", error);
    });
  }

  function unitOf(structure, atomIndex) {
    for (var i = 0; i < structure.units.length; i++) {
      var elements = structure.units[i].elements;
      for (var j = 0; j < elements.length; j++) {
        if (elements[j] === atomIndex) {
          return { unit: structure.units[i], index: j };
        }
      }
    }
    return null;
  }

  function locus(unit, indices) {
    return { unit: unit, indices: indices };
  }

  function contactLoci(structure, contact) {
    var first = unitOf(structure, contact[0]);
    var second = unitOf(structure, contact[1]);
    if (!first || !second) return null;
    var elements = [];
    if (first.unit === second.unit) {
      elements.push(locus(first.unit,
        new Int32Array([first.index, second.index])));
    } else {
      elements.push(locus(first.unit, new Int32Array([first.index])));
      elements.push(locus(second.unit, new Int32Array([second.index])));
    }
    return molstar.lib.structure.StructureElement.Loci(structure, elements);
  }

  /* 小结构整条 Loci：unit 内局部序号 0..n-1，不做任何全局序号换算 */
  function wholeStructureLoci(structure) {
    if (!structure) return null;
    var elements = [];
    for (var i = 0; i < structure.units.length; i++) {
      var unit = structure.units[i];
      if (!unit.elements.length) continue;
      elements.push(locus(unit, unitIndices(unit.elements.length)));
    }
    if (!elements.length) return null;
    return molstar.lib.structure.StructureElement.Loci(structure, elements);
  }

  function unitIndices(count) {
    var indices = new Int32Array(count);
    for (var i = 0; i < count; i++) indices[i] = i;
    return indices;
  }

  /* 小结构逐原子 entries（contact_pdb ≤2 原子；a↔b 顺序不影响连线） */
  function structureAtomEntries(structure) {
    var entries = [];
    for (var i = 0; i < structure.units.length; i++) {
      var unit = structure.units[i];
      for (var j = 0; j < unit.elements.length; j++) {
        entries.push(locus(unit, new Int32Array([j])));
      }
    }
    return entries;
  }

  /* camera manager 挂在 plugin 上（molstar.Viewer 实例不暴露 managers/setLayout） */
  function cameraManager() {
    var plugin = viewer && viewer.plugin;
    return plugin && plugin.managers ? plugin.managers.camera : null;
  }

  function focusLoci(loci) {
    var camera = cameraManager();
    if (!camera || !loci) return;
    try {
      camera.focusLoci(loci);
    } catch (error) {
      /* 相机聚焦失败不影响帧渲染，保留 Mol* 默认视角 */
    }
  }

  function focusMainOverlay(payload, structure) {
    if (!structure) return;
    focusLoci(payload.contact ? contactLoci(structure, payload.contact)
      : wholeStructureLoci(structure));
  }

  function atomLoci(structure, atomIndex) {
    var located = unitOf(structure, atomIndex);
    if (!located) return null;
    var elements = molstar.lib.structure.StructureElement;
    return elements.Loci(structure, [locus(located.unit,
      new Int32Array([located.index]))]);
  }

  function addDistance(first, second, text) {
    var options = {
      lineParams: { dashLength: DASH_LENGTH, linesSize: LINE_SIZE }
    };
    if (text) options.customText = text;
    var measurement = viewer.plugin.managers.structure.measurement;
    return measurement.addDistance(first, second, options).then(function (result) {
      if (result) contacts.push(result);
    });
  }

  /* 回退虚线：主结构 + payload.contact（旧行为） */
  function drawContact(structure, contact, text) {
    if (!structure || !contact) return Promise.resolve();
    var first = atomLoci(structure, contact[0]);
    var second = atomLoci(structure, contact[1]);
    if (!first || !second) return Promise.resolve();
    return addDistance(first, second, text);
  }

  /* 默认虚线：contact_pdb 的两个原子 + contact_label */
  function addContactPdb(payload) {
    var state = { loci: null, failed: false };
    if (!payload.contact_pdb) return Promise.resolve(state);
    return loadStructureCell(payload.contact_pdb).then(function (cell) {
      return clearDefaultReprs(cell).then(function () {
        var structure = structureDataOf(cell);
        var atoms = structure ? structureAtomEntries(structure) : [];
        if (atoms.length < 2) return state;
        var elements = molstar.lib.structure.StructureElement;
        var first = elements.Loci(structure, [atoms[0]]);
        var second = elements.Loci(structure, [atoms[1]]);
        return addDistance(first, second, contactText(payload)).then(function () {
          state.loci = lociOf(structure, [atoms[0], atoms[1]]);
          return state;
        });
      });
    }).catch(function (error) {
      console.warn("[dii] contact structure failed", error);
      state.failed = true;
      return state;
    });
  }

  function hexToInt(hex) {
    return parseInt(String(hex).replace("#", ""), 16);
  }

  function groupColor(index) {
    return GROUP_COLORS[index % GROUP_COLORS.length];
  }

  function groupLabel(group, index) {
    return group && group.label ? String(group.label) : "group " + (index + 1);
  }

  function bucketOf(buckets, unit) {
    for (var i = 0; i < buckets.length; i++) {
      if (buckets[i].unit === unit) return buckets[i];
    }
    var bucket = { unit: unit, indices: [] };
    buckets.push(bucket);
    return bucket;
  }

  /* 合并 entries 为一条 Loci：同一 unit 的原子必须合并成一条升序 indices */
  function lociOf(structure, entries) {
    var buckets = [];
    for (var i = 0; i < entries.length; i++) {
      var bucket = bucketOf(buckets, entries[i].unit);
      var indices = entries[i].indices;
      for (var j = 0; j < indices.length; j++) {
        if (bucket.indices.indexOf(indices[j]) < 0) {
          bucket.indices.push(indices[j]);
        }
      }
    }
    if (!buckets.length) return null;
    var elements = [];
    for (var k = 0; k < buckets.length; k++) {
      elements.push(locus(buckets[k].unit, sortedIndices(buckets[k].indices)));
    }
    return molstar.lib.structure.StructureElement.Loci(structure, elements);
  }

  /* 回退路径：主结构的 groups.atoms 与 contact 同坐标系（旧行为，保留） */
  function groupLoci(structure, atoms) {
    if (!atoms || !atoms.length) return null;
    var entries = [];
    for (var i = 0; i < atoms.length; i++) {
      var located = unitOf(structure, atoms[i]);
      if (!located) continue;
      entries.push(locus(located.unit, new Int32Array([located.index])));
    }
    return entries.length ? lociOf(structure, entries) : null;
  }

  /* Bundle.fromLoci 走 SortedArray，indices 必须升序 */
  function sortedIndices(indices) {
    indices.sort(function (left, right) { return left - right; });
    return new Int32Array(indices);
  }

  function groupParams(index) {
    return {
      type: GROUP_TYPE,
      typeParams: {
        sizeFactor: GROUP_SIZE_FACTOR,
        sizeAspectRatio: GROUP_SIZE_ASPECT
      },
      color: GROUP_COLOR_THEME,
      colorParams: { value: hexToInt(groupColor(index)) }
    };
  }

  /* addRepresentation 首参须是 cell/selector（裸 Loci 无法解析），因此先把 Loci
     转成 StructureElement.Bundle 建一个临时 component；tryCreateComponent 是
     async，必须等它给出的 selector，不能把 Promise 直接传给 addRepresentation。 */
  function groupComponent(loci, selector, group, index) {
    if (!loci) return Promise.resolve(null);
    var bundle = molstar.lib.structure.StructureElement.Bundle.fromLoci(loci);
    var spec = {
      type: { name: GROUP_BUNDLE_TYPE, params: bundle },
      nullIfEmpty: true,
      label: groupLabel(group, index)
    };
    return Promise.resolve(viewer.plugin.builders.structure.tryCreateComponent(
      selector, spec, GROUP_TAG_PREFIX + index));
  }

  function usableComponent(component) {
    return !!(component && (component.cell || component.ref));
  }

  function addGroupRepr(loci, selector, group, index) {
    return groupComponent(loci, selector, group, index)
      .then(function (component) {
        if (!usableComponent(component)) {
          console.warn("[dii] group component failed", index);
          return false;
        }
        var builders = viewer.plugin.builders.structure.representation;
        return builders.addRepresentation(
          component, groupParams(index), { tag: GROUP_TAG_PREFIX + index }
        ).then(function (result) {
          if (result) {
            groupReprs.push({ component: component, representation: result });
          }
          return !!result;
        });
      });
  }

  function runGroupSteps(steps) {
    var state = { failures: 0, added: 0 };
    var chain = Promise.resolve();
    for (var i = 0; i < steps.length; i++) {
      chain = chain.then(groupStep(steps[i], i, state));
    }
    return chain.then(function () {
      return { failed: state.failures > 0, added: state.added };
    });
  }

  function groupStep(step, index, state) {
    return function () {
      return Promise.resolve().then(step).then(function (ok) {
        if (ok) {
          state.added += 1;
          return;
        }
        state.failures += 1;
        console.warn("[dii] group representation failed", index);
      }).catch(function (error) {
        state.failures += 1;
        console.warn("[dii] group component failed", index, error);
      });
    };
  }

  /* 默认路径：一条 highlight PDB = 一个基团，整条小结构加球棍 */
  function addHighlightReprs(payload) {
    var highlights = payload.highlights || [];
    var steps = [];
    for (var i = 0; i < highlights.length; i++) {
      steps.push(highlightStep(highlights[i], i));
    }
    return runGroupSteps(steps).then(function (result) {
      highlightGroups = result.added;
      return result.failed;
    });
  }

  function highlightStep(highlight, index) {
    return function () {
      var pdb = highlight && highlight.pdb;
      if (!pdb) {
        console.warn("[dii] highlight group has no pdb", index);
        return false;
      }
      return loadStructureCell(pdb).then(function (cell) {
        var structure = structureDataOf(cell);
        if (!structure) throw new Error("empty highlight structure");
        return clearDefaultReprs(cell).then(function () {
          return addGroupRepr(wholeStructureLoci(structure), cell, highlight, index);
        });
      });
    };
  }

  /* 回退路径（highlights 缺失）：主结构 + groups.atoms + 旧 unitOf 定位 */
  function addGroupReprs(payload) {
    var groups = payload.groups || [];
    var structure = groups.length ? structureData() : null;
    var selector = groups.length ? structureSelector() : null;
    if (!structure || !selector) return Promise.resolve(false);
    var steps = [];
    for (var i = 0; i < groups.length; i++) {
      steps.push(groupReprStep(structure, selector, groups[i], i));
    }
    return runGroupSteps(steps).then(function (result) {
      return result.failed;
    });
  }

  function groupReprStep(structure, selector, group, index) {
    return function () {
      return addGroupRepr(groupLoci(structure, group.atoms), selector, group, index);
    };
  }

  function renderGroupLegend(groups) {
    if (!groupLegendEl) return;
    while (groupLegendEl.firstChild) {
      groupLegendEl.removeChild(groupLegendEl.firstChild);
    }
    for (var i = 0; i < groups.length; i++) {
      var item = document.createElement("span");
      item.className = GROUP_ITEM_CLASS;
      var dot = document.createElement("span");
      dot.className = GROUP_DOT_CLASS;
      dot.style.backgroundColor = groupColor(i);
      var label = document.createElement("span");
      label.className = GROUP_LABEL_CLASS;
      label.textContent = groupLabel(groups[i], i);
      item.appendChild(dot);
      item.appendChild(label);
      groupLegendEl.appendChild(item);
    }
  }

  function copyShallow(source) {
    var target = {};
    for (var key in source) {
      if (Object.prototype.hasOwnProperty.call(source, key)) {
        target[key] = source[key];
      }
    }
    return target;
  }

  function toggleSequenceClass(hide) {
    if (!host || !host.classList) return;
    if (hide) {
      host.classList.add(NO_SEQUENCE_CLASS);
    } else {
      host.classList.remove(NO_SEQUENCE_CLASS);
    }
  }

  /* regionState.top 是 bundle 内隐藏序列区的布局键（见 molstar.js 的
     layoutVisibilityClassName）；无该布局 API 时退回容器 CSS 类。 */
  function applySequence(visible) {
    var plugin = viewer && viewer.plugin;
    var layout = plugin && plugin.layout;
    var state = layout && layout.state;
    if (!state || !state.regionState ||
        typeof layout.updateProps !== "function") {
      toggleSequenceClass(!visible);
      return;
    }
    var regionState = copyShallow(state.regionState);
    regionState.top = visible ? REGION_FULL : REGION_HIDDEN;
    try {
      layout.updateProps({ regionState: regionState });
      toggleSequenceClass(false);
    } catch (error) {
      toggleSequenceClass(!visible);
    }
  }

  function readJson(response) {
    return response.json().then(function (data) {
      if (!response.ok) {
        var detail = data && (data.detail || data.error)
          ? (data.detail || data.error) : "HTTP " + response.status;
        throw new Error(String(detail));
      }
      return data;
    });
  }

  function pickNative(input) {
    setStatus("waiting for file selection…", false);
    fetch(PICK_ENDPOINT, {
      method: "POST",
      headers: { "Content-Type": "application/x-www-form-urlencoded" },
      body: PICK_BODY
    }).then(readJson).then(function (data) {
      var path = data && data.path ? String(data.path) : "";
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
    var input = document.getElementById(button.getAttribute(BROWSE_TARGET_ATTR));
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

  function initViewer() {
    if (!window.molstar || !molstar.Viewer) {
      setStatus("Mol* viewer is not available", true);
      return;
    }
    molstar.Viewer.create(host, VIEWER_OPTIONS).then(function (instance) {
      viewer = instance;
      ready = true;
      applySequence(!showSequence || showSequence.checked);
      loadFrame(0);
    }).catch(function (error) {
      setStatus("viewer init failed: " + error.message, true);
    });
  }

  function selectStructure(event) {
    event.preventDefault();
    var topology = topologyInput.value.trim();
    var trajectory = trajectoryInput.value.trim();
    if (!topology || !trajectory) {
      setStatus("topology and trajectory are required", true);
      return;
    }
    openButton.disabled = true;
    setStatus("opening structure…", false);
    var body = "topology=" + encodeURIComponent(topology) +
      "&trajectory=" + encodeURIComponent(trajectory);
    fetch("/api/structure/select", {
      method: "POST",
      headers: { "Content-Type": "application/x-www-form-urlencoded" },
      body: body
    }).then(function (response) {
      return response.json().then(function (data) {
        if (!response.ok) {
          throw new Error(data && data.detail ? String(data.detail)
            : "HTTP " + response.status);
        }
        return data;
      });
    }).then(function () {
      window.location.reload();
    }).catch(function (error) {
      openButton.disabled = false;
      setStatus(error.message, true);
    });
  }

  if (slider) {
    slider.addEventListener("input", function () {
      var frame = parseInt(slider.value, 10);
      if (isNaN(frame)) return;
      updateFrameLabel(frame);
      highlightFrame(frame);
      if (pendingTimer) clearTimeout(pendingTimer);
      pendingTimer = setTimeout(function () {
        pendingTimer = null;
        loadFrame(frame);
      }, DEBOUNCE_MS);
    });
  }

  if (showWater) {
    showWater.addEventListener("change", function () {
      clearCache();
      if (ready) loadFrame(latestFrame);
    });
  }

  if (showSequence) {
    showSequence.addEventListener("change", function () {
      applySequence(showSequence.checked);
    });
  }

  for (var i = 0; i < browseButtons.length; i++) {
    browseButtons[i].addEventListener("click", onBrowseClick);
  }

  if (form) form.addEventListener("submit", selectStructure);

  highlightFrame(0);
  if (cfg.has_structure) {
    initViewer();
  } else {
    setStatus("select topology and trajectory to load the structure", false);
  }
})();
