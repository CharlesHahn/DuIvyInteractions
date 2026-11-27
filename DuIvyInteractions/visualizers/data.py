# -*- coding: utf-8 -*-
"""Web 可视化数据层：数据源解析、总览、类型表格、下钻、降采样、最近数据源。

纯函数层（spec FR-009）：只依赖 core/io，零 web 依赖，全部可单元测试。
"""

import json
import os
import threading
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import h5py
import numpy as np

from ..core.datas import Group, Interaction
from ..io import EXPORTER_CLASSES, load_interactions
from ..io.interaction_exporter import InteractionExporter

MAX_SERIES_POINTS = 2000  # 下钻曲线降采样上限
MATRIX_FRAME_LIMIT = 200  # 帧视图最多列出的活跃对数
HIST_BINS = 24          # 指标分布直方图分箱数
FRAME_PAGE = 500        # 逐帧明细分页大小
OCCUPANCY_RANGE = (0.0, 1.0)  # 占有率直方图固定值域
STRING_TOP_N = 5        # 字符串指标展示的取值个数
VALUE_DIGITS = 4        # 指标取值保留小数位
ACTIVITY_BUCKETS = 120  # 每帧活跃数曲线的最大柱数
TABLE_PAGE_SIZE = 50    # 聚合表分页大小
SUMMARY_DIGITS = 4      # CSV 汇总数值小数位
BYTES_PER_KB = 1024
PT_METRIC = "pistacking_type"  # π-堆积类型指标名
PT_TYPES = ("P", "T", "N")     # π-堆积类型取值（并列时按此顺序优先）
METRIC_UNITS = {
    "distance": "Å", "offset": "Å", "dist_dw": "Å", "dist_wa": "Å",
    "dist_h_ow": "Å", "dist_d_ow": "Å", "dist_hw_a": "Å", "dist_ow_a": "Å",
    "angle": "°", "theta": "°", "omega": "°",
    "acc_angle": "°", "don_angle": "°",
    "angle_dh_ow": "°", "angle_ow_hw_a": "°", "angle_h_ow_a": "°",
}
WATER_GROUP_TYPE = "water"   # 水桥中的桥接水基团类型
METAL_GROUP_TYPE = "metal"   # 金属离子基团类型
METAL_COORDINATION_TYPE = "metal_coordination"  # 金属配位相互作用类型
DISTANCE_METRIC = "distance"  # 热点表平均距离所用指标名
DISTANCE_METRIC_BY_TYPE = {"water_bridge": "dist_d_ow"}  # 水桥无 distance，用供体–水氧
RECENT_LIMIT = 10       # 最近数据源记忆条数
RECENT_PATH = Path(os.environ.get("DII_CACHE_HOME", Path.home() / ".cache" / "duivyinteractions")) / "web_recent.json"
_RECENT_LOCK = threading.Lock()


def parse_datasource(path: str) -> List[Tuple[str, str, Interaction]]:
    """解析数据源（单个 h5 文件或输出目录）→ [(显示名, h5路径, Interaction)]。

    目录扫描全部 *.h5；无法解析的文件跳过并告警；同一文件内的全部 Interaction 均载入。
    """
    p = Path(path).expanduser().resolve()
    if p.is_file():
        candidates = [p]
    elif p.is_dir():
        candidates = sorted(p.glob("*.h5"))
    else:
        raise ValueError(f"path not found: {path}")
    items: List[Tuple[str, str, Interaction]] = []
    for h5 in candidates:
        for interaction in _load_interactions(h5):
            key = _unique_name(interaction.interaction_type, h5.name, items)
            items.append((key, str(h5), interaction))
    return items


def _load_interactions(h5: Path) -> List[Interaction]:
    """加载一个 h5 的全部 Interaction；无法解析时告警并返回空列表。"""
    try:
        loaded = load_interactions(str(h5))
    except (OSError, ValueError, KeyError, TypeError, AttributeError) as e:
        print(f"[WARN] viz: skip unreadable h5 {h5.name}: {e}")
        return []
    if not loaded:
        print(f"[WARN] viz: skip empty h5 {h5.name}")
    return loaded


def _unique_name(base: str, filename: str,
                 items: List[Tuple[str, str, Interaction]]) -> str:
    """生成唯一显示名（重复类型追加 @文件名，仍重复则追加序号）。"""
    key = base
    if any(key == name for name, _, _ in items):
        key = f"{base}@{filename}"
    index = 2
    while any(key == name for name, _, _ in items):
        key = f"{base}@{filename}~{index}"
        index += 1
    return key


def merge_batches(batches: Sequence[Sequence[Tuple[str, str, Interaction]]]
                  ) -> List[Tuple[str, str, Interaction]]:
    """按批次顺序合并 items；跨批次重名时用 _unique_name 追加 @文件名/序号。

    首批显示名已由 parse_datasource 保证批内唯一，原样透传。
    """
    merged: List[Tuple[str, str, Interaction]] = []
    for position, batch in enumerate(batches):
        if position == 0:
            merged.extend(batch)
            continue
        for _, h5_path, interaction in batch:
            key = _unique_name(interaction.interaction_type,
                               Path(h5_path).name, merged)
            merged.append((key, h5_path, interaction))
    return merged


def type_summary(items: List[Tuple[str, str, Interaction]]) -> List[Dict]:
    """导航用类型摘要：显示名、基础类型与对数。"""
    return [{"type": name, "data_type": it.interaction_type, "n_pairs": it.n_pairs}
            for name, _, it in items]


def _occupancy(interaction: Interaction) -> np.ndarray:
    """每对占有率数组（0 帧时全 0）。"""
    if interaction.n_frames == 0:
        return np.zeros(interaction.n_pairs, dtype=float)
    return np.asarray(interaction.occupancy(), dtype=float)


def overview(items: List[Tuple[str, str, Interaction]]) -> List[Dict]:
    """每类型卡片：类型名、对数、帧数、flags 名单与 occupancy 分布。"""
    cards = []
    for name, h5_path, it in items:
        summary = _summary(_occupancy(it), HIST_BINS, OCCUPANCY_RANGE)
        if it.n_frames == 0:
            summary["count"] = 0
        cards.append({"type": name, "data_type": it.interaction_type,
                      "h5_path": h5_path,
                      "n_pairs": it.n_pairs, "n_frames": it.n_frames,
                      "flags": sorted(it.flags.keys()), "occupancy": summary})
    return cards


def get_exporter(interaction: Interaction) -> InteractionExporter:
    """按 interaction_type 取导出器（pair 标签等复用 dii export 的单一来源）。"""
    cls = EXPORTER_CLASSES.get(interaction.interaction_type)
    if cls is None:
        raise ValueError(f"unknown interaction type: '{interaction.interaction_type}'")
    return cls()


def export_table(interaction: Interaction, search: str = "",
                 sort_by: str = "index", sort_dir: str = "asc",
                 min_occupancy: float = 0.0) -> Dict:
    """CSV 汇总表（列与 dii export summary 一致）。"""
    rows = [row for row in _build_rows(interaction, get_exporter(interaction))
            if row["occupancy"] > 0 and row["occupancy"] >= min_occupancy]
    rows = _sort_rows(_filter_rows(rows, search), sort_by, sort_dir)
    columns = _summary_columns(interaction)
    return {"columns": columns,
            "rows": [{column: _summary_value(row, column)
                      for column in columns} for row in rows],
            "total": len(rows)}


def _summary_columns(interaction: Interaction) -> List[str]:
    """CSV 列：pair_label/occupancy + 数值指标 avg/std + flags avg。"""
    columns = ["pair_label", "occupancy"]
    for name in numeric_metrics(interaction):
        columns += [f"avg_{name}", f"std_{name}"]
    columns += [f"avg_{name}" for name in interaction.flags]
    return columns


def _summary_value(row: Dict, column: str) -> str:
    """CSV 单元格：数值保留 SUMMARY_DIGITS 位，空值输出空串。"""
    if column == "pair_label":
        return row["pair_label"]
    value = row.get(column, "")
    if value == "" and f"flag_{column[4:]}" in row:
        value = row[f"flag_{column[4:]}"]
    if value == "" or value is None:
        return ""
    return f"{float(value):.{SUMMARY_DIGITS}f}"


def aggregate_table(interaction: Interaction, search: str = "",
                    sort_by: str = "occupancy", sort_dir: str = "desc",
                    page: int = 1, min_occupancy: float = 0.0) -> Dict:
    """聚合表：pair_label/occupancy/avg/std/flags，按页返回。"""
    rows = _table_rows(interaction, search, sort_by, sort_dir, min_occupancy)
    total = len(rows)
    pages = max(1, (total + TABLE_PAGE_SIZE - 1) // TABLE_PAGE_SIZE)
    page = min(max(1, page), pages)
    start = (page - 1) * TABLE_PAGE_SIZE
    return {"columns": _columns(interaction),
            "rows": rows[start:start + TABLE_PAGE_SIZE],
            "total": total, "pages": pages, "page": page,
            "units": _units(interaction)}


def _table_rows(interaction: Interaction, search: str, sort_by: str,
                sort_dir: str, min_occupancy: float) -> List[Dict]:
    """过滤（label + 占有率）并排序后的全部行。"""
    rows = [r for r in _build_rows(interaction, get_exporter(interaction))
            if r["occupancy"] >= min_occupancy]
    return _sort_rows(_filter_rows(rows, search), sort_by, sort_dir)


def _build_rows(interaction: Interaction, exporter: InteractionExporter) -> List[Dict]:
    """全 pair 的行数据（数值指标 avg/std、字符串指标 mode、flags 比例）。"""
    occ = _occupancy(interaction)
    return [_pair_row(interaction, i, float(occ[i]), exporter)
            for i in range(interaction.n_pairs)]


def _pair_row(interaction: Interaction, idx: int, occupancy: float,
              exporter: InteractionExporter) -> Dict:
    """单个 pair 的表格行。"""
    row: Dict = {"index": idx, "occupancy": occupancy,
                 "pair_label": exporter.get_pair_label(interaction, idx)}
    for name in numeric_metrics(interaction):
        active = _active_values(interaction.metrics[name][idx],
                                interaction.existence[idx])
        row[f"avg_{name}"] = float(active.mean()) if active.size else ""
        row[f"std_{name}"] = float(active.std()) if active.size else ""
    for name in _string_metrics(interaction):
        row[f"mode_{name}"] = _mode(interaction.metrics[name][idx],
                                    interaction.existence[idx])
    for name, arr in interaction.flags.items():
        active = arr[idx][interaction.existence[idx]]
        row[f"flag_{name}"] = float(np.mean(active)) if active.size else ""
    return row


def _mode(values: Any, existence: Any) -> str:
    """字符串指标在存在帧上的众数（无活跃帧返回空串）。"""
    active = np.asarray(values)[np.asarray(existence, dtype=bool)].astype(str)
    if active.size == 0:
        return ""
    uniq, counts = np.unique(active, return_counts=True)
    return str(uniq[int(np.argmax(counts))])


def _filter_rows(rows: List[Dict], search: str) -> List[Dict]:
    """按 pair_label 大小写不敏感搜索。"""
    if not search:
        return rows
    needle = search.lower()
    return [r for r in rows if needle in r["pair_label"].lower()]


SORT_ALIASES = {"label": "pair_label"}  # 页面排序名 → 行字段名


def _sort_rows(rows: List[Dict], sort_by: str, sort_dir: str) -> List[Dict]:
    """按列名排序（缺失/NaN 恒排末尾），未知列回退 index。"""
    sort_by = SORT_ALIASES.get(sort_by, sort_by)
    key = sort_by if rows and sort_by in rows[0] else "index"
    present = [r for r in rows if not _is_missing(r[key])]
    missing = [r for r in rows if _is_missing(r[key])]
    present.sort(key=lambda r: r[key], reverse=(sort_dir == "desc"))
    return present + missing


def _is_missing(value: Any) -> bool:
    """空串 / None / NaN 视为缺失值。"""
    if value is None or value == "":
        return True
    return isinstance(value, (float, np.floating)) and np.isnan(value)


def _columns(interaction: Interaction) -> List[str]:
    """表格列名顺序：index/pair_label/occupancy + 数值指标 + 字符串指标 + flags。"""
    numeric = numeric_metrics(interaction)
    cols = ["index", "pair_label", "occupancy"]
    cols += [f"avg_{name}" for name in numeric]
    cols += [f"std_{name}" for name in numeric]
    cols += [f"mode_{name}" for name in _string_metrics(interaction)]
    cols += [f"flag_{name}" for name in interaction.flags]
    return cols


def pair_drilldown(interaction: Interaction, idx: int) -> Dict:
    """单对下钻数据：元信息 + times + existence(0/1) + 各指标序列。"""
    if not 0 <= idx < interaction.n_pairs:
        raise ValueError(f"pair index out of range: {idx}")
    times = _times(interaction)
    edges = _bucket_edges(len(times)) if len(times) > MAX_SERIES_POINTS else None
    xs = _bucket_times(times, edges) if edges is not None else times
    data = {"idx": idx,
            "pair_label": get_exporter(interaction).get_pair_label(interaction, idx),
            "groups": _group_meta(interaction, idx), "times": xs,
            "downsampled": edges is not None,
            "existence": _existence(interaction, idx, edges),
            "series": _series(interaction, idx, xs, edges)}
    if PT_METRIC in interaction.metrics:
        data["types"] = _types(interaction, idx, edges)
    return data


def pair_types(interaction: Interaction, idx: int) -> List[str]:
    """逐帧 π-堆积类型（无该指标或越界时返回空列表）。"""
    if PT_METRIC not in interaction.metrics:
        return []
    if not 0 <= idx < interaction.n_pairs:
        return []
    return np.asarray(interaction.metrics[PT_METRIC][idx]).astype(str).tolist()


def _types(interaction: Interaction, idx: int,
           edges: Optional[np.ndarray]) -> Dict:
    """π-堆积类型序列（超长时取活跃帧桶内众数），与 times 等长。"""
    values = np.asarray(interaction.metrics[PT_METRIC][idx]).astype(str)
    if edges is None:
        return {"name": PT_METRIC, "values": values.tolist()}
    existence = np.asarray(interaction.existence[idx], dtype=bool)
    return {"name": PT_METRIC, "values": [
        _bucket_mode(values[start:end], existence[start:end])
        for start, end in zip(edges[:-1], edges[1:])]}


def _bucket_mode(values: np.ndarray, active: np.ndarray) -> str:
    """活跃帧上的众数（并列时按 PT_TYPES 顺序优先），无活跃帧返回空串。"""
    present = np.asarray(values)[np.asarray(active, dtype=bool)].astype(str)
    if present.size == 0:
        return ""
    uniq, counts = np.unique(present, return_counts=True)
    tied = set(uniq[counts == counts.max()].tolist())
    for name in PT_TYPES:
        if name in tied:
            return name
    return str(uniq[int(np.argmax(counts))])


def _bucket_edges(n_frames: int) -> np.ndarray:
    """超长帧的桶边界（桶数 = MAX_SERIES_POINTS）。"""
    return np.linspace(0, n_frames, MAX_SERIES_POINTS + 1).astype(int)


def _bucket_times(times: List[float], edges: np.ndarray) -> List[float]:
    """每个桶的中心时间。"""
    return [times[(start + end) // 2]
            for start, end in zip(edges[:-1], edges[1:])]


def _existence(interaction: Interaction, idx: int,
               edges: Optional[np.ndarray]) -> List[int]:
    """existence 序列（超长时桶内任一帧活跃即为 1）。"""
    row = interaction.existence[idx]
    if edges is None:
        return row.astype(int).tolist()
    return [int(row[start:end].any())
            for start, end in zip(edges[:-1], edges[1:])]


def _times(interaction: Interaction) -> List[float]:
    """全部帧时间（无帧时为空列表）。"""
    return np.asarray(interaction.times if interaction.n_frames else [],
                      dtype=float).tolist()


def _series(interaction: Interaction, idx: int, xs: List[float],
            edges: Optional[np.ndarray]) -> Dict[str, Dict]:
    """各数值指标序列（超长时桶内均值，字符串指标跳过）。"""
    series = {}
    for name, arr in interaction.metrics.items():
        if not _is_numeric(arr):
            continue
        values = np.asarray(arr[idx], dtype=float)
        if edges is None:
            series[name] = {"times": xs, "values": [float(v) for v in values]}
            continue
        series[name] = {"times": xs, "values": [
            _bucket_mean(values[start:end])
            for start, end in zip(edges[:-1], edges[1:])]}
    return series


def _bucket_mean(values: np.ndarray) -> float:
    """桶内有限值均值（无有限值返回 NaN）。"""
    finite = values[np.isfinite(values)]
    return float(finite.mean()) if finite.size else float("nan")


def pair_summary(interaction: Interaction, idx: int) -> Dict:
    """单对摘要：占有率与活跃帧上的指标均值/标准差。"""
    if not 0 <= idx < interaction.n_pairs:
        raise ValueError(f"pair index out of range: {idx}")
    values = {}
    for name in numeric_metrics(interaction):
        active = _active_values(interaction.metrics[name][idx],
                                interaction.existence[idx])
        values[name] = {"mean": float(active.mean()) if active.size else None,
                        "std": float(active.std()) if active.size else None,
                        "count": int(active.size),
                        "unit": METRIC_UNITS.get(name, "")}
    return {"occupancy": float(_occupancy(interaction)[idx]), "metrics": values,
            "string_metrics": _string_summary(interaction, idx),
            "flags": {name: {"active": int(arr[idx].sum()),
                             "frames": interaction.n_frames}
                      for name, arr in interaction.flags.items()}}


def _string_summary(interaction: Interaction, idx: int) -> Dict[str, List[Dict]]:
    """该对活跃帧上字符串指标的取值频次（最多 STRING_TOP_N 个）。"""
    summary = {}
    for name in _string_metrics(interaction):
        column = np.asarray(interaction.metrics[name][idx])
        values = column[interaction.existence[idx]]
        if not values.size:
            continue
        unique, counts = np.unique(values, return_counts=True)
        order = np.argsort(-counts)[:STRING_TOP_N]
        summary[name] = [{"value": str(unique[i]), "count": int(counts[i])}
                         for i in order]
    return summary


def _group_meta(interaction: Interaction, idx: int) -> List[Dict]:
    """基团元信息：残基、分子与原子级细节（名字/力场类型/元素/电荷/索引）。"""
    metas = []
    for g in interaction.groups[idx]:
        metas.append({
            "group_id": g.group_id, "group_type": g.group_type,
            "residue_name": g.residue_name, "residue_id": g.residue_id,
            "molecule": g.molecule,
            "atoms": [{"name": a.atom_name, "type": a.atom_type,
                       "element": a.atom_element, "charge": float(a.atom_charge),
                       "global_idx": int(a.atom_global_idx)}
                      for a in g.atoms],
        })
    return metas


def _recent_entries() -> List[str]:
    """读取最近数据源文件（不存在或损坏时视为空）。"""
    try:
        return json.loads(RECENT_PATH.read_text())[:RECENT_LIMIT]
    except (OSError, ValueError, TypeError, KeyError):
        return []


def load_recent() -> List[str]:
    """最近使用的数据源列表。"""
    with _RECENT_LOCK:
        return _recent_entries()


def add_recent(path: str) -> List[str]:
    """记录一次数据源使用（去重置顶，最多 RECENT_LIMIT 条）。"""
    with _RECENT_LOCK:
        entries = [p for p in _recent_entries() if p != path]
        entries.insert(0, path)
        entries = entries[:RECENT_LIMIT]
        try:
            RECENT_PATH.parent.mkdir(parents=True, exist_ok=True)
            RECENT_PATH.write_text(json.dumps(entries))
        except OSError as e:
            print(f"[WARN] viz: cannot persist recent datasource: {e}")
        return entries


def runs_from_mask(mask: Any) -> List[List[int]]:
    """布尔序列 → 连续 True 的 [start, end) 区段（矩阵行压缩）。"""
    arr = np.asarray(mask, dtype=np.int8).ravel()
    if arr.size == 0:
        return []
    edges = np.diff(np.concatenate(([0], arr, [0])))
    starts, ends = np.flatnonzero(edges == 1), np.flatnonzero(edges == -1)
    return [[int(s), int(e)] for s, e in zip(starts, ends)]


def numeric_metrics(interaction: Interaction) -> List[str]:
    """数值型指标名（跳过字符串指标，如 pistacking_type）。"""
    return [name for name, arr in interaction.metrics.items() if _is_numeric(arr)]


def _string_metrics(interaction: Interaction) -> List[str]:
    """字符串型指标名。"""
    return [name for name, arr in interaction.metrics.items() if not _is_numeric(arr)]


def matrix_view(interaction: Interaction, search: str = "",
                min_occupancy: float = 0.0, sort_by: str = "occupancy",
                sort_dir: str = "desc") -> Dict:
    """矩阵视图：行 = 过滤/排序后的 pair，existence 以 RLE 区段下发。"""
    return {
        "total_pairs": interaction.n_pairs,
        "n_frames": interaction.n_frames,
        "times": np.asarray(interaction.times, dtype=float).tolist(),
        "rows": _sort_matrix_rows(
            _matrix_rows(interaction, search, min_occupancy), sort_by, sort_dir),
    }


def _matrix_rows(interaction: Interaction, search: str,
                 min_occupancy: float) -> List[Dict]:
    """过滤后的矩阵行（idx + label + occupancy + existence 区段）。"""
    exporter = get_exporter(interaction)
    occ = _occupancy(interaction)
    rows = []
    for idx in range(interaction.n_pairs):
        label = exporter.get_pair_label(interaction, idx)
        if search and search.lower() not in label.lower():
            continue
        if float(occ[idx]) < min_occupancy:
            continue
        rows.append(_matrix_row(interaction, idx, label, float(occ[idx])))
    return rows


def _matrix_row(interaction: Interaction, idx: int, label: str,
                occupancy: float) -> Dict:
    """矩阵单行：idx/label/occupancy/runs；有 π-堆积类型时附加 typed。"""
    row: Dict = {"idx": idx, "label": label, "occupancy": occupancy,
                 "runs": runs_from_mask(interaction.existence[idx])}
    if PT_METRIC in interaction.metrics:
        row["typed"] = _typed_runs(interaction, idx)
    return row


def _typed_runs(interaction: Interaction, idx: int) -> Dict[str, List[List[int]]]:
    """该行 existence 区段按帧上类型切分（end 为开区间）。"""
    exists = np.asarray(interaction.existence[idx], dtype=bool)
    values = np.asarray(interaction.metrics[PT_METRIC][idx]).astype(str)
    names = list(PT_TYPES) + [v for v in np.unique(values) if v not in PT_TYPES]
    return {name: runs_from_mask(exists & (values == name)) for name in names}


def _sort_matrix_rows(rows: List[Dict], sort_by: str, sort_dir: str) -> List[Dict]:
    """矩阵行排序：occupancy / label / index，方向由 sort_dir 决定。"""
    if sort_by == "label":
        key = lambda row: row["label"]
    elif sort_by == "index":
        key = lambda row: row["idx"]
    else:
        key = lambda row: row["occupancy"]
    return sorted(rows, key=key, reverse=(sort_dir == "desc"))


def cell_values(interaction: Interaction, pair_idx: int, frame: int) -> Dict:
    """单 (pair, frame) 单元格：存在性、时刻、各指标值、单位与 flags。"""
    _check_bounds(interaction, pair_idx, frame)
    return {"pair": pair_idx, "frame": frame,
            "exists": bool(interaction.existence[pair_idx, frame]),
            "time": float(interaction.times[frame]),
            "values": _metric_values(interaction, pair_idx, frame),
            "units": _units(interaction),
            "flags": _flag_values(interaction, pair_idx, frame)}


def _units(interaction: Interaction) -> Dict[str, str]:
    """数值指标 → 单位（无单位为空串）。"""
    return {name: METRIC_UNITS.get(name, "") for name in numeric_metrics(interaction)}


def _flag_values(interaction: Interaction, pair_idx: int, frame: int) -> Dict:
    """该单元格上全部 flags 的布尔值。"""
    return {name: bool(arr[pair_idx, frame])
            for name, arr in interaction.flags.items()}


def _metric_values(interaction: Interaction, pair_idx: int,
                   frame: int) -> Dict:
    """该单元格全部指标值（NaN → None，字符串原样）。"""
    values = {}
    for name, arr in interaction.metrics.items():
        data = np.asarray(arr)
        value = data[pair_idx, frame]
        if not _is_numeric(data):
            values[name] = str(value)
        else:
            values[name] = None if np.isnan(float(value)) else round(float(value), VALUE_DIGITS)
    return values


def _check_bounds(interaction: Interaction, pair_idx: int, frame: int) -> None:
    """越界时抛 ValueError。"""
    if not 0 <= pair_idx < interaction.n_pairs:
        raise ValueError(f"pair index out of range: {pair_idx}")
    _check_frame(interaction, frame)


def _check_frame(interaction: Interaction, frame: int) -> None:
    """帧号越界 → ValueError。"""
    if not 0 <= frame < interaction.n_frames:
        raise ValueError(f"frame out of range: {frame}")


def frame_view(interaction: Interaction, frame: int,
               limit: int = MATRIX_FRAME_LIMIT) -> Dict:
    """某帧的活跃相互作用（按占有率降序，最多 limit 条）。"""
    _check_frame(interaction, frame)
    occ = _occupancy(interaction)
    active = np.flatnonzero(interaction.existence[:, frame])
    ordered = sorted(active.tolist(), key=lambda i: (-float(occ[i]), i))[:limit]
    return {"frame": frame, "time": float(interaction.times[frame]),
            "n_active": int(active.size), "shown": len(ordered),
            "metrics": numeric_metrics(interaction),
            "string_metrics": _string_metrics(interaction),
            "units": _units(interaction),
            "flags": sorted(interaction.flags),
            "pairs": _frame_pairs(interaction, ordered, frame, occ)}


def _frame_pairs(interaction: Interaction, ordered: List[int], frame: int,
                 occ: np.ndarray) -> List[Dict]:
    """帧视图每行：pair 标签、占有率、各指标值与 flags。"""
    if not ordered:
        return []
    exporter = get_exporter(interaction)
    return [{"idx": i, "label": exporter.get_pair_label(interaction, i),
             "occupancy": float(occ[i]),
             "values": _metric_values(interaction, i, frame),
             "flags": _flag_values(interaction, i, frame)}
            for i in ordered]


def interaction_meta(interaction: Interaction, h5_path: str = "") -> Dict:
    """单个 h5 的元数据：类型、形状、时间范围、指标、flags 与占有率分布。"""
    times = np.asarray(interaction.times, dtype=float)
    summary = _summary(_occupancy(interaction), HIST_BINS, OCCUPANCY_RANGE)
    if interaction.n_frames == 0:
        summary["count"] = 0
    return {"type": interaction.interaction_type, "h5_path": h5_path,
            "n_pairs": interaction.n_pairs, "n_frames": interaction.n_frames,
            "time_start": float(times[0]) if times.size else None,
            "time_end": float(times[-1]) if times.size else None,
            "metrics": list(interaction.metrics),
            "flags": sorted(interaction.flags),
            "occupancy": summary,
            "file": _file_info(h5_path)}


def _file_info(path: str) -> Dict:
    """h5 文件级信息：大小、修改时间、格式版本与 Interaction 数。"""
    if not path or not Path(path).is_file():
        return {}
    stat = Path(path).stat()
    info: Dict = {"size_kb": round(stat.st_size / BYTES_PER_KB, 1),
                  "mtime": datetime.fromtimestamp(stat.st_mtime).isoformat(
                      timespec="seconds")}
    try:
        with h5py.File(path, "r") as handle:
            info["format_version"] = str(handle.attrs.get("format_version", ""))
            info["n_interactions"] = int(handle.attrs.get("n_interactions", 0))
    except OSError:
        pass
    return info


def frame_activity(interaction: Interaction,
                   buckets: int = ACTIVITY_BUCKETS) -> List[int]:
    """每帧活跃对数（帧数超过 buckets 时按区间最大值聚合）。"""
    counts = np.asarray(interaction.existence, dtype=np.int64).sum(axis=0)
    if counts.size <= buckets:
        return counts.astype(int).tolist()
    edges = np.linspace(0, counts.size, buckets + 1).astype(int)
    return [int(counts[left:right].max())
            for left, right in zip(edges[:-1], edges[1:])]


def _is_numeric(array: Any) -> bool:
    """指标是否为数值型（字符串指标不能画曲线）。"""
    return np.asarray(array).dtype.kind not in ("U", "S", "O")


def metric_stats(interaction: Interaction, bins: int = HIST_BINS) -> List[Dict]:
    """各数值指标在存在帧上的分布（min/max/mean/std + 直方图）。"""
    stats = []
    for name, arr in interaction.metrics.items():
        if not _is_numeric(arr):
            continue
        stats.append({"name": name, "dtype": str(np.asarray(arr).dtype),
                      "unit": METRIC_UNITS.get(name, ""),
                      **_summary(_active_values(arr, interaction.existence), bins)})
    return stats


def _active_values(array: Any, existence: Any) -> np.ndarray:
    """存在帧上的指标值（丢弃 NaN）。"""
    values = np.asarray(array, dtype=float)[np.asarray(existence, dtype=bool)]
    return values[~np.isnan(values)]


def string_stats(interaction: Interaction, top: int = STRING_TOP_N) -> List[Dict]:
    """字符串指标在存在帧上的取值频次（最多 top 个）。"""
    stats = []
    for name in _string_metrics(interaction):
        values = np.asarray(interaction.metrics[name])[
            np.asarray(interaction.existence, dtype=bool)].astype(str)
        uniq, counts = np.unique(values, return_counts=True)
        order = np.argsort(-counts)[:top]
        stats.append({"name": name, "count": int(values.size),
                      "top": [{"value": str(uniq[i]), "count": int(counts[i])}
                              for i in order]})
    return stats


def _summary(values: np.ndarray, bins: int,
             bounds: Optional[Tuple[float, float]] = None) -> Dict:
    """一维数值的统计摘要与直方图（忽略非有限值，可指定固定分箱范围）。"""
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    counts, edges = np.histogram(finite, bins=bins, range=bounds)
    if finite.size == 0:
        return {"min": None, "max": None, "mean": None, "std": None,
                "count": 0, "hist": counts.tolist(), "edges": []}
    return {"min": float(finite.min()), "max": float(finite.max()),
            "mean": float(finite.mean()), "std": float(finite.std()),
            "count": int(finite.size), "hist": counts.tolist(),
            "edges": [round(float(edge), VALUE_DIGITS) for edge in edges]}


def pair_frames(interaction: Interaction, idx: int, start: int = 0,
                limit: int = FRAME_PAGE) -> Dict:
    """单对逐帧明细（分页）：帧号、时间、存在性、指标与 flags。"""
    if not 0 <= idx < interaction.n_pairs:
        raise ValueError(f"pair index out of range: {idx}")
    start = max(0, start)
    end = min(start + limit, interaction.n_frames)
    rows = [_frame_row(interaction, idx, frame) for frame in range(start, end)]
    return {"idx": idx, "start": start, "total": interaction.n_frames,
            "page_size": limit,
            "next_start": end if end < interaction.n_frames else None,
            "metrics": numeric_metrics(interaction),
            "string_metrics": _string_metrics(interaction),
            "units": _units(interaction),
            "flags": sorted(interaction.flags), "rows": rows}


def _frame_row(interaction: Interaction, idx: int, frame: int) -> Dict:
    """逐帧明细的一行。"""
    return {"frame": frame, "time": float(interaction.times[frame]),
            "exists": bool(interaction.existence[idx, frame]),
            "values": _metric_values(interaction, idx, frame),
            "flags": {name: bool(arr[idx, frame])
                      for name, arr in interaction.flags.items()}}


def residue_hotspots(items: Sequence[Tuple[str, str, Interaction]],
                     type_filter: Optional[str] = None,
                     include_water: bool = False) -> List[Dict[str, Any]]:
    """残基热点表：把相互作用对聚合到残基维度的数据行列表。

    残基身份严格取自 h5 的 residue_id/residue_name（不做 PDB 映射，label 为名字+编号）。
    occupancy 分母为被选中 Interaction 的最大帧数；frames_active 为跨对、跨
    interaction_type 的帧并集；types 明细以 interaction_type 为键。
    avg_distance 的指标按类型映射（水桥 dist_d_ow，其余 distance）。
    molecule 取该残基自身基团的 h5 molecule（多个值取字典序最小）；group_types
    为该残基自身基团的 group_type 去重排序；dominant_type 为 types 按
    (frames, pairs, type) 升序排序后的最后一个类型名。
    """
    selected = _hotspot_items(items, type_filter)
    n_frames = max((it.n_frames for _, _, it in selected), default=0)
    accumulator: Dict[Tuple[int, str], Dict[str, Any]] = {}
    for _, _, interaction in selected:
        for pair_idx in range(interaction.n_pairs):
            _add_pair_votes(accumulator, interaction, pair_idx,
                            include_water, n_frames)
    return _finalize_hotspots(accumulator, n_frames)


def _hotspot_items(items: Sequence[Tuple[str, str, Interaction]],
                   type_filter: Optional[str]
                   ) -> List[Tuple[str, str, Interaction]]:
    """按显示名过滤后的 items（type_filter 为 None 时全部保留）。"""
    if type_filter is None:
        return list(items)
    return [item for item in items if item[0] == type_filter]


def _group_key(group: Group) -> Tuple[int, str]:
    """残基身份键 (residue_id, residue_name)（严格取自 h5 原值）。"""
    return int(group.residue_id), str(group.residue_name)


def _pair_groups(pair: Sequence[Group], include_water: bool) -> List[Group]:
    """一对基团中参与聚合的基团（可排除桥接水）。"""
    return [group for group in pair
            if include_water or group.group_type != WATER_GROUP_TYPE]


def _residue_accumulator() -> Dict[str, Any]:
    """单个残基的累计容器。"""
    return {"frames": np.zeros(0, dtype=bool), "n_pairs": 0,
            "partners": set(), "types": {}, "dist_sum": 0.0, "dist_count": 0,
            "molecules": set(), "group_types": set()}


def _pad_frames(frames: Optional[np.ndarray], n_frames: int) -> np.ndarray:
    """布尔帧数组对齐到 n_frames（None 或偏短补 False，偏长截断）。"""
    padded = np.zeros(n_frames, dtype=bool)
    if frames is not None:
        size = min(frames.size, n_frames)
        padded[:size] = np.asarray(frames, dtype=bool)[:size]
    return padded


def _pair_distances(interaction: Interaction, pair_idx: int,
                    n_frames: int) -> np.ndarray:
    """一个基团对的逐帧距离（指标缺失或与 existence 形状不符时全 NaN）。"""
    values = np.full(n_frames, np.nan, dtype=float)
    array = _distance_array(interaction)
    if array is None:
        return values
    array = np.asarray(array, dtype=float)
    if array.shape != np.shape(interaction.existence):
        return values
    row = array[pair_idx]
    size = min(row.size, n_frames)
    values[:size] = row[:size]
    return values


def _distance_metric(interaction_type: str) -> str:
    """该类型用于热点平均距离的指标名（无类型映射时用 distance）。"""
    return DISTANCE_METRIC_BY_TYPE.get(interaction_type, DISTANCE_METRIC)


def _distance_array(interaction: Interaction) -> Optional[np.ndarray]:
    """该相互作用的距离指标数组（映射指标缺失时回落 distance）。"""
    name = _distance_metric(interaction.interaction_type)
    if name in interaction.metrics:
        return interaction.metrics[name]
    return interaction.metrics.get(DISTANCE_METRIC)


def _add_pair_votes(accumulator: Dict[Tuple[int, str], Dict[str, Any]],
                    interaction: Interaction, pair_idx: int,
                    include_water: bool, n_frames: int) -> None:
    """把一对基团的活跃帧、伙伴、类型与距离投票累加进它涉及的每个残基。"""
    groups = _pair_groups(interaction.groups[pair_idx], include_water)
    existence = _pad_frames(np.asarray(interaction.existence[pair_idx]), n_frames)
    distances = _pair_distances(interaction, pair_idx, n_frames)
    for key in dict.fromkeys(_group_key(group) for group in groups):
        record = accumulator.setdefault(key, _residue_accumulator())
        _vote_residue(record, key, groups, interaction.interaction_type,
                      existence, distances)


def _vote_residue(record: Dict[str, Any], key: Tuple[int, str],
                  groups: Sequence[Group], type_name: str,
                  existence: np.ndarray, distances: np.ndarray) -> None:
    """把一次基团对命中并入单残基累计容器（帧取并集，不求和）。"""
    record["frames"] = _pad_frames(record["frames"], existence.size) | existence
    record["n_pairs"] += 1
    record["partners"].update(_group_key(group) for group in groups
                              if _group_key(group) != key)
    _vote_group_meta(record, key, groups)
    _vote_type(record, type_name, existence)
    _vote_distance(record, existence, distances)


def _vote_group_meta(record: Dict[str, Any], key: Tuple[int, str],
                     groups: Sequence[Group]) -> None:
    """记录该残基自身基团的 molecule 与 group_type。"""
    for group in groups:
        if _group_key(group) == key:
            record["molecules"].add(str(group.molecule))
            record["group_types"].add(str(group.group_type))


def _vote_type(record: Dict[str, Any], type_name: str,
               existence: np.ndarray) -> None:
    """记录一次 interaction_type 命中（帧并集与对数）。"""
    entry = record["types"].get(type_name, {"frames": None, "pairs": 0})
    record["types"][type_name] = {
        "frames": _pad_frames(entry["frames"], existence.size) | existence,
        "pairs": entry["pairs"] + 1}


def _vote_distance(record: Dict[str, Any], existence: np.ndarray,
                   distances: np.ndarray) -> None:
    """累加活跃帧上的距离取值与有效值个数。"""
    active = distances[existence]
    record["dist_sum"] += float(np.nansum(active))
    record["dist_count"] += int(np.isfinite(active).sum())


def _finalize_hotspot(key: Tuple[int, str], record: Dict[str, Any],
                      n_frames: int) -> Dict[str, Any]:
    """单个残基的对外输出行。"""
    residue_id, residue_name = key
    frames_active = int(record["frames"].sum())
    types = [_hotspot_type_entry(name, entry, n_frames)
             for name, entry in sorted(record["types"].items())]
    return {"residue_id": residue_id, "residue_name": residue_name,
            "label": f"{residue_name}{residue_id}",
            "molecule": min(record["molecules"]) if record["molecules"] else "",
            "group_types": sorted(record["group_types"]),
            "occupancy": _hotspot_ratio(frames_active, n_frames),
            "frames_active": frames_active, "n_pairs": record["n_pairs"],
            "n_partners": len(record["partners"]), "n_types": len(types),
            "types": types, "dominant_type": _dominant_type(types),
            "avg_distance": _hotspot_distance(record)}


def _hotspot_type_entry(type_name: str, entry: Dict[str, Any],
                        n_frames: int) -> Dict[str, Any]:
    """单个 interaction_type 的明细行。"""
    frames = int(entry["frames"].sum())
    return {"type": type_name, "frames": frames,
            "occupancy": _hotspot_ratio(frames, n_frames),
            "pairs": entry["pairs"]}


def _dominant_type(types: Sequence[Dict[str, Any]]) -> Optional[str]:
    """主导类型：按 (frames, pairs, type) 升序排序后的最后一条（无类型返回 None）。"""
    if not types:
        return None
    return sorted(types, key=_dominant_key)[-1]["type"]


def _dominant_key(entry: Dict[str, Any]) -> Tuple[int, int, str]:
    """主导类型排序键：帧并集、对数、类型名。"""
    return entry["frames"], entry["pairs"], entry["type"]


def _hotspot_ratio(count: int, n_frames: int) -> float:
    """占比（帧数为 0 时返回 0.0）。"""
    return float(count) / float(n_frames) if n_frames else 0.0


def _hotspot_distance(record: Dict[str, Any]) -> Optional[float]:
    """活跃帧 distance 均值（无有效值返回 None）。"""
    if not record["dist_count"]:
        return None
    return record["dist_sum"] / record["dist_count"]


def _finalize_hotspots(accumulator: Dict[Tuple[int, str], Dict[str, Any]],
                       n_frames: int) -> List[Dict[str, Any]]:
    """全部残基输出行（按 occupancy 降序、n_pairs 降序、residue_id 升序）。"""
    rows = [_finalize_hotspot(key, record, n_frames)
            for key, record in accumulator.items()]
    return sorted(rows, key=_hotspot_sort_key)


def _hotspot_sort_key(row: Dict[str, Any]) -> Tuple[float, int, int, str]:
    """排序键：occupancy 降序、n_pairs 降序、residue_id 升序、residue_name 升序。"""
    return (-row["occupancy"], -row["n_pairs"], row["residue_id"],
            row["residue_name"])


def hotspot_segments(rows: Sequence[Dict[str, Any]],
                     exclude_unused: bool = True) -> List[Dict[str, Any]]:
    """把热点行按 h5 molecule 分组为图表分段（残基分布柱状图的数据）。

    exclude_unused=True 时丢弃 occupancy <= 0 的行；group_type 含 metal 但
    types 里没有 metal_coordination 的离子行始终丢弃。
    段内按 residue_id 升序；段间按段内最小 residue_id 升序（并列按 molecule 升序）。
    """
    kept = [row for row in rows if _segment_keep(row, exclude_unused)]
    segments = [_segment(molecule, group)
                for molecule, group in _molecule_groups(kept)]
    return sorted(segments, key=_segment_sort_key)


def _molecule_groups(rows: Sequence[Dict[str, Any]]
                     ) -> List[Tuple[str, List[Dict[str, Any]]]]:
    """按 molecule 聚合残基行（保持首次出现顺序）。"""
    grouped: Dict[str, List[Dict[str, Any]]] = {}
    for row in rows:
        grouped.setdefault(row["molecule"], []).append(row)
    return list(grouped.items())


def _segment_keep(row: Dict[str, Any], exclude_unused: bool) -> bool:
    """该残基行是否进入分段（零占用与未参与金属配位的离子被丢弃）。"""
    if exclude_unused and row["occupancy"] <= 0:
        return False
    return not _unused_ion(row)


def _unused_ion(row: Dict[str, Any]) -> bool:
    """金属基团但未参与 metal_coordination 的残基行。"""
    if METAL_GROUP_TYPE not in row["group_types"]:
        return False
    return METAL_COORDINATION_TYPE not in {entry["type"] for entry in row["types"]}


def _segment(molecule: str, rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    """一个 molecule 段：段内按 residue_id 升序，附残基数与最大占有率。"""
    residues = sorted(rows, key=lambda row: row["residue_id"])
    return {"molecule": molecule, "n_residues": len(residues),
            "max_occupancy": max(row["occupancy"] for row in residues),
            "residues": residues}


def _segment_sort_key(segment: Dict[str, Any]) -> Tuple[int, str]:
    """段排序键：段内最小 residue_id 升序、molecule 升序。"""
    return segment["residues"][0]["residue_id"], segment["molecule"]
