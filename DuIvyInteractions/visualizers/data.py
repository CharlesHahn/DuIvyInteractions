# -*- coding: utf-8 -*-
"""Web 可视化数据层：数据源解析、总览、类型表格、下钻、降采样、最近数据源。

纯函数层（spec FR-009）：只依赖 core/io，零 web 依赖，全部可单元测试。
"""

import json
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from ..core.datas import Interaction
from ..io import EXPORTER_CLASSES, load_interactions
from ..io.interaction_exporter import InteractionExporter

PAGE_SIZE = 50          # 类型表格分页大小
MAX_SERIES_POINTS = 2000  # 下钻曲线降采样上限
RECENT_LIMIT = 10       # 最近数据源记忆条数
RECENT_PATH = Path(os.environ.get("DII_CACHE_HOME", Path.home() / ".cache" / "duivyinteractions")) / "web_recent.json"


def parse_datasource(path: str) -> List[Tuple[str, str, Interaction]]:
    """解析数据源（单个 h5 文件或输出目录）→ [(类型名, h5路径, Interaction)]。

    目录扫描其中全部 *.h5；损坏/空数据文件跳过并告警（与 dii export 防护一致）。
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
        interaction = _load_single(h5)
        if interaction is not None:
            items.append((interaction.interaction_type, str(h5), interaction))
    return items


def _load_single(h5: Path) -> Optional[Interaction]:
    """加载单个 h5 的首个 Interaction；损坏/空数据返回 None 并告警。"""
    try:
        loaded = load_interactions(str(h5))
    except (OSError, ValueError) as e:
        print(f"[WARN] viz: skip corrupt h5 {h5.name}: {e}")
        return None
    if not loaded:
        print(f"[WARN] viz: skip empty h5 {h5.name}")
        return None
    return loaded[0]


def overview(items: List[Tuple[str, str, Interaction]]) -> List[Dict]:
    """每类型卡片：类型名、对数、帧数、flags 名单、occupancy 数组。"""
    cards = []
    for name, h5_path, it in items:
        cards.append({
            "type": name,
            "h5_path": h5_path,
            "n_pairs": it.n_pairs,
            "n_frames": it.n_frames,
            "flags": sorted(it.flags.keys()),
            "occupancies": it.occupancy().tolist(),
        })
    return cards


def get_exporter(interaction: Interaction) -> InteractionExporter:
    """按 interaction_type 取导出器（pair 标签等复用 dii export 的单一来源）。"""
    cls = EXPORTER_CLASSES.get(interaction.interaction_type)
    if cls is None:
        raise ValueError(f"unknown interaction type: '{interaction.interaction_type}'")
    return cls()


def type_table(interaction: Interaction, search: str = "",
               sort_by: str = "index", sort_dir: str = "asc",
               page: int = 1) -> Dict:
    """类型页表格：行 = 全部 pair，支持搜索 / 排序 / 分页（spec FR-004）。

    行字段：pair_label、occupancy、avg_<metric>/std_<metric>、flags 派生列。
    """
    exporter = get_exporter(interaction)
    rows = _build_rows(interaction, exporter)
    rows = _filter_rows(rows, search)
    rows = _sort_rows(rows, sort_by, sort_dir)
    total = len(rows)
    pages = max(1, (total + PAGE_SIZE - 1) // PAGE_SIZE)
    page = min(max(1, page), pages)
    return {
        "columns": _columns(interaction),
        "rows": rows[(page - 1) * PAGE_SIZE: page * PAGE_SIZE],
        "total": total, "pages": pages, "page": page,
    }


def _build_rows(interaction: Interaction, exporter: InteractionExporter) -> List[Dict]:
    """全 pair 的行数据（含 flags 派生列 = 活跃帧中 True 的比例）。"""
    occ = interaction.occupancy()
    rows = []
    for i in range(interaction.n_pairs):
        row: Dict = {
            "index": i,
            "pair_label": exporter.get_pair_label(interaction, i),
            "occupancy": float(occ[i]),
        }
        for name, arr in interaction.metrics.items():
            active = arr[i][interaction.existence[i]]
            row[f"avg_{name}"] = float(np.nanmean(active)) if active.size else ""
            row[f"std_{name}"] = float(np.nanstd(active)) if active.size else ""
        for fname, farr in interaction.flags.items():
            active = farr[i][interaction.existence[i]]
            row[f"flag_{fname}"] = float(np.mean(active)) if active.size else ""
        rows.append(row)
    return rows


def _filter_rows(rows: List[Dict], search: str) -> List[Dict]:
    """按 pair_label 大小写不敏感搜索。"""
    if not search:
        return rows
    needle = search.lower()
    return [r for r in rows if needle in r["pair_label"].lower()]


def _sort_rows(rows: List[Dict], sort_by: str, sort_dir: str) -> List[Dict]:
    """按列名排序（index/occupancy/avg_* 等），方向 asc/desc；未知列回退 index。"""
    if sort_by not in (rows[0] if rows else {}):
        sort_by = "index"
    return sorted(rows, key=lambda r: r[sort_by],
                  reverse=(sort_dir == "desc"))


def _columns(interaction: Interaction) -> List[str]:
    """表格列名顺序：index/pair_label/occupancy + 指标 + flags。"""
    cols = ["index", "pair_label", "occupancy"]
    cols += [f"avg_{n}" for n in interaction.metrics]
    cols += [f"std_{n}" for n in interaction.metrics]
    cols += [f"flag_{n}" for n in interaction.flags]
    return cols


def pair_drilldown(interaction: Interaction, idx: int) -> Dict:
    """单对下钻数据：元信息 + times + existence(0/1) + 各指标序列（可降采样）。"""
    exporter = get_exporter(interaction)
    times = np.asarray(interaction.times if interaction.n_frames else [],
                       dtype=float).tolist()
    ex = interaction.existence[idx].astype(int).tolist()
    series: Dict[str, List] = {}
    for name, arr in interaction.metrics.items():
        if np.asarray(arr).dtype.kind in ("U", "S", "O"):
            continue  # 字符串指标（如 pistacking_type）不进曲线
        values = np.asarray(arr[idx], dtype=float).tolist()
        xs, ys = _downsample(times, values, MAX_SERIES_POINTS)
        series[name] = {"times": xs, "values": ys}
    return {
        "type": interaction.interaction_type,
        "idx": idx,
        "pair_label": exporter.get_pair_label(interaction, idx),
        "groups": _group_meta(interaction, idx),
        "times": times,
        "existence": ex,
        "series": series,
        "flags": {n: interaction.flags[n][idx].astype(int).tolist()
                  for n in interaction.flags},
    }


def _group_meta(interaction: Interaction, idx: int) -> List[Dict]:
    """基团元信息（不含原子坐标，仅 residue/原子名列表）。"""
    metas = []
    for g in interaction.groups[idx]:
        metas.append({
            "group_id": g.group_id, "group_type": g.group_type,
            "residue_name": g.residue_name, "residue_id": g.residue_id,
            "molecule": g.molecule,
            "atom_names": [a.atom_name for a in g.atoms],
        })
    return metas


def _downsample(xs: List, ys: List, max_points: int) -> Tuple[List, List]:
    """超过 max_points 时均匀抽稀，保留首末两个点。"""
    n = len(ys)
    if n <= max_points:
        return xs, ys
    keep = np.unique(np.linspace(0, n - 1, max_points).astype(int))
    return [xs[k] for k in keep], [ys[k] for k in keep]


def load_recent() -> List[str]:
    """最近使用的数据源列表（文件损坏时静默重建）。"""
    try:
        return json.loads(RECENT_PATH.read_text())[:RECENT_LIMIT]
    except Exception:
        return []


def add_recent(path: str) -> List[str]:
    """记录一次数据源使用（去重置顶，最多 RECENT_LIMIT 条）。"""
    entries = [p for p in load_recent() if p != path]
    entries.insert(0, path)
    entries = entries[:RECENT_LIMIT]
    try:
        RECENT_PATH.parent.mkdir(parents=True, exist_ok=True)
        RECENT_PATH.write_text(json.dumps(entries))
    except OSError as e:
        print(f"[WARN] viz: cannot persist recent datasource: {e}")
    return entries
