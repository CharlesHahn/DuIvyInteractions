# -*- coding: utf-8 -*-
"""visualizers - Web 可视化纯函数层（数据组装 + 图表 JSON，零 web 依赖）。"""

from .data import (
    parse_datasource, overview, export_table, aggregate_table, pair_drilldown,
    pair_summary, load_recent, add_recent, get_exporter,
    runs_from_mask, numeric_metrics, matrix_view, cell_values, frame_view,
    interaction_meta, metric_stats, string_stats, pair_frames, frame_activity,
    type_summary,
)
from .figures import fig_pair

__all__ = [
    "parse_datasource", "overview", "export_table", "aggregate_table",
    "pair_drilldown", "pair_summary", "load_recent", "add_recent", "get_exporter",
    "runs_from_mask", "numeric_metrics", "matrix_view", "cell_values", "frame_view",
    "interaction_meta", "metric_stats", "string_stats", "pair_frames",
    "frame_activity", "type_summary", "fig_pair",
]
