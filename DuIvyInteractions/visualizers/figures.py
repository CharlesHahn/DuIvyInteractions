# -*- coding: utf-8 -*-
"""Web 可视化图表层：Plotly JSON 组装（纯 dict，不依赖 plotly.py）。

前端（plotly.js-basic）直接 Plotly.newPlot(el, figure_json)。
"""

from typing import Dict, List

import numpy as np

PLOT_BG = "rgba(0,0,0,0)"
FONT = {"size": 12}
MARGIN = {"l": 56, "r": 16, "t": 40, "b": 40}


def _layout(title: str, x_title: str, y_title: str,
            barmode: str = None) -> Dict:
    """统一布局（W1 全部图共用）。"""
    layout = {
        "title": {"text": title, "font": FONT},
        "xaxis": {"title": {"text": x_title}, "gridcolor": "#eee"},
        "yaxis": {"title": {"text": y_title}, "gridcolor": "#eee"},
        "margin": MARGIN,
        "paper_bgcolor": PLOT_BG,
        "plot_bgcolor": PLOT_BG,
        "font": FONT,
        "showlegend": True,
    }
    if barmode:
        layout["barmode"] = barmode
    return layout


def fig_existence_timeline(times: List[float], existence: List[int],
                           pair_label: str) -> Dict:
    """单对 existence 时间线（阶梯线，0/1）。"""
    trace = {
        "type": "scatter", "mode": "lines",
        "x": times, "y": existence,
        "line": {"shape": "hv", "width": 1.5},
        "name": "existence",
        "showlegend": False,
    }
    return {"data": [trace],
            "layout": _layout(f"Existence — {pair_label}", "Time (ps)", "Existence")}


def fig_metric_timeseries(times: List[float], values: List[float],
                          metric_name: str, pair_label: str) -> Dict:
    """单对单指标逐帧曲线（存在帧取值，缺失帧为 null 断点）。"""
    masked = _masked(times, values)
    trace = {
        "type": "scatter", "mode": "lines+markers",
        "x": masked["x"], "y": masked["y"],
        "line": {"width": 1.5}, "marker": {"size": 4},
        "name": metric_name,
        "showlegend": False,
    }
    return {"data": [trace],
            "layout": _layout(f"{metric_name} — {pair_label}",
                              "Time (ps)", metric_name)}


def fig_occupancy_hist(occupancies: List[float], type_name: str) -> Dict:
    """类型总览：全对 occupancy 直方图。"""
    trace = {
        "type": "histogram",
        "x": list(occupancies),
        "xbins": {"start": 0.0, "end": 1.0, "size": 0.05},
        "marker": {"color": "#4a7ebb"},
        "showlegend": False,
    }
    return {"data": [trace],
            "layout": _layout(f"Occupancy distribution — {type_name}",
                              "Occupancy", "Pairs")}


def _masked(times: List[float], values: List[float]) -> Dict[str, List]:
    """NaN 值置 None（plotly 断点语义），返回 {"x":..., "y":...}。"""
    ys = [None if (v is None or np.isnan(v)) else v for v in values]
    return {"x": list(times), "y": ys}
