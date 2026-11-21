# -*- coding: utf-8 -*-
"""Web 图表层：组装 Plotly JSON（前端用 plotly.js 渲染）。"""

from typing import Dict, List

import numpy as np

from .data import METRIC_UNITS

PLOT_BG = "rgba(0,0,0,0)"
FONT = {"size": 12}
EXISTENCE_ROW = 0.14       # existence 子图纵向占比
SUBPLOT_GAP = 0.05         # 子图间距（纵向占比）
CHART_COLORS = ("#4f46e5", "#0891b2", "#d97706", "#059669", "#db2777")
HEIGHT_BASE = 130          # 组合图基础高度（px）
HEIGHT_PER_METRIC = 120    # 每个指标子图追加高度（px）
PAIR_MARGIN = {"l": 72, "r": 20, "t": 24, "b": 44}
AXIS_FONT = {"size": 11}
EXISTENCE_RANGE = [0.0, 1.3]


def fig_pair(data: Dict) -> Dict:
    """下钻组合图：existence 带 + 每指标一行，共享时间轴并统一悬浮。"""
    metrics = list(data["series"])
    domains = _domains(len(metrics))
    traces = [_existence_trace(data)]
    traces += [_metric_trace(data["series"][name], name, i)
               for i, name in enumerate(metrics)]
    return {"data": traces,
            "layout": _pair_layout(metrics, domains)}


def _domains(n_metrics: int) -> List[List[float]]:
    """自上而下分配纵向区间：existence 一行 + 每指标一行（含间隔）。"""
    metric_h = (1 - EXISTENCE_ROW - SUBPLOT_GAP * n_metrics) / max(n_metrics, 1)
    heights = [EXISTENCE_ROW] + [metric_h] * n_metrics
    domains, top = [], 1.0
    for height in heights:
        domains.append([round(top - height, 4), round(top, 4)])
        top -= height + SUBPLOT_GAP
    return domains


def _pair_layout(metrics: List[str],
                 domains: List[List[float]]) -> Dict:
    """组合图布局：顶部 existence、下方每指标一行、各轴共享 x。"""
    layout = {
        "height": HEIGHT_BASE + HEIGHT_PER_METRIC * max(len(metrics), 1),
        "margin": PAIR_MARGIN,
        "paper_bgcolor": PLOT_BG, "plot_bgcolor": PLOT_BG, "font": FONT,
        "showlegend": False, "hovermode": "x unified", "hoversubplots": "axis",
        "xaxis": {"domain": [0, 1], "anchor": "y",
                  "showgrid": False, "zeroline": False,
                  "showticklabels": not metrics},
        "yaxis": {"domain": domains[0], "anchor": "x", "range": EXISTENCE_RANGE,
                  "tickvals": [0, 1],
                  "title": {"text": "exists", "font": AXIS_FONT},
                  "gridcolor": "#eef0f4", "zeroline": False},
    }
    layout.update(_metric_axes(metrics, domains))
    return layout


def _metric_axes(metrics: List[str], domains: List[List[float]]) -> Dict:
    """各指标子图的 x/y 轴定义（共享顶部 x 轴）。"""
    bottom = len(metrics) - 1
    axes = {}
    for i, name in enumerate(metrics):
        axes[f"yaxis{i + 2}"] = _metric_axis(name, domains[i + 1], i)
        axes[f"xaxis{i + 2}"] = _metric_xaxis(i, bottom)
    return axes


def _metric_axis(name: str, domain: List[float], i: int) -> Dict:
    """指标子图的 y 轴（锚定同名 x 轴，标题带单位）。"""
    unit = METRIC_UNITS.get(name)
    title = f"{name} ({unit})" if unit else name
    return {"domain": domain, "anchor": f"x{i + 2}",
            "title": {"text": title, "font": AXIS_FONT},
            "gridcolor": "#eef0f4", "zeroline": False}


def _metric_xaxis(i: int, bottom: int) -> Dict:
    """指标子图的 x 轴：匹配顶部轴缩放，仅最底部显示刻度。"""
    axis = {"domain": [0, 1], "anchor": f"y{i + 2}", "matches": "x",
            "showgrid": False, "zeroline": False,
            "showticklabels": i == bottom}
    if i == bottom:
        axis["title"] = {"text": "Time (ps)"}
    return axis


def _existence_trace(data: Dict) -> Dict:
    """existence 阶梯带（0/1 填充，顶部一行）。"""
    return {"type": "scatter", "mode": "lines",
            "x": data["times"], "y": data["existence"],
            "line": {"shape": "hv", "width": 1, "color": "#4f46e5"},
            "fill": "tozeroy", "fillcolor": "rgba(79,70,229,0.16)",
            "name": "existence", "showlegend": False,
            "hovertemplate": "existence %{y}<extra></extra>"}


def _metric_trace(series: Dict, name: str, i: int) -> Dict:
    """单个指标曲线（NaN → null 断点），第 i 行子图。"""
    masked = _masked(series["times"], series["values"])
    return {"type": "scatter", "mode": "lines",
            "x": masked["x"], "y": masked["y"],
            "line": {"width": 1.4, "color": CHART_COLORS[i % len(CHART_COLORS)]},
            "name": name, "showlegend": False,
            "xaxis": f"x{i + 2}", "yaxis": f"y{i + 2}",
            "hovertemplate": name + " %{y:.3f}<extra></extra>"}


def _masked(times: List[float], values: List[float]) -> Dict[str, List]:
    """NaN 值置 None（plotly 断点语义），返回 {"x":..., "y":...}。"""
    ys = [None if (v is None or np.isnan(v)) else v for v in values]
    return {"x": list(times), "y": ys}
