# -*- coding: utf-8 -*-
"""Web 图表层：组装 Plotly JSON（前端用 plotly.js 渲染）。"""

from typing import Dict, List

import numpy as np

from .data import METRIC_UNITS, PT_TYPES

PLOT_BG = "rgba(0,0,0,0)"
FONT = {"size": 12}
EXISTENCE_ROW = 0.14       # existence 子图纵向占比
TYPE_ROW = 0.1             # 类型带子图纵向占比
SUBPLOT_GAP = 0.05         # 子图间距（纵向占比）
# NPG 调色板（ggsci scale_color_npg，Nature 系列常用）
CHART_COLORS = ("#4dbbd5", "#e64b35", "#00a087", "#3c5488",
                "#f39b7f", "#8491b4", "#91d1c2", "#dc0000", "#7e6148")
HEIGHT_BASE = 130          # 组合图基础高度（px）
HEIGHT_PER_METRIC = 120    # 每个指标子图追加高度（px）
PAIR_MARGIN = {"l": 72, "r": 20, "t": 24, "b": 44}
AXIS_FONT = {"size": 11}
EXISTENCE_RANGE = [0.0, 1.3]
TYPE_RANGE = [0, 1.4]
# NPG 调色板（ggsci scale_color_npg，Nature 系列常用）
PT_COLORS: Dict[str, str] = {"P": "#3c5488", "T": "#e64b35", "N": "#8491b4"}


def fig_pair(data: Dict) -> Dict:
    """下钻组合图：existence 带 + 可选类型带 + 每指标一行，共享时间轴。"""
    metrics = list(data["series"])
    typed = "types" in data
    domains = _domains(len(metrics), typed)
    traces = [_existence_trace(data)]
    if typed:
        traces += _type_traces(data["times"], data["types"], data["existence"])
    axis_offset = 2 if typed else 1
    traces += [_metric_trace(data["series"][name], name, i, axis_offset)
               for i, name in enumerate(metrics)]
    return {"data": traces,
            "layout": _pair_layout(metrics, domains, typed)}


def _domains(n_metrics: int, typed: bool) -> List[List[float]]:
    """自上而下分配纵向区间：existence、可选类型带、每指标一行（含间隔）。"""
    fixed = EXISTENCE_ROW + (TYPE_ROW if typed else 0.0)
    gaps = n_metrics + (1 if typed else 0)
    metric_h = (1 - fixed - SUBPLOT_GAP * gaps) / max(n_metrics, 1)
    heights = [EXISTENCE_ROW] + ([TYPE_ROW] if typed else [])
    heights += [metric_h] * n_metrics
    domains, top = [], 1.0
    for height in heights:
        domains.append([round(top - height, 4), round(top, 4)])
        top -= height + SUBPLOT_GAP
    return domains


def _pair_layout(metrics: List[str], domains: List[List[float]],
                 typed: bool) -> Dict:
    """组合图布局：顶部 existence、可选类型带、下方每指标一行、各轴共享 x。"""
    layout = {
        "height": HEIGHT_BASE + HEIGHT_PER_METRIC * max(len(metrics), 1),
        "margin": PAIR_MARGIN,
        "paper_bgcolor": PLOT_BG, "plot_bgcolor": PLOT_BG, "font": FONT,
        "showlegend": typed, "hovermode": "x unified", "hoversubplots": "axis",
        "xaxis": {"domain": [0, 1], "anchor": "y",
                  "showgrid": False, "zeroline": False,
                  "showticklabels": not metrics},
        "yaxis": {"domain": domains[0], "anchor": "x", "range": EXISTENCE_RANGE,
                  "tickvals": [0, 1],
                  "title": {"text": "exists", "font": AXIS_FONT},
                  "gridcolor": "#eef0f4", "zeroline": False},
    }
    if typed:
        layout.update(_type_axis(domains[1]))
    layout.update(_metric_axes(metrics, domains, 2 if typed else 1))
    return layout


def _type_axis(domain: List[float]) -> Dict:
    """类型带子图的 x/y 轴（单行 0/1 带，隐藏 x 刻度）。"""
    return {"xaxis2": {"domain": [0, 1], "anchor": "y2", "matches": "x",
                       "showgrid": False, "zeroline": False,
                       "showticklabels": False},
            "yaxis2": {"domain": domain, "anchor": "x2", "range": TYPE_RANGE,
                       "tickvals": [1],
                       "title": {"text": "type", "font": AXIS_FONT},
                       "gridcolor": "#eef0f4", "zeroline": False}}


def _metric_axes(metrics: List[str], domains: List[List[float]],
                 offset: int) -> Dict:
    """各指标子图的 x/y 轴定义（共享顶部 x 轴）。"""
    bottom = len(metrics) - 1
    axes = {}
    for i, name in enumerate(metrics):
        axis = i + 1 + offset
        axes[f"yaxis{axis}"] = _metric_axis(name, domains[i + offset], axis)
        axes[f"xaxis{axis}"] = _metric_xaxis(axis, i == bottom)
    return axes


def _metric_axis(name: str, domain: List[float], axis: int) -> Dict:
    """指标子图的 y 轴（锚定同名 x 轴，标题带单位）。"""
    unit = METRIC_UNITS.get(name)
    title = f"{name} ({unit})" if unit else name
    return {"domain": domain, "anchor": f"x{axis}",
            "title": {"text": title, "font": AXIS_FONT},
            "gridcolor": "#eef0f4", "zeroline": False}


def _metric_xaxis(axis: int, show_ticks: bool) -> Dict:
    """指标子图的 x 轴：匹配顶部轴缩放，仅最底部显示刻度。"""
    xaxis = {"domain": [0, 1], "anchor": f"y{axis}", "matches": "x",
             "showgrid": False, "zeroline": False,
             "showticklabels": show_ticks}
    if show_ticks:
        xaxis["title"] = {"text": "Time (ps)"}
    return xaxis


def _existence_trace(data: Dict) -> Dict:
    """existence 阶梯带（0/1 填充，顶部一行）。"""
    return {"type": "scatter", "mode": "lines",
            "x": data["times"], "y": data["existence"],
            "line": {"shape": "hv", "width": 1, "color": "#3c5488"},
            "fill": "tozeroy", "fillcolor": "rgba(60,84,136,0.16)",
            "name": "existence", "showlegend": False,
            "hovertemplate": "existence %{y}<extra></extra>"}


def _type_traces(times: List[float], types: Dict,
                 existence: List[int]) -> List[Dict]:
    """类型带：每个出现类型一组方形点（仅活跃帧，y 恒为 1）。"""
    values = np.asarray(types["values"], dtype=str)
    active = np.asarray(existence, dtype=bool) & (values != "")
    xs = np.asarray(times, dtype=float)
    return [_type_trace(xs[active & (values == name)].tolist(), name)
            for name in PT_TYPES if bool(np.any(active & (values == name)))]


def _type_trace(xs: List[float], name: str) -> Dict:
    """单个类型的帧点带（颜色取自 PT_COLORS，图例可见）。"""
    return {"type": "scatter", "mode": "markers",
            "x": xs, "y": [1] * len(xs),
            "marker": {"symbol": "square", "size": 6,
                       "color": PT_COLORS[name]},
            "name": f"type {name}", "showlegend": True,
            "xaxis": "x2", "yaxis": "y2",
            "hovertemplate": f"type {name} %{{x}}<extra></extra>"}


def _metric_trace(series: Dict, name: str, i: int, offset: int) -> Dict:
    """单个指标曲线（NaN → null 断点），第 i 行子图。"""
    masked = _masked(series["times"], series["values"])
    axis = i + 1 + offset
    return {"type": "scatter", "mode": "lines",
            "x": masked["x"], "y": masked["y"],
            "line": {"width": 1.4, "color": CHART_COLORS[i % len(CHART_COLORS)]},
            "name": name, "showlegend": False,
            "xaxis": f"x{axis}", "yaxis": f"y{axis}",
            "hovertemplate": name + " %{y:.3f}<extra></extra>"}


def _masked(times: List[float], values: List[float]) -> Dict[str, List]:
    """NaN 值置 None（plotly 断点语义），返回 {"x":..., "y":...}。"""
    ys = [None if (v is None or np.isnan(v)) else v for v in values]
    return {"x": list(times), "y": ys}
