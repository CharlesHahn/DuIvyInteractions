# -*- coding: utf-8 -*-
"""Web 传输层：FastAPI 路由 + Jinja 模板（薄胶水，逻辑全在 visualizers/）。

create_app() 由 `dii web` 子命令调用；数据源状态存服务端内存（单用户本地）。
"""

from pathlib import Path
from typing import Dict, List

from fastapi import FastAPI, Form, Request
from fastapi.responses import HTMLResponse, JSONResponse, RedirectResponse, Response
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

from ..core.datas import Interaction
from ..visualizers import (
    add_recent, fig_existence_timeline, fig_metric_timeseries,
    fig_occupancy_hist, get_exporter, load_recent, overview,
    pair_drilldown, parse_datasource, type_table,
)

_WEB_DIR = Path(__file__).parent
TEMPLATES = Jinja2Templates(directory=str(_WEB_DIR / "templates"))  # Jinja2 3.1 内置 tojson 过滤器
MAX_PORT_TRIES = 10


def create_app() -> FastAPI:
    """创建 FastAPI 应用（数据源状态 = 服务端内存，单用户）。"""
    state: Dict = {"datasource": None, "items": []}
    app = FastAPI(title="DII Web", docs_url=None, redoc_url=None)
    app.mount("/static", StaticFiles(directory=str(_WEB_DIR / "static")), name="static")

    @app.get("/", response_class=HTMLResponse)
    def index(request: Request):
        return _render_index(request, state)

    @app.post("/select")
    def select(path: str = Form(...)):
        try:
            state["items"] = parse_datasource(path)
        except ValueError as e:
            print(f"[WARN] viz: {e}")
            return RedirectResponse("/", status_code=303)
        state["datasource"] = path if state["items"] else None
        if state["datasource"]:
            add_recent(path)
        return RedirectResponse("/", status_code=303)

    @app.get("/type/{name}", response_class=HTMLResponse)
    def type_page(name: str, request: Request, search: str = "",
                  sort_by: str = "index", sort_dir: str = "asc",
                  page: int = 1):
        return _render_type(request, state, name, search, sort_by, sort_dir, page)

    @app.get("/pair/{name}/{idx}", response_class=HTMLResponse)
    def pair_page(name: str, idx: int, request: Request):
        return _render_pair(request, state, name, idx)

    @app.get("/api/figures/{name}/{idx}.json")
    def figures_json(name: str, idx: int):
        it = _get_interaction(state, name)
        return JSONResponse(_pair_figures(pair_drilldown(it, idx)))

    @app.get("/download/csv/{name}")
    def download_csv(name: str, search: str = "", sort_by: str = "index",
                     sort_dir: str = "asc", page: int = 1):
        it = _get_interaction(state, name)
        csv = _to_csv(type_table(it, search, sort_by, sort_dir, page))
        return Response(csv, media_type="text/csv",
                        headers={"Content-Disposition":
                                 f"attachment; filename={name}_summary.csv"})

    return app


def _render_index(request: Request, state: Dict) -> HTMLResponse:
    """总览页：数据源表单 + 最近使用 + 类型卡片（含 occupancy 直方图）。"""
    cards = overview(state["items"])
    figs = [fig_occupancy_hist(c["occupancies"], c["type"])
            for c in cards if c["n_pairs"] > 0]
    return TEMPLATES.TemplateResponse(
        request, "index.html", { "title": "总览",
            "datasource": state["datasource"], "cards": cards,
            "figures": figs, "recent": load_recent(),
        })


def _render_type(request: Request, state: Dict, name: str, search: str,
                 sort_by: str, sort_dir: str, page: int) -> HTMLResponse:
    """类型页：搜索/排序/分页的表格视图。"""
    it = _get_interaction(state, name)
    table = type_table(it, search, sort_by, sort_dir, page)
    return TEMPLATES.TemplateResponse(request, "type.html", { "title": name, "type": name,
        "columns": table["columns"], "rows": table["rows"],
        "total": table["total"], "pages": table["pages"],
        "page": table["page"], "search": search,
        "sort_by": sort_by, "sort_dir": sort_dir,
    })


def _render_pair(request: Request, state: Dict, name: str, idx: int) -> HTMLResponse:
    """下钻页：元信息 + existence 时间线 + 各指标曲线。"""
    it = _get_interaction(state, name)
    idx = min(max(0, idx), it.n_pairs - 1)
    data = pair_drilldown(it, idx)
    return TEMPLATES.TemplateResponse(request, "pair.html", { "title": data["pair_label"],
        "type": name, "idx": idx, "n_pairs": it.n_pairs,
        "label": data["pair_label"], "meta": data,
        "metric_names": list(data["series"].keys()),
        "figures": _pair_figures(data),
    })


def _pair_figures(data: Dict) -> List[Dict]:
    """下钻页的 Plotly JSON 序列（existence + 各指标），带模板挂载的 div id。"""
    label = data["pair_label"]
    figs = [{"div": "figExistence",
             **fig_existence_timeline(data["times"], data["existence"], label)}]
    for name, s in data["series"].items():
        fig = fig_metric_timeseries(s["times"], s["values"], name, label)
        figs.append({"div": f"fig{len(figs)}", **fig})
    return figs


def _get_interaction(state: Dict, name: str) -> Interaction:
    """按类型名取当前数据源的 Interaction；未知类型报 404 语义的 ValueError。"""
    for tname, _, it in state["items"]:
        if tname == name:
            return it
    raise ValueError(f"interaction type not loaded: '{name}'")


def _to_csv(table: Dict) -> str:
    """表格数据 → CSV 字符串（含表头，行序与视图一致）。"""
    cols = table["columns"]
    lines = [",".join(cols)]
    for row in table["rows"]:
        lines.append(",".join(str(row[c]) for c in cols))
    return "\n".join(lines) + "\n"
