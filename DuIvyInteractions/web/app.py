# -*- coding: utf-8 -*-
"""Web 传输层：FastAPI 路由 + Jinja 模板。"""

import csv
import io
import math
from pathlib import Path
from typing import Any, Callable, Dict, Tuple, Union
from urllib.parse import parse_qs, quote

from fastapi import FastAPI, HTTPException, Request
from fastapi.exceptions import RequestValidationError
from fastapi.responses import HTMLResponse, JSONResponse, RedirectResponse, Response
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from starlette.concurrency import run_in_threadpool
from starlette.exceptions import HTTPException as StarletteHTTPException
from starlette.middleware.gzip import GZipMiddleware

from ..core.datas import Interaction
from ..visualizers import (
    add_recent, aggregate_table, cell_values, export_table, fig_pair,
    frame_activity, frame_view, interaction_meta, load_recent, matrix_view,
    metric_stats, overview, pair_drilldown, pair_frames, pair_summary,
    parse_datasource, string_stats, type_summary,
)
from .picker import pick_directory

_WEB_DIR = Path(__file__).parent
TEMPLATES = Jinja2Templates(directory=str(_WEB_DIR / "templates"))
FRAME_LIMIT_MAX = 2000  # 帧检查器单次最多返回的活跃对数


def create_app() -> FastAPI:
    """创建 FastAPI 应用（数据源状态 = 服务端内存，单用户）。"""
    state: Dict[str, Any] = {"source": None}
    app = FastAPI(title="DII Web", docs_url=None, redoc_url=None,
                  openapi_url=None)
    app.add_middleware(GZipMiddleware, minimum_size=1024)
    app.mount("/static", StaticFiles(directory=str(_WEB_DIR / "static")), name="static")
    _register_pages(app, state)
    _register_api(app, state)
    return app


def _register_pages(app: FastAPI, state: Dict[str, Any]) -> None:
    """页面路由 + HTML 错误页。"""
    _register_home(app, state)
    _register_detail(app, state)
    _register_error_page(app, state)


def _register_home(app: FastAPI, state: Dict[str, Any]) -> None:
    """总览与数据源选择。"""

    @app.get("/", response_class=HTMLResponse)
    def index(request: Request, error: str = "", warn: str = "") -> Response:
        return _render_index(request, state, error, warn)

    @app.post("/select")
    async def select(request: Request) -> Response:
        body = await request.body()
        return await run_in_threadpool(_select, state,
                                       _form_value(body, "path"))


def _register_detail(app: FastAPI, state: Dict[str, Any]) -> None:
    """类型页与单对下钻页。"""

    @app.get("/type/{name}", response_class=HTMLResponse)
    def type_page(request: Request, name: str, search: str = "",
                  min_occupancy: str = "0", sort_by: str = "occupancy",
                  sort_dir: str = "desc", page: int = 1) -> Response:
        sort_by, sort_dir = _parse_sort(sort_by, sort_dir)
        return _render_type(request, state, name, search,
                            _parse_occupancy(min_occupancy), sort_by,
                            sort_dir, page)

    @app.get("/pair/{name}/{idx}", response_class=HTMLResponse)
    def pair_page(request: Request, name: str, idx: int) -> Response:
        return _render_pair(request, state, name, idx)


def _register_error_page(app: FastAPI, state: Dict[str, Any]) -> None:
    """未知类型 / 非法参数 / 未捕获异常 → API 返回 JSON，页面返回 HTML 错误页。"""
    _register_http_error(app, state)
    _register_validation_error(app, state)
    _register_unexpected_error(app, state)


def _register_http_error(app: FastAPI, state: Dict[str, Any]) -> None:
    """HTTP 异常处理（404 等）。"""

    @app.exception_handler(StarletteHTTPException)
    async def http_error(request: Request,
                         exc: StarletteHTTPException) -> Response:
        if request.url.path.startswith("/api/"):
            return JSONResponse({"detail": exc.detail}, status_code=exc.status_code)
        if request.url.path.startswith("/download/"):
            return Response(str(exc.detail), media_type="text/plain",
                            status_code=exc.status_code)
        return _error_page(request, state, exc.status_code, str(exc.detail))


def _register_validation_error(app: FastAPI, state: Dict[str, Any]) -> None:
    """查询参数校验失败处理（422）。"""

    @app.exception_handler(RequestValidationError)
    async def invalid_query(request: Request,
                            exc: RequestValidationError) -> Response:
        message = "invalid parameters: " + ", ".join(
            str(item.get("loc")) for item in exc.errors())
        if request.url.path.startswith("/api/"):
            return JSONResponse({"detail": message}, status_code=422)
        return _error_page(request, state, 422, message)


def _register_unexpected_error(app: FastAPI, state: Dict[str, Any]) -> None:
    """未捕获异常处理（500）。"""

    @app.exception_handler(Exception)
    async def unexpected(request: Request, exc: Exception) -> Response:
        print(f"[ERROR] web: {exc!r}")
        if request.url.path.startswith("/api/"):
            return JSONResponse({"detail": "internal error"}, status_code=500)
        return _error_page(request, state, 500, "internal error")


def _error_page(request: Request, state: Dict[str, Any], status: int,
                message: str) -> HTMLResponse:
    """渲染错误页（保留侧栏与数据源上下文）。"""
    return TEMPLATES.TemplateResponse(
        request, "error.html",
        {"title": str(status), "status": status, "message": message,
         **_nav(state)},
        status_code=status)


def _register_api(app: FastAPI, state: Dict[str, Any]) -> None:
    """JSON 数据接口。"""
    _register_matrix_api(app, state)
    _register_cell_api(app, state)
    _register_pair_api(app, state)
    _register_csv_api(app, state)
    _register_pick_api(app)


def _register_matrix_api(app: FastAPI, state: Dict[str, Any]) -> None:
    """矩阵接口。"""

    @app.get("/api/matrix/{name}", response_model=None)
    def api_matrix(name: str, search: str = "", min_occupancy: float = 0.0,
                   sort_by: str = "occupancy",
                   sort_dir: str = "desc") -> Dict[str, Any]:
        sort_by, sort_dir = _parse_sort(sort_by, sort_dir)
        return _guard(lambda: matrix_view(
            _get_interaction(state, name), search,
            _clamp_occupancy(min_occupancy), sort_by, sort_dir))


def _register_cell_api(app: FastAPI, state: Dict[str, Any]) -> None:
    """单元格与帧检查。"""

    @app.get("/api/cell/{name}/{pair}/{frame}", response_model=None)
    def api_cell(name: str, pair: int, frame: int) -> Dict[str, Any]:
        return _guard(lambda: cell_values(_get_interaction(state, name), pair, frame))

    @app.get("/api/frame/{name}/{frame}", response_model=None)
    def api_frame(name: str, frame: int, limit: int = 200) -> Dict[str, Any]:
        limit = min(max(limit, 1), FRAME_LIMIT_MAX)
        data = _guard(lambda: frame_view(_get_interaction(state, name),
                                         frame, limit))
        data["limit_max"] = FRAME_LIMIT_MAX
        return data


def _register_pair_api(app: FastAPI, state: Dict[str, Any]) -> None:
    """逐帧明细。"""

    @app.get("/api/frames/{name}/{idx}", response_model=None)
    def api_frames(name: str, idx: int, start: int = 0) -> Dict[str, Any]:
        return _guard(lambda: pair_frames(_get_interaction(state, name), idx, start))


def _register_csv_api(app: FastAPI, state: Dict[str, Any]) -> None:
    """CSV 导出。"""

    @app.get("/download/csv/{name}", response_model=None)
    def download_csv(name: str, search: str = "", sort_by: str = "index",
                     sort_dir: str = "asc",
                     min_occupancy: str = "0") -> Response:
        sort_by, sort_dir = _parse_sort(sort_by, sort_dir)
        it = _get_interaction(state, name)
        table = _guard(lambda: export_table(
            it, search, sort_by, sort_dir, _parse_occupancy(min_occupancy)))
        filename = quote(f"{name}_summary.csv")
        return Response(_to_csv(table), media_type="text/csv",
                        headers={"Content-Disposition":
                                 f"attachment; filename={filename}; "
                                 f"filename*=UTF-8''{filename}"})


def _register_pick_api(app: FastAPI) -> None:
    """原生目录选择框。"""

    @app.post("/api/pick", response_model=None)
    def pick() -> Union[Dict[str, Any], JSONResponse]:
        try:
            return {"path": pick_directory()}
        except (RuntimeError, ValueError, OSError) as e:
            return JSONResponse({"path": "", "error": str(e)}, status_code=400)


SORT_FIELDS = ("occupancy", "label", "index")


def _parse_sort(sort_by: str, sort_dir: str) -> Tuple[str, str]:
    """排序字段白名单 + 方向归一（非法值回退 occupancy desc）。"""
    field = sort_by if sort_by in SORT_FIELDS else "occupancy"
    return field, ("asc" if sort_dir == "asc" else "desc")


def _clamp_occupancy(value: float) -> float:
    """把占有率下限限制在 [0, 1]；非有限值回落 0。"""
    if not math.isfinite(value):
        return 0.0
    return min(max(value, 0.0), 1.0)


def _parse_occupancy(raw: str) -> float:
    """解析占有率下限：空/非法值回落 0，并限制在 [0, 1]。"""
    try:
        return _clamp_occupancy(float(raw))
    except (TypeError, ValueError):
        return 0.0


def _guard(call: Callable[[], Dict[str, Any]]) -> Dict[str, Any]:
    """把越界/未知类型等 ValueError 转为 404（JSON 接口统一出口）。"""
    try:
        return call()
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))


def _form_value(body: bytes, key: str) -> str:
    """解析 URL-encoded 表单体的单值字段（缺失/解码失败返回空串）。"""
    try:
        values = parse_qs(body.decode("utf-8"))
    except UnicodeDecodeError:
        return ""
    return values.get(key, [""])[0]


def _select(state: Dict[str, Any], path: str) -> RedirectResponse:
    """解析数据源并跳回总览；失败或跳过文件时带回提示（不静默）。"""
    if not path.strip():
        return RedirectResponse(f"/?error={quote('no path provided')}",
                                status_code=303)
    expected = _source_files(path)
    try:
        items = parse_datasource(path)
    except ValueError as e:
        return RedirectResponse(f"/?error={quote(str(e))}", status_code=303)
    if not items:
        return RedirectResponse(f"/?error={quote(f'no usable .h5 found: {path}')}",
                                status_code=303)
    state["source"] = (path, items)
    add_recent(path)
    loaded = len({h5 for _, h5, _ in items})
    return _loaded_redirect(expected, loaded)


def _loaded_redirect(expected: int, loaded: int) -> RedirectResponse:
    """载入成功后的跳转；有文件被跳过时带上告警。"""
    skipped = max(0, expected - loaded)
    if skipped:
        message = quote(f"{skipped} file(s) skipped (unreadable or empty)")
        return RedirectResponse(f"/?warn={message}", status_code=303)
    return RedirectResponse("/", status_code=303)


def _source_files(path: str) -> int:
    """数据源中 h5 文件数（单文件 1、目录计数、不存在 0）。"""
    source = Path(path).expanduser()
    if source.is_file():
        return 1
    if source.is_dir():
        return len(list(source.glob("*.h5")))
    return 0


def _source(state: Dict[str, Any]) -> Any:
    """当前数据源 (路径, items)；未载入时为 ("", [])。"""
    source = state.get("source")
    return source if source else ("", [])


def _nav_source(path: str, items: Any, active_type: str = "") -> Dict[str, Any]:
    """导航上下文（侧栏类型列表 + 当前数据源）。"""
    return {"types": type_summary(items),
            "datasource": path or None,
            "ds_name": Path(path).name if path else "",
            "active_type": active_type}


def _nav(state: Dict[str, Any], active_type: str = "") -> Dict[str, Any]:
    """按当前状态生成导航上下文。"""
    path, items = _source(state)
    return _nav_source(path, items, active_type)


def _find_interaction(items: Any, name: str) -> Tuple[str, Interaction]:
    """在 items 中按显示名查找 (h5 路径, Interaction)；未知 → 404。"""
    for tname, h5_path, it in items:
        if tname == name:
            return h5_path, it
    raise HTTPException(status_code=404, detail=f"interaction type not loaded: '{name}'")


def _render_index(request: Request, state: Dict[str, Any], error: str,
                  warn: str) -> HTMLResponse:
    """总览页：唯一数据源入口 + 类型卡片。"""
    path, items = _source(state)
    return TEMPLATES.TemplateResponse(request, "index.html", {
        "title": "overview", "cards": overview(items),
        "recent": load_recent(), "error": error, "warning": warn,
        **_nav_source(path, items),
    })


def _render_type(request: Request, state: Dict[str, Any], name: str, search: str,
                 min_occupancy: float, sort_by: str, sort_dir: str,
                 page: int) -> HTMLResponse:
    """类型页：元数据 + 分布 + pairs × frames 矩阵 + 帧检查器 + 聚合表。"""
    path, items = _source(state)
    h5_path, it = _find_interaction(items, name)
    clamped = _clamp_occupancy(min_occupancy)
    table = aggregate_table(it, search, sort_by, sort_dir, page, clamped)
    return TEMPLATES.TemplateResponse(request, "type.html", {
        "title": name, "type": name, "search": search,
        "min_occupancy": clamped,
        "sort_by": sort_by, "sort_dir": sort_dir,
        "n_pairs": it.n_pairs, "n_frames": it.n_frames,
        "meta": interaction_meta(it, h5_path),
        "stats": metric_stats(it), "string_stats": string_stats(it),
        "activity": frame_activity(it), "table": table,
        **_nav_source(path, items, name),
    })


def _render_pair(request: Request, state: Dict[str, Any], name: str,
                 idx: int) -> HTMLResponse:
    """下钻页：基团元信息 + existence/指标对齐图 + 逐帧明细。"""
    path, items = _source(state)
    _, it = _find_interaction(items, name)
    if not 0 <= idx < it.n_pairs:
        raise HTTPException(status_code=404, detail=f"pair index out of range: {idx}")
    data = pair_drilldown(it, idx)
    return TEMPLATES.TemplateResponse(request, "pair.html", {
        "title": data["pair_label"], "type": name, "idx": idx,
        "n_pairs": it.n_pairs, "n_frames": it.n_frames,
        "pair_label": data["pair_label"], "meta": data["groups"],
        "summary": pair_summary(it, idx), "downsampled": data["downsampled"],
        "flags": sorted(it.flags), "figure": fig_pair(data),
        **_nav_source(path, items, name),
    })


def _get_interaction(state: Dict[str, Any], name: str) -> Interaction:
    """按类型名取当前数据源的 Interaction；未知类型 → 404。"""
    return _find_interaction(_source(state)[1], name)[1]


def _to_csv(table: Dict[str, Any]) -> str:
    """表格数据 → CSV 字符串（含表头，行序与视图一致）。"""
    buffer = io.StringIO()
    writer = csv.writer(buffer)
    writer.writerow(table["columns"])
    for row in table["rows"]:
        writer.writerow([row[column] for column in table["columns"]])
    return buffer.getvalue()
