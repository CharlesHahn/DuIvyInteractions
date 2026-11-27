# -*- coding: utf-8 -*-
"""Web 传输层：FastAPI 路由 + Jinja 模板。"""

import csv
import io
import math
import os
import secrets
import threading
import time
from contextlib import asynccontextmanager
from pathlib import Path
from typing import (
    Any, AsyncIterator, Callable, Dict, List, Optional, Sequence, Tuple)
from urllib.parse import parse_qs, quote

from fastapi import FastAPI, HTTPException, Request
from fastapi.exceptions import RequestValidationError
from fastapi.responses import (
    FileResponse, HTMLResponse, JSONResponse, RedirectResponse, Response)
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from starlette.concurrency import run_in_threadpool
from starlette.exceptions import HTTPException as StarletteHTTPException
from starlette.middleware.gzip import GZipMiddleware

from ..core.datas import Interaction
from ..io.h5 import RunMeta, RunRecord, load_run_meta
from ..visualizers import (
    add_recent, aggregate_table, cell_values, export_table, fig_pair,
    frame_activity, frame_view, get_exporter, interaction_meta, load_recent,
    matrix_view, metric_stats, overview, pair_drilldown, pair_frames,
    pair_summary, pair_types, parse_datasource, string_stats, type_summary,
)
# 聚合表首屏与 /api/rows 复用同一行组装/单位函数，避免两套格式化逻辑漂移。
from ..visualizers.data import (
    TABLE_PAGE_SIZE, _columns, _table_rows, _units, merge_batches)
from ..visualizers.structure import (
    alignment_errors, interaction_expressions, load_universe, topology_pdb,
)
from .fsbrowse import (
    BYTES_PER_KB, H5_SUFFIX, LIST_LIMIT, SIZE_DECIMALS, confine,
    list_directory,
)
from .jobs import (
    ANY_PAIRING, DEFAULT_FF, DEFAULT_STRATEGY, FORCE_FIELDS, INTERACTIONS,
    PAIRINGS, STATUS_CANCELLED, STRATEGIES, Job, JobManager, parse_run_form,
)
from .picker import pick_directory, pick_file

_WEB_DIR = Path(__file__).parent
TEMPLATES = Jinja2Templates(directory=str(_WEB_DIR / "templates"))
FRAME_LIMIT_MAX = 2000  # 帧检查器单次最多返回的活跃对数
ROWS_LIMIT_MIN = 1      # /api/rows 单次最少返回行数
ROWS_LIMIT_MAX = 200    # /api/rows 单次最多返回行数
ROWS_OFFSET_MIN = 0     # /api/rows 起始行号下限
ROWS_DEFAULT_LIMIT = TABLE_PAGE_SIZE  # 首屏与聚合表分页大小一致
CELL_DASH = "—"         # 聚合表空值占位（与页面渲染一致）
TOKEN_COOKIE = "dii_token"
SERVER_BROWSE_HINT = ("native picker unavailable on this host; "
                      "use the server folder browser")
AUTH_PAGE = (
    "<!doctype html><html><head><meta charset=\"utf-8\">"
    "<title>token required</title></head><body>"
    "<h1>token required</h1>"
    "<p>append ?token=... to the URL to open this page.</p>"
    "</body></html>"
)
AUTH_FAIL_LIMIT = 20  # 单个来源在窗口内允许的 token 失败次数
AUTH_FAIL_WINDOW = 60.0  # 失败计数的滑动窗口（秒）
TOO_MANY_FAILURES = "too many failed attempts"
FORBIDDEN_HOST = "forbidden host"
UNSUPPORTED_PATH = "unsupported path (expect a directory or .h5 file)"
STRUCTURE_NOT_SELECTED = "structure not selected"
STRUCTURE_PAIR_REQUIRED = "type and idx required"
# 3D 页自动载入：h5 的 run/ 组记录不可用时的 reason（前端据此提示并回退手选）
AUTO_NO_RUN_META = "no-run-meta"      # h5 没有 run/ 组
AUTO_FILE_MISSING = "missing"         # 记录的文件不存在
AUTO_FILE_CHANGED = "changed"         # 记录的文件大小与指纹不符
AUTO_DETAIL_NO_RUN_META = "no run metadata in this h5"
AUTO_DETAIL_FILE_MISSING = "recorded file is missing"
AUTO_DETAIL_FILE_CHANGED = "recorded file changed"
AUTO_SOURCE_RUN = "run"               # 成功响应里的来源标记
STRUCTURE_SUFFIXES = frozenset({
    ".tpr", ".xtc", ".trr", ".gro", ".pdb", ".dcd", ".nc", ".itp", ".psf",
    ".prmtop", ".parm7", ".mol2"})
PICK_DIR = "dir"
PICK_FILE = "file"
DATASET_KIND = "dataset"
STRUCTURE_KIND = "structure"
RUN_LOG_LINES = 300  # 任务页最多展示的日志尾部行数
ELAPSED_DECIMALS = 1  # 页面/JSON 中秒数保留的小数位
JOB_NOT_FOUND = "job not found"
CANNOT_CANCEL = "job is not running"
PAIR_VIEW_CURVES = "curves"
PAIR_VIEW_3D = "3d"
PAIR_VIEWS = (PAIR_VIEW_CURVES, PAIR_VIEW_3D)
# base.html / index.html 存在内联脚本且不在本次授权范围内，故允许 inline
# script；若将来把内联脚本外置，可收紧为 script-src 'self'。
SECURITY_HEADERS = {
    "X-Content-Type-Options": "nosniff",
    "X-Frame-Options": "DENY",
    "Referrer-Policy": "no-referrer",
    "Content-Security-Policy": (
        "default-src 'self'; script-src 'self' 'unsafe-inline'; "
        "style-src 'self' 'unsafe-inline'; img-src 'self' data:; "
        "connect-src 'self'"),
}
# Mol* viewer bundle 在加载期使用 new Function；仅 3D 页需要 'unsafe-eval'。
# Mol* 内置 WASM 以 base64 data: URL 经 fetch 载入 → connect-src 需放行 data:。
STRUCTURE_CSP = (
    "default-src 'self'; script-src 'self' 'unsafe-inline' 'unsafe-eval'; "
    "style-src 'self' 'unsafe-inline'; img-src 'self' data:; "
    "connect-src 'self' data:")
# 静态资源必须每次使用前向服务器校验：/static/* 由 StaticFiles 提供 ETag 与
# Last-Modified，未变更时回 304；若缺此头，浏览器启发式缓存会长期沿用旧 JS。
STATIC_PREFIX = "/static/"
STATIC_CACHE_CONTROL = "no-cache"
# 拓扑/轨迹按当前会话的结构动态生成，响应不得进入浏览器或代理缓存：
# 复用上一次会话的轨迹文件会让 Mol* 报 element count mismatch。
SESSION_CACHE_CONTROL = "no-store"
CACHE_CONTROL_HEADER = "Cache-Control"
_AUTH_FAILURES: Dict[str, List[float]] = {}
_AUTH_LOCK = threading.Lock()


def create_app(token: str = "", browse_mode: str = "picker", root: str = "",
               allowed_hosts: Optional[List[str]] = None) -> FastAPI:
    """创建 FastAPI 应用（数据源状态 = 服务端内存，单用户）。"""
    manager = JobManager(root=root)
    state: Dict[str, Any] = {"sources": [], "batch_items": {}, "items": [],
                             "structure": None, "structure_version": 0,
                             "jobs": manager}
    app = FastAPI(title="DII Web", docs_url=None, redoc_url=None,
                  openapi_url=None, lifespan=_job_lifespan(manager))
    app.state.jobs = manager
    app.add_middleware(GZipMiddleware, minimum_size=1024)
    app.mount("/static", StaticFiles(directory=str(_WEB_DIR / "static")), name="static")
    if token:
        _install_token_guard(app, token)
    if allowed_hosts:
        _install_host_guard(app, allowed_hosts)
    _install_security_headers(app)
    _register_pages(app, state, browse_mode, root)
    _register_api(app, state, browse_mode, root)
    return app


def _job_lifespan(manager: JobManager) -> Callable[[FastAPI], Any]:
    """应用生命周期：退出时终止仍在运行的计算任务。"""

    @asynccontextmanager
    async def lifespan(_app: FastAPI) -> AsyncIterator[None]:
        yield
        manager.shutdown()

    return lifespan


def _install_token_guard(app: FastAPI, token: str) -> None:
    """安装 access token 校验中间件（query 或 cookie 通过，query 通过后写 cookie）。"""

    @app.middleware("http")
    async def guard(request: Request, call_next: Callable) -> Response:
        client = _client_host(request)
        if _auth_blocked(client):
            return Response(TOO_MANY_FAILURES, media_type="text/plain",
                            status_code=429)
        query_token = request.query_params.get("token")
        cookie_token = request.cookies.get(TOKEN_COOKIE)
        if not (_token_ok(token, query_token) or _token_ok(token, cookie_token)):
            _record_auth_failure(client)
            return _unauthorized(request)
        _clear_auth_failures(client)
        response = await call_next(request)
        if _token_ok(token, query_token):
            response.set_cookie(TOKEN_COOKIE, token, httponly=True,
                                samesite="lax", path="/")
        return response


def _install_host_guard(app: FastAPI, allowed_hosts: List[str]) -> None:
    """安装 Host 头校验中间件（防 DNS rebinding）。"""
    allowed = frozenset(host.lower() for host in allowed_hosts)

    @app.middleware("http")
    async def host_guard(request: Request, call_next: Callable) -> Response:
        host = _host_name(request.headers.get("host", ""))
        if host in allowed:
            return await call_next(request)
        return Response(FORBIDDEN_HOST, media_type="text/plain", status_code=403)


def _install_security_headers(app: FastAPI) -> None:
    """安装安全响应头中间件（不覆盖响应已有的同名头）；/static/ 追加缓存校验头。"""

    @app.middleware("http")
    async def security_headers(request: Request, call_next: Callable) -> Response:
        response = await call_next(request)
        for name, value in SECURITY_HEADERS.items():
            response.headers.setdefault(name, value)
        if request.url.path.startswith(STATIC_PREFIX):
            response.headers.setdefault(CACHE_CONTROL_HEADER,
                                        STATIC_CACHE_CONTROL)
        return response


def _host_name(value: str) -> str:
    """提取 Host 头中的主机名（去掉端口与 IPv6 方括号）。"""
    host = (value or "").strip().lower()
    if host.startswith("["):
        end = host.find("]")
        return host[1:end] if end > 0 else host
    if host.count(":") == 1:
        name, _, port = host.partition(":")
        return name if port.isdigit() else host
    return host


def _client_host(request: Request) -> str:
    """客户端来源地址（缺失时为空串）。"""
    return request.client.host if request.client else ""


def _auth_blocked(client: str) -> bool:
    """来源是否已达到失败次数上限（同时清理过期记录）。"""
    with _AUTH_LOCK:
        failures = _recent_failures(client)
        if failures:
            _AUTH_FAILURES[client] = failures
        else:
            _AUTH_FAILURES.pop(client, None)
        return len(failures) >= AUTH_FAIL_LIMIT


def _record_auth_failure(client: str) -> None:
    """记录一次 token 失败（窗口内最多保留 AUTH_FAIL_LIMIT 条）。"""
    with _AUTH_LOCK:
        failures = _recent_failures(client)
        failures.append(time.monotonic())
        _AUTH_FAILURES[client] = failures[-AUTH_FAIL_LIMIT:]


def _recent_failures(client: str) -> List[float]:
    """窗口内的失败时间戳（调用方需持有 _AUTH_LOCK）。"""
    cutoff = time.monotonic() - AUTH_FAIL_WINDOW
    return [stamp for stamp in _AUTH_FAILURES.get(client, ()) if stamp > cutoff]


def _clear_auth_failures(client: str) -> None:
    """清空来源的失败计数。"""
    with _AUTH_LOCK:
        _AUTH_FAILURES.pop(client, None)


def _token_ok(token: str, supplied: Optional[str]) -> bool:
    """判断请求携带的 token 是否与配置一致。"""
    if not (token and supplied):
        return False
    return secrets.compare_digest(token.encode("utf-8"), supplied.encode("utf-8"))


def _unauthorized(request: Request) -> Response:
    """未携带有效 token 时的响应选择（API JSON / 下载纯文本 / 页面 HTML）。"""
    if request.url.path.startswith("/api/"):
        return JSONResponse({"detail": "token required"}, status_code=401)
    if request.url.path.startswith("/download/"):
        return Response("token required", media_type="text/plain",
                        status_code=401)
    return HTMLResponse(AUTH_PAGE, status_code=401)


def _register_pages(app: FastAPI, state: Dict[str, Any], browse_mode: str,
                    root: str) -> None:
    """页面路由 + HTML 错误页。"""
    _register_home(app, state, browse_mode, root)
    _register_detail(app, state, browse_mode)
    _register_structure_page(app, state)
    _register_run_page(app, state, browse_mode, root)
    _register_error_page(app, state)


def _register_home(app: FastAPI, state: Dict[str, Any], browse_mode: str,
                   root: str) -> None:
    """总览与数据源选择。"""

    @app.get("/", response_class=HTMLResponse)
    def index(request: Request, error: str = "", warn: str = "") -> Response:
        return _render_index(request, state, error, warn, browse_mode)

    @app.post("/select")
    async def select(request: Request) -> Response:
        body = await request.body()
        return await run_in_threadpool(_select, state,
                                       _form_value(body, "path"), root)

    @app.post("/select/remove")
    async def select_remove(request: Request) -> Response:
        body = await request.body()
        return await run_in_threadpool(_remove_batch, state,
                                       _form_value(body, "path"))

    @app.post("/select/clear")
    async def select_clear() -> Response:
        return await run_in_threadpool(_clear_batches, state)


def _register_detail(app: FastAPI, state: Dict[str, Any],
                     browse_mode: str) -> None:
    """类型页与单对下钻页（含 3D 视图开关）。"""

    @app.get("/type/{name}", response_class=HTMLResponse)
    def type_page(request: Request, name: str, search: str = "",
                  min_occupancy: str = "0", sort_by: str = "occupancy",
                  sort_dir: str = "desc", page: int = 1,
                  frame: int = 0) -> Response:
        sort_by, sort_dir = _parse_sort(sort_by, sort_dir)
        return _render_type(request, state, name, search,
                            _parse_occupancy(min_occupancy), sort_by,
                            sort_dir, page, frame)

    @app.get("/pair/{name}/{idx}", response_class=HTMLResponse)
    def pair_page(request: Request, name: str, idx: int,
                  view: str = PAIR_VIEW_CURVES,
                  frame: int = 0) -> Response:
        return _render_pair(request, state, name, idx, view, frame,
                            browse_mode)


def _register_structure_page(app: FastAPI, state: Dict[str, Any]) -> None:
    """旧 3D 独立页 → 单对页 3D 视图（302，保留调用方帧号）。"""

    @app.get("/structure/{name}/{idx}", response_class=HTMLResponse)
    def structure_page(request: Request, name: str, idx: int,
                       frame: int = 0) -> Response:
        requested = frame if "frame" in request.query_params else None
        return _structure_redirect(name, idx, requested)


def _register_run_page(app: FastAPI, state: Dict[str, Any], browse_mode: str,
                       root: str) -> None:
    """计算发起页、任务详情页与表单提交。"""

    @app.get("/run", response_class=HTMLResponse)
    def run_page(request: Request, error: str = "") -> Response:
        return _render_run(request, state, browse_mode, root, error)

    @app.get("/run/{job_id}", response_class=HTMLResponse)
    def job_page(request: Request, job_id: str) -> Response:
        return _render_job(request, state, job_id, browse_mode, root)

    @app.post("/run")
    async def start_run(request: Request) -> Response:
        body = await request.body()
        return await run_in_threadpool(_start_run, state["jobs"], root,
                                       _form_fields(body))


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


def _register_api(app: FastAPI, state: Dict[str, Any],
                  browse_mode: str, root: str) -> None:
    """JSON 数据接口。"""
    _register_matrix_api(app, state)
    _register_rows_api(app, state)
    _register_cell_api(app, state)
    _register_pair_api(app, state)
    _register_csv_api(app, state)
    _register_list_api(app, browse_mode, root)
    _register_pick_api(app, browse_mode)
    _register_structure_select(app, state, root)
    _register_structure_auto(app, state, root)
    _register_structure_export(app, state)
    _register_structure_interactions(app, state)
    _register_run_api(app, state)


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


def _register_rows_api(app: FastAPI, state: Dict[str, Any]) -> None:
    """聚合表行窗口接口（前端 load more 按需取数，不整页刷新）。"""

    @app.get("/api/rows/{name}", response_model=None)
    def api_rows(name: str, search: str = "", min_occupancy: str = "0",
                 sort_by: str = "occupancy", sort_dir: str = "desc",
                 offset: str = "", limit: str = "") -> Dict[str, Any]:
        sort_by, sort_dir = _parse_sort(sort_by, sort_dir)
        interaction = _get_interaction(state, name)
        return _guard(lambda: _rows_view(
            interaction, search, _parse_occupancy(min_occupancy), sort_by,
            sort_dir, _parse_int(offset, ROWS_OFFSET_MIN),
            _parse_int(limit, ROWS_DEFAULT_LIMIT)))


def _rows_view(interaction: Interaction, search: str, min_occupancy: float,
               sort_by: str, sort_dir: str, offset: int,
               limit: int) -> Dict[str, Any]:
    """聚合表行窗口：与 /type 聚合表共用行组装与单元格格式化。"""
    rows = _table_rows(interaction, search, sort_by, sort_dir, min_occupancy)
    offset, limit = _clamp_rows_window(offset, limit)
    columns = _columns(interaction)
    return {"offset": offset, "limit": limit, "total": len(rows),
            "columns": columns, "units": _units(interaction),
            "rows": [_rows_entry(row, columns)
                     for row in rows[offset:offset + limit]]}


def _clamp_rows_window(offset: int, limit: int) -> Tuple[int, int]:
    """行窗口钳制：offset 不为负，limit 落在 [ROWS_LIMIT_MIN, ROWS_LIMIT_MAX]。"""
    return (max(offset, ROWS_OFFSET_MIN),
            min(max(limit, ROWS_LIMIT_MIN), ROWS_LIMIT_MAX))


def _parse_int(raw: str, default: int) -> int:
    """整数查询参数：空/非法值回落 default（随后由行窗口钳制）。"""
    try:
        return int(raw)
    except (TypeError, ValueError):
        return default


def _rows_entry(row: Dict[str, Any], columns: List[str]) -> Dict[str, Any]:
    """单个 pair 的行窗口条目（values 与 columns 逐列对齐）。"""
    return {"idx": row["index"], "label": row["pair_label"],
            "values": [_cell_text(column, row[column]) for column in columns]}


def _attach_display_cells(table: Dict[str, Any]) -> Dict[str, Any]:
    """给聚合表首屏行补显示文本（与 /api/rows 使用同一格式化）。"""
    table["cells"] = [[_cell_text(column, row[column])
                       for column in table["columns"]]
                      for row in table["rows"]]
    return table


def _cell_text(column: str, value: Any) -> str:
    """聚合表单元格文本：占有率/flags 百分比，avg/std 保留三位小数。"""
    if column == "occupancy" or column.startswith("flag_"):
        return CELL_DASH if value == "" else f"{float(value) * 100:.1f}%"
    if column.startswith("avg_") or column.startswith("std_"):
        return CELL_DASH if value == "" else f"{float(value):.3f}"
    return "" if value is None else str(value)


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


def _register_list_api(app: FastAPI, browse_mode: str, root: str) -> None:
    """服务端目录浏览（仅 server 模式启用；kind 切换数据源/结构文件过滤）。"""

    @app.get("/api/list", response_model=None)
    def api_list(path: str = "", kind: str = DATASET_KIND) -> dict:
        if browse_mode != "server":
            raise HTTPException(status_code=404,
                                detail="server browse disabled")
        if kind not in (DATASET_KIND, STRUCTURE_KIND):
            raise HTTPException(status_code=400,
                                detail=f"unsupported kind: {kind}")
        return _guard(lambda: _list_entries(path, kind, root))


def _list_entries(path: str, kind: str, root: str) -> dict:
    """目录列表：structure 只列结构后缀，dataset 只列 *.h5。"""
    data = list_directory(path, root=root)
    if kind != STRUCTURE_KIND:
        return data
    try:
        files = _structure_files(Path(data["path"]), data["root"])
    except OSError:
        raise ValueError(f"cannot read directory: {path}") from None
    kept = files[:max(LIST_LIMIT - len(data["dirs"]), 0)]
    return {**data, "files": kept,
            "truncated": len(data["dirs"]) + len(files) > LIST_LIMIT}


def _structure_files(directory: Path, root: str) -> List[Dict[str, object]]:
    """目录下的结构文件条目（后缀白名单，跳过越界符号链接，按名排序）。"""
    base = Path(root) if root else None
    files: List[Dict[str, object]] = []
    for entry in sorted(directory.iterdir(), key=lambda item: item.name):
        if entry.suffix.lower() not in STRUCTURE_SUFFIXES or not entry.is_file():
            continue
        if base is not None and not _inside_root(entry, base):
            continue
        size = _size_kb(entry)
        if size is not None:
            files.append({"name": entry.name, "size_kb": size})
    return files


def _inside_root(entry: Path, base: Path) -> bool:
    """符号链接条目解析后是否仍在浏览边界内。"""
    try:
        resolved = entry.resolve()
    except OSError:
        return False
    return resolved == base or resolved.is_relative_to(base)


def _size_kb(entry: Path) -> Optional[float]:
    """文件大小（KB，保留一位小数）；读取失败返回 None。"""
    try:
        return round(entry.stat().st_size / BYTES_PER_KB, SIZE_DECIMALS)
    except OSError:
        return None


def _register_pick_api(app: FastAPI, browse_mode: str) -> None:
    """原生目录/文件选择框（仅 picker 模式可用）。"""

    @app.post("/api/pick", response_model=None)
    async def pick(request: Request) -> Response:
        if browse_mode != "picker":
            return JSONResponse({"path": "", "error": SERVER_BROWSE_HINT},
                                status_code=400)
        body = await request.body()
        mode = _param(request, body, "mode", PICK_DIR)
        kind = _param(request, body, "kind", STRUCTURE_KIND)
        try:
            path = await run_in_threadpool(_picked_path, mode, kind)
        except (RuntimeError, ValueError, OSError) as e:
            return JSONResponse({"path": "", "error": str(e)}, status_code=400)
        return {"path": path}


def _picked_path(mode: str, kind: str) -> str:
    """按 mode 打开原生对话框（file → 文件；dir → 目录；其它报错）。"""
    if mode == PICK_FILE:
        return pick_file(kind)
    if mode == PICK_DIR:
        return pick_directory()
    raise ValueError(f"unsupported mode: {mode}")


def _param(request: Request, body: bytes, key: str, default: str) -> str:
    """表单体优先、其次查询参数的字符串参数（缺省用 default）。"""
    return _form_value(body, key) or request.query_params.get(key, default)


def _register_structure_select(app: FastAPI, state: Dict[str, Any],
                               root: str) -> None:
    """选择结构并记录拓扑/轨迹路径。"""

    @app.post("/api/structure/select", response_model=None)
    async def select_structure(request: Request) -> Response:
        body = await request.body()
        return await run_in_threadpool(
            _select_structure, state, root,
            _form_value(body, "topology"), _form_value(body, "trajectory"))


def _register_structure_auto(app: FastAPI, state: Dict[str, Any],
                             root: str) -> None:
    """按 h5 的 run/ 记录自动载入结构（3D 页首次打开时前端尝试）。"""

    @app.post("/api/structure/auto", response_model=None)
    async def auto_structure(request: Request) -> Response:
        body = await request.body()
        return await run_in_threadpool(_auto_structure, state, root,
                                       _form_value(body, "type"))


def _register_structure_export(app: FastAPI, state: Dict[str, Any]) -> None:
    """拓扑 PDB（该对原子打 B-factor 标记）与原始轨迹文件。"""

    @app.get("/api/structure/topology", response_model=None)
    def api_structure_topology(type: str = "",
                               idx: Optional[int] = None) -> Response:
        return _topology_response(state, type, idx)

    @app.get("/api/structure/trajectory", response_model=None)
    def api_structure_trajectory() -> Response:
        structure = _structure_or_400(state)
        return FileResponse(structure["trajectory"],
                            media_type="application/octet-stream",
                            headers=_session_cache_headers())


def _register_structure_interactions(app: FastAPI,
                                     state: Dict[str, Any]) -> None:
    """相互作用元信息（pair_label/pt_types/选择表达式）。"""

    @app.get("/api/structure/interactions", response_model=None)
    def api_structure_interactions(type: str = "",
                                   idx: int = 0) -> Response:
        return _interactions_response(state, type, idx)


def _register_run_api(app: FastAPI, state: Dict[str, Any]) -> None:
    """任务轮询、取消与日志下载接口。"""

    @app.get("/api/run/{job_id}", response_model=None)
    def api_run(job_id: str, since: str = "0") -> Dict[str, Any]:
        job = _job_or_404(state["jobs"], job_id)
        return _job_payload(job, _parse_since(since))

    @app.post("/api/run/{job_id}/cancel", response_model=None)
    def api_cancel(job_id: str) -> Response:
        return _cancel_response(state["jobs"], job_id)

    @app.get("/download/run-log/{job_id}", response_model=None)
    def download_run_log(job_id: str) -> Response:
        return _run_log_response(_job_or_404(state["jobs"], job_id))


def _select_structure(state: Dict[str, Any], root: str, topology: str,
                      trajectory: str) -> Response:
    """打开并记录结构（会话版本号自增，供前端 URL 变体使用）。"""
    try:
        payload = _install_structure(state, _structure_path(root, topology),
                                     _structure_path(root, trajectory))
    except ValueError as error:
        return _structure_bad_request(error)
    return JSONResponse(payload)


def _install_structure(state: Dict[str, Any], top_path: Path,
                       traj_path: Path) -> Dict[str, Any]:
    """打开结构并写入会话状态（失败抛 ValueError，版本号自增）。"""
    universe = load_universe(str(top_path), str(traj_path))
    structure = {
        "topology": str(top_path), "trajectory": str(traj_path),
        "universe": universe, "n_atoms": int(universe.atoms.n_atoms),
        "n_frames": len(universe.trajectory),
        "has_box": universe.dimensions is not None,
        "lock": threading.Lock(),
    }
    state["structure_version"] = int(state.get("structure_version", 0)) + 1
    structure["version"] = state["structure_version"]
    state["structure"] = structure
    return {"ok": True, "n_atoms": structure["n_atoms"],
            "n_frames": structure["n_frames"],
            "has_box": structure["has_box"]}


def _structure_bad_request(error: ValueError) -> JSONResponse:
    """结构路径/打开失败 → 400（select 与 auto 共用）。"""
    return JSONResponse({"detail": str(error)}, status_code=400)


def _auto_structure(state: Dict[str, Any], root: str,
                    name: str) -> Response:
    """按 h5 的 run/ 记录自动载入结构（无记录/文件缺失或变更 → 409）。"""
    h5_path, _ = _find_interaction(_source(state)[1], name)
    record = load_run_meta(h5_path)
    if record is None:
        return _no_run_meta_conflict()
    try:
        top_path, traj_path = _recorded_paths(root, record)
    except ValueError as error:
        return _structure_bad_request(error)
    conflict = _fingerprint_conflict(record, top_path, traj_path)
    if conflict:
        return conflict
    try:
        payload = _install_structure(state, top_path, traj_path)
    except ValueError as error:
        return _structure_bad_request(error)
    return _auto_payload(payload, top_path, traj_path)


def _recorded_paths(root: str, record: RunRecord) -> Tuple[Path, Path]:
    """run/ 记录的拓扑/轨迹路径（同样受 root 限制）。"""
    return (_structure_path(root, record.meta.topology),
            _structure_path(root, record.meta.trajectory))


def _fingerprint_conflict(record: RunRecord, top_path: Path,
                          traj_path: Path) -> Optional[JSONResponse]:
    """记录文件缺失/大小不符 → 409 响应（一致为 None）。"""
    if not (top_path.is_file() and traj_path.is_file()):
        return _auto_conflict(AUTO_FILE_MISSING, AUTO_DETAIL_FILE_MISSING,
                              top_path, traj_path)
    if (_size_changed(top_path, record.topology_size)
            or _size_changed(traj_path, record.trajectory_size)):
        return _auto_conflict(AUTO_FILE_CHANGED, AUTO_DETAIL_FILE_CHANGED,
                              top_path, traj_path)
    return None


def _size_changed(path: Path, recorded: int) -> bool:
    """记录了大小时才比较（h5 未记录指纹时记为 -1，视为未记录）。"""
    if recorded <= 0:
        return False
    try:
        return path.stat().st_size != recorded
    except OSError:
        return True


def _auto_conflict(reason: str, detail: str, top_path: Path,
                   traj_path: Path) -> JSONResponse:
    """409：记录的结构不可用（返回路径供前端预填输入框）。"""
    return JSONResponse({"detail": detail, "reason": reason,
                         "topology": str(top_path),
                         "trajectory": str(traj_path)}, status_code=409)


def _no_run_meta_conflict() -> JSONResponse:
    """409：该 h5 没有 run/ 溯源组。"""
    return JSONResponse({"detail": AUTO_DETAIL_NO_RUN_META,
                         "reason": AUTO_NO_RUN_META}, status_code=409)


def _auto_payload(payload: Dict[str, Any], top_path: Path,
                  traj_path: Path) -> JSONResponse:
    """自动载入成功响应（附来源标记与记录路径）。"""
    return JSONResponse({**payload, "source": AUTO_SOURCE_RUN,
                         "topology": str(top_path),
                         "trajectory": str(traj_path)})


def _structure_path(root: str, path: str) -> Path:
    """拓扑/轨迹路径（非空，且在 root 之内）。"""
    if not path.strip():
        raise ValueError("structure path required")
    return confine(root, path)


def _structure_or_400(state: Dict[str, Any]) -> Dict[str, Any]:
    """已选择的结构状态（未选择 → 400）。"""
    structure = state.get("structure")
    if not structure:
        raise HTTPException(status_code=400, detail=STRUCTURE_NOT_SELECTED)
    return structure


def _session_cache_headers() -> Dict[str, str]:
    """会话相关响应的缓存头（拓扑/轨迹共用）。"""
    return {CACHE_CONTROL_HEADER: SESSION_CACHE_CONTROL}


def _locked_topology(structure: Dict[str, Any], interaction: Interaction,
                     idx: int) -> str:
    """加锁生成带该对 B-factor 标记的拓扑 PDB 文本。"""
    with structure["lock"]:
        return topology_pdb(structure["universe"], interaction, idx)


def _alignment_failure(structure: Dict[str, Any],
                       interaction: Interaction) -> Optional[JSONResponse]:
    """h5 与结构不一致时的 400 响应（一致为 None）。"""
    errors = alignment_errors(structure["universe"], interaction)
    if not errors:
        return None
    return JSONResponse({"detail": errors[0], "errors": errors},
                        status_code=400)


def _topology_response(state: Dict[str, Any], name: str,
                       idx: Optional[int]) -> Response:
    """该对打标后的拓扑 PDB（type/idx 缺失 → 400）。"""
    structure = _structure_or_400(state)
    if not name or idx is None:
        raise HTTPException(status_code=400, detail=STRUCTURE_PAIR_REQUIRED)
    interaction = _pair_interaction(state, name, idx)
    failure = _alignment_failure(structure, interaction)
    if failure:
        return failure
    return Response(_locked_topology(structure, interaction, idx),
                    media_type="chemical/x-pdb",
                    headers=_session_cache_headers())


def _interactions_response(state: Dict[str, Any], name: str,
                           idx: int) -> Response:
    """该对 payload（h5 与结构不一致 → 400 含 errors 列表）。"""
    structure = _structure_or_400(state)
    interaction = _pair_interaction(state, name, idx)
    failure = _alignment_failure(structure, interaction)
    if failure:
        return failure
    return JSONResponse(_interactions_payload(structure, interaction, idx))


def _interactions_payload(structure: Dict[str, Any], interaction: Interaction,
                          idx: int) -> Dict[str, Any]:
    """该对的标签/帧数/存在性/键型/选择表达式/原子数。"""
    return {
        "pair_label": get_exporter(interaction).get_pair_label(interaction,
                                                               idx),
        "n_frames": structure["n_frames"],
        "existence": [int(value) for value in interaction.existence[idx]],
        "pt_types": pair_types(interaction, idx),
        "residues": interaction_expressions(
            interaction, idx, structure["universe"])["residues"],
        "n_atoms": structure["n_atoms"],
    }


def _pair_interaction(state: Dict[str, Any], name: str,
                      idx: int) -> Interaction:
    """当前数据源中该类型第 idx 对（类型未知/对号越界 → 404）。"""
    _, interaction = _find_interaction(_source(state)[1], name)
    if not 0 <= idx < interaction.n_pairs:
        raise HTTPException(status_code=404,
                            detail=f"pair index out of range: {idx}")
    return interaction


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


def _pair_view(raw: str) -> str:
    """单对页视图白名单：curves / 3d；其它值回落 curves（页面不 422）。"""
    return raw if raw in PAIR_VIEWS else PAIR_VIEW_CURVES


def _structure_redirect(name: str, idx: int,
                        frame: Optional[int]) -> RedirectResponse:
    """旧 3D 独立页跳到单对页 3D 视图；调用方给了帧号则原样保留。"""
    url = f"/pair/{quote(name, safe='')}/{idx}?view={PAIR_VIEW_3D}"
    if frame is not None:
        url += f"&frame={frame}"
    return RedirectResponse(url, status_code=302)


def _clamp_frame(frame: int, n_frames: int) -> int:
    """把请求帧号钳制到 [0, n_frames-1]；无帧时恒为 0。"""
    if n_frames <= 0:
        return 0
    return min(max(frame, 0), n_frames - 1)


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


def _batch_key(path: str) -> str:
    """批次键：绝对解析路径（同一路径的不同写法视为同一批次）。"""
    return str(Path(path).expanduser().resolve())


def _add_batch(state: Dict[str, Any], path: str,
               items: List[Tuple[str, str, Interaction]]) -> None:
    """按路径追加/原位替换一个批次，并重建合并 items。"""
    key = _batch_key(path)
    if key not in state["batch_items"]:
        state["sources"].append({"path": key, "name": Path(key).name})
    state["batch_items"][key] = items
    state["items"] = _merge_sources(state)


def _merge_sources(state: Dict[str, Any]) -> List[Tuple[str, str, Interaction]]:
    """按批次顺序合并当前全部批次的 items（跨批次重名自动去重）。"""
    return merge_batches([state["batch_items"][source["path"]]
                          for source in state["sources"]])


def _remove_batch(state: Dict[str, Any], path: str) -> RedirectResponse:
    """移除一个批次并跳回总览；未载入的路径带回提示。"""
    key = _batch_key(path)
    if key not in state["batch_items"]:
        return _select_error(f"batch not loaded: {path}")
    del state["batch_items"][key]
    state["sources"] = [source for source in state["sources"]
                        if source["path"] != key]
    state["items"] = _merge_sources(state)
    return RedirectResponse("/", status_code=303)


def _clear_batches(state: Dict[str, Any]) -> RedirectResponse:
    """清空全部批次并跳回总览。"""
    state["sources"] = []
    state["batch_items"] = {}
    state["items"] = []
    return RedirectResponse("/", status_code=303)


def _select(state: Dict[str, Any], path: str, root: str = "") -> RedirectResponse:
    """追加一个数据源批次并跳回总览；同路径重选替换该批次。"""
    if not path.strip():
        return _select_error("no path provided")
    try:
        target = _select_target(path, root)
    except ValueError as e:
        return _select_error(str(e))
    try:
        items = parse_datasource(str(target))
    except ValueError as e:
        return _select_error(str(e))
    if not items:
        return _select_error(f"no usable .h5 found: {path}")
    _add_batch(state, str(target), items)
    add_recent(path)
    loaded = len({h5 for _, h5, _ in items})
    return _loaded_redirect(_source_files(str(target)), loaded)


def _select_target(path: str, root: str) -> Path:
    """校验数据源路径：必须在 root 内，且为目录或 *.h5 文件。"""
    target = confine(root, path)
    if not _selectable(target):
        raise ValueError(UNSUPPORTED_PATH)
    return target


def _selectable(target: Path) -> bool:
    """判断目标是否为可载入的数据源（目录或 *.h5 文件）。"""
    return target.is_dir() or (target.is_file()
                               and target.suffix.lower() == H5_SUFFIX)


def _select_error(message: str) -> RedirectResponse:
    """选择失败时跳回总览并带回错误提示。"""
    return RedirectResponse(f"/?error={quote(message)}", status_code=303)


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


def _source(state: Dict[str, Any]
            ) -> Tuple[str, List[Tuple[str, str, Interaction]]]:
    """当前数据源的组合标签与合并 items；未载入时为 ("", [])。"""
    sources = state.get("sources") or []
    if not sources:
        return "", []
    paths = [source["path"] for source in sources]
    return (paths[0] if len(paths) == 1 else " + ".join(paths)), state["items"]


def _batch_chunks(
        state: Dict[str, Any]
) -> List[Tuple[Dict[str, Any], List[Tuple[str, str, Interaction]]]]:
    """按批次切分合并 items（merge_batches 保持批次顺序与条数）。"""
    merged = state.get("items") or []
    chunks: List[Tuple[Dict[str, Any], List[Tuple[str, str, Interaction]]]] = []
    offset = 0
    for source in state["sources"]:
        size = len(state["batch_items"][source["path"]])
        chunks.append((source, merged[offset:offset + size]))
        offset += size
    return chunks


def _type_entry(name: str, interaction: Interaction) -> Dict[str, Any]:
    """侧栏类型条目：唯一链接名、基础类型名与对数。"""
    return {"name": name, "base": interaction.interaction_type,
            "n_pairs": interaction.n_pairs}


def _batch_summaries(state: Dict[str, Any]) -> List[Dict[str, Any]]:
    """导航用批次摘要：路径、名称、类型条数与类型明细。"""
    return [{"path": source["path"], "name": source["name"],
             "n_types": len(chunk),
             "types": [_type_entry(name, it) for name, _, it in chunk]}
            for source, chunk in _batch_chunks(state)]


def _batch_groups(state: Dict[str, Any]) -> List[Dict[str, Any]]:
    """总览分组：每个批次一张区块，卡片取该批次的合并 items。"""
    return [{"path": source["path"], "name": source["name"],
             "n_types": len(chunk), "cards": overview(chunk)}
            for source, chunk in _batch_chunks(state)]


def _ds_name(batches: Sequence[Dict[str, Any]]) -> str:
    """顶栏数据源名：单批次用批名，多批次为首批名 + 其余数量。"""
    if len(batches) == 1:
        return batches[0]["name"]
    return f"{batches[0]['name']} +{len(batches) - 1}"


def _nav_source(path: str, items: List[Tuple[str, str, Interaction]],
                active_type: str = "",
                batches: Sequence[Dict[str, Any]] = ()) -> Dict[str, Any]:
    """导航上下文（侧栏类型列表 + 当前数据源与批次）。"""
    return {"types": type_summary(items),
            "datasource": path or None,
            "ds_name": _ds_name(batches) if batches
            else (Path(path).name if path else ""),
            "batches": list(batches),
            "active_type": active_type}


def _nav(state: Dict[str, Any], active_type: str = "") -> Dict[str, Any]:
    """按当前状态生成导航上下文。"""
    path, items = _source(state)
    return _nav_source(path, items, active_type, _batch_summaries(state))


def _find_interaction(items: Any, name: str) -> Tuple[str, Interaction]:
    """在 items 中按显示名查找 (h5 路径, Interaction)；未知 → 404。"""
    for tname, h5_path, it in items:
        if tname == name:
            return h5_path, it
    raise HTTPException(status_code=404, detail=f"interaction type not loaded: '{name}'")


def _render_index(request: Request, state: Dict[str, Any], error: str,
                  warn: str, browse_mode: str) -> HTMLResponse:
    """总览页：按数据集分组的类型卡片；无数据时只展示两个入口。"""
    path, items = _source(state)
    return TEMPLATES.TemplateResponse(request, "index.html", {
        "title": "overview", "cards": overview(items),
        "groups": _batch_groups(state), "recent": load_recent(),
        "error": error, "warning": warn, "browse_mode": browse_mode,
        "add": bool(request.query_params.get("add")),
        "add_path": request.query_params.get("path", ""),
        **_nav_source(path, items, batches=_batch_summaries(state)),
    })


def _render_type(request: Request, state: Dict[str, Any], name: str, search: str,
                 min_occupancy: float, sort_by: str, sort_dir: str,
                 page: int, frame: int) -> HTMLResponse:
    """类型页：元数据 + 分布 + pairs × frames 矩阵 + 帧检查器 + 聚合表。"""
    path, items = _source(state)
    h5_path, it = _find_interaction(items, name)
    clamped = _clamp_occupancy(min_occupancy)
    table = _attach_display_cells(
        aggregate_table(it, search, sort_by, sort_dir, page, clamped))
    return TEMPLATES.TemplateResponse(request, "type.html", {
        "title": name, "type": name, "search": search,
        "min_occupancy": clamped,
        "sort_by": sort_by, "sort_dir": sort_dir,
        "n_pairs": it.n_pairs, "n_frames": it.n_frames,
        "frame": _clamp_frame(frame, it.n_frames),
        "meta": interaction_meta(it, h5_path),
        "stats": metric_stats(it), "string_stats": string_stats(it),
        "activity": frame_activity(it), "table": table,
        "rows_limit": ROWS_DEFAULT_LIMIT,
        "rows_offset": (table["page"] - 1) * TABLE_PAGE_SIZE,
        **_nav_source(path, items, name, _batch_summaries(state)),
    })


def _render_pair(request: Request, state: Dict[str, Any], name: str,
                 idx: int, view: str = PAIR_VIEW_CURVES, frame: int = 0,
                 browse_mode: str = "picker") -> HTMLResponse:
    """下钻页：基团元信息 + existence/指标对齐图 + 逐帧明细；

    view=3d 时同一页面改为 3D 区块（与 curves 视图共享帧号）。
    """
    path, items = _source(state)
    h5_path, it = _find_interaction(items, name)
    if not 0 <= idx < it.n_pairs:
        raise HTTPException(status_code=404, detail=f"pair index out of range: {idx}")
    view = _pair_view(view)
    data = pair_drilldown(it, idx)
    context: Dict[str, Any] = {
        "title": data["pair_label"], "type": name, "idx": idx,
        "n_pairs": it.n_pairs, "n_frames": it.n_frames,
        "pair_label": data["pair_label"], "meta": data["groups"],
        "summary": pair_summary(it, idx), "downsampled": data["downsampled"],
        "flags": sorted(it.flags), "figure": fig_pair(data),
        "view": view, "frame": _clamp_frame(frame, it.n_frames),
        "browse_mode": browse_mode,
        **_nav_source(path, items, name, _batch_summaries(state)),
    }
    if view == PAIR_VIEW_3D:
        context.update(_structure_context(state, it, idx, h5_path))
        context["frame"] = _clamp_frame(context["frame"], context["n_frames"])
    response = TEMPLATES.TemplateResponse(request, "pair.html", context)
    if view == PAIR_VIEW_3D:
        response.headers["Content-Security-Policy"] = STRUCTURE_CSP
    return response


def _get_interaction(state: Dict[str, Any], name: str) -> Interaction:
    """按类型名取当前数据源的 Interaction；未知类型 → 404。"""
    return _find_interaction(_source(state)[1], name)[1]


def _structure_mismatch(structure: Optional[Dict[str, Any]],
                        run_meta: Optional[Dict[str, str]]) -> bool:
    """当前结构与 h5 记录的输入是否不一致（无结构 → True，无记录 → False）。"""
    if structure is None:
        return True
    if run_meta is None:
        return False
    current = [structure.get(key, "") for key in ("topology", "trajectory")]
    recorded = [run_meta.get(key, "") for key in ("topology", "trajectory")]
    return any(_canonical_path(a) != _canonical_path(b)
               for a, b in zip(current, recorded))


def _canonical_path(path: str) -> str:
    """路径归一（abspath + realpath，兼容符号链接与相对写法差异）。"""
    return os.path.realpath(os.path.abspath(path))


def _auto_attempt(structure: Optional[Dict[str, Any]],
                  run_meta: Optional[Dict[str, str]]) -> bool:
    """3D 页是否需要自动载入（有记录且与当前结构不一致）。"""
    return run_meta is not None and _structure_mismatch(structure, run_meta)


def _structure_context(state: Dict[str, Any], interaction: Interaction,
                       idx: int, h5_path: str = "") -> Dict[str, Any]:
    """3D 视图上下文：结构状态 + existence 条 + 基团标签（pair/旧页共用）。"""
    structure, run_meta = state.get("structure"), _pair_run_meta(h5_path)
    return {
        "n_frames": structure["n_frames"] if structure
        else interaction.n_frames,
        "has_structure": structure is not None,
        "pair_label": get_exporter(interaction).get_pair_label(interaction,
                                                              idx),
        "existence": [int(value) for value in interaction.existence[idx]],
        "pt_types": pair_types(interaction, idx),
        "topology": structure["topology"] if structure else "",
        "trajectory": structure["trajectory"] if structure else "",
        "n_atoms": structure["n_atoms"] if structure else 0,
        "structure_version": structure["version"] if structure else 0,
        "run_meta": run_meta,
        "auto_attempt": _auto_attempt(structure, run_meta),
    }


def _pair_run_meta(h5_path: str) -> Optional[Dict[str, Any]]:
    """当前 pair 所属 h5 的 run/ 元信息（无 h5 / 无 run/ → None）。"""
    if not h5_path:
        return None
    try:
        record = load_run_meta(h5_path)
    except OSError:
        return None
    return _run_meta_payload(record.meta) if record else None


def _run_meta_payload(meta: RunMeta) -> Dict[str, Any]:
    """3D 页展示用的 run/ 字段（argv/interactions 转列表）。"""
    return {"topology": meta.topology, "trajectory": meta.trajectory,
            "ff": meta.ff, "strategy": meta.strategy,
            "created_at": meta.created_at, "dii_version": meta.dii_version,
            "argv": list(meta.argv),
            "interactions": list(meta.interactions)}


def _render_run(request: Request, state: Dict[str, Any], browse_mode: str,
                root: str, error: str) -> HTMLResponse:
    """计算发起页：参数表单默认值 + 最近任务。"""
    return TEMPLATES.TemplateResponse(request, "run.html", {
        "title": "new calculation",
        "browse_mode": browse_mode,
        "root": root,
        "ffs": list(FORCE_FIELDS),
        "strategies": list(STRATEGIES),
        "pairings": list(PAIRINGS),
        "interactions": list(INTERACTIONS),
        "form": _run_form_defaults(),
        "error": error,
        "recent_jobs": _recent_jobs(state["jobs"]),
        **_nav(state),
    })


def _run_form_defaults() -> Dict[str, Any]:
    """表单默认值（默认全选 8 类相互作用）。"""
    return {"topology": "", "trajectory": "", "output": "",
            "ff": DEFAULT_FF, "strategy": DEFAULT_STRATEGY,
            "interactions": list(INTERACTIONS), "select": "",
            "pairing": ANY_PAIRING, "mark_redundant": False}


def _render_job(request: Request, state: Dict[str, Any], job_id: str,
                browse_mode: str, root: str) -> HTMLResponse:
    """任务详情页：状态、日志尾部、产出目录与 h5 列表。"""
    job = _job_or_404(state["jobs"], job_id)
    return TEMPLATES.TemplateResponse(request, "job.html", {
        "title": "job", "browse_mode": browse_mode, "root": root,
        "job": _job_view(job), "recent_jobs": _recent_jobs(state["jobs"]),
        **_nav(state),
    })


def _job_view(job: Job) -> Dict[str, Any]:
    """任务页所需视图（日志只取尾部 RUN_LOG_LINES 行）。"""
    return {"id": job.id, "status": job.status,
            "summary": job.params.summary(),
            "params_preview": job.params.preview(),
            "lines": job.lines()[-RUN_LOG_LINES:],
            "lines_total": job.line_total(),
            "returncode": job.returncode,
            "elapsed": _seconds(job.elapsed()),
            "output_dir": job.output_dir,
            "h5_files": job.h5_files(),
            "error": job.error}


def _job_or_404(manager: JobManager, job_id: str) -> Job:
    """按 id 取任务；未知 → 404（页面走 HTML 错误页，API 走 JSON）。"""
    job = manager.get(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail=JOB_NOT_FOUND)
    return job


def _recent_jobs(manager: JobManager) -> List[Dict[str, Any]]:
    """最近任务摘要（新 → 旧）。"""
    return [{"id": job.id, "status": job.status,
             "summary": job.params.summary(),
             "elapsed": _seconds(job.elapsed()),
             "output": job.output_dir}
            for job in manager.recent()]


def _seconds(value: float) -> float:
    """秒数保留 ELAPSED_DECIMALS 位小数（显示用）。"""
    return round(value, ELAPSED_DECIMALS)


def _start_run(manager: JobManager, root: str,
               fields: Dict[str, List[str]]) -> RedirectResponse:
    """校验表单并启动任务；失败跳回表单页带回提示（与 /select 一致）。"""
    try:
        job = manager.start(parse_run_form(fields, root))
    except ValueError as error:
        return RedirectResponse(f"/run?error={quote(str(error))}",
                                status_code=303)
    return RedirectResponse(f"/run/{job.id}", status_code=303)


def _form_fields(body: bytes) -> Dict[str, List[str]]:
    """URL-encoded 表单体 → 字段多值映射（保留空值，支持重复键）。"""
    try:
        return parse_qs(body.decode("utf-8"), keep_blank_values=True)
    except UnicodeDecodeError:
        return {}


def _job_payload(job: Job, since: int) -> Dict[str, Any]:
    """轮询接口响应（lines 为 since 之后的日志增量）。"""
    lines, next_index = job.lines_since(since)
    return {"id": job.id, "status": job.status,
            "returncode": job.returncode,
            "elapsed": _seconds(job.elapsed()),
            "lines": lines, "next": next_index,
            "output_dir": job.output_dir, "h5_files": job.h5_files(),
            "error": job.error}


def _parse_since(raw: str) -> int:
    """日志游标：非整数或负数一律视为 0。"""
    try:
        value = int(raw)
    except (TypeError, ValueError):
        return 0
    return max(value, 0)


def _cancel_response(manager: JobManager, job_id: str) -> Response:
    """取消运行中任务（未知 → 404，已结束 → 400）。"""
    _job_or_404(manager, job_id)
    if not manager.cancel(job_id):
        raise HTTPException(status_code=400, detail=CANNOT_CANCEL)
    return JSONResponse({"status": STATUS_CANCELLED})


def _run_log_response(job: Job) -> Response:
    """任务日志全文（text/plain 附件）。"""
    body = "\n".join(job.lines())
    if body:
        body += "\n"
    filename = f"dii-run-{job.id}.log"
    return Response(body, media_type="text/plain",
                    headers={"Content-Disposition":
                             f'attachment; filename="{filename}"'})


def _to_csv(table: Dict[str, Any]) -> str:
    """表格数据 → CSV 字符串（含表头，行序与视图一致）。"""
    buffer = io.StringIO()
    writer = csv.writer(buffer)
    writer.writerow(table["columns"])
    for row in table["rows"]:
        writer.writerow([row[column] for column in table["columns"]])
    return buffer.getvalue()
