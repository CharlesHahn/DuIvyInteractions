# -*- coding: utf-8 -*-
"""Web 传输层：FastAPI 路由 + Jinja 模板。"""

import csv
import io
import math
import secrets
import threading
import time
from collections import OrderedDict
from pathlib import Path
from typing import Any, Callable, Dict, FrozenSet, List, Optional, Tuple
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
    frame_activity, frame_view, get_exporter, interaction_meta, load_recent,
    matrix_view, metric_stats, overview, pair_drilldown, pair_frames,
    pair_summary, parse_datasource, string_stats, type_summary,
)
from ..visualizers.structure import (
    FRAME_CACHE_SIZE, Highlight, alignment_errors, closest_contact,
    contact_label, contact_pdb, frame_view as structure_frame_view, group_atoms,
    group_label, highlights, involves_water, load_universe,
)
from .fsbrowse import (
    BYTES_PER_KB, H5_SUFFIX, LIST_LIMIT, SIZE_DECIMALS, confine,
    list_directory,
)
from .picker import pick_directory, pick_file

_WEB_DIR = Path(__file__).parent
TEMPLATES = Jinja2Templates(directory=str(_WEB_DIR / "templates"))
FRAME_LIMIT_MAX = 2000  # 帧检查器单次最多返回的活跃对数
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
STRUCTURE_SUFFIXES = frozenset({
    ".tpr", ".xtc", ".trr", ".gro", ".pdb", ".dcd", ".nc", ".itp", ".psf",
    ".prmtop", ".parm7", ".mol2"})
SOLVENT_FULL = "full"
SOLVENT_HIDE = "hide"
PICK_DIR = "dir"
PICK_FILE = "file"
DATASET_KIND = "dataset"
STRUCTURE_KIND = "structure"
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
STRUCTURE_CSP = (
    "default-src 'self'; script-src 'self' 'unsafe-inline' 'unsafe-eval'; "
    "style-src 'self' 'unsafe-inline'; img-src 'self' data:; "
    "connect-src 'self'")
_AUTH_FAILURES: Dict[str, List[float]] = {}
_AUTH_LOCK = threading.Lock()
_FRAME_CACHE: "OrderedDict[Tuple[Any, ...], Tuple[str, List[int]]]" = OrderedDict()
_FRAME_LOCK = threading.Lock()


def create_app(token: str = "", browse_mode: str = "picker", root: str = "",
               allowed_hosts: Optional[List[str]] = None) -> FastAPI:
    """创建 FastAPI 应用（数据源状态 = 服务端内存，单用户）。"""
    state: Dict[str, Any] = {"source": None, "structure": None}
    app = FastAPI(title="DII Web", docs_url=None, redoc_url=None,
                  openapi_url=None)
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
    """安装安全响应头中间件（不覆盖响应已有的同名头）。"""

    @app.middleware("http")
    async def security_headers(request: Request, call_next: Callable) -> Response:
        response = await call_next(request)
        for name, value in SECURITY_HEADERS.items():
            response.headers.setdefault(name, value)
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
    _register_detail(app, state)
    _register_structure_page(app, state, browse_mode)
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


def _register_structure_page(app: FastAPI, state: Dict[str, Any],
                             browse_mode: str) -> None:
    """3D 结构页（Mol* viewer + 帧滑块）。"""

    @app.get("/structure/{name}/{idx}", response_class=HTMLResponse)
    def structure_page(request: Request, name: str, idx: int) -> Response:
        return _render_structure(request, state, name, idx, browse_mode)


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
    _register_cell_api(app, state)
    _register_pair_api(app, state)
    _register_csv_api(app, state)
    _register_list_api(app, browse_mode, root)
    _register_pick_api(app, browse_mode)
    _register_structure_api(app, state, root)


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


def _register_structure_api(app: FastAPI, state: Dict[str, Any],
                            root: str) -> None:
    """3D 结构接口：选择结构 + 逐帧 PDB/接触对。"""

    @app.post("/api/structure/select", response_model=None)
    async def select_structure(request: Request) -> Response:
        body = await request.body()
        return await run_in_threadpool(
            _select_structure, state, root,
            _form_value(body, "topology"), _form_value(body, "trajectory"))

    @app.get("/api/structure/view", response_model=None)
    def api_structure_view(type: str = "", idx: int = 0, frame: int = 0,
                           solvent: str = SOLVENT_FULL) -> Response:
        return _structure_view(state, type, idx, frame, solvent)


def _select_structure(state: Dict[str, Any], root: str, topology: str,
                      trajectory: str) -> Response:
    """打开并记录结构；两个路径先受 root 限制。"""
    try:
        top_path = _structure_path(root, topology)
        traj_path = _structure_path(root, trajectory)
        universe = load_universe(str(top_path), str(traj_path))
    except ValueError as error:
        return JSONResponse({"detail": str(error)}, status_code=400)
    structure = {
        "topology": str(top_path), "trajectory": str(traj_path),
        "universe": universe, "n_atoms": int(universe.atoms.n_atoms),
        "n_frames": len(universe.trajectory),
        "has_box": universe.dimensions is not None,
        "lock": threading.Lock(),
    }
    state["structure"] = structure
    return JSONResponse({"ok": True, "n_atoms": structure["n_atoms"],
                         "n_frames": structure["n_frames"],
                         "has_box": structure["has_box"]})


def _structure_path(root: str, path: str) -> Path:
    """拓扑/轨迹路径（非空，且在 root 之内）。"""
    if not path.strip():
        raise ValueError("structure path required")
    return confine(root, path)


def _structure_view(state: Dict[str, Any], name: str, idx: int, frame: int,
                    solvent: str = SOLVENT_FULL) -> Response:
    """单帧视图响应：PDB 文本 + existence + 最近接触原子对 + 溶剂模式。"""
    structure = state.get("structure")
    if not structure:
        return JSONResponse({"detail": STRUCTURE_NOT_SELECTED},
                            status_code=400)
    _, interaction = _find_interaction(_source(state)[1], name)
    if not 0 <= idx < interaction.n_pairs:
        raise HTTPException(status_code=404,
                            detail=f"pair index out of range: {idx}")
    if not 0 <= frame < structure["n_frames"] or frame >= interaction.n_frames:
        raise HTTPException(status_code=404,
                            detail=f"frame out of range: {frame}")
    errors = alignment_errors(structure["universe"], interaction)
    if errors:
        return JSONResponse({"detail": errors[0], "errors": errors},
                            status_code=400)
    water = involves_water(interaction, idx)
    mode = _solvent_mode(solvent, water)
    active = bool(interaction.existence[idx, frame])
    contact, group_views, contact_view = _pair_geometry(
        structure, interaction, idx, frame, active)
    text, shown = _frame_text(structure, frame, mode,
                              _pair_atoms(interaction, idx))
    rows = _shown_row_map(shown, mode)
    return JSONResponse({
        "frame": frame, "n_frames": structure["n_frames"], "active": active,
        "existence": [int(value) for value in interaction.existence[idx]],
        "contact": _shown_contact(contact, rows),
        "groups": _shown_groups(interaction, idx, rows,
                                structure["universe"]),
        "highlights": [{"label": view.label, "pdb": view.pdb}
                       for view in group_views],
        "contact_pdb": contact_view,
        "pdb": text,
        "involves_water": water,
        "solvent": mode,
        "n_atoms_total": structure["n_atoms"],
        "n_atoms_shown": len(shown),
        "contact_label": _pair_contact_label(interaction, idx, contact,
                                             structure["universe"]),
    })


def _solvent_mode(solvent: str, water: bool) -> str:
    """实际生效的溶剂模式：仅当请求 hide 且该对与水无关时才隐藏水。"""
    return SOLVENT_HIDE if solvent == SOLVENT_HIDE and not water else SOLVENT_FULL


def _pair_atoms(interaction: Interaction, idx: int) -> FrozenSet[int]:
    """该对全部基团原子的全局索引（保证标注原子出现在输出中）。"""
    return frozenset(atom for group in interaction.groups[idx]
                     for atom in group.atom_indices)


def _pair_contact_label(interaction: Interaction, idx: int,
                        contact: Optional[Tuple[int, int]],
                        universe: Any) -> str:
    """接触对标签：全局索引先换成组内位置（多基团只看首末基团）。"""
    if contact is None:
        return contact_label(interaction, idx, None, universe)
    groups = interaction.groups[idx]
    return contact_label(interaction, idx, (
        _group_position(groups[0], contact[0]),
        _group_position(groups[-1], contact[1])), universe)


def _group_position(group: Any, global_index: int) -> int:
    """全局原子索引在基团 atoms 中的位置（缺失时为 -1）。"""
    for position, atom in enumerate(group.atoms):
        if atom.atom_global_idx == global_index:
            return position
    return -1


def _shown_groups(interaction: Interaction, idx: int,
                  rows: Optional[Dict[int, int]],
                  universe: Any) -> List[Dict[str, Any]]:
    """该对每个基团的标签与展示行号原子（编号取 universe，原子与 contact 同映射）。"""
    return [{"label": group_label(interaction, group, universe),
             "atoms": _shown_atoms(indices, rows)}
            for group, indices in zip(interaction.groups[idx],
                                      group_atoms(interaction, idx))]


def _shown_row_map(shown: List[int], mode: str) -> Optional[Dict[int, int]]:
    """hide 模式下全局原子索引 → 展示行号；full 模式索引即行号，返回 None。"""
    if mode == SOLVENT_FULL:
        return None
    return {index: row for row, index in enumerate(shown)}


def _shown_atoms(indices: List[int],
                 rows: Optional[Dict[int, int]]) -> List[int]:
    """全局原子索引 → 展示行号（rows 为 None 时恒等；展示中缺失的原子跳过）。"""
    if rows is None:
        return [int(index) for index in indices]
    return [rows[int(index)] for index in indices if int(index) in rows]


def _shown_contact(contact: Optional[Tuple[int, int]],
                   rows: Optional[Dict[int, int]]) -> Optional[List[int]]:
    """把全局索引接触对映射为 PDB 展示行号（rows 为 None 时索引即行号）。"""
    if not contact:
        return None
    if rows is None:
        return [int(contact[0]), int(contact[1])]
    if not all(int(index) in rows for index in contact):
        return None
    return [rows[int(index)] for index in contact]


def _frame_text(structure: Dict[str, Any], frame: int,
                mode: str = SOLVENT_FULL,
                keep_atoms: FrozenSet[int] = frozenset()
                ) -> Tuple[str, List[int]]:
    """单帧 PDB 文本与展示索引（LRU 缓存，键含拓扑/轨迹/帧/溶剂模式）。"""
    key = (structure["topology"], structure["trajectory"], frame, mode)
    if mode == SOLVENT_HIDE:
        key += (keep_atoms,)
    with _FRAME_LOCK:
        cached = _FRAME_CACHE.get(key)
        if cached is not None:
            _FRAME_CACHE.move_to_end(key)
            return cached
    with structure["lock"]:
        text, shown = structure_frame_view(structure["universe"], frame,
                                           hide_water=mode == SOLVENT_HIDE,
                                           keep_atoms=keep_atoms)
    with _FRAME_LOCK:
        _FRAME_CACHE[key] = (text, shown)
        _FRAME_CACHE.move_to_end(key)
        while len(_FRAME_CACHE) > FRAME_CACHE_SIZE:
            _FRAME_CACHE.popitem(last=False)
    return text, shown


def _pair_geometry(structure: Dict[str, Any], interaction: Interaction, idx: int,
                   frame: int, active: bool
                   ) -> Tuple[Optional[Tuple[int, int]], List[Highlight],
                              Optional[str]]:
    """同一帧内取最近接触对、各基团小结构与接触对小结构（只设置一次帧）。"""
    groups = interaction.groups[idx]
    universe = structure["universe"]
    with structure["lock"]:
        universe.trajectory[frame]
        contact = (closest_contact(universe, groups[0].atom_indices,
                                   groups[-1].atom_indices) if active else None)
        return (contact, highlights(interaction, idx, universe),
                contact_pdb(universe, contact))


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


def _select(state: Dict[str, Any], path: str, root: str = "") -> RedirectResponse:
    """解析数据源并跳回总览；失败或跳过文件时带回提示（不静默）。"""
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
    state["source"] = (path, items)
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
                  warn: str, browse_mode: str) -> HTMLResponse:
    """总览页：唯一数据源入口 + 类型卡片。"""
    path, items = _source(state)
    return TEMPLATES.TemplateResponse(request, "index.html", {
        "title": "overview", "cards": overview(items),
        "recent": load_recent(), "error": error, "warning": warn,
        "browse_mode": browse_mode,
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


def _render_structure(request: Request, state: Dict[str, Any], name: str,
                      idx: int, browse_mode: str) -> HTMLResponse:
    """3D 结构页：pair 信息 + 当前结构状态 + existence 条 + 浏览方式。"""
    path, items = _source(state)
    _, interaction = _find_interaction(items, name)
    if not 0 <= idx < interaction.n_pairs:
        raise HTTPException(status_code=404,
                            detail=f"pair index out of range: {idx}")
    structure = state.get("structure")
    response = TEMPLATES.TemplateResponse(request, "structure.html", {
        "title": "structure", "type": name, "idx": idx,
        "n_frames": structure["n_frames"] if structure
        else interaction.n_frames,
        "has_structure": structure is not None,
        "pair_label": get_exporter(interaction).get_pair_label(interaction,
                                                              idx),
        "existence": [int(value) for value in interaction.existence[idx]],
        "topology": structure["topology"] if structure else "",
        "trajectory": structure["trajectory"] if structure else "",
        "n_atoms": structure["n_atoms"] if structure else 0,
        "browse_mode": browse_mode,
        **_nav_source(path, items, name),
    })
    response.headers["Content-Security-Policy"] = STRUCTURE_CSP
    return response


def _to_csv(table: Dict[str, Any]) -> str:
    """表格数据 → CSV 字符串（含表头，行序与视图一致）。"""
    buffer = io.StringIO()
    writer = csv.writer(buffer)
    writer.writerow(table["columns"])
    for row in table["rows"]:
        writer.writerow([row[column] for column in table["columns"]])
    return buffer.getvalue()
