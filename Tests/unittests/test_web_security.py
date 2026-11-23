# -*- coding: utf-8 -*-
"""Web 安全加固测试：root 边界、/select 白名单、Host 校验、安全头、token 限流。"""

import os
import time
from pathlib import Path
from urllib.parse import unquote

import pytest

from DuIvyInteractions.DII import _browse_root, _print_web_banner
from DuIvyInteractions.web import app as web_app
from DuIvyInteractions.web.fsbrowse import ROOT_ERROR, confine, list_directory

TOKEN = "s3cr3t-token"
UNSUPPORTED = "unsupported path (expect a directory or .h5 file)"
FAIL_TEXT = "too many failed attempts"
SECURITY_HEADERS = {
    "x-content-type-options": "nosniff",
    "x-frame-options": "DENY",
    "referrer-policy": "no-referrer",
    "content-security-policy": (
        "default-src 'self'; script-src 'self' 'unsafe-inline'; "
        "style-src 'self' 'unsafe-inline'; img-src 'self' data:; "
        "connect-src 'self'"),
}


@pytest.fixture(autouse=True)
def clean_state(monkeypatch, tmp_path):
    """隔离 recent.json，并在每个用例前后清空限流计数。"""
    monkeypatch.setattr("DuIvyInteractions.visualizers.data.RECENT_PATH",
                        tmp_path / "recent.json")
    web_app._AUTH_FAILURES.clear()
    yield
    web_app._AUTH_FAILURES.clear()


def _client(base_url=None, **kwargs):
    """按关键字参数创建 TestClient。"""
    from fastapi.testclient import TestClient
    app = web_app.create_app(**kwargs)
    if base_url:
        return TestClient(app, base_url=base_url)
    return TestClient(app)


def _root(tmp_path: Path) -> Path:
    """构造含 sub/nested.h5 的受限根目录。"""
    root = tmp_path / "root"
    (root / "sub").mkdir(parents=True)
    (root / "sub" / "nested.h5").write_bytes(b"x")
    return root


def _symlink(link: Path, target: Path) -> None:
    """创建目录符号链接（平台不支持时跳过）。"""
    try:
        link.symlink_to(target, target_is_directory=True)
    except OSError:
        pytest.skip("symlinks not supported")


def _file_symlink(link: Path, target: Path) -> None:
    """创建文件符号链接（平台不支持时跳过）。"""
    try:
        link.symlink_to(target)
    except OSError:
        pytest.skip("symlinks not supported")


def _fake_source(monkeypatch) -> None:
    """把 parse_datasource 替换为伪实现（不依赖真实 h5 内容）。"""
    def fake(path: str):
        base = Path(path)
        if base.is_dir():
            return [("fake", str(base / "a.h5"), object())]
        return [("fake", str(base), object())]

    monkeypatch.setattr(web_app, "parse_datasource", fake)


def _assert_security_headers(response) -> None:
    """断言四个安全响应头均存在且值符合预期。"""
    for name, value in SECURITY_HEADERS.items():
        assert response.headers[name] == value


# ============================================================
# confine 纯函数
# ============================================================

class TestConfine:

    def test_absolute_path_inside_root(self, tmp_path):
        root = _root(tmp_path)
        assert confine(str(root), str(root / "sub")) == (root / "sub").resolve()

    def test_relative_path_inside_root(self, tmp_path):
        root = _root(tmp_path)
        assert confine(str(root), "sub") == (root / "sub").resolve()

    def test_root_itself_allowed(self, tmp_path):
        root = _root(tmp_path)
        assert confine(str(root), ".") == root.resolve()
        assert confine(str(root), "") == root.resolve()

    def test_parent_escape_rejected(self, tmp_path):
        root = _root(tmp_path)
        with pytest.raises(ValueError, match="path outside root"):
            confine(str(root), "..")

    def test_absolute_outside_rejected(self, tmp_path):
        root = _root(tmp_path)
        with pytest.raises(ValueError, match="path outside root"):
            confine(str(root), str(tmp_path / "elsewhere"))

    def test_symlink_outside_rejected(self, tmp_path):
        root = _root(tmp_path)
        outside = tmp_path / "outside"
        outside.mkdir()
        _symlink(root / "escape", outside)
        with pytest.raises(ValueError, match="path outside root"):
            confine(str(root), str(root / "escape"))

    def test_symlink_inside_allowed(self, tmp_path):
        root = _root(tmp_path)
        _symlink(root / "inside", root / "sub")
        assert confine(str(root), str(root / "inside")) == (root / "sub").resolve()

    def test_empty_root_is_unrestricted(self, tmp_path):
        assert confine("", str(tmp_path)) == tmp_path.resolve()
        assert confine("", "") == Path.home().resolve()

    def test_error_message_carries_no_path(self, tmp_path):
        root = _root(tmp_path)
        with pytest.raises(ValueError) as excinfo:
            confine(str(root), "..")
        assert str(excinfo.value) == ROOT_ERROR
        assert str(root) not in str(excinfo.value)


# ============================================================
# list_directory 的 root 边界
# ============================================================

class TestListDirectoryRoot:

    def test_parent_is_none_at_root(self, tmp_path):
        root = _root(tmp_path)
        data = list_directory(str(root), root=str(root))
        assert data["parent"] is None
        assert data["root"] == str(root.resolve())
        assert data["rel"] == "."

    def test_rel_and_parent_in_subdir(self, tmp_path):
        root = _root(tmp_path)
        data = list_directory(str(root / "sub"), root=str(root))
        assert data["rel"] == "sub"
        assert data["parent"] == str(root.resolve())

    def test_rel_is_posix_for_nested_dir(self, tmp_path):
        root = _root(tmp_path)
        nested = root / "sub" / "nested"
        nested.mkdir()
        data = list_directory(str(nested), root=str(root))
        assert data["rel"] == "sub/nested"

    def test_hidden_dirs_skipped(self, tmp_path):
        root = _root(tmp_path)
        (root / ".hidden").mkdir()
        assert ".hidden" not in list_directory(str(root), root=str(root))["dirs"]

    def test_dirs_sorted_case_insensitive(self, tmp_path):
        root = tmp_path / "ordered"
        (root / "beta").mkdir(parents=True)
        (root / "Alpha").mkdir()
        data = list_directory(str(root), root=str(root))
        assert data["dirs"] == ["Alpha", "beta"]

    def test_symlink_dir_outside_root_skipped(self, tmp_path):
        root = _root(tmp_path)
        outside = tmp_path / "outside"
        outside.mkdir()
        _symlink(root / "escape", outside)
        assert "escape" not in list_directory(str(root), root=str(root))["dirs"]

    def test_symlink_dir_inside_root_listed(self, tmp_path):
        root = _root(tmp_path)
        _symlink(root / "inside", root / "sub")
        assert "inside" in list_directory(str(root), root=str(root))["dirs"]

    def test_symlink_h5_outside_root_skipped(self, tmp_path):
        root = _root(tmp_path)
        outside = tmp_path / "outside"
        outside.mkdir()
        target = outside / "target.h5"
        target.write_bytes(b"x" * 2048)
        _file_symlink(root / "link.h5", target)
        assert list_directory("", root=str(root))["files"] == []
        inside = root / "inside.h5"
        inside.write_bytes(b"x")
        _file_symlink(root / "alias.h5", inside)
        names = [item["name"]
                 for item in list_directory("", root=str(root))["files"]]
        assert names == ["alias.h5", "inside.h5"]

    def test_absolute_outside_path_raises(self, tmp_path):
        root = _root(tmp_path)
        with pytest.raises(ValueError) as excinfo:
            list_directory(str(tmp_path), root=str(root))
        assert str(excinfo.value) == ROOT_ERROR
        assert str(root) not in str(excinfo.value)

    def test_parent_escape_path_raises(self, tmp_path):
        root = _root(tmp_path)
        with pytest.raises(ValueError, match="path outside root"):
            list_directory(str(root / ".."), root=str(root))

    def test_root_boundary_disabled_by_default(self, tmp_path):
        root = _root(tmp_path)
        data = list_directory(str(root), root="")
        assert data["parent"] == str(tmp_path.resolve())
        assert data["root"] == ""
        assert data["rel"] == "."

    def test_missing_path_keeps_prefix(self, tmp_path):
        root = _root(tmp_path)
        with pytest.raises(ValueError, match="path not found"):
            list_directory(str(root / "missing"), root=str(root))

    def test_file_path_keeps_prefix(self, tmp_path):
        root = _root(tmp_path)
        target = root / "single.h5"
        target.write_bytes(b"x")
        with pytest.raises(ValueError, match="not a directory"):
            list_directory(str(target), root=str(root))


# ============================================================
# /select 路径限制与类型白名单
# ============================================================

class TestSelect:

    def test_directory_allowed(self, monkeypatch, tmp_path):
        _fake_source(monkeypatch)
        root = _root(tmp_path)
        client = _client(root=str(root))
        resp = client.post("/select", data={"path": str(root)},
                           follow_redirects=False)
        assert resp.status_code == 303
        assert resp.headers["location"] == "/"

    def test_relative_directory_allowed(self, monkeypatch, tmp_path):
        _fake_source(monkeypatch)
        root = _root(tmp_path)
        client = _client(root=str(root))
        resp = client.post("/select", data={"path": "sub"},
                           follow_redirects=False)
        assert resp.status_code == 303
        assert resp.headers["location"] == "/"

    def test_h5_file_allowed(self, monkeypatch, tmp_path):
        _fake_source(monkeypatch)
        target = tmp_path / "data.h5"
        target.write_bytes(b"x")
        client = _client(root=str(tmp_path))
        resp = client.post("/select", data={"path": str(target)},
                           follow_redirects=False)
        assert resp.status_code == 303
        assert resp.headers["location"] == "/"

    def test_uppercase_h5_allowed(self, monkeypatch, tmp_path):
        _fake_source(monkeypatch)
        target = tmp_path / "DATA.H5"
        target.write_bytes(b"x")
        client = _client(root=str(tmp_path))
        resp = client.post("/select", data={"path": str(target)},
                           follow_redirects=False)
        assert resp.status_code == 303
        assert resp.headers["location"] == "/"

    def test_outside_path_rejected(self, monkeypatch, tmp_path):
        _fake_source(monkeypatch)
        root = _root(tmp_path)
        client = _client(root=str(root))
        resp = client.post("/select", data={"path": str(tmp_path)},
                           follow_redirects=False)
        assert resp.status_code == 303
        assert unquote(resp.headers["location"]) == f"/?error={ROOT_ERROR}"

    def test_parent_escape_rejected(self, monkeypatch, tmp_path):
        _fake_source(monkeypatch)
        root = _root(tmp_path)
        client = _client(root=str(root))
        resp = client.post("/select", data={"path": str(root / "..")},
                           follow_redirects=False)
        assert unquote(resp.headers["location"]) == f"/?error={ROOT_ERROR}"

    def test_non_h5_file_rejected(self, monkeypatch, tmp_path):
        _fake_source(monkeypatch)
        root = _root(tmp_path)
        target = root / "notes.txt"
        target.write_text("x", encoding="utf-8")
        client = _client(root=str(root))
        resp = client.post("/select", data={"path": str(target)},
                           follow_redirects=False)
        assert unquote(resp.headers["location"]) == f"/?error={UNSUPPORTED}"

    def test_empty_path_reported(self, tmp_path):
        client = _client(root=str(tmp_path))
        resp = client.post("/select", data={"path": "   "},
                           follow_redirects=False)
        assert unquote(resp.headers["location"]) == "/?error=no path provided"

    def test_default_root_keeps_external_directory(self, monkeypatch, tmp_path):
        _fake_source(monkeypatch)
        outside = tmp_path / "outside"
        outside.mkdir()
        (outside / "a.h5").write_bytes(b"x")
        client = _client()
        resp = client.post("/select", data={"path": str(outside)},
                           follow_redirects=False)
        assert resp.status_code == 303
        assert resp.headers["location"] == "/"

    def test_default_root_rejects_non_h5_file(self, tmp_path):
        target = tmp_path / "notes.txt"
        target.write_text("x", encoding="utf-8")
        client = _client()
        resp = client.post("/select", data={"path": str(target)},
                           follow_redirects=False)
        assert unquote(resp.headers["location"]) == f"/?error={UNSUPPORTED}"


# ============================================================
# /api/list
# ============================================================

class TestListApi:

    def test_outside_path_404(self, tmp_path):
        root = _root(tmp_path)
        client = _client(browse_mode="server", root=str(root))
        resp = client.get("/api/list", params={"path": str(tmp_path)})
        assert resp.status_code == 404
        assert resp.json()["detail"] == ROOT_ERROR

    def test_root_listing_reports_boundary(self, tmp_path):
        root = _root(tmp_path)
        client = _client(browse_mode="server", root=str(root))
        data = client.get("/api/list", params={"path": str(root)}).json()
        assert data["parent"] is None
        assert data["root"] == str(root.resolve())
        assert data["rel"] == "."

    def test_subdir_listing_reports_rel(self, tmp_path):
        root = _root(tmp_path)
        client = _client(browse_mode="server", root=str(root))
        data = client.get("/api/list",
                          params={"path": str(root / "sub")}).json()
        assert data["parent"] == str(root.resolve())
        assert data["rel"] == "sub"

    def test_default_root_unrestricted(self, tmp_path):
        (tmp_path / "a.h5").write_bytes(b"x")
        client = _client(browse_mode="server")
        resp = client.get("/api/list", params={"path": str(tmp_path)})
        assert resp.status_code == 200
        assert resp.json()["parent"] == str(tmp_path.parent)

    def test_missing_path_404_keeps_prefix(self, tmp_path):
        root = _root(tmp_path)
        client = _client(browse_mode="server", root=str(root))
        resp = client.get("/api/list", params={"path": str(root / "missing")})
        assert resp.status_code == 404
        assert "path not found" in resp.json()["detail"]


# ============================================================
# Host 头校验
# ============================================================

class TestHostGuard:

    def test_allowed_host_passes(self):
        client = _client(base_url="http://allowed.example",
                         allowed_hosts=["allowed.example"])
        assert client.get("/").status_code == 200

    def test_allowed_host_with_port_passes(self):
        client = _client(base_url="http://allowed.example:8000",
                         allowed_hosts=["allowed.example"])
        assert client.get("/").status_code == 200

    def test_ipv6_loopback_host_passes(self):
        client = _client(allowed_hosts=["::1", "localhost"])
        resp = client.get("/", headers={"host": "[::1]:8888"})
        assert resp.status_code == 200

    def test_malicious_host_forbidden(self):
        client = _client(base_url="http://evil.example:8000",
                         allowed_hosts=["allowed.example"])
        resp = client.get("/")
        assert resp.status_code == 403
        assert resp.text == "forbidden host"
        assert resp.headers["content-type"].startswith("text/plain")

    def test_allowlist_disabled_by_default(self):
        client = _client(base_url="http://evil.example")
        assert client.get("/").status_code == 200

    def test_host_checked_before_token(self):
        client = _client(base_url="http://evil.example",
                         allowed_hosts=["allowed.example"], token=TOKEN)
        assert client.get("/").status_code == 403


class TestHostName:

    @pytest.mark.parametrize("raw, expected", [
        ("example.com", "example.com"),
        ("example.com:8000", "example.com"),
        ("127.0.0.1:5000", "127.0.0.1"),
        ("[::1]:8888", "::1"),
        ("[::1]", "::1"),
        ("::1", "::1"),
        ("  EVIL.example:80  ", "evil.example"),
        ("", ""),
    ])
    def test_extracts_host(self, raw, expected):
        assert web_app._host_name(raw) == expected


# ============================================================
# 安全响应头
# ============================================================

class TestSecurityHeaders:

    def test_page_has_headers(self):
        _assert_security_headers(_client().get("/"))

    def test_api_has_headers(self, tmp_path):
        resp = _client(browse_mode="server").get(
            "/api/list", params={"path": str(tmp_path)})
        _assert_security_headers(resp)

    def test_unauthorized_has_headers(self):
        resp = _client(token=TOKEN).get("/")
        assert resp.status_code == 401
        _assert_security_headers(resp)

    def test_forbidden_host_has_headers(self):
        client = _client(base_url="http://evil.example",
                         allowed_hosts=["allowed.example"])
        resp = client.get("/")
        assert resp.status_code == 403
        _assert_security_headers(resp)

    def test_csp_allows_inline_scripts(self):
        csp = _client().get("/").headers["content-security-policy"]
        assert "script-src 'self' 'unsafe-inline'" in csp


# ============================================================
# token 失败限流
# ============================================================

class TestAuthRateLimit:

    def test_limit_plus_one_is_429(self):
        client = _client(token=TOKEN)
        limit = web_app.AUTH_FAIL_LIMIT
        for _ in range(limit):
            assert client.get("/?token=wrong").status_code == 401
        resp = client.get("/?token=wrong")
        assert resp.status_code == 429
        assert resp.text == FAIL_TEXT
        assert resp.headers["content-type"].startswith("text/plain")

    def test_expired_failures_do_not_block(self):
        client = _client(token=TOKEN)
        for _ in range(web_app.AUTH_FAIL_LIMIT):
            client.get("/?token=wrong")
        host = next(iter(web_app._AUTH_FAILURES))
        stale = time.monotonic() - web_app.AUTH_FAIL_WINDOW - 1
        web_app._AUTH_FAILURES[host] = [stale] * web_app.AUTH_FAIL_LIMIT
        assert client.get("/?token=wrong").status_code == 401

    def test_success_clears_counter(self):
        failing = _client(token=TOKEN)
        successful = _client(token=TOKEN)
        for _ in range(5):
            assert failing.get("/?token=wrong").status_code == 401
        assert successful.get(f"/?token={TOKEN}").status_code == 200
        for _ in range(web_app.AUTH_FAIL_LIMIT):
            assert failing.get("/?token=wrong").status_code == 401

    def test_failures_stay_bounded(self):
        client = _client(token=TOKEN)
        for _ in range(web_app.AUTH_FAIL_LIMIT + 5):
            client.get("/?token=wrong")
        host = next(iter(web_app._AUTH_FAILURES))
        assert len(web_app._AUTH_FAILURES[host]) == web_app.AUTH_FAIL_LIMIT


# ============================================================
# dii web 的默认浏览根策略（本地不限制 / 远程强制 cwd）
# ============================================================

class TestBrowseRootPolicy:

    def test_explicit_root_wins_in_both_modes(self):
        expected = os.path.abspath("/tmp/x")
        assert _browse_root("/tmp/x", True) == expected
        assert _browse_root("/tmp/x", False) == expected

    def test_default_local_is_unrestricted(self):
        assert _browse_root("", True) == ""

    def test_default_remote_uses_cwd(self):
        assert _browse_root("", False) == os.getcwd()

    def test_banner_unrestricted(self, capsys):
        _print_web_banner(["http://127.0.0.1:8888"], "", "127.0.0.1", True, "")
        out = capsys.readouterr().out
        assert "unrestricted (loopback only)" in out
        assert "restricted to" not in out

    def test_banner_restricted(self, capsys):
        urls = ["http://127.0.0.1:8888"]
        _print_web_banner(urls, "", "127.0.0.1", True, "/some/dir")
        out = capsys.readouterr().out
        assert "restricted to: /some/dir" in out
        assert "unrestricted" not in out
