# -*- coding: utf-8 -*-
"""Web 鉴权与目录浏览测试：URL 推导、picker、fsbrowse、token 中间件、API。"""

import ipaddress
from pathlib import Path

import pytest

from DuIvyInteractions.DII import _is_loopback_host, _lan_urls, _local_ipv4s
from DuIvyInteractions.web import fsbrowse
from DuIvyInteractions.web.picker import picker_available

TOKEN = "s3cr3t-token"
PORT = 8000


@pytest.fixture(autouse=True)
def isolate_recent(monkeypatch, tmp_path):
    monkeypatch.setattr("DuIvyInteractions.visualizers.data.RECENT_PATH",
                        tmp_path / "recent.json")


def _client(**kwargs):
    from fastapi.testclient import TestClient
    from DuIvyInteractions.web.app import create_app
    return TestClient(create_app(**kwargs))


def _tree(root: Path) -> dict:
    data = fsbrowse.list_directory(str(root))
    return data


# ============================================================
# URL 推导
# ============================================================

class TestLoopbackHost:

    def test_accepts_local_names(self):
        assert _is_loopback_host("localhost")
        assert _is_loopback_host("127.0.0.1")
        assert _is_loopback_host("::1")

    def test_rejects_network_and_wildcard(self):
        assert not _is_loopback_host("0.0.0.0")
        assert not _is_loopback_host("192.168.1.5")
        assert not _is_loopback_host("")


class TestLanUrls:

    def test_wildcard_ipv4_lists_loopback_and_lan(self):
        urls = _lan_urls("0.0.0.0", PORT, ["192.168.1.5"])
        assert urls == ["http://127.0.0.1:8000",
                        "http://192.168.1.5:8000"]

    def test_wildcard_ipv6_lists_loopback_and_lan(self):
        urls = _lan_urls("::", PORT, ["10.0.0.7"])
        assert urls == ["http://127.0.0.1:8000", "http://10.0.0.7:8000"]

    def test_loopback_ip_lists_only_loopback(self):
        assert _lan_urls("127.0.0.1", PORT, ["192.168.1.5"]) == \
            ["http://127.0.0.1:8000"]

    def test_localhost_lists_only_loopback(self):
        assert _lan_urls("localhost", PORT, ["192.168.1.5"]) == \
            ["http://127.0.0.1:8000"]

    def test_ipv6_loopback_lists_only_loopback(self):
        assert _lan_urls("::1", PORT, []) == ["http://127.0.0.1:8000"]

    def test_specific_lan_host_lists_itself(self):
        assert _lan_urls("10.0.0.7", 9000, ["192.168.1.5"]) == \
            ["http://10.0.0.7:9000"]

    def test_urls_deduplicated_preserving_order(self):
        urls = _lan_urls("0.0.0.0", 8888,
                         ["127.0.0.1", "192.168.1.5", "192.168.1.5"])
        assert urls == ["http://127.0.0.1:8888",
                        "http://192.168.1.5:8888"]


class TestLocalIpv4s:

    def test_returns_list_of_strings(self):
        ips = _local_ipv4s()
        assert isinstance(ips, list)
        assert all(isinstance(ip, str) for ip in ips)

    def test_excludes_loopback_and_unspecified(self):
        for ip in _local_ipv4s():
            addr = ipaddress.IPv4Address(ip)
            assert not addr.is_loopback
            assert not addr.is_unspecified

    def test_unique_and_stable_order(self):
        ips = _local_ipv4s()
        assert ips == sorted(set(ips))


# ============================================================
# 原生目录对话框可用性
# ============================================================

class TestPickerAvailable:

    def test_returns_bool(self):
        assert isinstance(picker_available(), bool)


# ============================================================
# 目录浏览纯函数
# ============================================================

class TestListDirectory:

    def test_lists_subdirs_and_h5_files(self, tmp_path):
        root = tmp_path / "data"
        (root / "beta").mkdir(parents=True)
        (root / "Alpha").mkdir()
        (root / "sub" / "nested").mkdir(parents=True)
        (root / "a.h5").write_bytes(b"x" * 2048)
        (root / "B.H5").write_text("x", encoding="utf-8")
        (root / "note.txt").write_text("x", encoding="utf-8")
        data = _tree(root)
        assert data["path"] == str(root.resolve())
        assert data["parent"] == str(tmp_path.resolve())
        assert data["dirs"] == ["Alpha", "beta", "sub"]
        assert [f["name"] for f in data["files"]] == ["B.H5", "a.h5"]
        assert data["truncated"] is False

    def test_h5_named_directory_is_not_listed_as_file(self, tmp_path):
        root = tmp_path / "data"
        (root / "fake.h5").mkdir(parents=True)
        assert _tree(root)["files"] == []

    def test_empty_dir_returns_no_entries(self, tmp_path):
        root = tmp_path / "empty"
        root.mkdir()
        data = _tree(root)
        assert data["dirs"] == [] and data["files"] == []
        assert data["truncated"] is False

    def test_empty_path_uses_home(self):
        assert fsbrowse.list_directory("")["path"] == str(Path.home().resolve())

    def test_missing_path_raises(self, tmp_path):
        with pytest.raises(ValueError, match="path not found"):
            fsbrowse.list_directory(str(tmp_path / "nope"))

    def test_file_path_raises_not_a_directory(self, tmp_path):
        target = tmp_path / "single.h5"
        target.write_bytes(b"x")
        with pytest.raises(ValueError, match="not a directory"):
            fsbrowse.list_directory(str(target))

    def test_only_h5_files_listed(self, tmp_path):
        root = tmp_path / "h5only"
        root.mkdir()
        for name in ["one.h5", "two.H5", "three.csv", "four.hdf5"]:
            (root / name).write_text("x", encoding="utf-8")
        names = [f["name"] for f in _tree(root)["files"]]
        assert names == ["one.h5", "two.H5"]

    def test_size_kb_reported(self, tmp_path):
        root = tmp_path / "sized"
        root.mkdir()
        (root / "a.h5").write_bytes(b"x" * 2048)
        (root / "b.h5").write_bytes(b"x" * 512)
        sizes = {f["name"]: f["size_kb"] for f in _tree(root)["files"]}
        assert sizes == {"a.h5": 2.0, "b.h5": 0.5}

    def test_unreadable_dir_raises_value_error(self, tmp_path, monkeypatch):
        root = tmp_path / "locked"
        root.mkdir()

        def deny_read(path):
            raise PermissionError(13, "Permission denied")

        monkeypatch.setattr(Path, "iterdir", deny_read)
        with pytest.raises(ValueError, match="cannot read directory"):
            fsbrowse.list_directory(str(root))

    def test_unreadable_size_skips_file(self, tmp_path, monkeypatch):
        root = tmp_path / "sizeskip"
        root.mkdir()
        (root / "ok.h5").write_bytes(b"x")
        (root / "bad.h5").write_bytes(b"x")
        monkeypatch.setattr(fsbrowse, "_size_kb", lambda entry: (
            None if entry.name == "bad.h5" else 0.0))
        names = [f["name"] for f in _tree(root)["files"]]
        assert names == ["ok.h5"]

    def test_limit_keeps_dirs_before_files(self, tmp_path):
        root = tmp_path / "limited"
        for name in ["d1", "d2"]:
            (root / name).mkdir(parents=True)
        for name in ["a.h5", "b.h5"]:
            (root / name).write_bytes(b"x")
        data = fsbrowse.list_directory(str(root), limit=2)
        assert data["dirs"] == ["d1", "d2"]
        assert data["files"] == []
        assert data["truncated"] is True

    def test_limit_cuts_files(self, tmp_path):
        root = tmp_path / "cut"
        for name in ["d1", "d2"]:
            (root / name).mkdir(parents=True)
        for name in ["a.h5", "b.h5"]:
            (root / name).write_bytes(b"x")
        data = fsbrowse.list_directory(str(root), limit=3)
        assert data["dirs"] == ["d1", "d2"]
        assert [f["name"] for f in data["files"]] == ["a.h5"]
        assert data["truncated"] is True

    def test_limit_large_keeps_all_entries(self, tmp_path):
        root = tmp_path / "all"
        (root / "d1").mkdir(parents=True)
        (root / "a.h5").write_bytes(b"x")
        data = fsbrowse.list_directory(str(root), limit=10)
        assert data["dirs"] == ["d1"]
        assert [f["name"] for f in data["files"]] == ["a.h5"]
        assert data["truncated"] is False

    def test_limit_zero_returns_empty(self, tmp_path):
        root = tmp_path / "zero"
        (root / "d1").mkdir(parents=True)
        (root / "a.h5").write_bytes(b"x")
        data = fsbrowse.list_directory(str(root), limit=0)
        assert data["dirs"] == [] and data["files"] == []
        assert data["truncated"] is True

    def test_symlinked_dir_listed_but_not_recursed(self, tmp_path):
        root = tmp_path / "links"
        real = tmp_path / "real"
        (real / "inner").mkdir(parents=True)
        root.mkdir()
        try:
            (root / "linked").symlink_to(real, target_is_directory=True)
        except OSError:
            pytest.skip("symlinks not supported")
        data = _tree(root)
        assert "linked" in data["dirs"]
        assert "inner" not in data["dirs"]


# ============================================================
# token 中间件
# ============================================================

class TestTokenGuard:

    def test_page_without_token_401_html(self):
        resp = _client(token=TOKEN).get("/")
        assert resp.status_code == 401
        assert resp.headers["content-type"].startswith("text/html")

    def test_api_without_token_401_json(self):
        resp = _client(token=TOKEN).get("/api/matrix/hydrogen_bond")
        assert resp.status_code == 401
        assert resp.json() == {"detail": "token required"}

    def test_download_without_token_401_text(self):
        resp = _client(token=TOKEN).get("/download/csv/hydrogen_bond")
        assert resp.status_code == 401
        assert resp.headers["content-type"].startswith("text/plain")
        assert resp.text == "token required"

    def test_static_without_token_401(self):
        resp = _client(token=TOKEN).get("/static/app.css")
        assert resp.status_code == 401

    def test_static_with_query_token_passes(self):
        resp = _client(token=TOKEN).get(
            f"/static/app.css?token={TOKEN}")
        assert resp.status_code == 200

    def test_query_token_grants_access(self):
        assert _client(token=TOKEN).get("/?token=" + TOKEN).status_code == 200

    def test_query_token_sets_httponly_cookie(self):
        resp = _client(token=TOKEN).get("/?token=" + TOKEN)
        cookie = resp.headers["set-cookie"].lower()
        assert "dii_token=" + TOKEN.lower() in cookie
        assert "httponly" in cookie
        assert "samesite=lax" in cookie
        assert "path=/" in cookie

    def test_cookie_grants_access(self):
        resp = _client(token=TOKEN).get("/", cookies={"dii_token": TOKEN})
        assert resp.status_code == 200

    def test_api_passes_with_cookie(self, tmp_path):
        resp = _client(token=TOKEN, browse_mode="server").get(
            "/api/list", params={"path": str(tmp_path)},
            cookies={"dii_token": TOKEN})
        assert resp.status_code == 200

    def test_wrong_query_token_rejected(self):
        assert _client(token=TOKEN).get("/?token=wrong").status_code == 401

    def test_wrong_cookie_token_rejected(self):
        resp = _client(token=TOKEN).get("/", cookies={"dii_token": "wrong"})
        assert resp.status_code == 401

    def test_401_page_does_not_echo_token(self):
        resp = _client(token=TOKEN).get("/")
        assert TOKEN not in resp.text

    def test_default_app_has_no_auth(self):
        client = _client()
        assert client.get("/").status_code == 200
        assert client.get("/api/matrix/hydrogen_bond").status_code == 404

    def test_default_app_pick_unchanged(self, monkeypatch, tmp_path):
        monkeypatch.setattr("DuIvyInteractions.web.app.pick_directory",
                            lambda: str(tmp_path))
        client = _client()
        resp = client.post("/api/pick")
        assert resp.status_code == 200 and resp.json()["path"] == str(tmp_path)


# ============================================================
# /api/list 与 /api/pick 模式
# ============================================================

class TestBrowseModeApi:

    def test_list_disabled_in_picker_mode(self, tmp_path):
        resp = _client().get("/api/list", params={"path": str(tmp_path)})
        assert resp.status_code == 404
        assert resp.json()["detail"] == "server browse disabled"

    def test_list_enabled_in_server_mode(self, tmp_path):
        (tmp_path / "sub").mkdir()
        (tmp_path / "a.h5").write_bytes(b"x")
        data = _client(browse_mode="server").get(
            "/api/list", params={"path": str(tmp_path)}).json()
        assert data["dirs"] == ["sub"]
        assert data["files"][0]["name"] == "a.h5"
        assert data["truncated"] is False

    def test_list_missing_path_404_json(self, tmp_path):
        resp = _client(browse_mode="server").get(
            "/api/list", params={"path": str(tmp_path / "nope")})
        assert resp.status_code == 404
        assert "path not found" in resp.json()["detail"]

    def test_list_file_path_404_json(self, tmp_path):
        target = tmp_path / "single.h5"
        target.write_bytes(b"x")
        resp = _client(browse_mode="server").get(
            "/api/list", params={"path": str(target)})
        assert resp.status_code == 404
        assert "not a directory" in resp.json()["detail"]

    def test_pick_returns_400_in_server_mode(self):
        resp = _client(browse_mode="server").post("/api/pick")
        assert resp.status_code == 400
        data = resp.json()
        assert data["path"] == ""
        assert "native picker unavailable" in data["error"]
