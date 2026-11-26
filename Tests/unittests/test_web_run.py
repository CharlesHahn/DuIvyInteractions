# -*- coding: utf-8 -*-
"""W3 计算发起测试：表单归一、CLI 参数、JobManager、JSON 路由（不跑真实 dii）。"""

import shlex
import subprocess
import threading
import time
from pathlib import Path
from urllib.parse import unquote

import pytest

from DuIvyInteractions.web import app as web_app
from DuIvyInteractions.web import jobs
from DuIvyInteractions.web.jobs import RunParams, parse_run_form

TOKEN = "s3cr3t-token"
FORM = {"topology": "md.tpr", "trajectory": "md.xtc", "output": "out",
        "ff": "amber", "strategy": "two_pass",
        "interactions": list(jobs.INTERACTIONS)}
WAIT_TIMEOUT = 3.0
POLL_INTERVAL = 0.005


# ============================================================
# 测试替身
# ============================================================

class FakeStream:
    """模拟 Popen.stdout：逐行产出，可选在 EOF 前阻塞等待 release。"""

    def __init__(self, lines, release):
        self._lines = list(lines)
        self._release = release

    def __iter__(self):
        for line in self._lines:
            yield line
        self._release.wait(WAIT_TIMEOUT)


class FakeProcess:
    """模拟 subprocess.Popen：立即退出或阻塞到 terminate/kill。"""

    def __init__(self, lines=(), returncode=0, hold=False, terminate_works=True):
        self.release = threading.Event()
        if not hold:
            self.release.set()
        self.stdout = FakeStream(lines, self.release)
        self.returncode = None
        self.terminated = False
        self.killed = False
        self._exit_code = returncode
        self._terminate_works = terminate_works

    def terminate(self):
        self.terminated = True
        self._exit_code = -15
        if self._terminate_works:
            self.release.set()

    def kill(self):
        self.killed = True
        self._exit_code = -9
        self.release.set()

    def wait(self, timeout=None):
        if not self.release.wait(timeout):
            raise subprocess.TimeoutExpired("fake", timeout)
        if self.returncode is None:
            self.returncode = self._exit_code
        return self.returncode


def _spawn_of(fake):
    """返回始终给出同一假进程的 spawn 可注入函数。"""

    def spawn(argv, **kwargs):
        return fake

    return spawn


def _params(**overrides):
    """直接构造 RunParams（argv/summary 测试用）。"""
    values = {"topology": "/data/md.tpr", "trajectory": "/data/md.xtc",
              "output": "/data/out", "ff": "amber",
              "strategy": "two_pass",
              "interactions": ("hydrogen_bond", "pi_stacking")}
    values.update(overrides)
    return RunParams(**values)


def _form(**overrides):
    """合法表单字段（逐个用例覆盖单字段）。"""
    values = {"topology": ["md.tpr"], "trajectory": ["md.xtc"],
              "output": ["out"], "ff": ["amber"],
              "strategy": ["two_pass"],
              "interactions": list(jobs.INTERACTIONS)}
    values.update(overrides)
    return values


def _wait_until(predicate, timeout=WAIT_TIMEOUT):
    """轮询等待断言条件成立。"""
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if predicate():
            return
        time.sleep(POLL_INTERVAL)
    raise AssertionError("condition not reached in time")


def _wait_status(job, expected, timeout=WAIT_TIMEOUT):
    """等待任务进入指定终态（expected 可为字符串或集合）。"""
    statuses = {expected} if isinstance(expected, str) else set(expected)

    def reached():
        return job.status in statuses

    _wait_until(reached, timeout)


def _manager(spawn, now=None):
    """构造注入了假 spawn 的 JobManager。"""
    kwargs = {"spawn": spawn}
    if now is not None:
        kwargs["now"] = now
    return jobs.JobManager(**kwargs)


def _client(tmp_path, spawn, **kwargs):
    """TestClient（把 manager 的 spawn 换成假进程）。"""
    from fastapi.testclient import TestClient
    app = web_app.create_app(root=str(tmp_path), **kwargs)
    app.state.jobs._spawn = spawn
    return TestClient(app)


def _posted_job(client):
    """提交表单并返回 (响应, 任务 id)。"""
    response = client.post("/run", data=FORM, follow_redirects=False)
    return response, response.headers["location"].rsplit("/", 1)[-1]


# ============================================================
# parse_run_form
# ============================================================

class TestParseRunForm:

    def test_normalizes_paths_and_defaults(self, tmp_path):
        params = parse_run_form(_form(), str(tmp_path))
        assert params.topology == str(jobs.confine(str(tmp_path), "md.tpr"))
        assert params.trajectory == str(jobs.confine(str(tmp_path), "md.xtc"))
        assert params.output == str(jobs.confine(str(tmp_path), "out"))
        assert params.ff == "amber"
        assert params.strategy == "two_pass"
        assert params.pairing == "any"
        assert params.select == ()
        assert params.mark_redundant is False
        assert params.interactions == jobs.INTERACTIONS

    def test_missing_topology_rejected(self, tmp_path):
        with pytest.raises(ValueError, match="topology is required"):
            parse_run_form(_form(topology=[]), str(tmp_path))

    def test_missing_trajectory_rejected(self, tmp_path):
        with pytest.raises(ValueError, match="trajectory is required"):
            parse_run_form(_form(trajectory=["   "]), str(tmp_path))

    def test_missing_output_rejected(self, tmp_path):
        with pytest.raises(ValueError, match="output is required"):
            parse_run_form(_form(output=[]), str(tmp_path))

    def test_missing_interactions_rejected(self, tmp_path):
        with pytest.raises(ValueError, match="select at least one"):
            parse_run_form(_form(interactions=[]), str(tmp_path))

    def test_invalid_ff_rejected(self, tmp_path):
        with pytest.raises(ValueError, match="unsupported ff"):
            parse_run_form(_form(ff=["nope"]), str(tmp_path))

    def test_invalid_strategy_rejected(self, tmp_path):
        with pytest.raises(ValueError, match="unsupported strategy"):
            parse_run_form(_form(strategy=["nope"]), str(tmp_path))

    def test_invalid_pairing_rejected(self, tmp_path):
        with pytest.raises(ValueError, match="unsupported pairing"):
            parse_run_form(_form(pairing=["nope"]), str(tmp_path))

    def test_mark_redundant_requires_hydrophobic(self, tmp_path):
        fields = _form(interactions=["pi_stacking"], mark_redundant=["on"])
        with pytest.raises(ValueError, match="hydrophobic"):
            parse_run_form(fields, str(tmp_path))

    def test_mark_redundant_allowed_with_hydrophobic(self, tmp_path):
        params = parse_run_form(_form(mark_redundant=["on"]), str(tmp_path))
        assert params.mark_redundant is True

    def test_select_multiline_normalized(self, tmp_path):
        fields = _form(select=["  segid A\n\n segid B  \n\tsegid C\n"])
        assert parse_run_form(fields, str(tmp_path)).select == \
            ("segid A", "segid B", "segid C")

    def test_duplicate_interactions_deduplicated_in_order(self, tmp_path):
        fields = _form(interactions=["pi_stacking,hydrogen_bond",
                                     "pi_stacking", "HYDROPHOBIC"])
        params = parse_run_form(fields, str(tmp_path))
        assert params.interactions == \
            ("pi_stacking", "hydrogen_bond", "hydrophobic")

    def test_root_escape_rejected(self, tmp_path):
        with pytest.raises(ValueError, match="outside root"):
            parse_run_form(_form(topology=["../md.tpr"]), str(tmp_path))

    def test_empty_root_keeps_absolute_path(self):
        target = str(Path("out").resolve())
        assert parse_run_form(_form(output=["out"]), "").output == target


# ============================================================
# RunParams.argv / preview / summary
# ============================================================

class TestRunParams:

    def test_argv_basic_shape(self):
        assert _params().argv() == [
            jobs.sys.executable, "-m", "DuIvyInteractions.DII", "run",
            "-t", "/data/md.tpr", "-f", "/data/md.xtc", "-o", "/data/out",
            "--ff", "amber", "--strategy", "two_pass",
            "--interactions", "hydrogen_bond,pi_stacking"]

    def test_argv_uses_one_select_flag_per_selection(self):
        argv = _params(select=("segid A", "segid B")).argv()
        assert argv[-4:] == ["--select", "segid A", "--select", "segid B"]

    def test_argv_pairing_omitted_for_any(self):
        assert "--pairing" not in _params(pairing="any").argv()

    def test_argv_pairing_between_appended(self):
        assert _params(pairing="between").argv()[-2:] == \
            ["--pairing", "between"]

    def test_argv_mark_redundant_only_when_set(self):
        assert "--mark-redundant" not in _params().argv()
        assert _params(mark_redundant=True).argv()[-1] == "--mark-redundant"

    def test_preview_is_shell_quoted(self):
        params = _params(topology="/data/my dir/md.tpr")
        preview = params.preview()
        assert "'/data/my dir/md.tpr'" in preview
        assert shlex.split(preview) == params.argv()

    def test_summary_plural(self):
        assert _params().summary() == "amber / two_pass · 2 types"

    def test_summary_singular(self):
        assert _params(interactions=("hydrogen_bond",)).summary() == \
            "amber / two_pass · 1 type"


# ============================================================
# JobManager
# ============================================================

class TestJobManager:

    def test_start_marks_running(self):
        fake = FakeProcess(hold=True)
        manager = _manager(_spawn_of(fake))
        job = manager.start(_params())
        assert job.status == jobs.STATUS_RUNNING
        assert manager.get(job.id) is job
        assert manager.recent() == [job]
        manager.shutdown()

    def test_lines_since_incremental(self):
        fake = FakeProcess(lines=["[INFO] a", "[WARN] b", "error: c"],
                           hold=True)
        manager = _manager(_spawn_of(fake))
        job = manager.start(_params())
        _wait_until(lambda: len(job.lines()) == 3)
        assert job.lines() == ["[INFO] a", "[WARN] b", "error: c"]
        assert job.lines_since(0) == (job.lines(), 3)
        assert job.lines_since(2) == (["error: c"], 3)
        fake.release.set()
        _wait_status(job, jobs.STATUS_DONE)

    def test_lines_since_after_buffer_truncation(self, monkeypatch):
        monkeypatch.setattr(jobs, "LOG_LIMIT", 3)
        fake = FakeProcess(lines=["l1", "l2", "l3", "l4", "l5"], hold=True)
        manager = _manager(_spawn_of(fake))
        job = manager.start(_params())
        _wait_until(lambda: job.lines() == ["l3", "l4", "l5"])
        assert job.line_total() == 5
        assert job.lines_since(0) == (["l3", "l4", "l5"], 5)
        assert job.lines_since(4) == (["l5"], 5)
        assert job.lines_since(job.line_total()) == ([], 5)
        fake.release.set()
        _wait_status(job, jobs.STATUS_DONE)
        assert job.line_total() == 5

    def test_exit_zero_done_and_scans_h5(self, tmp_path):
        output = tmp_path / "out"
        output.mkdir()
        (output / "b.h5").write_bytes(b"x")
        (output / "a.h5").write_bytes(b"x")
        (output / "note.txt").write_text("x")
        fake = FakeProcess(returncode=0)
        manager = _manager(_spawn_of(fake))
        job = manager.start(_params(output=str(output)))
        _wait_status(job, jobs.STATUS_DONE)
        assert job.returncode == 0
        assert job.error is None
        assert job.h5_files() == ["a.h5", "b.h5"]

    def test_nonzero_exit_failed(self):
        fake = FakeProcess(returncode=3)
        manager = _manager(_spawn_of(fake))
        job = manager.start(_params())
        _wait_status(job, jobs.STATUS_FAILED)
        assert job.returncode == 3
        assert "3" in job.error

    def test_cancel_marks_cancelled(self):
        fake = FakeProcess(hold=True)
        manager = _manager(_spawn_of(fake))
        job = manager.start(_params())
        assert manager.cancel(job.id) is True
        assert job.status == jobs.STATUS_CANCELLED
        assert fake.terminated and not fake.killed
        _wait_until(lambda: job.returncode is not None)

    def test_cancel_finished_job_returns_false(self):
        manager = _manager(_spawn_of(FakeProcess(returncode=0)))
        job = manager.start(_params())
        _wait_status(job, jobs.STATUS_DONE)
        assert manager.cancel(job.id) is False
        assert job.status == jobs.STATUS_DONE

    def test_cancel_unknown_returns_false(self):
        assert _manager(_spawn_of(FakeProcess())).cancel("nope") is False

    def test_cancel_escalates_to_kill_after_grace(self, monkeypatch):
        monkeypatch.setattr(jobs, "CANCEL_GRACE_SECONDS", 0.05)
        fake = FakeProcess(hold=True, terminate_works=False)
        manager = _manager(_spawn_of(fake))
        job = manager.start(_params())
        assert manager.cancel(job.id) is True
        assert fake.terminated and fake.killed
        assert job.status == jobs.STATUS_CANCELLED

    def test_start_while_running_rejected(self):
        fake = FakeProcess(hold=True)
        manager = _manager(_spawn_of(fake))
        job = manager.start(_params())
        with pytest.raises(ValueError, match="already in progress"):
            manager.start(_params())
        manager.cancel(job.id)

    def test_recent_keeps_newest_and_caps_at_limit(self):
        manager = _manager(lambda argv, **kwargs: FakeProcess(returncode=0))
        ids = []
        for _ in range(jobs.MAX_JOBS + 2):
            job = manager.start(_params())
            ids.append(job.id)
            _wait_status(job, jobs.STATUS_DONE)
        recent = manager.recent()
        assert len(recent) == jobs.MAX_JOBS
        assert [job.id for job in recent] == \
            list(reversed(ids[-jobs.MAX_JOBS:]))
        assert manager.get(ids[0]) is None
        assert len(manager.recent(2)) == 2

    def test_h5_files_empty_while_running(self, tmp_path):
        output = tmp_path / "out"
        output.mkdir()
        (output / "a.h5").write_bytes(b"x")
        manager = _manager(_spawn_of(FakeProcess(hold=True)))
        job = manager.start(_params(output=str(output)))
        assert job.h5_files() == []
        manager.shutdown()

    def test_shutdown_cancels_running_job(self):
        fake = FakeProcess(hold=True)
        manager = _manager(_spawn_of(fake))
        job = manager.start(_params())
        manager.shutdown()
        assert job.status == jobs.STATUS_CANCELLED
        assert fake.terminated

    def test_spawn_failure_marks_failed(self):
        def spawn(argv, **kwargs):
            raise OSError("no exec")

        manager = _manager(spawn)
        job = manager.start(_params())
        assert job.status == jobs.STATUS_FAILED
        assert "no exec" in job.error
        assert job.returncode is None


# ============================================================
# 路由（JSON API + 表单提交）
# ============================================================

class TestRunRoutes:

    def test_post_run_redirects_to_job(self, tmp_path):
        client = _client(tmp_path, _spawn_of(FakeProcess(returncode=0)))
        response, job_id = _posted_job(client)
        assert response.status_code == 303
        assert response.headers["location"] == f"/run/{job_id}"
        job = client.app.state.jobs.get(job_id)
        assert job is not None
        assert job.params.topology == str(jobs.confine(str(tmp_path), "md.tpr"))
        _wait_status(job, jobs.STATUS_DONE)

    def test_post_run_invalid_form_redirects_with_error(self, tmp_path):
        client = _client(tmp_path, _spawn_of(FakeProcess(hold=True)))
        form = dict(FORM)
        form.pop("interactions")
        response = client.post("/run", data=form, follow_redirects=False)
        assert response.status_code == 303
        location = response.headers["location"]
        assert location.startswith("/run?error=")
        assert "select at least one interaction" in unquote(location)

    def test_post_run_conflict_redirects_with_error(self, tmp_path):
        client = _client(tmp_path, _spawn_of(FakeProcess(hold=True)))
        _posted_job(client)
        second = client.post("/run", data=FORM, follow_redirects=False)
        assert second.status_code == 303
        assert "already in progress" in unquote(second.headers["location"])
        client.app.state.jobs.shutdown()

    def test_api_run_first_round(self, tmp_path):
        lines = ["[INFO] one", "[WARN] two", "error: three"]
        client = _client(tmp_path, _spawn_of(
            FakeProcess(lines=lines, hold=True)))
        _, job_id = _posted_job(client)
        _wait_until(lambda: len(client.app.state.jobs.get(job_id).lines()) == 3)
        first = client.get(f"/api/run/{job_id}").json()
        assert first["id"] == job_id
        assert first["status"] == "running"
        assert first["returncode"] is None
        assert first["lines"] == lines
        assert first["next"] == 3
        assert first["h5_files"] == []
        assert first["error"] is None
        assert isinstance(first["elapsed"], float)
        assert first["output_dir"] == str(tmp_path / "out")
        client.app.state.jobs.shutdown()

    def test_api_run_incremental_then_final(self, tmp_path):
        lines = ["[INFO] one", "[WARN] two", "error: three"]
        fake = FakeProcess(lines=lines, hold=True)
        client = _client(tmp_path, _spawn_of(fake))
        _, job_id = _posted_job(client)
        _wait_until(lambda: len(client.app.state.jobs.get(job_id).lines()) == 3)
        incremental = client.get(f"/api/run/{job_id}",
                                  params={"since": 2}).json()
        assert incremental["lines"] == ["error: three"]
        assert incremental["next"] == 3
        fake.release.set()
        _wait_status(client.app.state.jobs.get(job_id), jobs.STATUS_DONE)
        final = client.get(f"/api/run/{job_id}").json()
        assert final["status"] == "done"
        assert final["returncode"] == 0

    def test_job_view_lines_total_matches_cursor(self, tmp_path):
        lines = [f"[INFO] line {index}" for index in range(350)]
        client = _client(tmp_path, _spawn_of(FakeProcess(lines=lines)))
        _, job_id = _posted_job(client)
        job = client.app.state.jobs.get(job_id)
        _wait_status(job, jobs.STATUS_DONE)
        view = web_app._job_view(job)
        assert view["lines_total"] == 350
        assert view["lines"] == lines[-web_app.RUN_LOG_LINES:]
        replay = client.get(f"/api/run/{job_id}",
                            params={"since": view["lines_total"]}).json()
        assert replay["lines"] == []
        assert replay["next"] == job.line_total() == 350

    def test_api_run_invalid_or_out_of_range_since_reads_from_zero(self, tmp_path):
        lines = ["[INFO] one", "[WARN] two"]
        client = _client(tmp_path, _spawn_of(
            FakeProcess(lines=lines, hold=True)))
        _, job_id = _posted_job(client)
        _wait_until(lambda: len(client.app.state.jobs.get(job_id).lines()) == 2)
        for since in ("abc", "-1", "999"):
            body = client.get(f"/api/run/{job_id}",
                              params={"since": since}).json()
            assert body["lines"] == lines
            assert body["next"] == 2
        client.app.state.jobs.shutdown()

    def test_api_run_unknown_404_json(self, tmp_path):
        client = _client(tmp_path, _spawn_of(FakeProcess()))
        response = client.get("/api/run/missing")
        assert response.status_code == 404
        assert response.json() == {"detail": "job not found"}

    def test_cancel_api_returns_cancelled(self, tmp_path):
        fake = FakeProcess(hold=True)
        client = _client(tmp_path, _spawn_of(fake))
        _, job_id = _posted_job(client)
        response = client.post(f"/api/run/{job_id}/cancel")
        assert response.status_code == 200
        assert response.json() == {"status": "cancelled"}
        assert client.app.state.jobs.get(job_id).status == "cancelled"
        assert fake.terminated

    def test_cancel_api_finished_job_400(self, tmp_path):
        client = _client(tmp_path, _spawn_of(FakeProcess(returncode=0)))
        _, job_id = _posted_job(client)
        _wait_status(client.app.state.jobs.get(job_id), jobs.STATUS_DONE)
        response = client.post(f"/api/run/{job_id}/cancel")
        assert response.status_code == 400
        assert response.json()["detail"] == web_app.CANNOT_CANCEL

    def test_cancel_api_unknown_404(self, tmp_path):
        client = _client(tmp_path, _spawn_of(FakeProcess()))
        assert client.post("/api/run/missing/cancel").status_code == 404

    def test_download_run_log(self, tmp_path):
        fake = FakeProcess(lines=["[INFO] start", "error: boom"],
                           returncode=1)
        client = _client(tmp_path, _spawn_of(fake))
        _, job_id = _posted_job(client)
        _wait_status(client.app.state.jobs.get(job_id), jobs.STATUS_FAILED)
        response = client.get(f"/download/run-log/{job_id}")
        assert response.status_code == 200
        assert response.text == "[INFO] start\nerror: boom\n"
        assert response.headers["content-disposition"] == \
            f'attachment; filename="dii-run-{job_id}.log"'
        assert response.headers["content-type"].startswith("text/plain")

    def test_download_run_log_unknown_404(self, tmp_path):
        client = _client(tmp_path, _spawn_of(FakeProcess()))
        assert client.get("/download/run-log/missing").status_code == 404

    def test_job_page_unknown_404(self, tmp_path):
        client = _client(tmp_path, _spawn_of(FakeProcess()))
        assert client.get("/run/missing").status_code == 404

    def test_token_required_for_run_endpoints(self, tmp_path):
        client = _client(tmp_path, _spawn_of(FakeProcess()), token=TOKEN)
        api = client.get("/api/run/any")
        assert api.status_code == 401
        assert api.json() == {"detail": "token required"}
        assert client.post("/run", data=FORM).status_code == 401
        assert client.get("/download/run-log/any").status_code == 401

    def test_lifespan_shutdown_cancels_running_job(self, tmp_path):
        from fastapi.testclient import TestClient
        fake = FakeProcess(hold=True)
        app = web_app.create_app(root=str(tmp_path))
        app.state.jobs._spawn = _spawn_of(fake)
        with TestClient(app) as client:
            _, job_id = _posted_job(client)
            job = app.state.jobs.get(job_id)
            assert job.status == "running"
        assert job.status == "cancelled"
        assert fake.terminated
