# -*- coding: utf-8 -*-
"""计算任务层：表单归一、CLI 参数生成、后台 subprocess 与日志缓冲。"""

import secrets
import shlex
import subprocess
import sys
import threading
import time
from collections import deque
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Deque, Dict, List, Optional, Tuple

from .fsbrowse import H5_SUFFIX, confine

LOG_LIMIT = 2000  # 每个任务在内存中保留的日志行上限
MAX_JOBS = 5  # 最近任务保留数
CANCEL_GRACE_SECONDS = 3.0  # terminate 后等待进程退出的宽限期（秒）

FORCE_FIELDS = ("amber", "gromos", "charmm", "opls")
STRATEGIES = ("two_pass", "per_frame", "per_tuple")
PAIRINGS = ("any", "between")
INTERACTIONS = ("halogen_bond", "hydrogen_bond", "hydrophobic",
                "metal_coordination", "pi_cation", "pi_stacking",
                "salt_bridge", "water_bridge")
DEFAULT_FF = "amber"
DEFAULT_STRATEGY = "two_pass"
ANY_PAIRING = "any"
HYDROPHOBIC = "hydrophobic"
TRUTHY_VALUES = frozenset({"1", "true", "on", "yes"})

STATUS_PENDING = "pending"
STATUS_RUNNING = "running"
STATUS_DONE = "done"
STATUS_FAILED = "failed"
STATUS_CANCELLED = "cancelled"
ACTIVE_STATUSES = frozenset({STATUS_PENDING, STATUS_RUNNING})

CLI_MODULE = "DuIvyInteractions.DII"
JOB_ID_BYTES = 4
NO_INTERACTION = "select at least one interaction"
RUN_IN_PROGRESS = "a run is already in progress"
MARK_REDUNDANT_NEEDS_HYDROPHOBIC = (
    "mark_redundant requires the hydrophobic interaction")
SPAWN_ERROR = "cannot start process: {}"
EXIT_ERROR = "process exited with status {}"
LINE_SEPARATORS = "\r\n"


@dataclass(frozen=True)
class RunParams:
    """一次 `dii run` 的规范化参数。"""

    topology: str
    trajectory: str
    output: str
    ff: str
    strategy: str
    interactions: Tuple[str, ...]
    select: Tuple[str, ...] = ()
    pairing: str = ANY_PAIRING
    mark_redundant: bool = False

    def argv(self) -> List[str]:
        """完整 CLI 参数列表（模块入口 + 全部显式选项）。"""
        parts = [sys.executable, "-m", CLI_MODULE, "run",
                 "-t", self.topology, "-f", self.trajectory, "-o", self.output,
                 "--ff", self.ff, "--strategy", self.strategy,
                 "--interactions", ",".join(self.interactions)]
        parts.extend(_select_argv(self.select))
        if self.pairing != ANY_PAIRING:
            parts.extend(["--pairing", self.pairing])
        if self.mark_redundant:
            parts.append("--mark-redundant")
        return parts

    def preview(self) -> str:
        """shell 引号安全的单行命令（供页面显示）。"""
        return " ".join(shlex.quote(part) for part in self.argv())

    def summary(self) -> str:
        """形如 "amber / two_pass · 8 types" 的一行摘要。"""
        count = len(self.interactions)
        return f"{self.ff} / {self.strategy} · {count} {_type_noun(count)}"


def parse_run_form(fields: Dict[str, List[str]], root: str) -> RunParams:
    """校验并归一提交的表单字段为 RunParams（非法输入 raise ValueError）。"""
    topology = _path_field(fields, "topology", root)
    trajectory = _path_field(fields, "trajectory", root)
    output = _path_field(fields, "output", root)
    ff = _choice(fields, "ff", FORCE_FIELDS, DEFAULT_FF)
    strategy = _choice(fields, "strategy", STRATEGIES, DEFAULT_STRATEGY)
    pairing = _choice(fields, "pairing", PAIRINGS, ANY_PAIRING)
    interactions = _interactions(fields)
    mark_redundant = _flag(fields, "mark_redundant")
    if mark_redundant and HYDROPHOBIC not in interactions:
        raise ValueError(MARK_REDUNDANT_NEEDS_HYDROPHOBIC)
    return RunParams(topology=topology, trajectory=trajectory, output=output,
                     ff=ff, strategy=strategy, interactions=interactions,
                     select=_selects(fields), pairing=pairing,
                     mark_redundant=mark_redundant)


class Job:
    """一次后台计算的进度、日志与产出。"""

    def __init__(self, job_id: str, params: RunParams,
                 now: Callable[[], float] = time.monotonic) -> None:
        self.id = job_id
        self.params = params
        self.status = STATUS_PENDING
        self.started_at = now()
        self.ended_at: Optional[float] = None
        self.returncode: Optional[int] = None
        self.output_dir = params.output
        self.error: Optional[str] = None
        self.process: Optional[Any] = None
        self._now = now
        self._lines: Deque[str] = deque(maxlen=LOG_LIMIT)
        self._total = 0
        self._lock = threading.Lock()

    def lines(self) -> List[str]:
        """当前保留的全部日志行（快照）。"""
        with self._lock:
            return list(self._lines)

    def lines_since(self, index: int) -> Tuple[List[str], int]:
        """返回 index 之后的日志增量与下一个游标（越界索引视为 0）。"""
        with self._lock:
            return self._slice_locked(index), self._total

    def line_total(self) -> int:
        """累计写入的日志行总数（与 lines_since 的游标同口径）。"""
        with self._lock:
            return self._total

    def elapsed(self) -> float:
        """从创建/启动到结束（或当前）的秒数。"""
        end = self.ended_at if self.ended_at is not None else self._now()
        return max(end - self.started_at, 0.0)

    def h5_files(self) -> List[str]:
        """任务结束后产出目录下的 *.h5 文件名（按名排序）。"""
        if self.status in ACTIVE_STATUSES:
            return []
        return _scan_h5(self.output_dir)

    def append_line(self, line: str) -> None:
        """追加一行日志并推进总游标。"""
        with self._lock:
            self._lines.append(line)
            self._total += 1

    def mark_running(self) -> None:
        """标记任务进入运行态。"""
        with self._lock:
            self.status = STATUS_RUNNING

    def fail_start(self, error: str) -> None:
        """标记任务启动失败（进程未创建）。"""
        with self._lock:
            self.status = STATUS_FAILED
            self.error = error
            self.ended_at = self._now()

    def try_cancel(self) -> bool:
        """把运行中任务切到 cancelled；非运行态返回 False。"""
        with self._lock:
            if self.status not in ACTIVE_STATUSES or self.process is None:
                return False
            self.status = STATUS_CANCELLED
            self.ended_at = self._now()
            return True

    def finish(self, returncode: int) -> None:
        """记录进程退出码与终态（已取消的任务保持 cancelled）。"""
        with self._lock:
            self.returncode = returncode
            if self.ended_at is None:
                self.ended_at = self._now()
            if self.status == STATUS_CANCELLED:
                return
            self.status = STATUS_DONE if returncode == 0 else STATUS_FAILED
            if returncode != 0:
                self.error = EXIT_ERROR.format(returncode)

    def _slice_locked(self, index: int) -> List[str]:
        """返回 index 之后的保留日志（调用方需持有 _lock）。"""
        retained_start = self._total - len(self._lines)
        if index < 0 or index > self._total:
            index = 0
        return list(self._lines)[max(index, retained_start) - retained_start:]


class JobManager:
    """进程内任务表：同时只跑一个计算，保留最近 MAX_JOBS 个任务。"""

    def __init__(self, root: str = "",
                 spawn: Optional[Callable[..., Any]] = None,
                 now: Callable[[], float] = time.monotonic) -> None:
        self.root = root
        self._spawn = spawn or subprocess.Popen
        self._now = now
        self._lock = threading.Lock()
        self._jobs: Dict[str, Job] = {}

    def start(self, params: RunParams) -> Job:
        """创建并启动任务（已有运行中任务 → ValueError）。"""
        with self._lock:
            if any(job.status in ACTIVE_STATUSES
                   for job in self._jobs.values()):
                raise ValueError(RUN_IN_PROGRESS)
            job = Job(secrets.token_hex(JOB_ID_BYTES), params, self._now)
            self._jobs[job.id] = job
            self._trim_locked()
        return self._launch(job)

    def get(self, job_id: str) -> Optional[Job]:
        """按 id 取任务（未知返回 None）。"""
        with self._lock:
            return self._jobs.get(job_id)

    def recent(self, limit: int = MAX_JOBS) -> List[Job]:
        """最近任务（新 → 旧，最多 limit 个）。"""
        with self._lock:
            jobs = list(self._jobs.values())
        return list(reversed(jobs))[:max(limit, 0)]

    def cancel(self, job_id: str) -> bool:
        """终止运行中任务（terminate，宽限超时后 kill）；不可取消返回 False。"""
        job = self.get(job_id)
        if job is None or not job.try_cancel():
            return False
        _terminate(job.process)
        return True

    def shutdown(self) -> None:
        """终止全部运行中任务（应用退出时调用）。"""
        with self._lock:
            active = [job.id for job in self._jobs.values()
                      if job.status in ACTIVE_STATUSES]
        for job_id in active:
            self.cancel(job_id)

    def _launch(self, job: Job) -> Job:
        """spawn 子进程并启动日志读取线程。"""
        try:
            job.process = self._spawn(job.params.argv(), stdout=subprocess.PIPE,
                                      stderr=subprocess.STDOUT, text=True,
                                      bufsize=1, cwd=None)
        except OSError as error:
            job.fail_start(SPAWN_ERROR.format(error))
            return job
        job.mark_running()
        threading.Thread(target=self._watch, args=(job,), daemon=True).start()
        return job

    def _watch(self, job: Job) -> None:
        """逐行读日志，EOF 后等待退出码并写终态。"""
        process = job.process
        try:
            for line in process.stdout or ():
                job.append_line(line.rstrip(LINE_SEPARATORS))
        finally:
            job.finish(process.wait())

    def _trim_locked(self) -> None:
        """超出 MAX_JOBS 时淘汰最旧的终态任务（调用方需持有 _lock）。"""
        while len(self._jobs) > MAX_JOBS:
            stale = next((key for key, job in self._jobs.items()
                          if job.status not in ACTIVE_STATUSES), None)
            if stale is None:
                return
            del self._jobs[stale]


def _path_field(fields: Dict[str, List[str]], key: str, root: str) -> str:
    """必填路径字段：非空且经 root 边界归一。"""
    raw = _first(fields, key)
    if not raw:
        raise ValueError(f"{key} is required")
    return str(confine(root, raw))


def _choice(fields: Dict[str, List[str]], key: str, allowed: Tuple[str, ...],
            default: str) -> str:
    """白名单字段（缺失用默认值，大小写归一）。"""
    value = (_first(fields, key) or default).lower()
    if value not in allowed:
        raise ValueError(f"unsupported {key}: {value}")
    return value


def _interactions(fields: Dict[str, List[str]]) -> Tuple[str, ...]:
    """相互作用类型：支持重复键/逗号分隔，去重保序，空则报错。"""
    names: List[str] = []
    for raw in fields.get("interactions", []):
        names.extend(_split_names(raw))
    unique = tuple(dict.fromkeys(names))
    if not unique:
        raise ValueError(NO_INTERACTION)
    return unique


def _split_names(raw: str) -> List[str]:
    """单值按逗号拆分并归一（小写、去空白、丢空项）。"""
    return [part.strip().lower() for part in (raw or "").split(",")
            if part.strip()]


def _selects(fields: Dict[str, List[str]]) -> Tuple[str, ...]:
    """选择表达式：多行文本按行拆分、去空白、丢空行。"""
    lines: List[str] = []
    for raw in fields.get("select", []):
        lines.extend(part.strip() for part in (raw or "").splitlines())
    return tuple(line for line in lines if line)


def _flag(fields: Dict[str, List[str]], key: str) -> bool:
    """复选框字段（任一真值即视为勾选）。"""
    return any(value.strip().lower() in TRUTHY_VALUES
               for value in fields.get(key, []))


def _first(fields: Dict[str, List[str]], key: str) -> str:
    """字段首个值（缺失/空值列表返回空串）。"""
    values = fields.get(key) or [""]
    return (values[0] or "").strip()


def _select_argv(selects: Tuple[str, ...]) -> List[str]:
    """每个选择表达式生成一对 --select 参数。"""
    parts: List[str] = []
    for selection in selects:
        parts.extend(["--select", selection])
    return parts


def _type_noun(count: int) -> str:
    """类型计数的单复数。"""
    return "type" if count == 1 else "types"


def _scan_h5(directory: str) -> List[str]:
    """目录下 *.h5 文件名（按名排序；目录不可读返回空列表）。"""
    try:
        entries = sorted(Path(directory).iterdir(), key=lambda item: item.name)
    except OSError:
        return []
    return [entry.name for entry in entries
            if entry.is_file() and entry.suffix.lower() == H5_SUFFIX]


def _terminate(process: Optional[Any]) -> None:
    """terminate 子进程，宽限期内未退出则 kill。"""
    if process is None:
        return
    try:
        process.terminate()
        process.wait(timeout=CANCEL_GRACE_SECONDS)
    except subprocess.TimeoutExpired:
        process.kill()
        process.wait()
    except OSError:
        pass
