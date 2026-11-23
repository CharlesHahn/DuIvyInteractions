# -*- coding: utf-8 -*-
"""原生目录选择对话框：子进程运行 tkinter，返回所选目录路径。"""

import os
import subprocess
import sys

DIALOG_TIMEOUT = 300  # 秒：用户长期不操作时避免请求永久挂起

# 子进程脚本：选择目录；取消时 stdout 为空串
_DIALOG_SCRIPT = r"""
import sys
import tkinter as tk
from tkinter import filedialog

root = tk.Tk()
root.withdraw()
root.attributes("-topmost", True)
path = filedialog.askdirectory(title="Select a folder containing *.h5", mustexist=True)
root.destroy()
sys.stdout.write(path or "")
"""


def picker_available() -> bool:
    """判断当前主机是否能弹出原生目录对话框。"""
    try:
        import tkinter  # noqa: F401
        if sys.platform.startswith("linux"):
            return bool(os.environ.get("DISPLAY")
                        or os.environ.get("WAYLAND_DISPLAY"))
        return True
    except Exception:
        return False


def pick_directory() -> str:
    """弹原生目录对话框并返回所选路径；用户取消返回空串。"""
    proc = _run_dialog()
    if proc.returncode != 0:
        raise RuntimeError(_last_error_line(proc.stderr))
    return proc.stdout.strip()


def _run_dialog() -> subprocess.CompletedProcess:
    """在子进程里执行目录对话框（UTF-8 管道，避免中文路径乱码）。"""
    env = dict(os.environ, PYTHONIOENCODING="utf-8")
    try:
        return subprocess.run(
            [sys.executable, "-c", _DIALOG_SCRIPT],
            capture_output=True, encoding="utf-8", env=env,
            timeout=DIALOG_TIMEOUT,
        )
    except subprocess.TimeoutExpired:
        raise RuntimeError("dialog timed out (no selection within timeout)")


def _last_error_line(stderr: str) -> str:
    """子进程 stderr 的最后一行作为错误信息（无 GUI / 无 tkinter 时）。"""
    lines = [line.strip() for line in (stderr or "").splitlines() if line.strip()]
    return lines[-1] if lines else "dialog process failed"
