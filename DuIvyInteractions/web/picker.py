# -*- coding: utf-8 -*-
"""原生目录/文件选择对话框：子进程运行 tkinter，返回所选路径。"""

import os
import subprocess
import sys

DIALOG_TIMEOUT = 300  # 秒：用户长期不操作时避免请求永久挂起

# kind → tkinter filetypes；未知 kind 直接报错（不静默退回全部文件）
FILE_TYPES = {
    "structure": [
        ("MD files",
         "*.tpr *.xtc *.trr *.gro *.pdb *.dcd *.nc *.itp *.psf *.prmtop "
         "*.parm7 *.mol2"),
        ("All files", "*"),
    ],
}

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

# 子进程脚本：选择文件；filetypes 由调用方按 kind 注入（repr 后为合法字面量）
_FILE_DIALOG_SCRIPT = r"""
import sys
import tkinter as tk
from tkinter import filedialog

root = tk.Tk()
root.withdraw()
root.attributes("-topmost", True)
path = filedialog.askopenfilename(title="Select a topology or trajectory",
                                  filetypes={filetypes})
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


def pick_file(kind: str = "structure") -> str:
    """弹原生文件对话框并返回所选路径；用户取消返回空串。"""
    proc = _run_dialog(_file_script(kind))
    if proc.returncode != 0:
        raise RuntimeError(_last_error_line(proc.stderr))
    return proc.stdout.strip()


def _file_script(kind: str) -> str:
    """把 kind 对应的 filetypes 注入文件对话框脚本。"""
    if kind not in FILE_TYPES:
        raise ValueError(f"unsupported file kind: {kind}")
    return _FILE_DIALOG_SCRIPT.format(filetypes=repr(FILE_TYPES[kind]))


def _run_dialog(script: str = _DIALOG_SCRIPT) -> subprocess.CompletedProcess:
    """在子进程里执行选择对话框（UTF-8 管道，避免中文路径乱码）。"""
    env = dict(os.environ, PYTHONIOENCODING="utf-8")
    try:
        return subprocess.run(
            [sys.executable, "-c", script],
            capture_output=True, encoding="utf-8", env=env,
            timeout=DIALOG_TIMEOUT,
        )
    except subprocess.TimeoutExpired:
        raise RuntimeError("dialog timed out (no selection within timeout)")


def _last_error_line(stderr: str) -> str:
    """子进程 stderr 的最后一行作为错误信息（无 GUI / 无 tkinter 时）。"""
    lines = [line.strip() for line in (stderr or "").splitlines() if line.strip()]
    return lines[-1] if lines else "dialog process failed"
