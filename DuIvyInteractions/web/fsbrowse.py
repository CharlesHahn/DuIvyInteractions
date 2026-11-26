# -*- coding: utf-8 -*-
"""服务端目录浏览的纯函数层：列出直接子目录与 *.h5 文件。"""

from pathlib import Path
from typing import Dict, List, Optional

LIST_LIMIT = 500  # 单次返回的条目上限
H5_SUFFIX = ".h5"
BYTES_PER_KB = 1024
SIZE_DECIMALS = 1
ROOT_ERROR = "path outside root"
HIDDEN_PREFIX = "."


def confine(root: str, path: str) -> Path:
    """把 path 限制在 root 内并解析为绝对路径（root 为空时不做限制）。"""
    raw = (path or "").strip()
    if not root:
        return Path(raw or Path.home()).expanduser().resolve()
    base = _absolute(root)
    target = Path(raw).expanduser() if raw else base
    resolved = (target if target.is_absolute() else base / target).resolve()
    if resolved != base and not resolved.is_relative_to(base):
        raise ValueError(ROOT_ERROR)
    return resolved


def list_directory(path: str, limit: int = LIST_LIMIT, root: str = "") -> dict:
    """列出目录的直接子目录与 *.h5 文件，条目总数受 limit 限制。"""
    directory = _checked(confine(root, path), path)
    base = _base(root)
    try:
        dirs, files = _subdirs(directory, base), _h5_files(directory, base)
    except OSError:
        raise ValueError(f"cannot read directory: {path}") from None
    kept_dirs = dirs[:max(limit, 0)]
    kept_files = files[:max(limit - len(kept_dirs), 0)]
    return {
        "path": str(directory),
        "parent": _parent(directory, base),
        "root": str(base) if base is not None else "",
        "rel": _rel(directory, base),
        "dirs": kept_dirs,
        "files": kept_files,
        "truncated": len(dirs) + len(files) > limit,
    }


def _absolute(root: str) -> Path:
    """root 的绝对路径。"""
    return Path(root).expanduser().resolve()


def _base(root: str) -> Optional[Path]:
    """浏览边界；root 为空时返回 None（不限制）。"""
    return _absolute(root) if root else None


def _checked(directory: Path, raw: str) -> Path:
    """校验路径存在且为目录（错误文案保留原始 path）。"""
    if not directory.exists():
        raise ValueError(f"path not found: {raw}")
    if not directory.is_dir():
        raise ValueError(f"not a directory: {raw}")
    return directory


def _parent(directory: Path, base: Optional[Path]) -> Optional[str]:
    """父目录路径；到达 root 边界或文件系统根目录时返回 None。"""
    if base is not None:
        return None if directory == base else str(directory.parent)
    return str(directory.parent) if directory.parent != directory else None


def _rel(directory: Path, base: Optional[Path]) -> str:
    """相对浏览边界的 POSIX 路径（边界自身为 "."；无边界时恒为 "."）。"""
    return directory.relative_to(base).as_posix() if base is not None else "."


def _subdirs(directory: Path, base: Optional[Path]) -> List[str]:
    """直接子目录名（跳过隐藏目录与越界符号链接，按小写名排序）。"""
    names = [entry.name for entry in directory.iterdir()
             if _visible_dir(entry) and _inside_base(entry, base)]
    return sorted(names, key=str.lower)


def _visible_dir(entry: Path) -> bool:
    """是否为非隐藏子目录（跟随符号链接判断）。"""
    return not entry.name.startswith(HIDDEN_PREFIX) and entry.is_dir()


def _inside_base(entry: Path, base: Optional[Path]) -> bool:
    """符号链接条目解析后是否仍在浏览边界内（无边界时直接放行）。"""
    if base is None:
        return True
    try:
        resolved = entry.resolve()
    except OSError:
        return False
    return resolved == base or resolved.is_relative_to(base)


def _h5_files(directory: Path,
              base: Optional[Path]) -> List[Dict[str, object]]:
    """直接 *.h5 文件（跳过越界符号链接，按名字排序；大小读取失败的文件跳过）。"""
    files: List[Dict[str, object]] = []
    for entry in sorted(directory.iterdir(), key=lambda item: item.name):
        if entry.suffix.lower() != H5_SUFFIX or not entry.is_file():
            continue
        if not _inside_base(entry, base):
            continue
        size_kb = _size_kb(entry)
        if size_kb is not None:
            files.append({"name": entry.name, "size_kb": size_kb})
    return files


def _size_kb(entry: Path) -> Optional[float]:
    """文件大小（KB，保留一位小数）；读取失败返回 None。"""
    try:
        return round(entry.stat().st_size / BYTES_PER_KB, SIZE_DECIMALS)
    except OSError:
        return None
