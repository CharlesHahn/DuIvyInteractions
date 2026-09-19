# -*- coding: utf-8 -*-
"""load_interactions 错误路径单元测试：损坏/截断/缺失/版本不匹配等。

补充 test_io_h5.py 未覆盖的 load 侧错误行为（该文件已覆盖 roundtrip
与 path 校验，本文件专注损坏与格式错误文件）。
"""

import shutil
from pathlib import Path

import h5py
import numpy as np
import pytest

from DuIvyInteractions.io.h5 import save_interactions, load_interactions
from DuIvyInteractions.core.datas import Interaction


TMP_DIR = Path(__file__).parent.parent.parent / "test_temp" / "h5_errors_test"


@pytest.fixture(scope="module", autouse=True)
def _tmp_cleanup():
    if TMP_DIR.exists():
        shutil.rmtree(TMP_DIR, ignore_errors=True)
    TMP_DIR.mkdir(parents=True)
    yield
    shutil.rmtree(TMP_DIR, ignore_errors=True)


def _make_interaction() -> Interaction:
    """空对盐桥 Interaction（合法：groups=[] 对应 existence (0, 5)）。"""
    return Interaction(
        interaction_type="salt_bridge", groups=[],
        existence=np.zeros((0, 5), dtype=bool),
        metrics={"distance": np.zeros((0, 5))},
        times=np.arange(5, dtype=float) * 10.0,
    )


class TestCorruptedFiles:
    """损坏/非 h5 文件 → OSError。"""

    def test_garbage_text_raises_oserror(self):
        p = TMP_DIR / "garbage.h5"
        p.write_text("this is not h5", encoding="utf-8")
        with pytest.raises(OSError):
            load_interactions(str(p))

    def test_empty_file_raises_oserror(self):
        p = TMP_DIR / "empty_bytes.h5"
        p.write_bytes(b"")
        with pytest.raises(OSError):
            load_interactions(str(p))

    def test_truncated_h5_raises_oserror(self):
        """合法 h5 截断一半 → h5py 打开/读取失败（OSError）。"""
        p = TMP_DIR / "truncated.h5"
        save_interactions([_make_interaction()], str(p))
        data = p.read_bytes()
        p.write_bytes(data[: len(data) // 2])
        with pytest.raises(OSError):
            load_interactions(str(p))


class TestMissingPath:
    """不存在路径 / 目录 → OSError（FileNotFoundError/PermissionError 是其子类）。"""

    def test_missing_file_raises_filenotfound(self):
        with pytest.raises(FileNotFoundError):
            load_interactions(str(TMP_DIR / "nofile.h5"))

    def test_directory_raises_oserror(self):
        with pytest.raises(OSError):
            load_interactions(str(TMP_DIR))


class TestFormatVersion:
    """格式版本不匹配 → ValueError。"""

    def test_version_mismatch_raises(self):
        p = TMP_DIR / "version.h5"
        save_interactions([_make_interaction()], str(p))
        with h5py.File(p, "a") as f:
            f.attrs["format_version"] = "9.9"
        with pytest.raises(ValueError, match="Unsupported format version"):
            load_interactions(str(p))

    def test_empty_h5_no_attrs_raises(self):
        """h5 文件但无 format_version 属性 → unknown 版本 → ValueError。"""
        p = TMP_DIR / "no_attrs.h5"
        with h5py.File(p, "w"):
            pass
        with pytest.raises(ValueError, match="Unsupported format version"):
            load_interactions(str(p))

    def test_missing_n_interactions_raises(self):
        """有版本号但无 n_interactions 属性 → KeyError。"""
        p = TMP_DIR / "no_n.h5"
        with h5py.File(p, "w") as f:
            f.attrs["format_version"] = "1.0"
        with pytest.raises(KeyError):
            load_interactions(str(p))


class TestLoadPathValidation:
    """load 侧的路径参数校验（save 侧已有 test_io_h5 覆盖）。"""

    def test_non_str_path_raises_typeerror(self):
        with pytest.raises(TypeError):
            load_interactions(123)

    def test_empty_string_raises_valueerror(self):
        with pytest.raises(ValueError):
            load_interactions("")

    def test_pathlib_accepted(self):
        """Path 对象合法：不存在的 Path 应抛 FileNotFoundError 而非 TypeError。"""
        with pytest.raises(FileNotFoundError):
            load_interactions(Path(TMP_DIR) / "nofile.h5")
