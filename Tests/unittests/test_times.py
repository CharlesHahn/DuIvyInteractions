# -*- coding: utf-8 -*-
"""Interaction.times 字段的单元测试。

覆盖：
1. __post_init__ 校验（shape 不匹配、类型转换）
2. 检测器返回的 times 正确性（shape、单调性、值）
"""

import pytest
import numpy as np
from pathlib import Path

from DuIvyInteractions.core.datas import Interaction, Group, AtomData
from DuIvyInteractions.io.h5 import save_interactions, load_interactions


# ============================================================
# fixtures
# ============================================================

TEMP_DIR = Path(__file__).parent.parent.parent / "test_temp"


@pytest.fixture(autouse=True)
def setup_temp_dir():
    TEMP_DIR.mkdir(exist_ok=True)
    yield
    import shutil
    if TEMP_DIR.exists():
        shutil.rmtree(TEMP_DIR, ignore_errors=True)


def _make_group(gid: int = 0) -> Group:
    """构造最小 Group。"""
    return Group(
        group_id=gid, group_type="H_donor",
        molecule="test", residue_name="RES", residue_id=1,
        atoms=[AtomData(gid * 2, 0, "N", "N", "N", -0.5, 14.0),
               AtomData(gid * 2 + 1, 1, "H", "H", "H", 0.3, 1.0)],
    )


def _make_interaction(n_pairs: int = 1, n_frames: int = 5,
                      times=None) -> Interaction:
    """构造最小 Interaction。"""
    groups = [(_make_group(i), _make_group(i + 100))
              for i in range(n_pairs)]
    existence = np.zeros((n_pairs, n_frames), dtype=bool)
    metrics = {"distance": np.full((n_pairs, n_frames), 3.0)}
    if times is None:
        times = np.arange(n_frames, dtype=float) * 10.0
    return Interaction(
        interaction_type="test", groups=groups,
        existence=existence, metrics=metrics, times=times,
    )


# ============================================================
# 1. __post_init__ 校验
# ============================================================

class TestTimesValidation:
    """Interaction.__post_init__ 对 times 的校验。"""

    def test_correct_shape_accepted(self):
        """正确 shape 不报错。"""
        it = _make_interaction(n_pairs=2, n_frames=10)
        assert it.times.shape == (10,)
        assert it.times.dtype == np.float64

    def test_wrong_shape_raises(self):
        """times shape != n_frames 应报错。"""
        with pytest.raises(ValueError, match="times shape"):
            _make_interaction(n_frames=5, times=np.array([0.0, 1.0]))

    def test_scalar_raises(self):
        """传标量应报错。"""
        with pytest.raises(ValueError, match="times shape"):
            _make_interaction(n_frames=5, times=42.0)

    def test_list_converted_to_float64(self):
        """传 list 应自动转为 float64 ndarray。"""
        it = _make_interaction(n_frames=3, times=[0.0, 10.0, 20.0])
        assert isinstance(it.times, np.ndarray)
        assert it.times.dtype == np.float64
        np.testing.assert_array_equal(it.times, [0.0, 10.0, 20.0])

    def test_int_array_converted_to_float64(self):
        """传 int ndarray 应自动转为 float64。"""
        it = _make_interaction(n_frames=3, times=np.array([0, 10, 20]))
        assert it.times.dtype == np.float64

    def test_empty_interaction(self):
        """n_pairs=0, n_frames=5 → times shape (5,)。"""
        it = _make_interaction(n_pairs=0, n_frames=5)
        assert it.times.shape == (5,)

    def test_zero_frames(self):
        """n_frames=0 → times shape (0,)。"""
        groups = []
        existence = np.empty((0, 0), dtype=bool)
        metrics = {"distance": np.empty((0, 0))}
        it = Interaction(
            interaction_type="test", groups=groups,
            existence=existence, metrics=metrics,
            times=np.array([], dtype=float),
        )
        assert it.times.shape == (0,)


# ============================================================
# 2. h5 往返中的 times
# ============================================================

class TestTimesH5Roundtrip:
    """times 在 h5 保存/加载中无损往返。"""

    def test_roundtrip_preserves_times(self):
        times_orig = np.array([0.0, 100.0, 200.0, 300.0, 400.0])
        it = _make_interaction(n_frames=5, times=times_orig)

        path = str(TEMP_DIR / "times_roundtrip.h5")
        save_interactions([it], path)
        loaded = load_interactions(path)[0]

        np.testing.assert_array_equal(loaded.times, times_orig)
        assert loaded.times.dtype == np.float64


# ============================================================
# 3. 检测器 times 正确性（用真实数据）
# ============================================================

TPR_FILE = Path(__file__).parent.parent / "test_MD_case" / "md.tpr"
XTC_FILE = Path(__file__).parent.parent / "test_MD_case" / "md1ns.xtc"


@pytest.fixture(scope="module")
def saltbridge_result():
    """运行盐桥 TwoPass 检测器，返回 Interaction 列表。"""
    import MDAnalysis as mda
    from DuIvyInteractions.system_readers import GmxTprReader
    from DuIvyInteractions.group_identifiers import AmberFFGroupIdentifier
    from DuIvyInteractions.interaction_detectors import SaltBridgeDetectorTwoPass

    sd = GmxTprReader().read(str(TPR_FILE))
    groups = AmberFFGroupIdentifier().identify(sd)
    filtered = [g for g in groups
                if g.group_type in ("charged_positive", "charged_negative")]
    u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    return SaltBridgeDetectorTwoPass().detect(filtered, trajectory=u.trajectory)


@pytest.fixture(scope="module")
def expected_n_frames():
    import MDAnalysis as mda
    u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    return u.trajectory.n_frames


class TestDetectorTimes:
    """检测器返回的 times 正确性。"""

    def test_times_shape(self, saltbridge_result, expected_n_frames):
        """times shape 应为 (n_frames,)。"""
        it = saltbridge_result[0]
        assert it.times.shape == (expected_n_frames,)

    def test_times_dtype(self, saltbridge_result):
        """times dtype 应为 float64。"""
        it = saltbridge_result[0]
        assert it.times.dtype == np.float64

    def test_times_monotonically_increasing(self, saltbridge_result):
        """times 应单调递增（物理时间）。"""
        it = saltbridge_result[0]
        diffs = np.diff(it.times)
        assert np.all(diffs > 0), f"非单调递增: diffs={diffs[diffs <= 0]}"

    def test_times_first_frame_is_zero(self, saltbridge_result):
        """第一帧时间应为 0。"""
        it = saltbridge_result[0]
        assert it.times[0] == 0.0

    def test_times_unit_ps(self, saltbridge_result, expected_n_frames):
        """时间单位应为 ps（101 帧 xtc 的最后一帧 = 1000 ps）。"""
        it = saltbridge_result[0]
        assert it.times[-1] == pytest.approx(1000.0, abs=1.0)

    def test_times_consistent_across_pairs(self, saltbridge_result):
        """同一 Interaction 内 times 只有一份，不随 pair 变化。"""
        it = saltbridge_result[0]
        assert it.times.ndim == 1
        assert it.times.shape[0] == it.n_frames
