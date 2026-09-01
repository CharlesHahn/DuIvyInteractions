# -*- coding: utf-8 -*-
"""盐桥 TwoPass 检测器集成测试。

验证策略三（两轮遍历 + 稀疏存储）的两个接口：
1. run_pass1() — 稀疏结果
2. run_pass2() — Pass2 补全全帧 metric
3. detect() — Pass1 + Pass2 便捷方法
"""

import pytest
from pathlib import Path

import numpy as np
import MDAnalysis as mda

from DuIvyInteractions.input_readers import GmxTprReader
from DuIvyInteractions.group_identifiers import AmberFFGroupIdentifier
from DuIvyInteractions.core.datas import InteractionSparse
from DuIvyInteractions.interaction_detectors import (
    SaltBridgeDetectorTwoPass, SaltBridgeDetectorPerFrame)


TPR_FILE = Path(__file__).parent.parent / "test_MD_case" / "md.tpr"
XTC_FILE = Path(__file__).parent.parent / "test_MD_case" / "md1ns.xtc"


@pytest.fixture(scope="module")
def system_data():
    reader = GmxTprReader()
    return reader.read(str(TPR_FILE))


@pytest.fixture(scope="module")
def groups(system_data):
    return AmberFFGroupIdentifier().identify(system_data)


@pytest.fixture(scope="module")
def filtered_groups(groups):
    return [g for g in groups
            if g.group_type in ("charged_positive", "charged_negative")]


@pytest.fixture(scope="module")
def detector():
    return SaltBridgeDetectorTwoPass()


@pytest.fixture(scope="module")
def sparse_result(detector, filtered_groups):
    """Pass1 only：稀疏结果。"""
    u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    return detector.run_pass1(filtered_groups, u.trajectory), detector


@pytest.fixture(scope="module")
def salt_bridges(detector, filtered_groups):
    """Pass1 + Pass2：完整结果。"""
    u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    return detector.detect(filtered_groups, trajectory=u.trajectory)


@pytest.fixture(scope="module")
def n_frames():
    u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    return u.trajectory.n_frames


# ============================================================
# Pass1 only：稀疏结果
# ============================================================

class TestPass1Only:
    """验证 run_pass1 返回的稀疏结果。"""

    def test_returns_interaction_sparse(self, sparse_result):
        sparse, _ = sparse_result
        assert isinstance(sparse, InteractionSparse)

    def test_non_empty(self, sparse_result):
        sparse, _ = sparse_result
        assert sparse.n_pairs > 0

    def test_pair_count(self, sparse_result):
        """应发现 47 个 unique pair。"""
        sparse, _ = sparse_result
        assert sparse.n_pairs == 47

    def test_each_pair_has_groups(self, sparse_result):
        """每个 pair 应有 groups 字段（正电, 负电）。"""
        sparse, _ = sparse_result
        for group_ids, data in sparse.data.items():
            assert "groups" in data
            g1, g2 = data["groups"]
            assert g1.group_type == "charged_positive"
            assert g2.group_type == "charged_negative"

    def test_each_pair_has_frames(self, sparse_result):
        """每个 pair 应有 frames 列表（非空）。"""
        sparse, _ = sparse_result
        for group_ids, data in sparse.data.items():
            assert "frames" in data
            assert len(data["frames"]) > 0

    def test_each_pair_has_metrics(self, sparse_result):
        """每个 pair 应有 metrics 字典，包含 distance。"""
        sparse, _ = sparse_result
        for group_ids, data in sparse.data.items():
            assert "metrics" in data
            assert "distance" in data["metrics"]
            assert len(data["metrics"]["distance"]) == len(data["frames"])

    def test_sparse_only_stores_active_frames(self, sparse_result):
        """稀疏结果只存储 active 的帧，metrics 长度应等于 frames 长度。"""
        sparse, _ = sparse_result
        for group_ids, data in sparse.data.items():
            assert len(data["frames"]) == len(data["metrics"]["distance"])

    def test_total_active_frame_pairs(self, sparse_result):
        """全部 active frame-pairs 应为 2473。"""
        sparse, _ = sparse_result
        total = sum(len(d["frames"]) for d in sparse.data.values())
        assert total == 2473


# ============================================================
# run_pass2：Pass2 补全
# ============================================================

class TestRunPass2:
    """验证 run_pass2 补全全帧 metric。"""

    def test_returns_interaction(self, sparse_result):
        sparse, det = sparse_result
        u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
        result = det.run_pass2(sparse, u.trajectory)
        assert len(result) > 0
        assert result[0].interaction_type == "salt_bridge"

    def test_n_pairs(self, sparse_result):
        sparse, det = sparse_result
        u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
        result = det.run_pass2(sparse, u.trajectory)
        assert result[0].n_pairs == 47

    def test_no_nan(self, sparse_result):
        """Pass2 应补全全部帧，无 NaN。"""
        sparse, det = sparse_result
        u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
        result = det.run_pass2(sparse, u.trajectory)
        assert np.sum(np.isnan(result[0].metrics["distance"])) == 0

    def test_total_values(self, sparse_result):
        """全部 47×101 个值应为真实 metric。"""
        sparse, det = sparse_result
        u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
        result = det.run_pass2(sparse, u.trajectory)
        assert result[0].metrics["distance"].shape == (47, 101)

    def test_active_distance_range(self, sparse_result):
        """活跃帧距离 ≤ 5.5Å。"""
        sparse, det = sparse_result
        u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
        result = det.run_pass2(sparse, u.trajectory)
        it = result[0]
        for i in range(it.n_pairs):
            active = it.metrics["distance"][i][it.existence[i]]
            if len(active) > 0:
                assert np.all(active <= 5.5)


# ============================================================
# detect：Pass1 + Pass2 完整流程
# ============================================================

class TestDetect:
    """验证 detect() = run_pass1 + run_pass2。"""

    def test_has_results(self, salt_bridges):
        assert len(salt_bridges) > 0

    def test_interaction_type(self, salt_bridges):
        assert salt_bridges[0].interaction_type == "salt_bridge"

    def test_n_pairs(self, salt_bridges):
        assert salt_bridges[0].n_pairs == 47

    def test_all_opposite_charge(self, salt_bridges):
        for g1, g2 in salt_bridges[0].groups:
            assert g1.group_type == "charged_positive"
            assert g2.group_type == "charged_negative"

    def test_no_nan(self, salt_bridges):
        """detect() 应补全全部帧。"""
        assert np.sum(np.isnan(salt_bridges[0].metrics["distance"])) == 0

    def test_occupancy_range(self, salt_bridges):
        occ = salt_bridges[0].occupancy()
        assert np.all(occ >= 0) and np.all(occ <= 1)

    def test_residue_coverage(self, salt_bridges):
        names = set()
        for g1, g2 in salt_bridges[0].groups:
            names.update([g1.residue_name, g2.residue_name])
        assert {"ARG", "LYS", "ASP", "GLU"}.issubset(names)


# ============================================================
# tuple_filter 测试
# ============================================================

class TestTupleFilter:
    """验证 tuple_filter 参数。"""

    def test_tuple_filter_reduces_pairs(self, filtered_groups):
        """只检测蛋白间的盐桥应减少 pair 数。"""
        u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
        det = SaltBridgeDetectorTwoPass()
        MOL_RBD = "seg_0_RBD_pro"
        MOL_KRAS = "seg_2_KRAS_pro"
        res = det.detect(
            filtered_groups, trajectory=u.trajectory,
            tuple_filter=lambda gt: {gt[0].molecule, gt[1].molecule} == {MOL_RBD, MOL_KRAS})
        assert len(res) > 0
        assert res[0].n_pairs < 47

    def test_tuple_filter_all_inter_protein(self, filtered_groups):
        """过滤后的 pair 应全部来自不同蛋白。"""
        u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
        det = SaltBridgeDetectorTwoPass()
        MOL_RBD = "seg_0_RBD_pro"
        MOL_KRAS = "seg_2_KRAS_pro"
        res = det.detect(
            filtered_groups, trajectory=u.trajectory,
            tuple_filter=lambda gt: {gt[0].molecule, gt[1].molecule} == {MOL_RBD, MOL_KRAS})
        for g1, g2 in res[0].groups:
            assert g1.molecule != g2.molecule


# ============================================================
# 检测器元信息
# ============================================================

class TestDetectorMeta:

    def test_name(self):
        assert SaltBridgeDetectorTwoPass().name == "salt_bridge"

    def test_required_group_types(self):
        types = SaltBridgeDetectorTwoPass().required_group_types
        assert "charged_positive" in types
        assert "charged_negative" in types

    def test_metric_names(self):
        assert "distance" in SaltBridgeDetectorTwoPass().metric_names


# ============================================================
# 与 PerFrame 交叉验证
# ============================================================

class TestVsPerFrame:

    def test_pair_count_matches(self, salt_bridges):
        u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
        per_frame = SaltBridgeDetectorPerFrame().detect(
            [g for g in AmberFFGroupIdentifier().identify(
                GmxTprReader().read(str(TPR_FILE)))
             if g.group_type in ("charged_positive", "charged_negative")],
            trajectory=u.trajectory)
        assert salt_bridges[0].n_pairs == per_frame[0].n_pairs

    def test_pair_set_matches(self, salt_bridges):
        u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
        per_frame = SaltBridgeDetectorPerFrame().detect(
            [g for g in AmberFFGroupIdentifier().identify(
                GmxTprReader().read(str(TPR_FILE)))
             if g.group_type in ("charged_positive", "charged_negative")],
            trajectory=u.trajectory)
        two = {(g1.group_id, g2.group_id) for g1, g2 in salt_bridges[0].groups}
        per = {(g1.group_id, g2.group_id) for g1, g2 in per_frame[0].groups}
        assert two == per
