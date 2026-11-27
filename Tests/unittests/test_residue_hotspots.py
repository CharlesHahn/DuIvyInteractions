# -*- coding: utf-8 -*-
"""残基热点聚合单元测试（visualizers.residue_hotspots）。

合成数据直接构造 Interaction/Group；真实数据对照 Tests/interaction_h5data 的 h5 原始数组。
"""

import time
from pathlib import Path

import h5py
import numpy as np
import pytest

from DuIvyInteractions.core.datas import AtomData, Group, Interaction
from DuIvyInteractions.visualizers import (
    hotspot_segments, parse_datasource, residue_hotspots)
from DuIvyInteractions.visualizers.data import (
    DISTANCE_METRIC_BY_TYPE, _dominant_type)

FRAMES = 4
H5_DATA_DIR = Path(__file__).parent.parent / "interaction_h5data"
WATER_GROUP_TYPE = "water"


def _atom(gidx: int) -> AtomData:
    """占位原子（热点聚合不读原子）。"""
    return AtomData(atom_global_idx=gidx, atom_idx_in_residue=0, atom_name="C",
                    atom_type="", atom_element="C", atom_charge=0.0,
                    atom_mass=0.0)


def _group(gid: int, gtype: str, res: str, rid: int,
           molecule: str = "test") -> Group:
    """按 (group_type, residue_name, residue_id[, molecule]) 造基团。"""
    return Group(group_id=gid, group_type=gtype, molecule=molecule,
                 residue_name=res, residue_id=rid, atoms=[_atom(gid)])


def _interaction(itype: str, specs, existence, distances=None) -> Interaction:
    """specs 每对为 [(group_type, residue_name, residue_id), ...]。"""
    groups = [tuple(_group(10 * i + pos, *spec) for pos, spec in enumerate(pair))
              for i, pair in enumerate(specs)]
    existence = np.asarray(existence, dtype=bool)
    metrics = {}
    if distances is not None:
        metrics["distance"] = np.asarray(distances, dtype=float)
    return Interaction(interaction_type=itype, groups=groups,
                       existence=existence, metrics=metrics,
                       times=np.arange(existence.shape[1], dtype=float))


def _item(itype: str, specs, existence, distances=None, display=None):
    """(显示名, h5 路径, Interaction) 三元组。"""
    return (display or itype, f"/tmp/{itype}.h5",
            _interaction(itype, specs, existence, distances))


def _item_with_metrics(itype: str, specs, existence, metrics, display=None):
    """(显示名, h5 路径, Interaction) 三元组（自定 metrics）。"""
    interaction = _interaction(itype, specs, existence)
    interaction.metrics = metrics
    return (display or itype, f"/tmp/{itype}.h5", interaction)


def _row(rows, residue_id: int, residue_name: str):
    """按残基身份取一行。"""
    return next(r for r in rows if r["residue_id"] == residue_id
                and r["residue_name"] == residue_name)


LYS_ASP = [[("H_donor", "LYS", 10), ("H_acceptor", "ASP", 100)]]
LYS_CHARGED = [[("charged_positive", "LYS", 10),
                ("charged_negative", "ASP", 100)]]
WATER_PAIR = [[("H_donor", "LYS", 10), (WATER_GROUP_TYPE, "SOL", 500),
               ("H_acceptor", "ASP", 100)]]


# ============================================================
# 合成数据：单类型单对
# ============================================================

class TestSinglePair:

    def test_frames_and_occupancy(self):
        rows = residue_hotspots([_item("hydrogen_bond", LYS_ASP,
                                       [[True, True, False, False]])])
        lys = _row(rows, 10, "LYS")
        assert lys["frames_active"] == 2
        assert lys["occupancy"] == pytest.approx(0.5)
        assert lys["n_pairs"] == 1
        assert lys["n_partners"] == 1
        assert lys["n_types"] == 1
        assert lys["label"] == "LYS10"

    def test_type_detail_and_average_distance(self):
        rows = residue_hotspots([_item("hydrogen_bond", LYS_ASP,
                                       [[True, True, False, False]],
                                       [[3.0, 2.0, 9.0, 9.0]])])
        lys = _row(rows, 10, "LYS")
        assert lys["types"] == [{"type": "hydrogen_bond", "frames": 2,
                                 "occupancy": 0.5, "pairs": 1}]
        assert lys["avg_distance"] == pytest.approx(2.5)

    def test_acceptor_symmetric(self):
        rows = residue_hotspots([_item("hydrogen_bond", LYS_ASP,
                                       [[False, True, True, True]])])
        asp = _row(rows, 100, "ASP")
        assert asp["frames_active"] == 3
        assert asp["occupancy"] == pytest.approx(0.75)
        assert asp["n_partners"] == 1


# ============================================================
# 合成数据：并集口径
# ============================================================

class TestUnionAcrossPairs:

    def test_frames_active_is_union_not_sum(self):
        specs = [[("H_donor", "LYS", 10), ("H_acceptor", "ASP", 100)],
                 [("H_donor", "LYS", 10), ("H_acceptor", "GLU", 200)]]
        existence = [[True, True, False, False],
                     [False, True, True, False]]
        lys = _row(residue_hotspots([_item("hydrogen_bond", specs, existence)]),
                   10, "LYS")
        assert lys["frames_active"] == 3
        assert lys["occupancy"] == pytest.approx(0.75)
        assert lys["n_pairs"] == 2
        assert lys["n_partners"] == 2


class TestCrossTypeUnion:

    def test_same_frame_two_types_counted_once(self):
        hb = _item("hydrogen_bond", LYS_ASP, [[True, True, False, False]])
        sb = _item("salt_bridge", LYS_CHARGED, [[True, False, True, False]])
        lys = _row(residue_hotspots([hb, sb]), 10, "LYS")
        assert lys["frames_active"] == 3
        assert lys["occupancy"] == pytest.approx(0.75)
        assert lys["n_types"] == 2


# ============================================================
# 合成数据：伙伴与类型明细
# ============================================================

class TestPartners:

    def test_three_distinct_partners(self):
        specs = [[("H_donor", "LYS", 10), ("H_acceptor", "ASP", 100)],
                 [("H_donor", "LYS", 10), ("H_acceptor", "GLU", 200)],
                 [("H_donor", "LYS", 10), ("H_acceptor", "SER", 300)]]
        existence = np.zeros((3, FRAMES), dtype=bool)
        lys = _row(residue_hotspots([_item("hydrogen_bond", specs, existence)]),
                   10, "LYS")
        assert lys["n_partners"] == 3

    def test_same_partner_twice_counted_once(self):
        specs = [[("H_donor", "LYS", 10), ("H_acceptor", "ASP", 100)],
                 [("H_donor", "LYS", 10), ("H_acceptor", "ASP", 100)]]
        existence = np.zeros((2, FRAMES), dtype=bool)
        lys = _row(residue_hotspots([_item("hydrogen_bond", specs, existence)]),
                   10, "LYS")
        assert lys["n_pairs"] == 2
        assert lys["n_partners"] == 1


class TestTypeDetail:

    def test_each_type_frames_and_occupancy(self):
        hb = _item("hydrogen_bond", LYS_ASP, [[True, True, False, False]])
        sb = _item("salt_bridge", LYS_CHARGED, [[True, False, True, False]])
        lys = _row(residue_hotspots([hb, sb]), 10, "LYS")
        by_type = {entry["type"]: entry for entry in lys["types"]}
        assert by_type["hydrogen_bond"]["frames"] == 2
        assert by_type["hydrogen_bond"]["occupancy"] == pytest.approx(0.5)
        assert by_type["salt_bridge"]["frames"] == 2
        assert by_type["salt_bridge"]["occupancy"] == pytest.approx(0.5)
        assert [entry["type"] for entry in lys["types"]] == \
            ["hydrogen_bond", "salt_bridge"]


# ============================================================
# 合成数据：过滤、水、边界
# ============================================================

class TestTypeFilter:

    def test_only_selected_display_name(self):
        hb = _item("hydrogen_bond", LYS_ASP, [[True, True, False, False]])
        sb = _item("salt_bridge", LYS_CHARGED, [[True, False, False, False]],
                   display="salt_bridge@lig.h5")
        lys = _row(residue_hotspots([hb, sb], type_filter="salt_bridge@lig.h5"),
                   10, "LYS")
        assert lys["n_types"] == 1
        assert lys["types"][0]["type"] == "salt_bridge"
        assert lys["frames_active"] == 1

    def test_unknown_filter_returns_empty(self):
        hb = _item("hydrogen_bond", LYS_ASP, [[True, False, False, False]])
        assert residue_hotspots([hb], type_filter="pi_stacking") == []


class TestWater:

    def test_default_excludes_water(self):
        rows = residue_hotspots([_item("water_bridge", WATER_PAIR,
                                       [[True, False, False, False]])])
        assert all(row["residue_name"] != "SOL" for row in rows)
        assert _row(rows, 10, "LYS")["n_partners"] == 1

    def test_include_water_adds_water_residue(self):
        rows = residue_hotspots([_item("water_bridge", WATER_PAIR,
                                       [[True, False, False, False]])],
                                include_water=True)
        sol = _row(rows, 500, "SOL")
        assert sol["frames_active"] == 1
        assert sol["n_partners"] == 2


class TestZeroFrames:

    def test_no_frames_occupancy_zero(self):
        rows = residue_hotspots([_item("salt_bridge", LYS_CHARGED,
                                       np.zeros((1, 0), dtype=bool))])
        lys = _row(rows, 10, "LYS")
        assert lys["frames_active"] == 0
        assert lys["occupancy"] == 0.0
        assert lys["types"] == [{"type": "salt_bridge", "frames": 0,
                                 "occupancy": 0.0, "pairs": 1}]


class TestAverageDistance:

    def test_only_active_frames(self):
        rows = residue_hotspots([_item("hydrogen_bond", LYS_ASP,
                                       [[True, False, True, False]],
                                       [[2.0, 99.0, 4.0, 99.0]])])
        assert _row(rows, 10, "LYS")["avg_distance"] == pytest.approx(3.0)

    def test_missing_metric_is_none(self):
        rows = residue_hotspots([_item("hydrogen_bond", LYS_ASP,
                                       [[True, True, False, False]])])
        assert _row(rows, 10, "LYS")["avg_distance"] is None

    def test_shape_mismatch_skipped(self):
        item = _item("hydrogen_bond", LYS_ASP, [[True, True, False, False]],
                     [[2.0, 2.0, 2.0, 2.0]])
        item[2].metrics["distance"] = np.zeros((1, FRAMES + 1))
        assert _row(residue_hotspots([item]), 10, "LYS")["avg_distance"] is None


class TestAverageDistanceByType:

    def test_water_bridge_uses_d_ow(self):
        item = _item_with_metrics(
            "water_bridge", WATER_PAIR, [[True, True, False, False]],
            {"dist_d_ow": np.array([[2.8, 3.2, 9.0, 9.0]]),
             "dist_h_ow": np.array([[1.0, 1.0, 9.0, 9.0]])})
        rows = residue_hotspots([item])
        assert _row(rows, 10, "LYS")["avg_distance"] == pytest.approx(3.0)
        assert _row(rows, 100, "ASP")["avg_distance"] == pytest.approx(3.0)

    def test_water_bridge_without_d_ow_is_none(self):
        item = _item_with_metrics(
            "water_bridge", WATER_PAIR, [[True, True, False, False]],
            {"dist_h_ow": np.array([[1.0, 1.0, 9.0, 9.0]])})
        assert _row(residue_hotspots([item]), 10, "LYS")["avg_distance"] is None

    def test_water_bridge_falls_back_to_distance(self):
        item = _item_with_metrics(
            "water_bridge", WATER_PAIR, [[True, True, False, False]],
            {"distance": np.array([[3.4, 3.6, 9.0, 9.0]])})
        assert _row(residue_hotspots([item]), 10, "LYS")["avg_distance"] == \
            pytest.approx(3.5)

    def test_hydrogen_bond_still_uses_distance(self):
        item = _item_with_metrics(
            "hydrogen_bond", LYS_ASP, [[True, True, False, False]],
            {"distance": np.array([[2.0, 2.4, 9.0, 9.0]]),
             "dist_d_ow": np.array([[5.0, 5.0, 9.0, 9.0]])})
        assert _row(residue_hotspots([item]), 10, "LYS")["avg_distance"] == \
            pytest.approx(2.2)

    def test_metric_map_only_water_bridge(self):
        assert DISTANCE_METRIC_BY_TYPE == {"water_bridge": "dist_d_ow"}


class TestFrameDenominator:

    def test_denominator_is_max_frames_across_items(self):
        short = _item("hydrogen_bond", LYS_ASP, [[True, True, False, False]])
        long_specs = [[("charged_positive", "ARG", 40),
                       ("charged_negative", "GLU", 400)]]
        long_item = _item("salt_bridge", long_specs, np.zeros((1, 8), dtype=bool))
        lys = _row(residue_hotspots([short, long_item]), 10, "LYS")
        assert lys["frames_active"] == 2
        assert lys["occupancy"] == pytest.approx(2 / 8)


class TestDuplicateResidueInPair:

    def test_same_residue_twice_counts_one_pair(self):
        specs = [[("hydrophobic", "LEU", 50), ("hydrophobic", "LEU", 50)]]
        rows = residue_hotspots([_item("hydrophobic", specs,
                                       [[True, True, False, False]])])
        leu = _row(rows, 50, "LEU")
        assert leu["n_pairs"] == 1
        assert leu["frames_active"] == 2
        assert leu["n_partners"] == 0


# ============================================================
# 合成数据：排序与空输入
# ============================================================

class TestOrdering:

    def test_occupancy_then_pairs_then_residue_id(self):
        specs = [[("H_donor", "LYS", 30), ("H_acceptor", "ASP", 130)],
                 [("H_donor", "LYS", 30), ("H_acceptor", "GLU", 230)],
                 [("H_donor", "LYS", 20), ("H_acceptor", "ASP", 120)],
                 [("H_donor", "LYS", 25), ("H_acceptor", "ASP", 125)],
                 [("H_donor", "ARG", 40), ("H_acceptor", "ASP", 140)]]
        existence = [[True, True, True, True],
                     [True, True, True, True],
                     [True, True, True, True],
                     [True, True, True, True],
                     [True, True, False, False]]
        rows = residue_hotspots([_item("hydrogen_bond", specs, existence)])
        assert [row["label"] for row in rows] == [
            "LYS30", "LYS20", "LYS25", "ASP120", "ASP125", "ASP130",
            "GLU230", "ARG40", "ASP140"]


class TestEmptyInput:

    def test_empty_items_returns_empty_list(self):
        assert residue_hotspots([]) == []


class TestAnomalousResidueData:

    def test_same_id_two_names_kept_apart(self):
        specs = [[("H_donor", "LYS", 10), ("H_acceptor", "ASP", 100)],
                 [("H_donor", "ARG", 10), ("H_acceptor", "GLU", 200)]]
        existence = [[True, False, False, False],
                     [False, True, False, False]]
        rows = residue_hotspots([_item("hydrogen_bond", specs, existence)])
        same_id = [row for row in rows if row["residue_id"] == 10]
        assert {row["residue_name"] for row in same_id} == {"LYS", "ARG"}
        assert all(row["frames_active"] == 1 and row["n_pairs"] == 1
                   for row in same_id)
        assert {row["label"] for row in same_id} == {"LYS10", "ARG10"}


# ============================================================
# 残基元信息：molecule / group_types / types[].pairs / dominant_type
# ============================================================

class TestRowMetadata:

    def test_molecule_and_group_types_passed_through(self):
        specs = [[("H_donor", "LYS", 10, "seg_0_RBD_pro"),
                  ("H_acceptor", "ASP", 100, "seg_2_KRAS_pro")]]
        rows = residue_hotspots([_item("hydrogen_bond", specs,
                                       [[True, False, False, False]])])
        lys = _row(rows, 10, "LYS")
        assert lys["molecule"] == "seg_0_RBD_pro"
        assert lys["group_types"] == ["H_donor"]
        asp = _row(rows, 100, "ASP")
        assert asp["molecule"] == "seg_2_KRAS_pro"
        assert asp["group_types"] == ["H_acceptor"]

    def test_multiple_molecules_take_smallest(self):
        specs = [[("H_donor", "LYS", 10, "seg_2"),
                  ("H_acceptor", "ASP", 100, "seg_0")],
                 [("H_donor", "LYS", 10, "seg_1"),
                  ("H_acceptor", "GLU", 200, "seg_3")]]
        rows = residue_hotspots([_item("hydrogen_bond", specs,
                                       np.zeros((2, FRAMES), dtype=bool))])
        assert _row(rows, 10, "LYS")["molecule"] == "seg_1"

    def test_group_types_union_of_own_groups(self):
        hb = _item("hydrogen_bond",
                   [[("H_donor", "LYS", 10, "seg_0"),
                     ("H_acceptor", "ASP", 100, "seg_0")]],
                   np.zeros((1, FRAMES), dtype=bool))
        sb = _item("salt_bridge",
                   [[("charged_positive", "LYS", 10, "seg_0"),
                     ("charged_negative", "ASP", 100, "seg_0")]],
                   np.zeros((1, FRAMES), dtype=bool))
        lys = _row(residue_hotspots([hb, sb]), 10, "LYS")
        assert lys["group_types"] == ["H_donor", "charged_positive"]


class TestTypePairs:

    def test_pairs_per_type(self):
        specs = [[("H_donor", "LYS", 10), ("H_acceptor", "ASP", 100)],
                 [("H_donor", "LYS", 10), ("H_acceptor", "GLU", 200)]]
        existence = [[True, False, False, False],
                     [False, True, False, False]]
        lys = _row(residue_hotspots([_item("hydrogen_bond", specs, existence)]),
                   10, "LYS")
        assert lys["types"][0]["pairs"] == 2
        assert lys["n_pairs"] == 2

    def test_self_pair_counts_once(self):
        specs = [[("hydrophobic", "LEU", 50), ("hydrophobic", "LEU", 50)]]
        rows = residue_hotspots([_item("hydrophobic", specs,
                                       [[True, True, False, False]])])
        leu = _row(rows, 50, "LEU")
        assert leu["types"][0]["pairs"] == 1
        assert leu["n_pairs"] == 1


class TestDominantType:

    def test_more_frames_wins(self):
        hb = _item("hydrogen_bond", LYS_ASP, [[True, True, False, False]])
        sb = _item("salt_bridge", LYS_CHARGED, [[True, False, False, False]])
        assert _row(residue_hotspots([hb, sb]), 10, "LYS")["dominant_type"] \
            == "hydrogen_bond"

    def test_pairs_break_frame_tie(self):
        hb = _item("hydrogen_bond", LYS_ASP, [[True, True, False, False]])
        sb = _item("salt_bridge",
                   [[("charged_positive", "LYS", 10),
                     ("charged_negative", "ASP", 100)],
                    [("charged_positive", "LYS", 10),
                     ("charged_negative", "GLU", 200)]],
                   [[True, False, False, False],
                    [False, True, False, False]])
        lys = _row(residue_hotspots([hb, sb]), 10, "LYS")
        assert lys["dominant_type"] == "salt_bridge"

    def test_name_breaks_tie(self):
        hb = _item("hydrogen_bond", LYS_ASP, [[True, True, False, False]])
        sb = _item("salt_bridge", LYS_CHARGED, [[True, True, False, False]])
        assert _row(residue_hotspots([hb, sb]), 10, "LYS")["dominant_type"] \
            == "salt_bridge"

    def test_no_types_returns_none(self):
        assert _dominant_type([]) is None


# ============================================================
# 分段（hotspot_segments）
# ============================================================

class TestHotspotSegments:

    def test_segment_and_inner_ordering(self):
        specs = [[("H_donor", "LYS", 30, "seg_1"),
                  ("H_acceptor", "ASP", 130, "seg_1")],
                 [("H_donor", "LYS", 20, "seg_0"),
                  ("H_acceptor", "ASP", 120, "seg_0")],
                 [("H_donor", "LYS", 25, "seg_0"),
                  ("H_acceptor", "ASP", 125, "seg_0")]]
        existence = [[True, False, False, False]] * 3
        rows = residue_hotspots([_item("hydrogen_bond", specs, existence)])
        segments = hotspot_segments(rows)
        assert [seg["molecule"] for seg in segments] == ["seg_0", "seg_1"]
        assert [r["residue_id"] for r in segments[0]["residues"]] == \
            [20, 25, 120, 125]
        assert [r["residue_id"] for r in segments[1]["residues"]] == [30, 130]

    def test_zero_occupancy_dropped_by_default(self):
        specs = [[("H_donor", "LYS", 10), ("H_acceptor", "ASP", 100)],
                 [("H_donor", "ARG", 40), ("H_acceptor", "GLU", 200)]]
        existence = [[True, False, False, False],
                     [False, False, False, False]]
        rows = residue_hotspots([_item("hydrogen_bond", specs, existence)])
        labels = [r["label"] for seg in hotspot_segments(rows)
                  for r in seg["residues"]]
        assert "ARG40" not in labels and "GLU200" not in labels
        kept = [r["label"] for seg in hotspot_segments(rows, exclude_unused=False)
                for r in seg["residues"]]
        assert {"ARG40", "GLU200"} <= set(kept)

    def test_unused_ion_dropped_used_ion_kept(self):
        used = _item("metal_coordination",
                     [[("metal", "MG", 311, "seg_4_Mg"),
                       ("metal_binding", "GLU", 15, "seg_2_KRAS_pro")]],
                     [[True, False, False, False]])
        stray = _item("hydrogen_bond",
                      [[("metal", "NA", 900, "seg_6_NA"),
                        ("H_acceptor", "ASP", 100, "seg_0")]],
                      [[True, False, False, False]])
        rows = residue_hotspots([used, stray])
        labels = [r["label"] for seg in hotspot_segments(rows)
                  for r in seg["residues"]]
        assert "MG311" in labels
        assert "NA900" not in labels
        kept = [r["label"] for seg in hotspot_segments(rows, exclude_unused=False)
                for r in seg["residues"]]
        assert "NA900" not in kept      # 离子规则与 exclude_unused 无关

    def test_counts_and_max_occupancy(self):
        specs = [[("H_donor", "LYS", 10, "seg_0"),
                  ("H_acceptor", "ASP", 100, "seg_0")],
                 [("H_donor", "LYS", 12, "seg_0"),
                  ("H_acceptor", "ASP", 102, "seg_0")]]
        existence = [[True, True, True, True], [True, True, False, False]]
        rows = residue_hotspots([_item("hydrogen_bond", specs, existence)])
        segment = hotspot_segments(rows)[0]
        assert segment["n_residues"] == 4
        assert segment["max_occupancy"] == pytest.approx(1.0)

    def test_residues_are_original_rows(self):
        rows = residue_hotspots([_item("hydrogen_bond", LYS_ASP,
                                       [[True, False, False, False]])])
        residues = [r for seg in hotspot_segments(rows) for r in seg["residues"]]
        assert all(any(r is row for row in rows) for r in residues)

    def test_empty_rows_returns_empty_list(self):
        assert hotspot_segments([]) == []


# ============================================================
# 真实数据（Tests/interaction_h5data/）
# ============================================================

def _real_items(filename: str):
    """真实 h5 → items（数据缺失时跳过）。"""
    path = H5_DATA_DIR / filename
    if not path.exists():
        pytest.skip(f"real h5 test data not available: {filename}")
    return parse_datasource(str(path))


def _real_all_items():
    """全部 8 类真实 h5 → items（数据缺失时跳过）。"""
    if not (H5_DATA_DIR / "water_bridge.h5").exists():
        pytest.skip("real h5 test data not available")
    return parse_datasource(str(H5_DATA_DIR))


def _bruteforce_frames(interaction: Interaction, residue_id: int,
                       residue_name: str, include_water: bool = False) -> int:
    """参考实现：逐对 OR 出该残基的活跃帧掩码（不复用被测函数）。"""
    mask = np.zeros(interaction.n_frames, dtype=bool)
    for pair_idx, pair in enumerate(interaction.groups):
        for group in pair:
            if not include_water and group.group_type == WATER_GROUP_TYPE:
                continue
            if (int(group.residue_id) == residue_id
                    and str(group.residue_name) == residue_name):
                mask |= np.asarray(interaction.existence[pair_idx], dtype=bool)
                break
    return int(mask.sum())


def _raw_group_keys(filename: str):
    """直接读 h5 原始数组 → [(residue_id, residue_name)]。"""
    with h5py.File(H5_DATA_DIR / filename, "r") as handle:
        group = handle["interaction_0/groups"]
        ids = group["residue_id"][:]
        names = [n.decode() if isinstance(n, bytes) else str(n)
                 for n in group["residue_name"][:]]
    return [(int(i), n) for i, n in zip(ids, names)]


class TestRealSaltBridge:

    def test_charged_residue_names_present(self):
        rows = residue_hotspots(_real_items("salt_bridge.h5"))
        assert {"ARG", "LYS", "ASP", "GLU"} <= {r["residue_name"] for r in rows}

    def test_top1_frames_match_bruteforce(self):
        items = _real_items("salt_bridge.h5")
        top = residue_hotspots(items)[0]
        assert top["frames_active"] == _bruteforce_frames(
            items[0][2], top["residue_id"], top["residue_name"])


class TestRealAllTypes:

    def test_types_cap_and_occupancy_range(self):
        rows = residue_hotspots(_real_all_items())
        assert rows
        assert max(row["n_types"] for row in rows) <= 8
        assert all(0.0 <= row["occupancy"] <= 1.0 for row in rows)


class TestRealWaterBridge:

    def test_default_excludes_water_residue(self):
        rows = residue_hotspots(_real_items("water_bridge.h5"))
        assert rows
        assert all(row["residue_name"] != "SOL" for row in rows)

    def test_include_water_contains_water_residue(self):
        rows = residue_hotspots(_real_items("water_bridge.h5"),
                                include_water=True)
        assert any(row["residue_name"] == "SOL" for row in rows)


class TestRealResidueIdentity:

    def test_residue_id_equals_h5_raw_value(self):
        rows = residue_hotspots(_real_items("salt_bridge.h5"))
        sample_id, _ = _raw_group_keys("salt_bridge.h5")[0]
        row = next(r for r in rows if r["residue_id"] == sample_id)
        assert row["residue_id"] == sample_id

    def test_residue_name_and_label_equal_h5_raw_value(self):
        rows = residue_hotspots(_real_items("salt_bridge.h5"))
        sample_id, sample_name = _raw_group_keys("salt_bridge.h5")[0]
        row = next(r for r in rows if r["residue_id"] == sample_id
                   and r["residue_name"] == sample_name)
        assert row["residue_name"] == sample_name
        assert row["label"] == f"{sample_name}{sample_id}"

    def test_identity_set_equals_h5_groups(self):
        rows = residue_hotspots(_real_items("salt_bridge.h5"))
        assert {(r["residue_id"], r["residue_name"]) for r in rows} == \
            set(_raw_group_keys("salt_bridge.h5"))


class TestRealPerformance:

    def test_all_types_aggregate_under_two_seconds(self):
        items = _real_all_items()
        start = time.perf_counter()
        rows = residue_hotspots(items)
        elapsed = time.perf_counter() - start
        print(f"\n[hotspots] 8 类聚合 {len(rows)} 个残基，耗时 {elapsed:.3f}s")
        assert rows


class TestRealSegments:

    def test_segments_keep_metal_and_drop_stray_ions(self):
        rows = residue_hotspots(_real_all_items())
        segments = hotspot_segments(rows)
        print(f"\n[segments] 共 {len(segments)} 段: " + ", ".join(
            f"{seg['molecule']}={seg['n_residues']}/{seg['max_occupancy']:.3f}"
            for seg in segments))
        assert segments
        residues = [r for seg in segments for r in seg["residues"]]
        assert "MG311" in {r["label"] for r in residues}
        stray = [r for r in residues if "metal" in r["group_types"]
                 and "metal_coordination" not in {t["type"] for t in r["types"]}]
        assert stray == []
