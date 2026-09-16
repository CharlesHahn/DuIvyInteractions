# -*- coding: utf-8 -*-
"""导出器单元测试。

覆盖：
1. to_xvg_metric — 指标折线图
2. to_xvg_count — 每帧活跃数量
3. to_xpm_existence — existence 热力图
4. save 方法 — 文件写入
5. 8 个子类的 get_pair_label — 标签格式
"""

import pytest
import numpy as np
from pathlib import Path

from DuIvyInteractions.core.datas import Interaction, Group, AtomData
from DuIvyInteractions.io.h5 import load_interactions
from DuIvyInteractions.io.hydrogen_bond_exporter import HydrogenBondExporter
from DuIvyInteractions.io.saltbridge_exporter import SaltBridgeExporter
from DuIvyInteractions.io.pi_stacking_exporter import PiStackingExporter
from DuIvyInteractions.io.pi_cation_exporter import PiCationExporter
from DuIvyInteractions.io.halogen_bond_exporter import HalogenBondExporter
from DuIvyInteractions.io.hydrophobic_exporter import HydrophobicExporter
from DuIvyInteractions.io.metal_coordination_exporter import MetalCoordinationExporter
from DuIvyInteractions.io.water_bridge_exporter import WaterBridgeExporter


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


def _make_atom(gidx: int, name: str, element: str, charge: float = 0.0) -> AtomData:
    return AtomData(
        atom_global_idx=gidx, atom_idx_in_residue=0,
        atom_name=name, atom_type="", atom_element=element,
        atom_charge=charge, atom_mass=0.0,
    )


def _make_group(gid: int, gtype: str, atoms: list,
                res_name: str = "RES", res_id: int = 1) -> Group:
    return Group(
        group_id=gid, group_type=gtype, molecule="test",
        residue_name=res_name, residue_id=res_id, atoms=atoms,
    )


@pytest.fixture
def saltbridge_interaction():
    """构造一个简单的盐桥 Interaction（3 pairs, 5 frames）。"""
    pos1 = _make_group(0, "charged_positive",
                        [_make_atom(10, "NZ", "N", 0.5),
                         _make_atom(11, "HZ1", "H", 0.3)],
                        "LYS", 10)
    pos2 = _make_group(1, "charged_positive",
                        [_make_atom(20, "NE", "N", 0.4),
                         _make_atom(21, "HE", "H", 0.3)],
                        "ARG", 20)
    pos3 = _make_group(2, "charged_positive",
                        [_make_atom(30, "NZ", "N", 0.5),
                         _make_atom(31, "HZ1", "H", 0.3)],
                        "LYS", 30)
    neg1 = _make_group(3, "charged_negative",
                        [_make_atom(40, "OD1", "O", -0.6),
                         _make_atom(41, "OD2", "O", -0.6)],
                        "ASP", 40)
    neg2 = _make_group(4, "charged_negative",
                        [_make_atom(50, "OE1", "O", -0.6),
                         _make_atom(51, "OE2", "O", -0.6)],
                        "GLU", 50)

    existence = np.array([
        [True,  True,  False, True,  True],   # pair0: 4/5
        [False, True,  True,  True,  False],   # pair1: 3/5
        [True,  True,  True,  True,  True],    # pair2: 5/5
    ])
    metrics = {
        "distance": np.array([
            [3.5, 3.2, 6.0, 3.8, 3.1],
            [5.0, 4.5, 4.8, 4.2, 5.5],
            [3.0, 3.1, 3.2, 3.3, 3.4],
        ]),
    }
    times = np.array([0.0, 10.0, 20.0, 30.0, 40.0])

    return Interaction(
        interaction_type="salt_bridge",
        groups=[(pos1, neg1), (pos2, neg2), (pos3, neg1)],
        existence=existence, metrics=metrics, times=times,
    )


@pytest.fixture
def saltbridge_exporter():
    return SaltBridgeExporter()


# ============================================================
# 1. to_xvg_metric
# ============================================================

class TestToXvgMetric:

    def test_returns_xvg(self, saltbridge_interaction, saltbridge_exporter):
        xvg = saltbridge_exporter.to_xvg_metric(
            saltbridge_interaction, "distance")
        assert hasattr(xvg, 'data_columns')

    def test_title(self, saltbridge_interaction, saltbridge_exporter):
        xvg = saltbridge_exporter.to_xvg_metric(
            saltbridge_interaction, "distance")
        assert "Salt Bridge" in xvg.title
        assert "Charge Center Distance" in xvg.title

    def test_xlabel_default(self, saltbridge_interaction, saltbridge_exporter):
        xvg = saltbridge_exporter.to_xvg_metric(
            saltbridge_interaction, "distance")
        assert xvg.xlabel == "Time (ps)"

    def test_ylabel_from_metric_labels(self, saltbridge_interaction, saltbridge_exporter):
        xvg = saltbridge_exporter.to_xvg_metric(
            saltbridge_interaction, "distance")
        assert "Charge Center Distance" in xvg.ylabel

    def test_column_count(self, saltbridge_interaction, saltbridge_exporter):
        xvg = saltbridge_exporter.to_xvg_metric(
            saltbridge_interaction, "distance")
        # 1 (times) + 3 (pairs) = 4
        assert xvg.column_num == 4

    def test_row_count(self, saltbridge_interaction, saltbridge_exporter):
        xvg = saltbridge_exporter.to_xvg_metric(
            saltbridge_interaction, "distance")
        assert xvg.row_num == 5

    def test_first_column_is_times(self, saltbridge_interaction, saltbridge_exporter):
        xvg = saltbridge_exporter.to_xvg_metric(
            saltbridge_interaction, "distance")
        np.testing.assert_array_equal(
            xvg.data_columns[0], [0.0, 10.0, 20.0, 30.0, 40.0])

    def test_metric_data_correct(self, saltbridge_interaction, saltbridge_exporter):
        xvg = saltbridge_exporter.to_xvg_metric(
            saltbridge_interaction, "distance")
        # pair0 distance
        np.testing.assert_array_almost_equal(
            xvg.data_columns[1], [3.5, 3.2, 6.0, 3.8, 3.1])

    def test_legends_count(self, saltbridge_interaction, saltbridge_exporter):
        xvg = saltbridge_exporter.to_xvg_metric(
            saltbridge_interaction, "distance")
        assert len(xvg.legends) == 3

    def test_pair_indices_filter(self, saltbridge_interaction, saltbridge_exporter):
        xvg = saltbridge_exporter.to_xvg_metric(
            saltbridge_interaction, "distance", pair_indices=[0, 2])
        assert xvg.column_num == 3  # 1 times + 2 pairs
        assert len(xvg.legends) == 2

    def test_invalid_metric_raises(self, saltbridge_interaction, saltbridge_exporter):
        with pytest.raises(ValueError, match="not found"):
            saltbridge_exporter.to_xvg_metric(
                saltbridge_interaction, "nonexistent")

    def test_empty_pairs(self, saltbridge_exporter):
        """n_pairs=0 时只产生时间列，无 metric 列。"""
        it = Interaction(
            interaction_type="test", groups=[], existence=np.empty((0, 5)),
            metrics={"distance": np.empty((0, 5))},
            times=np.arange(5, dtype=float),
        )
        xvg = saltbridge_exporter.to_xvg_metric(it, "distance")
        assert xvg.column_num == 1  # 只有时间列
        assert len(xvg.legends) == 0


# ============================================================
# 2. to_xvg_count
# ============================================================

class TestToXvgCount:

    def test_returns_xvg(self, saltbridge_interaction, saltbridge_exporter):
        xvg = saltbridge_exporter.to_xvg_count(saltbridge_interaction)
        assert hasattr(xvg, 'data_columns')

    def test_column_count(self, saltbridge_interaction, saltbridge_exporter):
        xvg = saltbridge_exporter.to_xvg_count(saltbridge_interaction)
        assert xvg.column_num == 2

    def test_first_column_is_times(self, saltbridge_interaction, saltbridge_exporter):
        xvg = saltbridge_exporter.to_xvg_count(saltbridge_interaction)
        np.testing.assert_array_equal(
            xvg.data_columns[0], [0.0, 10.0, 20.0, 30.0, 40.0])

    def test_count_correct(self, saltbridge_interaction, saltbridge_exporter):
        xvg = saltbridge_exporter.to_xvg_count(saltbridge_interaction)
        # existence: [[T,T,F,T,T], [F,T,T,T,F], [T,T,T,T,T]]
        # count:     [2, 3, 2, 3, 2]
        np.testing.assert_array_equal(
            xvg.data_columns[1], [2, 3, 2, 3, 2])

    def test_row_count(self, saltbridge_interaction, saltbridge_exporter):
        xvg = saltbridge_exporter.to_xvg_count(saltbridge_interaction)
        assert xvg.row_num == 5

    def test_title(self, saltbridge_interaction, saltbridge_exporter):
        xvg = saltbridge_exporter.to_xvg_count(saltbridge_interaction)
        assert "Salt Bridge" in xvg.title
        assert "Count" in xvg.title


# ============================================================
# 3. to_xpm_existence
# ============================================================

class TestToXpmExistence:

    def test_returns_xpm(self, saltbridge_interaction, saltbridge_exporter):
        xpm = saltbridge_exporter.to_xpm_existence(saltbridge_interaction)
        assert hasattr(xpm, 'value_matrix')

    def test_width_equals_frames(self, saltbridge_interaction, saltbridge_exporter):
        xpm = saltbridge_exporter.to_xpm_existence(saltbridge_interaction)
        assert xpm.width == 5

    def test_height_equals_pairs(self, saltbridge_interaction, saltbridge_exporter):
        xpm = saltbridge_exporter.to_xpm_existence(saltbridge_interaction)
        assert xpm.height == 3

    def test_value_matrix_is_01(self, saltbridge_interaction, saltbridge_exporter):
        xpm = saltbridge_exporter.to_xpm_existence(saltbridge_interaction)
        flat = [v for row in xpm.value_matrix for v in row]
        assert set(flat) == {0, 1}

    def test_value_matrix_matches_existence(self, saltbridge_interaction, saltbridge_exporter):
        xpm = saltbridge_exporter.to_xpm_existence(saltbridge_interaction)
        expected = saltbridge_interaction.existence.astype(int).tolist()
        assert xpm.value_matrix == expected

    def test_xaxis_is_times(self, saltbridge_interaction, saltbridge_exporter):
        xpm = saltbridge_exporter.to_xpm_existence(saltbridge_interaction)
        np.testing.assert_array_equal(
            xpm.xaxis, [0.0, 10.0, 20.0, 30.0, 40.0])

    def test_yaxis_is_sequential(self, saltbridge_interaction, saltbridge_exporter):
        xpm = saltbridge_exporter.to_xpm_existence(saltbridge_interaction)
        assert xpm.yaxis == [0, 1, 2]

    def test_legend_has_mapping(self, saltbridge_interaction, saltbridge_exporter):
        xpm = saltbridge_exporter.to_xpm_existence(saltbridge_interaction)
        assert "0:" in xpm.legend
        assert "1:" in xpm.legend
        assert "2:" in xpm.legend

    def test_legend_indices_match_yaxis(self, saltbridge_interaction, saltbridge_exporter):
        xpm = saltbridge_exporter.to_xpm_existence(saltbridge_interaction)
        entries = xpm.legend.split(" ")
        for i, entry in enumerate(entries):
            idx = int(entry.split(":")[0])
            assert idx == i

    def test_type_is_discrete(self, saltbridge_interaction, saltbridge_exporter):
        xpm = saltbridge_exporter.to_xpm_existence(saltbridge_interaction)
        assert xpm.type == "Discrete"


# ============================================================
# 4. save 方法
# ============================================================

class TestSaveMethods:

    def test_save_xvg_creates_file(self, saltbridge_interaction, saltbridge_exporter):
        path = str(TEMP_DIR / "test.xvg")
        saltbridge_exporter.save_xvg(
            saltbridge_interaction, "distance", path)
        assert Path(path).exists()

    def test_save_xvg_count_creates_file(self, saltbridge_interaction, saltbridge_exporter):
        path = str(TEMP_DIR / "test_count.xvg")
        saltbridge_exporter.save_xvg_count(saltbridge_interaction, path)
        assert Path(path).exists()

    def test_save_xpm_creates_file(self, saltbridge_interaction, saltbridge_exporter):
        path = str(TEMP_DIR / "test.xpm")
        saltbridge_exporter.save_xpm(saltbridge_interaction, path)
        assert Path(path).exists()

    def test_xvg_roundtrip(self, saltbridge_interaction, saltbridge_exporter):
        """保存后用 DuIvyTools XVG 读回验证。"""
        from DuIvyTools.DuIvyTools.FileParser.xvgParser import XVG
        path = str(TEMP_DIR / "roundtrip.xvg")
        saltbridge_exporter.save_xvg(
            saltbridge_interaction, "distance", path)
        loaded = XVG(path)
        assert loaded.column_num == 4
        assert loaded.row_num == 5

    def test_xpm_roundtrip(self, saltbridge_interaction, saltbridge_exporter):
        """保存后用 DuIvyTools XPM 读回验证。"""
        from DuIvyTools.DuIvyTools.FileParser.xpmParser import XPM
        path = str(TEMP_DIR / "roundtrip.xpm")
        saltbridge_exporter.save_xpm(saltbridge_interaction, path)
        loaded = XPM(path)
        assert loaded.width == 5
        assert loaded.height == 3


# ============================================================
# 5. get_pair_label — 8 个子类
# ============================================================

class TestGetPairLabel:

    def _make_hbond_interaction(self):
        donor = _make_group(0, "H_donor",
                            [_make_atom(100, "NE", "N"),
                             _make_atom(101, "HE", "H")],
                            "ARG", 10)
        acceptor = _make_group(1, "H_acceptor",
                               [_make_atom(200, "OD1", "O")],
                               "ASP", 20)
        return Interaction(
            interaction_type="hydrogen_bond",
            groups=[(donor, acceptor)],
            existence=np.array([[True]]),
            metrics={"distance": np.array([[3.0]]),
                     "angle": np.array([[160.0]])},
            times=np.array([0.0]),
        )

    def test_hydrogen_bond(self):
        it = self._make_hbond_interaction()
        label = HydrogenBondExporter().get_pair_label(it, 0)
        assert "ARG10" in label
        assert "NE(100)" in label
        assert "HE(101)" in label
        assert "ASP20" in label
        assert "OD1(200)" in label
        assert "···" in label

    def test_saltbridge(self):
        pos = _make_group(0, "charged_positive",
                          [_make_atom(10, "NZ", "N"),
                           _make_atom(11, "HZ1", "H")],
                          "LYS", 10)
        neg = _make_group(1, "charged_negative",
                          [_make_atom(20, "OD1", "O"),
                           _make_atom(21, "OD2", "O")],
                          "ASP", 20)
        it = Interaction(
            interaction_type="salt_bridge",
            groups=[(pos, neg)],
            existence=np.array([[True]]),
            metrics={"distance": np.array([[3.0]])},
            times=np.array([0.0]),
        )
        label = SaltBridgeExporter().get_pair_label(it, 0)
        assert "LYS10(10-11)" in label
        assert "ASP20(20-21)" in label

    def test_pi_stacking(self):
        ring1 = _make_group(0, "aromatic_ring",
                            [_make_atom(10, "CG", "C"),
                             _make_atom(11, "CD1", "C")],
                            "PHE", 10)
        ring2 = _make_group(1, "aromatic_ring",
                            [_make_atom(20, "CG", "C"),
                             _make_atom(21, "CD1", "C")],
                            "PHE", 20)
        it = Interaction(
            interaction_type="pi_stacking",
            groups=[(ring1, ring2)],
            existence=np.array([[True]]),
            metrics={"distance": np.array([[4.0]]),
                     "angle": np.array([[10.0]]),
                     "offset": np.array([[1.0]]),
                     "pistacking_type": np.array([["P"]])},
            times=np.array([0.0]),
        )
        label = PiStackingExporter().get_pair_label(it, 0)
        assert "PHE10(10-11)" in label
        assert "PHE20(20-21)" in label

    def test_pi_cation(self):
        ring = _make_group(0, "aromatic_ring",
                           [_make_atom(10, "CG", "C")],
                           "PHE", 10)
        cation = _make_group(1, "charged_positive",
                             [_make_atom(20, "NZ", "N")],
                             "LYS", 20)
        it = Interaction(
            interaction_type="pi_cation",
            groups=[(ring, cation)],
            existence=np.array([[True]]),
            metrics={"distance": np.array([[5.0]]),
                     "offset": np.array([[1.0]])},
            times=np.array([0.0]),
        )
        label = PiCationExporter().get_pair_label(it, 0)
        assert "PHE10(10-10)" in label
        assert "LYS20(20-20)" in label

    def test_halogen_bond(self):
        donor = _make_group(0, "halogen_donor",
                            [_make_atom(100, "C1", "C"),
                             _make_atom(101, "CL1", "Cl")],
                            "D927", 1)
        acceptor = _make_group(1, "halogen_acceptor",
                               [_make_atom(200, "OE1", "O"),
                                _make_atom(201, "OE2", "O")],
                               "GLU", 20)
        it = Interaction(
            interaction_type="halogen_bond",
            groups=[(donor, acceptor)],
            existence=np.array([[True]]),
            metrics={"distance": np.array([[3.5]]),
                     "don_angle": np.array([[170.0]]),
                     "acc_angle": np.array([[120.0]])},
            times=np.array([0.0]),
        )
        label = HalogenBondExporter().get_pair_label(it, 0)
        assert "C1(100)" in label
        assert "CL1(101)" in label
        assert "OE1(200)" in label

    def test_hydrophobic(self):
        g1 = _make_group(0, "hydrophobic",
                         [_make_atom(10, "CB", "C")],
                         "LEU", 10)
        g2 = _make_group(1, "hydrophobic",
                         [_make_atom(20, "CB", "C")],
                         "ILE", 20)
        it = Interaction(
            interaction_type="hydrophobic",
            groups=[(g1, g2)],
            existence=np.array([[True]]),
            metrics={"distance": np.array([[3.8]])},
            times=np.array([0.0]),
        )
        label = HydrophobicExporter().get_pair_label(it, 0)
        assert "LEU10:CB(10)" in label
        assert "ILE20:CB(20)" in label

    def test_metal_coordination(self):
        metal = _make_group(0, "metal",
                            [_make_atom(1, "MG", "Mg")],
                            "MG", 1)
        binding = _make_group(1, "metal_binding",
                              [_make_atom(100, "OD1", "O")],
                              "ASP", 10)
        it = Interaction(
            interaction_type="metal_coordination",
            groups=[(metal, binding)],
            existence=np.array([[True]]),
            metrics={"distance": np.array([[2.0]])},
            times=np.array([0.0]),
        )
        label = MetalCoordinationExporter().get_pair_label(it, 0)
        assert "MG1:MG(1)" in label
        assert "ASP10:OD1(100)" in label

    def test_water_bridge(self):
        donor = _make_group(0, "H_donor",
                            [_make_atom(100, "NE", "N"),
                             _make_atom(101, "HE", "H")],
                            "ARG", 10)
        water = _make_group(1, "water",
                            [_make_atom(500, "OW", "O"),
                             _make_atom(501, "HW1", "H"),
                             _make_atom(502, "HW2", "H")],
                            "SOL", 50)
        acceptor = _make_group(2, "H_acceptor",
                               [_make_atom(200, "OD1", "O")],
                               "ASP", 20)
        it = Interaction(
            interaction_type="water_bridge",
            groups=[(donor, water, acceptor)],
            existence=np.array([[True]]),
            metrics={"dist_dw": np.array([[2.8]]),
                     "dist_wa": np.array([[2.9]]),
                     "theta": np.array([[150.0]]),
                     "omega": np.array([[100.0]])},
            times=np.array([0.0]),
        )
        label = WaterBridgeExporter().get_pair_label(it, 0)
        assert "ARG10:NE(100)-HE(101)" in label
        assert "SOL50:OW(500)" in label
        assert "ASP20:OD1(200)" in label
        assert label.count("···") == 2


# ============================================================
# 6. 用真实数据验证（盐桥 h5）
# ============================================================

class TestRealData:

    @pytest.fixture(scope="class")
    def saltbridge_h5(self):
        h5_path = Path(__file__).parent.parent / "interaction_h5data" / "salt_bridge.h5"
        if not h5_path.exists():
            pytest.skip("salt_bridge.h5 not found")
        return load_interactions(str(h5_path))[0]

    def test_xvg_metric_save(self, saltbridge_h5):
        exporter = SaltBridgeExporter()
        path = str(TEMP_DIR / "real_distance.xvg")
        exporter.save_xvg(saltbridge_h5, "distance", path)
        assert Path(path).exists()

    def test_xvg_count_save(self, saltbridge_h5):
        exporter = SaltBridgeExporter()
        path = str(TEMP_DIR / "real_count.xvg")
        exporter.save_xvg_count(saltbridge_h5, path)
        assert Path(path).exists()

    def test_xpm_save(self, saltbridge_h5):
        exporter = SaltBridgeExporter()
        path = str(TEMP_DIR / "real_existence.xpm")
        exporter.save_xpm(saltbridge_h5, path)
        assert Path(path).exists()

    def test_xpm_roundtrip_alignment(self, saltbridge_h5):
        """验证真实数据的 XPM 对齐。"""
        from DuIvyTools.DuIvyTools.FileParser.xpmParser import XPM
        exporter = SaltBridgeExporter()
        xpm = exporter.to_xpm_existence(saltbridge_h5)

        assert xpm.width == saltbridge_h5.n_frames
        assert xpm.height == saltbridge_h5.n_pairs
        assert len(xpm.xaxis) == xpm.width
        assert len(xpm.yaxis) == xpm.height
        assert len(xpm.value_matrix) == xpm.height
        assert all(len(row) == xpm.width for row in xpm.value_matrix)

    def test_xvg_count_values(self, saltbridge_h5):
        """验证 count 数值等于 existence 按帧求和。"""
        exporter = SaltBridgeExporter()
        xvg = exporter.to_xvg_count(saltbridge_h5)
        expected = np.sum(saltbridge_h5.existence, axis=0).astype(int).tolist()
        np.testing.assert_array_equal(xvg.data_columns[1], expected)


# ============================================================
# 7. to_csv_summary
# ============================================================

class TestCsvSummary:

    def test_creates_file(self, saltbridge_interaction, saltbridge_exporter):
        path = str(TEMP_DIR / "summary.csv")
        saltbridge_exporter.to_csv_summary(saltbridge_interaction, path)
        assert Path(path).exists()

    def test_row_count(self, saltbridge_interaction, saltbridge_exporter):
        import csv as csv_mod
        path = str(TEMP_DIR / "summary.csv")
        saltbridge_exporter.to_csv_summary(saltbridge_interaction, path)
        with open(path) as f:
            reader = csv_mod.reader(f)
            rows = list(reader)
        # n_pairs + 1 header
        assert len(rows) == saltbridge_interaction.n_pairs + 1

    def test_headers(self, saltbridge_interaction, saltbridge_exporter):
        import csv as csv_mod
        path = str(TEMP_DIR / "summary.csv")
        saltbridge_exporter.to_csv_summary(saltbridge_interaction, path)
        with open(path) as f:
            reader = csv_mod.reader(f)
            headers = next(reader)
        assert headers[0] == "pair_label"
        assert headers[1] == "occupancy"
        assert "avg_distance" in headers
        assert "std_distance" in headers

    def test_occupancy_values(self, saltbridge_interaction, saltbridge_exporter):
        import csv as csv_mod
        path = str(TEMP_DIR / "summary.csv")
        saltbridge_exporter.to_csv_summary(saltbridge_interaction, path)
        with open(path) as f:
            reader = csv_mod.reader(f)
            next(reader)
            rows = list(reader)
        # pair0: [T,T,F,T,T] → 4/5 = 0.8
        assert float(rows[0][1]) == pytest.approx(0.8000, abs=0.001)
        # pair2: [T,T,T,T,T] → 5/5 = 1.0
        assert float(rows[2][1]) == pytest.approx(1.0000, abs=0.001)

    def test_avg_only_active_frames(self, saltbridge_interaction, saltbridge_exporter):
        """avg_distance 仅在活跃帧上计算。"""
        import csv as csv_mod
        path = str(TEMP_DIR / "summary.csv")
        saltbridge_exporter.to_csv_summary(saltbridge_interaction, path)
        with open(path) as f:
            reader = csv_mod.reader(f)
            next(reader)
            rows = list(reader)
        # pair0: active frames = [0,1,3,4], distance = [3.5, 3.2, 3.8, 3.1]
        # mean = 3.4
        expected_avg = np.mean([3.5, 3.2, 3.8, 3.1])
        assert float(rows[0][2]) == pytest.approx(expected_avg, abs=0.01)

    def test_std_only_active_frames(self, saltbridge_interaction, saltbridge_exporter):
        """std_distance 仅在活跃帧上计算。"""
        import csv as csv_mod
        path = str(TEMP_DIR / "summary.csv")
        saltbridge_exporter.to_csv_summary(saltbridge_interaction, path)
        with open(path) as f:
            reader = csv_mod.reader(f)
            next(reader)
            rows = list(reader)
        # pair0: active distance = [3.5, 3.2, 3.8, 3.1], std
        expected_std = np.std([3.5, 3.2, 3.8, 3.1])
        assert float(rows[0][3]) == pytest.approx(expected_std, abs=0.01)

    def test_pair_label_format(self, saltbridge_interaction, saltbridge_exporter):
        import csv as csv_mod
        path = str(TEMP_DIR / "summary.csv")
        saltbridge_exporter.to_csv_summary(saltbridge_interaction, path)
        with open(path) as f:
            reader = csv_mod.reader(f)
            next(reader)
            rows = list(reader)
        # pair0: LYS10(10-11)···ASP40(40-41)
        assert "LYS10" in rows[0][0]
        assert "ASP40" in rows[0][0]
        assert "···" in rows[0][0]

    def test_real_data(self):
        """真实数据验证。"""
        import csv as csv_mod
        h5_path = Path(__file__).parent.parent / "interaction_h5data" / "salt_bridge.h5"
        if not h5_path.exists():
            pytest.skip("salt_bridge.h5 not found")
        it = load_interactions(str(h5_path))[0]
        exporter = SaltBridgeExporter()
        path = str(TEMP_DIR / "real_summary.csv")
        exporter.to_csv_summary(it, path)
        with open(path) as f:
            reader = csv_mod.reader(f)
            rows = list(reader)
        assert len(rows) == it.n_pairs + 1
        # 第一个 pair 占位率应为 1.0
        assert float(rows[1][1]) == pytest.approx(1.0, abs=0.001)


# ============================================================
# 8. XPM 手动构建（对抗性测试：原 refresh 重映射 bug 回归防护）
# ============================================================

class TestDiscreteXpmBuild:
    """对抗性测试：XPM 手动构建（原 refresh 重映射 bug 回归防护）。"""

    def _mk_saltbridge_it(self, existence, times):
        pos = _make_group(0, "charged_positive",
                          [_make_atom(10, "NZ", "N")], "LYS", 10)
        neg = _make_group(1, "charged_negative",
                          [_make_atom(20, "OD1", "O")], "ASP", 20)
        n_pairs, n_frames = existence.shape
        dist = np.where(existence, 3.0, 9.0).astype(float)
        return Interaction("salt_bridge", [(pos, neg)] * n_pairs,
                           existence=existence,
                           metrics={"distance": dist},
                           times=times)

    def test_all_active_existence_not_remapped(self):
        """全活跃矩阵保持 [1,1,1]（原 bug：被 refresh 重映射为全 0）。"""
        it = self._mk_saltbridge_it(
            existence=np.array([[True, True, True]]),
            times=np.array([0.0, 10.0, 20.0]))
        xpm = SaltBridgeExporter().to_xpm_existence(it)
        assert np.array(xpm.value_matrix).tolist() == [[1, 1, 1]]
        assert xpm.colors == ["#FFFFFF", "#38A7D0"]
        assert xpm.notes == ["No", "Yes"]

    def test_stacking_all_active_types_not_shifted(self):
        """[T,P] 全活跃保持 [[1,2]]（原 bug：P 被重映射成 1 显示为 T）。"""
        ring1 = _make_group(0, "aromatic_ring",
                            [_make_atom(10, "CG", "C")], "PHE", 10)
        ring2 = _make_group(1, "aromatic_ring",
                            [_make_atom(20, "CG", "C")], "TYR", 20)
        it = Interaction("pi_stacking", [(ring1, ring2)],
                         existence=np.array([[True, True]]),
                         metrics={"distance": np.array([[4.0, 4.1]]),
                                  "pistacking_type": np.array([["T", "P"]])},
                         times=np.array([0.0, 10.0]))
        xpm = PiStackingExporter().to_xpm_stacking_type(it)
        assert np.array(xpm.value_matrix).tolist() == [[1, 2]]
        assert xpm.notes == ["None", "T-shaped", "Parallel"]
        assert xpm.colors == ["#FFFFFF", "#F67088", "#38A7D0"]

    def test_partial_existence_unchanged(self):
        """部分活跃（值集合完整）输出与旧版一致（回归）。"""
        it = self._mk_saltbridge_it(
            existence=np.array([[True, False, True]]),
            times=np.array([0.0, 10.0, 20.0]))
        xpm = SaltBridgeExporter().to_xpm_existence(it)
        assert np.array(xpm.value_matrix).tolist() == [[1, 0, 1]]

    def test_roundtrip_all_active(self):
        """全活跃 XPM 保存后读回一致。"""
        from DuIvyTools.DuIvyTools.FileParser.xpmParser import XPM
        it = self._mk_saltbridge_it(
            existence=np.array([[True, True, True]]),
            times=np.array([0.0, 10.0, 20.0]))
        xpm = SaltBridgeExporter().to_xpm_existence(it)
        path = str(TEMP_DIR / "all_active.xpm")
        xpm.save(path)
        loaded = XPM(path)
        assert np.array(loaded.value_matrix).tolist() == [[1, 1, 1]]
        assert loaded.colors == ["#FFFFFF", "#38A7D0"]
        assert loaded.notes == ["No", "Yes"]

    def test_negative_index_rejected(self):
        """负索引拒绝（原 refresh 静默错乱）。"""
        with pytest.raises(ValueError):
            SaltBridgeExporter()._build_discrete_xpm(
                value_matrix=np.array([[-1, 1]]),
                colors=["#FFFFFF", "#38A7D0"], notes=["No", "Yes"],
                title="t", legend="l", xlabel="x", ylabel="y",
                times=np.array([0.0, 1.0]))

    def test_float_index_rejected(self):
        """float dtype 拒绝。"""
        with pytest.raises(TypeError):
            SaltBridgeExporter()._build_discrete_xpm(
                value_matrix=np.array([[1.0, 0.0]]),
                colors=["#FFFFFF", "#38A7D0"], notes=["No", "Yes"],
                title="t", legend="l", xlabel="x", ylabel="y",
                times=np.array([0.0, 1.0]))

    def test_index_out_of_range_rejected(self):
        """越界正索引拒绝。"""
        with pytest.raises(ValueError):
            SaltBridgeExporter()._build_discrete_xpm(
                value_matrix=np.array([[0, 2]]),
                colors=["#FFFFFF", "#38A7D0"], notes=["No", "Yes"],
                title="t", legend="l", xlabel="x", ylabel="y",
                times=np.array([0.0, 1.0]))

    def test_too_many_colors_rejected(self):
        """超过字符表容量（82）拒绝，不静默截断。"""
        with pytest.raises(ValueError):
            SaltBridgeExporter()._build_discrete_xpm(
                value_matrix=np.zeros((1, 1), dtype=int),
                colors=[f"#{i:06x}" for i in range(90)],
                notes=[str(i) for i in range(90)],
                title="t", legend="l", xlabel="x", ylabel="y",
                times=np.array([0.0]))

    def test_colors_notes_length_mismatch_rejected(self):
        """colors 与 notes 长度不一致拒绝。"""
        with pytest.raises(ValueError):
            SaltBridgeExporter()._build_discrete_xpm(
                value_matrix=np.array([[0, 1]]),
                colors=["#FFFFFF", "#38A7D0"], notes=["No"],
                title="t", legend="l", xlabel="x", ylabel="y",
                times=np.array([0.0, 1.0]))


# ============================================================
# 9. DII export 对抗性测试（空数据/多类型/重复类型）
# ============================================================

class TestDIIExportAdversarial:
    """dii export 边界场景（n_pairs=0、空帧、多类型、重复类型）。"""

    def _save_and_export(self, interactions, filename):
        """保存 h5 并运行 _run_export，返回输出目录。"""
        from DuIvyInteractions.io.h5 import save_interactions
        from DuIvyInteractions.DII import _run_export
        import argparse
        h5_path = TEMP_DIR / filename
        out_dir = TEMP_DIR / (filename.replace(".h5", "_out"))
        save_interactions(interactions, str(h5_path))
        args = argparse.Namespace(input=str(h5_path), output=str(out_dir))
        _run_export(args)
        return out_dir

    def _make_saltbridge(self, n_pairs=2, n_frames=5, exist=True):
        """构造盐桥 Interaction。"""
        from DuIvyInteractions.io.saltbridge_exporter import SaltBridgeExporter
        pos = _make_group(0, "charged_positive",
                          [_make_atom(10, "NZ", "N"), _make_atom(11, "HZ1", "H")],
                          "LYS", 10)
        neg = _make_group(1, "charged_negative",
                          [_make_atom(20, "OD1", "O"), _make_atom(21, "OD2", "O")],
                          "ASP", 20)
        existence = np.full((n_pairs, n_frames), exist, dtype=bool)
        metrics = {"distance": np.full((n_pairs, n_frames), 3.5)}
        groups = [(pos, neg)] * n_pairs
        return Interaction(
            interaction_type="salt_bridge", groups=groups,
            existence=existence, metrics=metrics,
            times=np.arange(n_frames, dtype=float) * 10.0)

    def test_n_pairs_zero_skipped(self):
        """n_pairs=0 不崩溃、不产生文件。"""
        it = self._make_saltbridge(n_pairs=0)
        out = self._save_and_export([it], "empty_pairs.h5")
        files = list(out.glob("*")) if out.exists() else []
        assert files == []

    def test_empty_frames_overview(self):
        """空帧概览不 IndexError（先有 0 帧 → 跳过导出）。"""
        it = self._make_saltbridge(n_pairs=0, n_frames=0)
        out = self._save_and_export([it], "empty_frames.h5")
        files = list(out.glob("*")) if out.exists() else []
        assert files == []

    def test_multiple_types_all_exported(self):
        """多类型 h5：所有类型都导出，不丢弃。"""
        from DuIvyInteractions.io.pi_stacking_exporter import PiStackingExporter
        sb = self._make_saltbridge(n_pairs=1)
        # π-stacking Interaction
        ring = _make_group(0, "aromatic_ring", [_make_atom(10, "CG", "C")], "PHE", 10)
        ring2 = _make_group(1, "aromatic_ring", [_make_atom(20, "CG", "C")], "PHE", 20)
        pi = Interaction(
            interaction_type="pi_stacking", groups=[(ring, ring2)],
            existence=np.ones((1, 5), dtype=bool),
            metrics={"distance": np.full((1, 5), 4.0),
                     "angle": np.full((1, 5), 10.0),
                     "offset": np.full((1, 5), 1.0),
                     "pistacking_type": np.full((1, 5), "P", dtype="U1")},
            times=np.arange(5, dtype=float) * 10.0)
        out = self._save_and_export([sb, pi], "multi.h5")
        names = [p.name for p in out.glob("*")]
        assert any("salt_bridge" in n for n in names)
        assert any("pi_stacking" in n for n in names)

    def test_same_type_repeated_gets_index(self):
        """同类型重复：第二个带 _2 序号，不覆盖。"""
        sb1 = self._make_saltbridge(n_pairs=1)
        sb2 = self._make_saltbridge(n_pairs=1)
        out = self._save_and_export([sb1, sb2], "repeat.h5")
        names = sorted(p.name for p in out.glob("*"))
        assert any("_2_" in n for n in names)
        # 两个 summary 都在
        summaries = [n for n in names if n.endswith("_summary.csv")]
        assert len(summaries) == 2
