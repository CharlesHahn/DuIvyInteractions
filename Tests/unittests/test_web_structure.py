# -*- coding: utf-8 -*-
"""3D 结构可视化测试：structure 纯函数层 + /api/structure/* 路由。

不依赖真实轨迹：用 tmp_path 构造最小 PDB 与 h5。
"""

import collections
import re
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pytest

from DuIvyInteractions.core.datas import AtomData, Group, Interaction
from DuIvyInteractions.io.h5 import load_interactions, save_interactions
from DuIvyInteractions.visualizers import get_exporter
from DuIvyInteractions.visualizers.structure import (
    SAMPLE_ATOMS, STRUCTURE_CACHE_SIZE, WATER_RESIDUES, alignment_errors,
    contact_label, group_atoms, group_label, interaction_expressions,
    involves_water, load_universe, residue_chains, topology_cif,
)
from DuIvyInteractions.web.picker import FILE_TYPES, pick_file

TOKEN = "s3cr3t-token"
BOX = (10.0, 10.0, 10.0)
ATOMS = [("N1", (0.0, 0.0, 0.0), "N"),
         ("H1", (0.0, 0.0, 1.0), "H"),
         ("O1", (0.0, 0.0, 2.8), "O"),
         ("O2", (0.0, 0.0, 4.0), "O")]
WATER_ATOMS = [("OW", (9.0, 9.0, 9.0), "O", "SOL", 2),
               ("HW1", (9.4, 9.0, 9.0), "H", "SOL", 2),
               ("N1", (0.0, 0.0, 0.0), "N", "ALA", 1),
               ("H1", (0.0, 0.0, 1.0), "H", "ALA", 1),
               ("O1", (0.0, 0.0, 2.8), "O", "ALA", 1),
               ("O2", (0.0, 0.0, 4.0), "O", "ALA", 1)]
# 水桥三组：SOL2 在前（组 1 的原子索引 0/1）、ALA1 居中（组 0 的 2/3）、LYS3 在后
BRIDGE_ATOMS = [("OW", (9.0, 9.0, 9.0), "O", "SOL", 2),
                ("HW1", (9.4, 9.0, 9.0), "H", "SOL", 2),
                ("N1", (0.0, 0.0, 0.0), "N", "ALA", 1),
                ("H1", (0.0, 0.0, 1.0), "H", "ALA", 1),
                ("O1", (0.0, 0.0, 2.8), "O", "LYS", 3),
                ("O2", (0.0, 0.0, 4.0), "O", "LYS", 3)]
# 蛋白 / 多原子配体 / 水 / 单原子金属离子 / 单原子卤素离子
MIXED_ATOMS = [("N1", (0.0, 0.0, 0.0), "N", "ALA", 1),
               ("CA", (0.0, 0.0, 1.0), "C", "ALA", 1),
               ("C1", (1.0, 0.0, 0.0), "C", "MOL", 5),
               ("C2", (2.0, 0.0, 0.0), "C", "MOL", 5),
               ("OW", (3.0, 0.0, 0.0), "O", "SOL", 7),
               ("HW1", (3.4, 0.0, 0.0), "H", "SOL", 7),
               ("HW2", (3.0, 0.4, 0.0), "H", "SOL", 7),
               ("NA", (4.0, 0.0, 0.0), "NA", "NA", 9),
               ("CL", (5.0, 0.0, 0.0), "CL", "CL", 10)]
# 两个 segment（PDB chainID 列不同），各 1 残基且 resid 相同 → chain 字母必须不同
SEGMENT_ATOMS = [("N1", (0.0, 0.0, 0.0), "N", "ALA", 51, "A"),
                 ("N1", (1.0, 0.0, 0.0), "N", "GLY", 51, "B")]
# h5 的 residue_id 是 0 基索引，universe 的 resid 是 1 基 auth 编号（差 1）。
# 蛋白段与配体段分属两个 segment（PDB chainID A / B）。
NUMBERED_ATOMS = [("CZ", (0.0, 0.0, 0.0), "C", "TYR", 51, "A"),
                  ("CE1", (4.0, 0.0, 2.0), "C", "TYR", 51, "A"),
                  ("N16", (0.0, 0.0, 1.2), "N", "MOL", 143, "B"),
                  ("N17", (0.0, 0.0, 5.0), "N", "MOL", 143, "B")]
D927_DIR = Path(__file__).resolve().parents[1] / "test_MD_case_amber"
D927_TOPOLOGY = D927_DIR / "md.tpr"
D927_TRAJECTORY = D927_DIR / "md1ns.xtc"
H5_DATA_DIR = Path(__file__).resolve().parents[1] / "interaction_h5data"


def _cryst1(box) -> str:
    return ("CRYST1%9.3f%9.3f%9.3f  90.00  90.00  90.00 P 1           1"
            % box)


def _pdb_atom_lines(atoms) -> list:
    return [f"ATOM  {serial:5d} {name:>4s} ALA A   1    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}"
            for serial, (name, (x, y, z), element)
            in enumerate(atoms, start=1)]


def _chain_atom_lines(atoms) -> list:
    """atoms: (name, xyz, element, resname, resid, chain) → ATOM 行。"""
    return [f"ATOM  {serial:5d} {name:>4s} {resname:>3s} {chain}{resid:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}"
            for serial, (name, (x, y, z), element, resname, resid, chain)
            in enumerate(atoms, start=1)]


def _residue_atom_lines(atoms) -> list:
    """(name, xyz, element, resname, resid) → chain A 的 ATOM 行。"""
    return _chain_atom_lines([tuple(atom) + ("A",) for atom in atoms])


def _pdb_text(atoms, box=BOX) -> str:
    lines = [_cryst1(box)] if box else []
    return "\n".join(lines + _pdb_atom_lines(atoms) + ["END"]) + "\n"


def _write_pdb(path: Path, atoms=None, box=BOX) -> Path:
    path.write_text(_pdb_text(atoms or ATOMS, box), encoding="utf-8")
    return path


def _chain_pdb_text(atoms, box=BOX) -> str:
    lines = [_cryst1(box)] if box else []
    return "\n".join(lines + _chain_atom_lines(atoms) + ["END"]) + "\n"


def _write_chain_pdb(path: Path, atoms, box=BOX) -> Path:
    path.write_text(_chain_pdb_text(atoms, box), encoding="utf-8")
    return path


def _water_pdb_text(atoms, box=BOX) -> str:
    lines = [_cryst1(box)] if box else []
    return "\n".join(lines + _residue_atom_lines(atoms) + ["END"]) + "\n"


def _write_water_pdb(path: Path, atoms=None, box=BOX) -> Path:
    path.write_text(_water_pdb_text(atoms or WATER_ATOMS, box),
                    encoding="utf-8")
    return path


def _write_multi_pdb(path: Path, atoms=None, frames: int = 2,
                     box=BOX) -> Path:
    head = [_cryst1(box)] if box else []
    body = _pdb_atom_lines(atoms or ATOMS)
    models = []
    for frame in range(frames):
        models += [f"MODEL     {frame + 1:4d}"] + body + ["ENDMDL"]
    path.write_text("\n".join(head + models + ["END"]) + "\n",
                    encoding="utf-8")
    return path


def _atom(index: int, name: str, atype: str = "") -> AtomData:
    return AtomData(atom_global_idx=index, atom_idx_in_residue=0,
                    atom_name=name, atom_type=atype, atom_element="C",
                    atom_charge=0.0, atom_mass=0.0)


def _group(group_id: int, group_type: str, atoms,
           residue_name: str = "LYS", residue_id: int = 1,
           molecule: str = "test") -> Group:
    return Group(group_id=group_id, group_type=group_type, molecule=molecule,
                 residue_name=residue_name, residue_id=residue_id,
                 atoms=list(atoms))


def _single(interaction_type: str, groups) -> Interaction:
    existence = np.zeros((1, 1), dtype=bool)
    existence[0, 0] = True
    return Interaction(interaction_type=interaction_type, groups=list(groups),
                       existence=existence,
                       metrics={"distance": np.full((1, 1), 3.0)},
                       times=np.zeros(1))


def _protein_interaction() -> Interaction:
    groups = [(_group(0, "H_donor", [_atom(2, "N1", "N"),
                                     _atom(3, "H1", "H")], "ALA", 1),
               _group(1, "H_acceptor", [_atom(4, "O1", "O"),
                                        _atom(5, "O2", "O")], "ALA", 1))]
    return _single("hydrogen_bond", groups)


def _water_bridge_interaction() -> Interaction:
    groups = [(_group(0, "H_donor", [_atom(2, "N1", "N"),
                                     _atom(3, "H1", "H")], "ALA", 1),
               _group(1, "water", [_atom(0, "OW", "O"),
                                   _atom(1, "HW1", "H")], "SOL", 2),
               _group(2, "H_acceptor", [_atom(4, "O1", "O"),
                                        _atom(5, "O2", "O")], "LYS", 3))]
    return _single("water_bridge", groups)


def _numbered_interaction() -> Interaction:
    """h5 残基编号（50/142）与 NUMBERED_ATOMS 的 auth 编号（51/143）相差 1。"""
    groups = [(_group(0, "aromatic_ring", [_atom(0, "CZ", "C"),
                                           _atom(1, "CE1", "C")], "TYR", 50),
               _group(1, "aromatic_ring", [_atom(2, "N16", "N"),
                                           _atom(3, "N17", "N")], "MOL", 142))]
    return _single("pi_stacking", groups)


def _interaction(active: bool = True) -> Interaction:
    existence = np.zeros((1, 1), dtype=bool)
    existence[0, 0] = active
    groups = [(_group(0, "H_donor", [_atom(0, "N1", "N"),
                                     _atom(1, "H1", "H")]),
               _group(1, "H_acceptor", [_atom(2, "O1", "O"),
                                        _atom(3, "O2", "O")]))]
    return Interaction(interaction_type="hydrogen_bond", groups=groups,
                       existence=existence,
                       metrics={"distance": np.full((1, 1), 3.0)},
                       times=np.zeros(1))


def _segment_interaction() -> Interaction:
    """两 segment 各 1 残基、resid 同为 51（SEGMENT_ATOMS 顺序）。"""
    groups = [(_group(0, "aromatic_ring", [_atom(0, "N1", "N")], "ALA", 51),
               _group(1, "aromatic_ring", [_atom(1, "N1", "N")], "GLY", 51))]
    return _single("pi_stacking", groups)


def _residue_selectors(universe, residues) -> list:
    """用真实残基构造最小 interaction，读回各组的 residue_selector。"""
    groups = tuple(_group(index, "aromatic_ring",
                          [_atom(int(res.atoms[0].index),
                                 str(res.atoms[0].name))],
                          str(res.resname), int(res.resid) - 1)
                   for index, res in enumerate(residues))
    rows = interaction_expressions(_single("pi_stacking", [groups]),
                                   0, universe)["residues"]
    return [row["residue_selector"] for row in rows]


@pytest.fixture(autouse=True)
def isolate_recent(monkeypatch, tmp_path):
    monkeypatch.setattr("DuIvyInteractions.visualizers.data.RECENT_PATH",
                        tmp_path / "recent.json")


@pytest.fixture
def pdb_file(tmp_path):
    return _write_pdb(tmp_path / "mini.pdb")


@pytest.fixture
def universe(pdb_file):
    return load_universe(str(pdb_file), str(pdb_file))


@pytest.fixture
def h5_dir(tmp_path):
    directory = tmp_path / "out"
    directory.mkdir()
    save_interactions([_interaction()], str(directory / "hydrogen_bond.h5"))
    return directory


def _client(**kwargs):
    from fastapi.testclient import TestClient
    from DuIvyInteractions.web.app import create_app
    return TestClient(create_app(**kwargs))


@pytest.fixture
def client(h5_dir):
    client = _client()
    client.post("/select", data={"path": str(h5_dir)})
    return client


def _open_structure(client, pdb_file):
    return client.post("/api/structure/select",
                       data={"topology": str(pdb_file),
                             "trajectory": str(pdb_file)})


def _atom_lines(text: str) -> list:
    """_atom_site 数据行解析为 dict（列名去掉 _atom_site. 前缀）。"""
    columns = [line.split(".", 1)[1] for line in text.splitlines()
               if line.startswith("_atom_site.")]
    return [dict(zip(columns, line.split()))
            for line in text.splitlines() if line.startswith("ATOM ")]


def _chain_rows(text: str) -> dict:
    """(label_comp_id, auth_seq_id) → label_asym_id（拓扑 mmCIF 解析）。"""
    return {(row["label_comp_id"], int(row["auth_seq_id"])):
            row["label_asym_id"] for row in _atom_lines(text)}


def _b_factor_rows(text: str) -> list:
    """(label_atom_id, B_iso_or_equiv) 列表（按行序）。"""
    return [(row["label_atom_id"], float(row["B_iso_or_equiv"]))
            for row in _atom_lines(text)]


def _b_band(expression: str) -> tuple:
    """从 'b > 99.75 and b < 100.25' 取 (下界, 上界)。"""
    low, high = re.findall(r"-?\d+\.\d+", expression)
    return float(low), float(high)


@pytest.fixture
def water_pdb(tmp_path):
    return _write_water_pdb(tmp_path / "water.pdb")


@pytest.fixture
def water_universe(water_pdb):
    return load_universe(str(water_pdb), str(water_pdb))


@pytest.fixture
def water_client(tmp_path):
    directory = tmp_path / "water_out"
    directory.mkdir()
    save_interactions([_protein_interaction()],
                      str(directory / "hydrogen_bond.h5"))
    client = _client()
    client.post("/select", data={"path": str(directory)})
    return client


@pytest.fixture
def water_bridge_client(tmp_path):
    directory = tmp_path / "water_bridge_out"
    directory.mkdir()
    save_interactions([_water_bridge_interaction()],
                      str(directory / "water_bridge.h5"))
    client = _client()
    client.post("/select", data={"path": str(directory)})
    return client


@pytest.fixture
def bridge_pdb(tmp_path):
    return _write_water_pdb(tmp_path / "bridge.pdb", BRIDGE_ATOMS)


@pytest.fixture
def bridge_universe(bridge_pdb):
    return load_universe(str(bridge_pdb), str(bridge_pdb))


@pytest.fixture
def numbered_pdb(tmp_path):
    return _write_chain_pdb(tmp_path / "numbered.pdb", NUMBERED_ATOMS)


@pytest.fixture
def numbered_universe(numbered_pdb):
    return load_universe(str(numbered_pdb), str(numbered_pdb))


@pytest.fixture
def numbered_client(tmp_path):
    directory = tmp_path / "numbered_out"
    directory.mkdir()
    save_interactions([_numbered_interaction()],
                      str(directory / "pi_stacking.h5"))
    client = _client()
    client.post("/select", data={"path": str(directory)})
    return client


@pytest.fixture
def mixed_pdb(tmp_path):
    return _write_water_pdb(tmp_path / "mixed.pdb", MIXED_ATOMS)


@pytest.fixture
def mixed_universe(mixed_pdb):
    return load_universe(str(mixed_pdb), str(mixed_pdb))


@pytest.fixture
def multi_pdb(tmp_path):
    return _write_multi_pdb(tmp_path / "multi.pdb", frames=2)


@pytest.fixture
def multi_universe(multi_pdb):
    return load_universe(str(multi_pdb), str(multi_pdb))


@pytest.fixture
def segment_pdb(tmp_path):
    return _write_chain_pdb(tmp_path / "segments.pdb", SEGMENT_ATOMS)


@pytest.fixture
def segment_universe(segment_pdb):
    return load_universe(str(segment_pdb), str(segment_pdb))


@pytest.fixture(scope="module")
def d927_universe():
    """真实 KRAS-RBD D927 体系（8 个 segment），只做冒烟级断言。"""
    return load_universe(str(D927_TOPOLOGY), str(D927_TRAJECTORY))


@pytest.fixture(scope="module")
def d927_topology_text(d927_universe):
    return topology_cif(None, None, d927_universe)


# ============================================================
# structure 纯函数层
# ============================================================

class TestConstants:

    def test_cache_and_sample_sizes(self):
        assert STRUCTURE_CACHE_SIZE == 2
        assert SAMPLE_ATOMS == 5

    def test_water_residue_constant_covers_common_names(self):
        assert {"SOL", "HOH", "WAT", "TIP3", "TIP4", "HO4", "HO5",
                "H2O"} <= WATER_RESIDUES


class TestLoadUniverse:

    def test_opens_pdb_topology_and_trajectory(self, pdb_file):
        universe = load_universe(str(pdb_file), str(pdb_file))
        assert universe.atoms.n_atoms == len(ATOMS)
        assert len(universe.trajectory) == 1

    def test_multi_model_pdb_keeps_every_frame(self, multi_pdb):
        universe = load_universe(str(multi_pdb), str(multi_pdb))
        assert universe.atoms.n_atoms == len(ATOMS)
        assert len(universe.trajectory) == 2

    def test_bad_file_raises_value_error(self, tmp_path):
        bad = tmp_path / "bad.txt"
        bad.write_text("not a structure")
        with pytest.raises(ValueError, match="cannot open structure"):
            load_universe(str(bad), str(bad))

    def test_missing_file_raises_value_error(self, tmp_path):
        with pytest.raises(ValueError, match="cannot open structure"):
            load_universe(str(tmp_path / "nope.tpr"), str(tmp_path / "nope.xtc"))


class TestAlignmentErrors:

    def test_consistent_interaction_has_no_errors(self, universe):
        assert alignment_errors(universe, _interaction()) == []

    def test_too_few_atoms_reported(self, universe):
        interaction = _interaction()
        interaction.groups[0][1].atoms.append(_atom(9, "X", "C"))
        errors = alignment_errors(universe, interaction)
        assert errors and "needs 10" in errors[0]

    def test_name_mismatch_reported(self, universe):
        interaction = _interaction()
        interaction.groups[0][0].atoms[0].atom_name = "ZZ"
        errors = alignment_errors(universe, interaction)
        assert errors and "name" in errors[0]

    def test_type_mismatch_reported(self, universe):
        interaction = _interaction()
        interaction.groups[0][0].atoms[0].atom_type = "C"
        errors = alignment_errors(universe, interaction)
        assert errors and "type" in errors[0]

    def test_only_first_sample_atoms_checked(self, tmp_path):
        atoms = [(f"C{i}", (float(i), 0.0, 0.0), "C") for i in range(10)]
        path = _write_pdb(tmp_path / "many.pdb", atoms)
        universe = load_universe(str(path), str(path))
        groups = [(_group(0, "H_donor",
                          [_atom(i, f"C{i}", "C") for i in range(5)] +
                          [_atom(5, "WRONG", "C")]),
                   _group(1, "H_acceptor", [_atom(6, "C6", "C")]))]
        interaction = Interaction(
            interaction_type="hydrogen_bond", groups=groups,
            existence=np.ones((1, 1), dtype=bool),
            metrics={"distance": np.full((1, 1), 3.0)},
            times=np.zeros(1))
        assert alignment_errors(universe, interaction) == []


class TestTopologyCif:
    """第 0 帧全原子拓扑 mmCIF：行数/chain 规则/编号/实体分类。"""

    def test_one_atom_line_per_atom(self, universe):
        rows = _atom_lines(topology_cif(None, None, universe))
        assert len(rows) == universe.atoms.n_atoms
        assert [row["label_atom_id"] for row in rows] == [
            str(name) for name in universe.atoms.names]

    def test_cif_block_header_and_columns(self, universe):
        lines = topology_cif(None, None, universe).splitlines()
        assert lines[0] == "data_dii"
        assert lines[2] == "loop_"
        columns = [line for line in lines if line.startswith("_atom_site.")]
        assert len(columns) == 20
        assert "_atom_site.label_alt_id" in columns
        assert "_atom_site.auth_seq_id" in columns
        assert "_atom_site.B_iso_or_equiv" in columns
        widths = {len(row) for row in _atom_lines("\n".join(lines))}
        assert widths == {20}

    def test_chain_rules_for_protein_water_and_ions(self, mixed_universe):
        rows = _chain_rows(topology_cif(None, None, mixed_universe))
        assert rows[("ALA", 1)] == "A"
        assert rows[("MOL", 5)] == "A"  # 同一 segment → 同一 chain 字母
        assert rows[("SOL", 7)] == "W"
        assert rows[("NA", 9)] == "I"
        assert rows[("CL", 10)] == "I"

    def test_segments_get_distinct_chain_letters(self, segment_universe):
        rows = _chain_rows(topology_cif(None, None, segment_universe))
        assert rows[("ALA", 51)] == "A"
        assert rows[("GLY", 51)] == "B"
        rows = interaction_expressions(_segment_interaction(), 0,
                                       segment_universe)["residues"]
        selectors = [row["residue_selector"] for row in rows]
        assert selectors == ["chain A and resi 51", "chain B and resi 51"]

    def test_too_many_segments_raise_value_error(self):
        universe = mda.Universe.empty(
            25, n_residues=25, n_segments=25,
            atom_resindex=np.arange(25), residue_segindex=np.arange(25))
        universe.add_TopologyAttr("resnames", ["ALA"] * 25)
        with pytest.raises(ValueError,
                           match="too many segments for chain letters"):
            residue_chains(universe)

    def test_resseq_is_auth_numbering(self, numbered_universe):
        rows = _atom_lines(topology_cif(None, None, numbered_universe))
        assert {int(row["auth_seq_id"]) for row in rows
                if row["label_comp_id"] == "TYR"} == {51}
        assert {int(row["auth_seq_id"]) for row in rows
                if row["label_comp_id"] == "MOL"} == {143}

    def test_long_water_number_is_not_truncated(self, d927_topology_text):
        """D2 关键：水 resid 15689 原样写入，且恰好 3 个原子行。"""
        rows = [row for row in _atom_lines(d927_topology_text)
                if row["label_comp_id"] == "SOL"
                and row["auth_seq_id"] == "15689"]
        assert len(rows) == 3
        assert [row["label_atom_id"] for row in rows] == ["OW", "HW1", "HW2"]

    def test_chain_sequence_pairs_are_unique(self, d927_universe,
                                             d927_topology_text):
        """全结构 (auth_asym_id, auth_seq_id) 唯一：0 冲突。"""
        rows = _atom_lines(d927_topology_text)
        pairs = {(row["auth_asym_id"], row["auth_seq_id"]) for row in rows}
        assert len(pairs) == len(d927_universe.atoms.residues)

    def test_cif_has_no_conect_and_ends_with_hash(self, universe):
        text = topology_cif(None, None, universe)
        assert "CONECT" not in text
        assert text.rstrip().endswith("#")

    def test_coordinates_are_wrapped_for_the_box(self, tmp_path):
        path = _write_pdb(tmp_path / "wrapped.pdb",
                          [("C1", (12.0, 2.0, 3.0), "C")])
        universe = load_universe(str(path), str(path))
        row = _atom_lines(topology_cif(None, None, universe))[0]
        assert float(row["Cartn_x"]) == pytest.approx(2.0)
        assert float(row["Cartn_y"]) == pytest.approx(2.0)

    def test_without_box_keeps_coordinates(self, tmp_path):
        path = _write_pdb(tmp_path / "nobox.pdb",
                          [("C1", (12.0, 2.0, 3.0), "C")], box=None)
        universe = load_universe(str(path), str(path))
        row = _atom_lines(topology_cif(None, None, universe))[0]
        assert float(row["Cartn_x"]) == pytest.approx(12.0)

    def test_cif_lists_expected_chains(self, mixed_universe):
        rows = _atom_lines(topology_cif(None, None, mixed_universe))
        assert len(rows) == mixed_universe.atoms.n_atoms
        assert sorted({row["label_asym_id"] for row in rows}) == [
            "A", "I", "W"]

    def test_entity_per_chain_and_class(self, mixed_universe):
        """蛋白/配体同链也分属不同实体；水与离子按类别分实体。"""
        rows = _atom_lines(topology_cif(None, None, mixed_universe))
        entities = {}
        for row in rows:
            entities.setdefault((row["label_asym_id"],
                                 row["label_comp_id"]),
                                row["label_entity_id"])
        assert entities[("A", "ALA")] != entities[("A", "MOL")]
        assert entities[("W", "SOL")] != entities[("A", "ALA")]
        assert entities[("I", "NA")] == entities[("I", "CL")]

    def test_tagged_atoms_get_unique_b_factors(self, numbered_universe):
        interaction = _numbered_interaction()
        rows = _b_factor_rows(topology_cif(interaction, 0, numbered_universe))
        tagged = [value for _, value in rows if value != 0.0]
        assert len(tagged) == len(
            [atom for group in interaction.groups[0] for atom in group.atoms])
        assert tagged == [100.0, 101.0, 102.0, 103.0]

    def test_untagged_atoms_keep_zero_b_factor(self, water_universe):
        rows = dict(_b_factor_rows(
            topology_cif(_protein_interaction(), 0, water_universe)))
        assert rows["OW"] == 0.0 and rows["HW1"] == 0.0
        assert [rows[name] for name in ("N1", "H1", "O1", "O2")] == [
            100.0, 101.0, 102.0, 103.0]

    def test_tag_order_matches_expression_order(self, numbered_universe):
        interaction = _numbered_interaction()
        tagged = [name for name, value
                  in _b_factor_rows(topology_cif(interaction, 0,
                                                 numbered_universe))
                  if value != 0.0]
        residues = interaction_expressions(
            interaction, 0, numbered_universe)["residues"]
        assert tagged == [row["name"] for group in residues
                          for row in group["atoms"]]

    def test_atom_id_is_global_index_plus_one(self, numbered_universe):
        rows = _atom_lines(topology_cif(None, None, numbered_universe))
        assert [int(row["id"]) for row in rows] == [1, 2, 3, 4]

    def test_large_atom_ids_are_not_truncated(self, d927_topology_text):
        rows = _atom_lines(d927_topology_text)
        assert rows[-1]["id"] == "116383"


class TestD927SegmentChains:
    """真实 D927（8 个 segment）的 chain 分配与蛋白残基选择器唯一性。"""

    SEGMENT_CHAINS = {"seg_0_RBD_pro": "A", "seg_1_D927": "B",
                      "seg_2_KRAS_pro": "C", "seg_3_GNP_neg": "D",
                      "seg_4_Mg": "I", "seg_5_SOL": "W",
                      "seg_6_NA": "I", "seg_7_CL": "I"}

    def test_one_chain_letter_per_segment(self, d927_universe):
        chains = residue_chains(d927_universe)
        for segment in d927_universe.atoms.segments:
            letters = {chains[int(index)]
                       for index in segment.residues.resindices}
            assert letters == {self.SEGMENT_CHAINS[segment.segid]}
        protein = {self.SEGMENT_CHAINS[str(segment.segid)]: len(segment.atoms)
                   for segment in d927_universe.atoms.segments
                   if str(segment.segid).endswith("_pro")}
        atom_chains = collections.Counter(
            chains[int(index)] for index in d927_universe.atoms.resindices)
        assert protein == {"A": atom_chains["A"], "C": atom_chains["C"]}

    def test_topology_chain_atom_counts_match_segments(self, d927_universe,
                                                       d927_topology_text):
        segment_atoms = {segment.segid: len(segment.atoms)
                         for segment in d927_universe.atoms.segments}
        counts = collections.Counter(row["label_asym_id"] for row in
                                     _atom_lines(d927_topology_text))
        assert counts["A"] == segment_atoms["seg_0_RBD_pro"] == 2355
        assert counts["B"] == segment_atoms["seg_1_D927"] == 53
        assert counts["C"] == segment_atoms["seg_2_KRAS_pro"] == 2648
        assert counts["D"] == segment_atoms["seg_3_GNP_neg"] == 45
        assert counts["W"] == segment_atoms["seg_5_SOL"] == 111063
        assert counts["I"] == (segment_atoms["seg_4_Mg"]
                               + segment_atoms["seg_6_NA"]
                               + segment_atoms["seg_7_CL"]) == 219
        assert sum(counts.values()) == d927_universe.atoms.n_atoms == 116383

    def test_protein_selectors_use_segment_letters(self, d927_universe):
        residues = d927_universe.atoms.residues
        rbd = next(res for res in residues
                   if res.segid == "seg_0_RBD_pro" and int(res.resid) == 51)
        kras = [res for res in residues if res.segid == "seg_2_KRAS_pro"]
        # 该体系 KRAS 编号为 144-310、与 RBD 1-142 不重叠，故没有 resid=51；
        # 同号冲突由合成两 segment 用例（TestTopologyPdb）覆盖。
        assert all(int(res.resid) != 51 for res in kras)
        selectors = _residue_selectors(d927_universe, [rbd, kras[0]])
        assert selectors[0] == "chain A and resi 51"
        assert selectors[1] == f"chain C and resi {int(kras[0].resid)}"


class TestInteractionExpressions:
    """Mol* 选择表达式（每组残基级 or 连接 + 每原子 b 区间）。"""

    ATOM_BANDS = ["b > 99.75 and b < 100.25",
                  "b > 100.75 and b < 101.25",
                  "b > 101.75 and b < 102.25",
                  "b > 102.75 and b < 103.25"]

    def test_labels_and_residue_expressions(self, numbered_universe):
        residues = interaction_expressions(
            _numbered_interaction(), 0, numbered_universe)["residues"]
        assert [row["label"] for row in residues] == [
            "TYR50 (PDB 51) · aromatic_ring",
            "MOL142 (PDB 143) · aromatic_ring"]
        assert [row["expression"] for row in residues] == [
            f"({self.ATOM_BANDS[0]}) or ({self.ATOM_BANDS[1]})",
            f"({self.ATOM_BANDS[2]}) or ({self.ATOM_BANDS[3]})"]

    def test_one_expression_per_atom(self, numbered_universe):
        residues = interaction_expressions(
            _numbered_interaction(), 0, numbered_universe)["residues"]
        assert [row["name"] for row in residues[0]["atoms"]] == ["CZ", "CE1"]
        assert [row["expression"] for row in residues[0]["atoms"]] == (
            self.ATOM_BANDS[:2])
        assert [row["name"] for row in residues[1]["atoms"]] == ["N16", "N17"]
        assert [row["expression"] for row in residues[1]["atoms"]] == (
            self.ATOM_BANDS[2:])

    def test_expressions_never_fall_back_to_chain_or_resi(
            self, numbered_universe):
        residues = interaction_expressions(
            _numbered_interaction(), 0, numbered_universe)["residues"]
        expressions = [row["expression"] for row in residues] + [
            atom["expression"] for group in residues
            for atom in group["atoms"]]
        assert all("chain" not in text and "resi" not in text
                   for text in expressions)

    def test_residue_selector_uses_topology_chain_rules(self, numbered_universe):
        residues = interaction_expressions(
            _numbered_interaction(), 0, numbered_universe)["residues"]
        assert all(set(row) == {"label", "residue_selector", "is_water",
                                "expression", "group_expression", "indices",
                                "residue_indices", "atoms"}
                   for row in residues)
        assert [row["residue_selector"] for row in residues] == [
            "chain A and resi 51", "chain B and resi 143"]
        chains = _chain_rows(topology_cif(None, None, numbered_universe))
        assert residues[0]["residue_selector"] == (
            f"chain {chains[('TYR', 51)]} and resi 51")
        assert residues[1]["residue_selector"] == (
            f"chain {chains[('MOL', 143)]} and resi 143")

    def test_water_group_expression_keeps_full_number(self, d927_universe):
        """D2 关键：CIF 后水组表达式是 chain W and resi 15689（无 b 区间）。"""
        interaction = load_interactions(
            str(H5_DATA_DIR / "water_bridge.h5"))[0]
        residues = interaction_expressions(
            interaction, 0, d927_universe)["residues"]
        water = [row for row in residues if row["is_water"]]
        assert len(water) == 1
        assert water[0]["group_expression"] == "chain W and resi 15689"
        assert "b >" not in water[0]["group_expression"]

    def test_water_bridge_groups_all_have_residue_selectors(self,
                                                            bridge_universe):
        residues = interaction_expressions(
            _water_bridge_interaction(), 0, bridge_universe)["residues"]
        selectors = [row["residue_selector"] for row in residues]
        assert all(selectors)
        assert selectors == ["chain A and resi 1", "chain W and resi 2",
                             "chain A and resi 3"]

    def test_is_water_flags_only_water_groups(self, bridge_universe,
                                              numbered_universe):
        bridge = interaction_expressions(
            _water_bridge_interaction(), 0, bridge_universe)["residues"]
        assert [row["is_water"] for row in bridge] == [False, True, False]
        stack = interaction_expressions(
            _numbered_interaction(), 0, numbered_universe)["residues"]
        assert [row["is_water"] for row in stack] == [False, False]

    def test_water_bridge_groups_have_disjoint_bands(self, bridge_universe):
        residues = interaction_expressions(
            _water_bridge_interaction(), 0, bridge_universe)["residues"]
        assert [row["label"] for row in residues] == [
            "ALA1 · H_donor", "SOL2 · water", "LYS3 · H_acceptor"]
        assert [row["name"] for row in residues[1]["atoms"]] == ["OW", "HW1"]
        for row in residues:
            assert row["expression"] == " or ".join(
                f"({atom['expression']})" for atom in row["atoms"])
            assert "resi" not in row["expression"]
        bands = [_b_band(row["expression"])
                 for group in residues for row in group["atoms"]]
        assert bands == [(99.75, 100.25), (100.75, 101.25),
                         (101.75, 102.25), (102.75, 103.25),
                         (103.75, 104.25), (104.75, 105.25)]
        for index, (low, high) in enumerate(bands):
            others = bands[:index] + bands[index + 1:]
            assert all(high < other_low or low > other_high
                       for other_low, other_high in others)


class TestStructureSelect:

    def test_select_reports_structure_size(self, client, pdb_file):
        response = _open_structure(client, pdb_file)
        assert response.status_code == 200
        data = response.json()
        assert data["ok"] is True
        assert data["n_atoms"] == len(ATOMS)
        assert data["n_frames"] == 1
        assert data["has_box"] is True

    def test_select_bad_path_returns_400(self, client, tmp_path):
        response = client.post("/api/structure/select",
                               data={"topology": str(tmp_path / "x.tpr"),
                                     "trajectory": str(tmp_path / "x.xtc")})
        assert response.status_code == 400
        assert "cannot open structure" in response.json()["detail"]

    def test_select_empty_path_returns_400(self, client):
        response = client.post("/api/structure/select",
                               data={"topology": "", "trajectory": ""})
        assert response.status_code == 400

    def test_select_outside_root_returns_400(self, tmp_path, h5_dir):
        root = tmp_path / "root"
        root.mkdir()
        outside = _write_pdb(tmp_path / "outside.pdb")
        client = _client(root=str(root))
        client.post("/select", data={"path": str(h5_dir)})
        response = client.post("/api/structure/select",
                               data={"topology": str(outside),
                                     "trajectory": str(outside)})
        assert response.status_code == 400
        assert response.json()["detail"] == "path outside root"


class TestTopologyRoute:

    def test_returns_full_topology_cif(self, client, pdb_file):
        _open_structure(client, pdb_file)
        response = client.get("/api/structure/topology",
                              params={"type": "hydrogen_bond", "idx": 0})
        assert response.status_code == 200
        assert response.headers["content-type"] == "chemical/x-cif"
        assert response.text.startswith("data_dii")
        rows = _b_factor_rows(response.text)
        assert len(rows) == len(ATOMS)
        assert [value for _, value in rows] == [100.0, 101.0, 102.0, 103.0]
        assert response.text.rstrip().endswith("#")

    def test_mmcif_is_parseable_by_biopython(self):
        """标准 mmCIF：Biopython 可解析；水 15689 恰好 3 个原子（无截断）。"""
        try:
            from Bio.PDB import MMCIFParser
        except ImportError:  # pragma: no cover - 环境缺依赖时显式跳过
            pytest.skip("Biopython unavailable: mmCIF parse check skipped")
        from io import StringIO
        client = _client()
        client.post("/select", data={"path": str(H5_DATA_DIR)})
        opened = client.post("/api/structure/select", data={
            "topology": str(D927_TOPOLOGY),
            "trajectory": str(D927_TRAJECTORY)})
        assert opened.status_code == 200
        topology = client.get("/api/structure/topology",
                              params={"type": "pi_stacking", "idx": 0})
        assert topology.status_code == 200
        structure = MMCIFParser(QUIET=True).get_structure(
            "dii", StringIO(topology.text))
        atoms = list(structure.get_atoms())
        assert len(atoms) == opened.json()["n_atoms"] == 116383
        water = [residue for residue in structure.get_residues()
                 if residue.resname.strip() == "SOL"
                 and residue.id[1] == 15689]
        assert len(water) == 1
        assert [atom.get_name() for atom in water[0].get_atoms()] == [
            "OW", "HW1", "HW2"]

    def test_tags_only_pair_atoms_in_full_structure(self, water_client,
                                                    water_pdb):
        _open_structure(water_client, water_pdb)
        response = water_client.get("/api/structure/topology",
                                    params={"type": "hydrogen_bond",
                                            "idx": 0})
        rows = dict(_b_factor_rows(response.text))
        assert len(rows) == len(WATER_ATOMS)
        assert rows["OW"] == 0.0 and rows["HW1"] == 0.0
        assert [rows[name] for name in ("N1", "H1", "O1", "O2")] == [
            100.0, 101.0, 102.0, 103.0]

    def test_missing_pair_parameters_return_400(self, client, pdb_file):
        _open_structure(client, pdb_file)
        assert client.get("/api/structure/topology").status_code == 400
        response = client.get("/api/structure/topology",
                              params={"type": "hydrogen_bond"})
        assert response.status_code == 400
        assert response.json()["detail"] == "type and idx required"

    def test_alignment_errors_return_400(self, client, pdb_file, monkeypatch):
        _open_structure(client, pdb_file)
        monkeypatch.setattr("DuIvyInteractions.web.app.alignment_errors",
                            lambda universe, interaction: ["atom 0: mismatch"])
        response = client.get("/api/structure/topology",
                              params={"type": "hydrogen_bond", "idx": 0})
        assert response.status_code == 400
        assert response.json()["errors"] == ["atom 0: mismatch"]

    def test_without_selection_returns_400(self, client):
        response = client.get("/api/structure/topology")
        assert response.status_code == 400
        assert response.json()["detail"] == "structure not selected"

    def test_served_with_token(self, h5_dir, pdb_file):
        client = _client(token=TOKEN)
        client.post("/select", params={"token": TOKEN},
                    data={"path": str(h5_dir)})
        selected = client.post("/api/structure/select",
                               params={"token": TOKEN},
                               data={"topology": str(pdb_file),
                                     "trajectory": str(pdb_file)})
        assert selected.status_code == 200
        response = client.get("/api/structure/topology",
                              params={"token": TOKEN,
                                      "type": "hydrogen_bond", "idx": 0})
        assert response.status_code == 200
        assert "ATOM" in response.text

    def test_requires_token(self):
        client = _client(token=TOKEN)
        response = client.get("/api/structure/topology")
        assert response.status_code == 401


class TestTrajectoryRoute:

    def test_returns_raw_trajectory_bytes(self, client, pdb_file):
        _open_structure(client, pdb_file)
        response = client.get("/api/structure/trajectory")
        assert response.status_code == 200
        assert response.headers["content-type"] == "application/octet-stream"
        assert response.content == pdb_file.read_bytes()
        assert len(response.content) == pdb_file.stat().st_size
        assert response.content[:4] == b"CRYS"

    def test_without_selection_returns_400(self, client):
        response = client.get("/api/structure/trajectory")
        assert response.status_code == 400
        assert response.json()["detail"] == "structure not selected"


class TestInteractionsRoute:

    def _payload(self, client, pdb_file, type_name="pi_stacking", idx=0):
        _open_structure(client, pdb_file)
        return client.get("/api/structure/interactions",
                          params={"type": type_name, "idx": idx})

    def test_payload_shape(self, numbered_client, numbered_pdb):
        response = self._payload(numbered_client, numbered_pdb)
        assert response.status_code == 200
        data = response.json()
        assert set(data) == {"pair_label", "n_frames", "existence",
                             "pt_types", "residues", "n_atoms"}
        assert "contacts" not in data
        assert data["pt_types"] == []
        assert data["n_frames"] == 1
        assert data["existence"] == [1]
        assert data["n_atoms"] == len(NUMBERED_ATOMS)

    def test_payload_carries_expressions(self, numbered_client, numbered_pdb):
        data = self._payload(numbered_client, numbered_pdb).json()
        assert [row["expression"] for row in data["residues"]] == [
            "(b > 99.75 and b < 100.25) or (b > 100.75 and b < 101.25)",
            "(b > 101.75 and b < 102.25) or (b > 102.75 and b < 103.25)"]
        assert len(data["residues"][0]["atoms"]) == 2

    def test_alignment_errors_return_400(self, numbered_client, numbered_pdb,
                                         monkeypatch):
        monkeypatch.setattr("DuIvyInteractions.web.app.alignment_errors",
                            lambda universe, interaction: ["atom 0: mismatch"])
        response = self._payload(numbered_client, numbered_pdb)
        assert response.status_code == 400
        assert response.json() == {"detail": "atom 0: mismatch",
                                   "errors": ["atom 0: mismatch"]}

    def test_unknown_type_returns_404(self, numbered_client, numbered_pdb):
        response = self._payload(numbered_client, numbered_pdb, type_name="nope")
        assert response.status_code == 404

    def test_pair_out_of_range_returns_404(self, numbered_client, numbered_pdb):
        response = self._payload(numbered_client, numbered_pdb, idx=9)
        assert response.status_code == 404

    def test_without_selection_returns_400(self, numbered_client):
        response = numbered_client.get("/api/structure/interactions",
                                       params={"type": "pi_stacking",
                                               "idx": 0})
        assert response.status_code == 400
        assert response.json()["detail"] == "structure not selected"


class TestPairLabel:
    """pair_label 必须与导出器给表格/CSV 的名字逐字一致。"""

    def test_pair_label_matches_exporter_for_synthetic_pair(
            self, numbered_client, numbered_pdb):
        _open_structure(numbered_client, numbered_pdb)
        data = numbered_client.get(
            "/api/structure/interactions",
            params={"type": "pi_stacking", "idx": 0}).json()
        interaction = _numbered_interaction()
        assert data["pair_label"] == get_exporter(
            interaction).get_pair_label(interaction, 0)

    def test_pair_label_matches_exporter_for_real_h5(self):
        client = _client()
        client.post("/select", data={"path": str(H5_DATA_DIR)})
        client.post("/api/structure/select",
                    data={"topology": str(D927_TOPOLOGY),
                          "trajectory": str(D927_TRAJECTORY)})
        data = client.get("/api/structure/interactions",
                          params={"type": "pi_stacking", "idx": 0}).json()
        interaction = load_interactions(
            str(H5_DATA_DIR / "pi_stacking.h5"))[0]
        assert data["pair_label"] == get_exporter(
            interaction).get_pair_label(interaction, 0)

    def test_contact_route_is_removed(self):
        # 拼接旧路径，避免"旧路由已删除"的检索残留（代码里不应再出现该字面量）
        removed_path = "/api/structure/" + "contact"
        response = _client().get(removed_path,
                                 params={"type": "pi_stacking", "idx": 0})
        assert response.status_code == 404


class TestStructurePage:

    def test_page_renders_without_structure(self, client):
        response = client.get("/structure/hydrogen_bond/0")
        assert response.status_code == 200
        html = response.text
        assert "molstar.js" in html and "structure.js" in html
        assert 'id="molstarViewer"' in html and 'id="structureConfig"' in html
        assert "structure not selected" in html

    def test_page_renders_with_structure(self, client, pdb_file):
        _open_structure(client, pdb_file)
        html = client.get("/structure/hydrogen_bond/0").text
        assert str(pdb_file) in html
        assert 'id="existenceBar"' in html
        assert "structure-status" in html

    def test_page_links_molstar_stylesheet(self, client):
        html = client.get("/structure/hydrogen_bond/0").text
        assert "/static/vendor/molstar.css" in html

    def test_molstar_stylesheet_served(self, client):
        response = client.get("/static/vendor/molstar.css",
                              headers={"Accept-Encoding": "identity"})
        assert response.status_code == 200
        assert int(response.headers["content-length"]) > 50000

    def test_page_unknown_type_returns_404(self, client):
        assert client.get("/structure/nope/0").status_code == 404

    def test_page_pair_out_of_range_returns_404(self, client):
        assert client.get("/structure/hydrogen_bond/9").status_code == 404

    def test_page_allows_molstar_eval(self, client):
        csp = client.get("/structure/hydrogen_bond/0").headers[
            "content-security-policy"]
        assert "'unsafe-eval'" in csp

    def test_other_pages_keep_strict_csp(self, client):
        csp = client.get("/").headers["content-security-policy"]
        assert "'unsafe-eval'" not in csp

    def test_pair_page_embeds_3d_viewer(self, client):
        html = client.get("/pair/hydrogen_bond/0").text
        assert 'id="molstarViewer"' in html
        assert "/static/structure.js" in html


# ============================================================
# 水相关判定 / 基团标签 / 接触标签
# ============================================================

class TestInvolvesWater:

    def test_water_bridge_type_is_water_related(self):
        assert involves_water(_water_bridge_interaction(), 0) is True

    def test_pi_stacking_is_not_water_related(self):
        interaction = _single("pi_stacking", [
            (_group(0, "aromatic_ring", [_atom(0, "CZ", "C")]),
             _group(1, "aromatic_ring", [_atom(1, "CZ", "C")]))])
        assert involves_water(interaction, 0) is False

    def test_hydrogen_bond_with_water_group_is_water_related(self):
        interaction = _single("hydrogen_bond", [
            (_group(0, "H_donor", [_atom(0, "N1", "N")]),
             _group(1, "water", [_atom(1, "OW", "O")], "SOL", 2))])
        assert involves_water(interaction, 0) is True


class TestGroupLabel:

    @staticmethod
    def _pair():
        return _single("pi_stacking", [
            (_group(0, "aromatic_ring", [_atom(0, "CZ", "C")], "TYR", 50),
             _group(1, "aromatic_ring", [_atom(2, "N16", "N")], "MOL", 142,
                    molecule="D927"))])

    def test_protein_group_is_residue_plus_group_type(self):
        interaction = self._pair()
        assert group_label(interaction, interaction.groups[0][0]) == (
            "TYR50 · aromatic_ring")

    def test_ligand_group_omits_molecule_suffix(self):
        interaction = self._pair()
        ligand = interaction.groups[0][1]
        assert group_label(interaction, ligand) == "MOL142 · aromatic_ring"
        assert "molecule" not in group_label(interaction, ligand)

    def test_water_group_keeps_water_group_type(self):
        interaction = _water_bridge_interaction()
        assert group_label(interaction, interaction.groups[0][1]) == (
            "SOL2 · water")

    def test_universe_numbering_appears_in_parentheses(self, numbered_universe):
        interaction = _numbered_interaction()
        assert group_label(interaction, interaction.groups[0][0],
                           numbered_universe) == (
            "TYR50 (PDB 51) · aromatic_ring")
        assert group_label(interaction, interaction.groups[0][1],
                           numbered_universe) == (
            "MOL142 (PDB 143) · aromatic_ring")

    def test_equal_numbering_omits_pdb_parenthesis(self, universe):
        interaction = _interaction()
        assert group_label(interaction, interaction.groups[0][0],
                           universe) == "LYS1 · H_donor"

    def test_long_residue_number_is_kept_verbatim(self, water_universe):
        """mmCIF 不再截断：水 resid 15689 在标签与选择器里都是 15689。"""
        interaction = _water_bridge_interaction()
        water_universe.atoms.residues.resids = [15689, 1]
        assert group_label(interaction, interaction.groups[0][1],
                           water_universe) == "SOL2 (PDB 15689) · water"

    def test_residue_selector_keeps_five_digits(self, water_universe):
        interaction = _water_bridge_interaction()
        water_universe.atoms.residues.resids = [15689, 1]
        residues = interaction_expressions(
            interaction, 0, water_universe)["residues"]
        assert residues[1]["residue_selector"] == "chain W and resi 15689"

    def test_without_universe_falls_back_to_h5_number(self):
        interaction = _numbered_interaction()
        assert group_label(interaction, interaction.groups[0][0]) == (
            "TYR50 · aromatic_ring")


class TestGroupAtoms:

    def test_returns_one_index_list_per_group(self):
        interaction = _water_bridge_interaction()
        assert group_atoms(interaction, 0) == [[2, 3], [0, 1], [4, 5]]

    def test_lengths_and_order_follow_group_atoms(self):
        first = [_atom(7, "CZ", "C"), _atom(3, "CE1", "C"),
                 _atom(9, "CD2", "C")]
        interaction = _single("pi_stacking", [
            (_group(0, "aromatic_ring", first, "TYR", 50),
             _group(1, "aromatic_ring", [_atom(1, "N16", "N")], "MOL", 142))])
        indices = group_atoms(interaction, 0)
        assert [len(row) for row in indices] == [3, 1]
        assert indices == [[7, 3, 9], [1]]
        assert group_atoms(interaction, 0) == indices


class TestContactLabel:

    @staticmethod
    def _pair():
        return _single("hydrogen_bond", [
            (_group(0, "H_donor", [_atom(0, "CZ", "C"),
                                   _atom(1, "OH", "O")], "TYR", 50),
             _group(1, "H_acceptor", [_atom(2, "N16", "N"),
                                      _atom(3, "N17", "N")], "MOL", 142))])

    def test_atom_pair_is_group_label_plus_atom_name(self):
        assert contact_label(self._pair(), 0, (0, 0)) == (
            "TYR50 · H_donor CZ ┄ MOL142 · H_acceptor N16")

    def test_without_contact_stops_at_group_label(self):
        assert contact_label(self._pair(), 0, None) == (
            "TYR50 · H_donor ┄ MOL142 · H_acceptor")

    def test_aromatic_ring_format_matches_view_spec(self):
        interaction = _single("pi_stacking", [
            (_group(0, "aromatic_ring", [_atom(0, "CZ", "C")], "TYR", 50),
             _group(1, "aromatic_ring", [_atom(2, "N16", "N")], "MOL",
                    142))])
        assert contact_label(interaction, 0, (0, 0)) == (
            "TYR50 · aromatic_ring CZ ┄ MOL142 · aromatic_ring N16")

    def test_water_bridge_labels_first_and_last_group(self):
        interaction = _water_bridge_interaction()
        assert contact_label(interaction, 0, (1, 0)) == (
            "ALA1 · H_donor H1 ┄ LYS3 · H_acceptor O1")

    def test_water_bridge_without_contact_keeps_group_types(self):
        interaction = _water_bridge_interaction()
        assert contact_label(interaction, 0, None) == (
            "ALA1 · H_donor ┄ LYS3 · H_acceptor")

    def test_atom_index_out_of_range_stops_at_group_label(self):
        assert contact_label(self._pair(), 0, (0, 9)) == (
            "TYR50 · H_donor CZ ┄ MOL142 · H_acceptor")

    def test_universe_numbering_is_used_when_given(self, numbered_universe):
        interaction = _numbered_interaction()
        assert contact_label(interaction, 0, (0, 0),
                             numbered_universe) == (
            "TYR50 (PDB 51) · aromatic_ring CZ ┄ "
            "MOL142 (PDB 143) · aromatic_ring N16")

    def test_without_contact_uses_universe_numbering(self, numbered_universe):
        interaction = _numbered_interaction()
        assert contact_label(interaction, 0, None,
                             numbered_universe) == (
            "TYR50 (PDB 51) · aromatic_ring ┄ MOL142 (PDB 143) · "
            "aromatic_ring")


# ============================================================
# /api/list 的 kind 过滤 与 /api/pick 的 file 模式
# ============================================================

class TestListKind:

    @pytest.fixture
    def server_client(self, tmp_path):
        for name in ("md.tpr", "md.xtc", "data.h5", "notes.txt"):
            (tmp_path / name).write_bytes(b"x")
        return _client(browse_mode="server", root=str(tmp_path))

    def test_structure_kind_lists_tpr_and_xtc_not_h5(self, server_client):
        data = server_client.get("/api/list",
                                 params={"path": "",
                                         "kind": "structure"}).json()
        assert [item["name"] for item in data["files"]] == [
            "md.tpr", "md.xtc"]
        assert data["truncated"] is False

    def test_structure_kind_skips_other_suffixes(self, server_client):
        data = server_client.get("/api/list",
                                 params={"path": "",
                                         "kind": "structure"}).json()
        names = {item["name"] for item in data["files"]}
        assert "notes.txt" not in names and "data.h5" not in names

    def test_dataset_kind_lists_only_h5(self, server_client):
        data = server_client.get("/api/list", params={"path": ""}).json()
        assert [item["name"] for item in data["files"]] == ["data.h5"]

    def test_unknown_kind_returns_400(self, server_client):
        resp = server_client.get("/api/list",
                                 params={"path": "", "kind": "nope"})
        assert resp.status_code == 400
        assert resp.json()["detail"] == "unsupported kind: nope"

    def test_structure_kind_still_disabled_in_picker_mode(self, tmp_path):
        (tmp_path / "md.tpr").write_bytes(b"x")
        resp = _client().get("/api/list",
                             params={"path": str(tmp_path),
                                     "kind": "structure"})
        assert resp.status_code == 404
        assert resp.json()["detail"] == "server browse disabled"


class TestPickFile:

    def test_file_mode_returns_selected_path(self, monkeypatch):
        monkeypatch.setattr(
            "DuIvyInteractions.web.app.pick_file",
            lambda kind="structure": "/data/md.tpr")
        resp = _client().post("/api/pick",
                              params={"mode": "file", "kind": "structure"})
        assert resp.status_code == 200
        assert resp.json() == {"path": "/data/md.tpr"}

    def test_kind_is_forwarded_to_pick_file(self, monkeypatch):
        seen = {}

        def fake(kind="structure"):
            seen["kind"] = kind
            return ""

        monkeypatch.setattr("DuIvyInteractions.web.app.pick_file", fake)
        _client().post("/api/pick", params={"mode": "file", "kind": "structure"})
        assert seen["kind"] == "structure"

    def test_dir_mode_is_the_default(self, monkeypatch, tmp_path):
        monkeypatch.setattr("DuIvyInteractions.web.app.pick_directory",
                            lambda: str(tmp_path))
        assert _client().post("/api/pick").json() == {"path": str(tmp_path)}

    def test_file_mode_reports_dialog_errors(self, monkeypatch):
        def boom(kind="structure"):
            raise RuntimeError("no display name and no $DISPLAY")

        monkeypatch.setattr("DuIvyInteractions.web.app.pick_file", boom)
        resp = _client().post("/api/pick", params={"mode": "file"})
        assert resp.status_code == 400
        assert "DISPLAY" in resp.json()["error"]

    def test_file_mode_disabled_in_server_mode(self):
        resp = _client(browse_mode="server").post(
            "/api/pick", params={"mode": "file"})
        assert resp.status_code == 400
        assert "native picker unavailable" in resp.json()["error"]

    def test_unknown_mode_returns_400(self):
        resp = _client().post("/api/pick", params={"mode": "nope"})
        assert resp.status_code == 400
        assert resp.json()["error"] == "unsupported mode: nope"


class TestPickerFileDialog:

    def test_structure_filetypes_cover_md_suffixes(self):
        pattern = FILE_TYPES["structure"][0][1]
        for suffix in ("*.tpr", "*.xtc", "*.pdb", "*.mol2"):
            assert suffix in pattern

    def test_unknown_kind_raises_value_error(self):
        with pytest.raises(ValueError, match="unsupported file kind"):
            pick_file("nope")


class TestStructurePageBrowseMode:

    def test_page_renders_in_server_browse_mode(self, h5_dir):
        client = _client(browse_mode="server")
        client.post("/select", data={"path": str(h5_dir)})
        assert client.get("/structure/hydrogen_bond/0").status_code == 200
