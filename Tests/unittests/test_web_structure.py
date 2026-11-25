# -*- coding: utf-8 -*-
"""3D 结构可视化测试：structure 纯函数层 + /api/structure/* 路由。

不依赖真实轨迹：用 tmp_path 构造最小 PDB 与 h5。
"""

import io
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pytest

from DuIvyInteractions.core.datas import AtomData, Group, Interaction
from DuIvyInteractions.io.h5 import save_interactions
from DuIvyInteractions.visualizers.structure import (
    FRAME_CACHE_SIZE, SAMPLE_ATOMS, STRUCTURE_CACHE_SIZE, WATER_RESIDUES,
    alignment_errors, closest_contact, contact_label, frame_pdb, frame_view,
    group_atoms, group_label, involves_water, load_universe,
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
# h5 的 residue_id 是 0 基索引，universe 的 resid 是 1 基 auth 编号（差 1）。
NUMBERED_ATOMS = [("CZ", (0.0, 0.0, 0.0), "C", "TYR", 51),
                  ("CE1", (4.0, 0.0, 2.0), "C", "TYR", 51),
                  ("N16", (0.0, 0.0, 1.2), "N", "MOL", 143),
                  ("N17", (0.0, 0.0, 5.0), "N", "MOL", 143)]


def _pdb_text(atoms, box=BOX) -> str:
    lines = []
    if box:
        lines.append("CRYST1%9.3f%9.3f%9.3f  90.00  90.00  90.00 P 1           1"
                     % box)
    for serial, (name, (x, y, z), element) in enumerate(atoms, start=1):
        lines.append(
            f"ATOM  {serial:5d} {name:>4s} ALA A   1    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}")
    lines.append("END")
    return "\n".join(lines) + "\n"


def _write_pdb(path: Path, atoms=None, box=BOX) -> Path:
    path.write_text(_pdb_text(atoms or ATOMS, box), encoding="utf-8")
    return path


def _water_pdb_text(atoms, box=BOX) -> str:
    lines = []
    if box:
        lines.append("CRYST1%9.3f%9.3f%9.3f  90.00  90.00  90.00 P 1           1"
                     % box)
    for serial, (name, (x, y, z), element, resname, resid) in enumerate(
            atoms, start=1):
        lines.append(
            f"ATOM  {serial:5d} {name:>4s} {resname:>3s} A{resid:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}")
    lines.append("END")
    return "\n".join(lines) + "\n"


def _write_water_pdb(path: Path, atoms=None) -> Path:
    path.write_text(_water_pdb_text(atoms or WATER_ATOMS), encoding="utf-8")
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


def _atom_lines(text: str):
    return [line for line in text.splitlines() if line.startswith("ATOM")]


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
def numbered_pdb(tmp_path):
    return _write_water_pdb(tmp_path / "numbered.pdb", NUMBERED_ATOMS)


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


# ============================================================
# structure 纯函数层
# ============================================================

class TestConstants:

    def test_cache_and_sample_sizes(self):
        assert STRUCTURE_CACHE_SIZE == 2
        assert FRAME_CACHE_SIZE == 8
        assert SAMPLE_ATOMS == 5


class TestLoadUniverse:

    def test_opens_pdb_topology_and_trajectory(self, pdb_file):
        universe = load_universe(str(pdb_file), str(pdb_file))
        assert universe.atoms.n_atoms == len(ATOMS)
        assert len(universe.trajectory) == 1

    def test_bad_file_raises_value_error(self, tmp_path):
        bad = tmp_path / "bad.txt"
        bad.write_text("not a structure")
        with pytest.raises(ValueError, match="cannot open structure"):
            load_universe(str(bad), str(bad))

    def test_missing_file_raises_value_error(self, tmp_path):
        with pytest.raises(ValueError, match="cannot open structure"):
            load_universe(str(tmp_path / "nope.tpr"), str(tmp_path / "nope.xtc"))


class TestFramePdb:

    def test_atom_count_and_order(self, universe):
        text = frame_pdb(universe, 0)
        names = [line[12:16].strip() for line in text.splitlines()
                 if line.startswith("ATOM")]
        assert len(names) == universe.atoms.n_atoms
        assert names == list(universe.atoms.names)

    def test_out_of_range_raises(self, universe):
        with pytest.raises(ValueError, match="frame out of range: 1"):
            frame_pdb(universe, 1)

    def test_box_wraps_coordinates(self, tmp_path):
        path = _write_pdb(tmp_path / "wrapped.pdb",
                          [("C1", (12.0, 2.0, 3.0), "C")])
        universe = load_universe(str(path), str(path))
        line = next(line for line in frame_pdb(universe, 0).splitlines()
                    if line.startswith("ATOM"))
        assert float(line[30:38]) == pytest.approx(2.0)
        assert float(line[38:46]) == pytest.approx(2.0)

    def test_without_box_keeps_coordinates(self, tmp_path):
        path = _write_pdb(tmp_path / "nobox.pdb",
                          [("C1", (12.0, 2.0, 3.0), "C")], box=None)
        universe = load_universe(str(path), str(path))
        assert universe.dimensions is None
        line = next(line for line in frame_pdb(universe, 0).splitlines()
                    if line.startswith("ATOM"))
        assert float(line[30:38]) == pytest.approx(12.0)

    def test_original_universe_is_not_modified(self, universe):
        before = universe.atoms.positions.copy()
        frame_pdb(universe, 0)
        assert np.array_equal(universe.atoms.positions, before)


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


class TestClosestContact:

    def test_returns_nearest_pair(self, universe):
        assert closest_contact(universe, [0, 1], [2, 3]) == (1, 2)

    def test_empty_group_returns_none(self, universe):
        assert closest_contact(universe, [], [2, 3]) is None
        assert closest_contact(universe, [0], []) is None

    def test_uses_minimum_image(self, tmp_path):
        atoms = [("A1", (9.5, 0.0, 0.0), "C"),
                 ("B1", (0.5, 0.0, 0.0), "O"),
                 ("C1", (5.0, 0.0, 0.0), "N")]
        path = _write_pdb(tmp_path / "pbc.pdb", atoms)
        universe = load_universe(str(path), str(path))
        assert closest_contact(universe, [0], [1, 2]) == (0, 1)


# ============================================================
# 路由
# ============================================================

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


class TestStructureView:

    def test_view_without_selection_returns_400(self, client):
        response = client.get("/api/structure/view",
                              params={"type": "hydrogen_bond", "idx": 0,
                                      "frame": 0})
        assert response.status_code == 400
        assert response.json()["detail"] == "structure not selected"

    def test_view_returns_frame_payload(self, client, pdb_file):
        _open_structure(client, pdb_file)
        response = client.get("/api/structure/view",
                              params={"type": "hydrogen_bond", "idx": 0,
                                      "frame": 0})
        assert response.status_code == 200
        data = response.json()
        assert data["frame"] == 0 and data["n_frames"] == 1
        assert data["active"] is True
        assert data["existence"] == [1]
        assert data["contact"] == [1, 2]
        assert data["pdb"].count("ATOM") == len(ATOMS)

    def test_view_inactive_frame_has_no_contact(self, tmp_path, pdb_file):
        directory = tmp_path / "inactive"
        directory.mkdir()
        save_interactions([_interaction(active=False)],
                          str(directory / "hydrogen_bond.h5"))
        client = _client()
        client.post("/select", data={"path": str(directory)})
        _open_structure(client, pdb_file)
        data = client.get("/api/structure/view",
                          params={"type": "hydrogen_bond", "idx": 0,
                                  "frame": 0}).json()
        assert data["active"] is False and data["contact"] is None

    def test_view_frame_out_of_range_returns_404(self, client, pdb_file):
        _open_structure(client, pdb_file)
        response = client.get("/api/structure/view",
                              params={"type": "hydrogen_bond", "idx": 0,
                                      "frame": 5})
        assert response.status_code == 404

    def test_view_unknown_type_returns_404(self, client, pdb_file):
        _open_structure(client, pdb_file)
        response = client.get("/api/structure/view",
                              params={"type": "nope", "idx": 0, "frame": 0})
        assert response.status_code == 404

    def test_view_pair_out_of_range_returns_404(self, client, pdb_file):
        _open_structure(client, pdb_file)
        response = client.get("/api/structure/view",
                              params={"type": "hydrogen_bond", "idx": 9,
                                      "frame": 0})
        assert response.status_code == 404

    def test_view_requires_token(self, h5_dir, pdb_file):
        client = _client(token=TOKEN)
        response = client.get("/api/structure/view",
                              params={"type": "hydrogen_bond", "idx": 0,
                                      "frame": 0})
        assert response.status_code == 401

    def test_alignment_errors_are_returned(self, client, pdb_file, monkeypatch):
        _open_structure(client, pdb_file)
        monkeypatch.setattr("DuIvyInteractions.web.app.alignment_errors",
                            lambda universe, interaction: ["atom 0: mismatch"])
        response = client.get("/api/structure/view",
                              params={"type": "hydrogen_bond", "idx": 0,
                                      "frame": 0})
        assert response.status_code == 400
        data = response.json()
        assert data["errors"] == ["atom 0: mismatch"]
        assert data["detail"] == "atom 0: mismatch"


class TestStructurePage:

    def test_page_renders_without_structure(self, client):
        response = client.get("/structure/hydrogen_bond/0")
        assert response.status_code == 200
        html = response.text
        assert "molstar.js" in html and "structure.js" in html
        assert 'id="frameSlider"' in html and 'id="structureConfig"' in html
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

    def test_pair_page_links_to_3d_view(self, client):
        html = client.get("/pair/hydrogen_bond/0").text
        assert 'href="/structure/hydrogen_bond/0"' in html
        assert "3D view" in html


# ============================================================
# frame_view / 水相关判定 / 接触标签
# ============================================================

class TestFrameView:

    def test_hide_water_drops_water_atoms(self, water_universe):
        text, shown = frame_view(water_universe, 0, hide_water=True)
        assert shown == [2, 3, 4, 5]
        assert text.count("ATOM") == len(shown)
        assert "SOL" not in text

    def test_shown_maps_each_pdb_line_to_global_atom(self, water_universe):
        text, shown = frame_view(water_universe, 0, hide_water=True)
        names = [line[12:16].strip() for line in text.splitlines()
                 if line.startswith("ATOM")]
        assert names == [water_universe.atoms[index].name for index in shown]

    def test_keep_atoms_forces_water_atoms(self, water_universe):
        text, shown = frame_view(water_universe, 0, hide_water=True,
                                 keep_atoms=frozenset({0}))
        assert shown == [0, 2, 3, 4, 5]
        names = [line[12:16].strip() for line in text.splitlines()
                 if line.startswith("ATOM")]
        assert names == ["OW", "N1", "H1", "O1", "O2"]

    def test_without_hide_returns_all_atoms(self, water_universe):
        text, shown = frame_view(water_universe, 0)
        assert shown == list(range(water_universe.atoms.n_atoms))
        assert text == frame_pdb(water_universe, 0)

    def test_out_of_range_raises(self, water_universe):
        with pytest.raises(ValueError, match="frame out of range: 1"):
            frame_view(water_universe, 1, hide_water=True)

    def test_original_universe_positions_unchanged(self, water_universe):
        before = water_universe.atoms.positions.copy()
        frame_view(water_universe, 0, hide_water=True)
        assert np.array_equal(water_universe.atoms.positions, before)

    def test_water_residue_constant_covers_common_names(self):
        assert {"SOL", "HOH", "WAT", "TIP3", "TIP4", "HO4", "HO5",
                "H2O"} <= WATER_RESIDUES


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

    def test_long_residue_number_is_truncated_to_four_columns(
            self, water_universe):
        interaction = _water_bridge_interaction()
        water_universe.atoms.residues.resids = [15689, 1]
        assert group_label(interaction, interaction.groups[0][1],
                           water_universe) == "SOL2 (PDB 5689) · water"

    def test_pdb_resid_truncates_to_four_digits(self):
        from DuIvyInteractions.visualizers.structure import _pdb_resid
        assert _pdb_resid(15689) == 5689
        assert _pdb_resid(51) == 51

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
# /api/structure/view 的 solvent 模式
# ============================================================

class TestStructureViewSolvent:

    def test_hide_reports_fewer_atoms_than_total(self, water_client, water_pdb):
        _open_structure(water_client, water_pdb)
        data = water_client.get("/api/structure/view",
                                params={"type": "hydrogen_bond", "idx": 0,
                                        "frame": 0,
                                        "solvent": "hide"}).json()
        assert data["solvent"] == "hide"
        assert data["involves_water"] is False
        assert data["n_atoms_total"] == len(WATER_ATOMS)
        assert data["n_atoms_shown"] == len(WATER_ATOMS) - 2
        assert "SOL" not in data["pdb"]

    def test_hide_contact_uses_display_indices(self, water_client, water_pdb):
        _open_structure(water_client, water_pdb)
        data = water_client.get("/api/structure/view",
                                params={"type": "hydrogen_bond", "idx": 0,
                                        "frame": 0,
                                        "solvent": "hide"}).json()
        assert data["contact"] == [1, 2]
        lines = [line for line in data["pdb"].splitlines()
                 if line.startswith("ATOM")]
        assert [lines[row][12:16].strip() for row in data["contact"]] == [
            "H1", "O1"]
        assert data["contact_label"] == (
            "ALA1 · H_donor H1 ┄ ALA1 · H_acceptor O1")

    def test_full_solvent_keeps_all_atoms(self, water_client, water_pdb):
        _open_structure(water_client, water_pdb)
        data = water_client.get("/api/structure/view",
                                params={"type": "hydrogen_bond", "idx": 0,
                                        "frame": 0}).json()
        assert data["solvent"] == "full"
        assert data["n_atoms_shown"] == data["n_atoms_total"]
        assert data["contact"] == [3, 4]

    def test_unknown_solvent_falls_back_to_full(self, water_client, water_pdb):
        _open_structure(water_client, water_pdb)
        data = water_client.get("/api/structure/view",
                                params={"type": "hydrogen_bond", "idx": 0,
                                        "frame": 0,
                                        "solvent": "weird"}).json()
        assert data["solvent"] == "full"

    def test_water_bridge_ignores_hide(self, water_bridge_client, water_pdb):
        _open_structure(water_bridge_client, water_pdb)
        data = water_bridge_client.get(
            "/api/structure/view",
            params={"type": "water_bridge", "idx": 0, "frame": 0,
                    "solvent": "hide"}).json()
        assert data["involves_water"] is True
        assert data["solvent"] == "full"
        assert data["n_atoms_shown"] == data["n_atoms_total"]
        assert "SOL" in data["pdb"]

    def test_frame_cache_keeps_solvent_modes_apart(self, water_client,
                                                   water_pdb):
        _open_structure(water_client, water_pdb)
        params = {"type": "hydrogen_bond", "idx": 0, "frame": 0}
        hide = water_client.get("/api/structure/view",
                                params={**params, "solvent": "hide"}).json()
        full = water_client.get("/api/structure/view",
                                params={**params, "solvent": "full"}).json()
        again = water_client.get("/api/structure/view",
                                 params={**params, "solvent": "hide"}).json()
        assert hide["pdb"] != full["pdb"]
        assert len(hide["pdb"]) < len(full["pdb"])
        assert again["pdb"] == hide["pdb"]
        assert again["n_atoms_shown"] == hide["n_atoms_shown"]


# ============================================================
# /api/structure/view 的 groups（两端基团原子 + 标签）
# ============================================================

class TestStructureViewGroups:

    def test_groups_follow_pair_groups_in_full_mode(self, client, pdb_file):
        _open_structure(client, pdb_file)
        data = client.get("/api/structure/view",
                          params={"type": "hydrogen_bond", "idx": 0,
                                  "frame": 0}).json()
        # h5 与 PDB 编号同为 1 → 标签只写 h5 编号（前缀取 h5 残基名）
        assert data["groups"] == [
            {"label": "LYS1 · H_donor", "atoms": [0, 1]},
            {"label": "LYS1 · H_acceptor", "atoms": [2, 3]}]
        assert data["n_atoms_shown"] == len(ATOMS)
        assert data["contact_label"] == (
            "LYS1 · H_donor H1 ┄ LYS1 · H_acceptor O1")

    def test_group_atoms_use_display_rows_in_hide_mode(self, water_client,
                                                       water_pdb):
        _open_structure(water_client, water_pdb)
        data = water_client.get("/api/structure/view",
                                params={"type": "hydrogen_bond", "idx": 0,
                                        "frame": 0,
                                        "solvent": "hide"}).json()
        assert data["n_atoms_shown"] == len(WATER_ATOMS) - 2
        assert data["groups"] == [
            {"label": "ALA1 · H_donor", "atoms": [0, 1]},
            {"label": "ALA1 · H_acceptor", "atoms": [2, 3]}]
        lines = [line for line in data["pdb"].splitlines()
                 if line.startswith("ATOM")]
        names = [lines[row][12:16].strip()
                 for group in data["groups"] for row in group["atoms"]]
        assert names == ["N1", "H1", "O1", "O2"]

    def test_every_group_atom_is_inside_shown_pdb(self, water_client,
                                                  water_pdb):
        _open_structure(water_client, water_pdb)
        data = water_client.get("/api/structure/view",
                                params={"type": "hydrogen_bond", "idx": 0,
                                        "frame": 0,
                                        "solvent": "hide"}).json()
        assert data["groups"][0]["atoms"] == [0, 1]  # 全局 [2, 3] 已映射
        for group in data["groups"]:
            assert all(0 <= atom < data["n_atoms_shown"]
                       for atom in group["atoms"])

    def test_water_bridge_reports_three_groups(self, water_bridge_client,
                                               water_pdb):
        _open_structure(water_bridge_client, water_pdb)
        data = water_bridge_client.get(
            "/api/structure/view",
            params={"type": "water_bridge", "idx": 0, "frame": 0,
                    "solvent": "hide"}).json()
        assert data["involves_water"] is True
        # h5 末组 LYS3 的原子在 PDB 中属 resid 1 → 补 (PDB 1)
        assert [group["label"] for group in data["groups"]] == [
            "ALA1 · H_donor", "SOL2 · water", "LYS3 (PDB 1) · H_acceptor"]
        assert data["groups"][1]["atoms"] == [0, 1]

    def test_view_labels_show_both_numberings(self, numbered_client,
                                              numbered_pdb):
        _open_structure(numbered_client, numbered_pdb)
        data = numbered_client.get("/api/structure/view",
                                   params={"type": "pi_stacking", "idx": 0,
                                           "frame": 0,
                                           "solvent": "hide"}).json()
        assert [group["label"] for group in data["groups"]] == [
            "TYR50 (PDB 51) · aromatic_ring",
            "MOL142 (PDB 143) · aromatic_ring"]
        assert data["contact_label"] == (
            "TYR50 (PDB 51) · aromatic_ring CZ ┄ "
            "MOL142 (PDB 143) · aromatic_ring N16")

    def test_missing_display_atom_is_skipped(self):
        from DuIvyInteractions.web.app import _shown_atoms
        assert _shown_atoms([0, 9], {0: 5}) == [5]
        assert _shown_atoms([0, 9], None) == [0, 9]


# ============================================================
# /api/structure/view 的 highlights（每组一个独立小结构）
# ============================================================

class TestStructureViewHighlights:

    def test_one_highlight_per_group_with_h5_atom_names(self, client,
                                                        pdb_file):
        _open_structure(client, pdb_file)
        data = client.get("/api/structure/view",
                          params={"type": "hydrogen_bond", "idx": 0,
                                  "frame": 0}).json()
        assert [view["label"] for view in data["highlights"]] == [
            "LYS1 · H_donor", "LYS1 · H_acceptor"]
        names = [[line[12:16].strip() for line in _atom_lines(view["pdb"])]
                 for view in data["highlights"]]
        assert names == [["N1", "H1"], ["O1", "O2"]]

    def test_contact_pdb_holds_the_closest_pair_in_two_chains(self, client,
                                                              pdb_file,
                                                              universe):
        _open_structure(client, pdb_file)
        data = client.get("/api/structure/view",
                          params={"type": "hydrogen_bond", "idx": 0,
                                  "frame": 0}).json()
        lines = _atom_lines(data["contact_pdb"])
        assert len(lines) == 2
        expected = closest_contact(universe, [0, 1], [2, 3])
        assert [line[12:16].strip() for line in lines] == [
            universe.atoms[index].name for index in expected]
        assert lines[0][21] == "A" and lines[1][21] == "B"

    def test_highlight_uses_universe_auth_residue_numbers(self, numbered_client,
                                                          numbered_pdb):
        _open_structure(numbered_client, numbered_pdb)
        data = numbered_client.get("/api/structure/view",
                                   params={"type": "pi_stacking", "idx": 0,
                                           "frame": 0,
                                           "solvent": "hide"}).json()
        assert [view["label"] for view in data["highlights"]] == [
            "TYR50 (PDB 51) · aromatic_ring",
            "MOL142 (PDB 143) · aromatic_ring"]
        first = _atom_lines(data["highlights"][0]["pdb"])
        second = _atom_lines(data["highlights"][1]["pdb"])
        assert {line[17:20].strip() for line in first} == {"TYR"}
        assert {line[22:26].strip() for line in first} == {"51"}
        assert {line[21] for line in first} == {"A"}
        assert {line[17:20].strip() for line in second} == {"MOL"}
        assert {line[22:26].strip() for line in second} == {"143"}
        assert {line[21] for line in second} == {"B"}

    def test_highlight_coordinates_match_main_pdb(self, water_client,
                                                  water_pdb):
        _open_structure(water_client, water_pdb)
        data = water_client.get("/api/structure/view",
                                params={"type": "hydrogen_bond", "idx": 0,
                                        "frame": 0,
                                        "solvent": "hide"}).json()
        main = {line[12:16].strip(): [float(line[30:38]),
                                      float(line[38:46]),
                                      float(line[46:54])]
                for line in _atom_lines(data["pdb"])}
        checked = 0
        for view in data["highlights"]:
            for line in _atom_lines(view["pdb"]):
                coords = [float(line[30:38]), float(line[38:46]),
                          float(line[46:54])]
                assert coords == pytest.approx(main[line[12:16].strip()])
                checked += 1
        assert checked == 4

    def test_water_bridge_has_three_highlights(self, water_bridge_client,
                                               water_pdb):
        _open_structure(water_bridge_client, water_pdb)
        data = water_bridge_client.get(
            "/api/structure/view",
            params={"type": "water_bridge", "idx": 0, "frame": 0,
                    "solvent": "hide"}).json()
        assert [view["label"] for view in data["highlights"]] == [
            "ALA1 · H_donor", "SOL2 · water", "LYS3 (PDB 1) · H_acceptor"]
        assert [len(_atom_lines(view["pdb"]))
                for view in data["highlights"]] == [2, 2, 2]
        assert _atom_lines(data["highlights"][1]["pdb"])[0][12:16].strip() == (
            "OW")

    def test_inactive_frame_has_highlights_but_no_contact_pdb(self, tmp_path,
                                                              pdb_file):
        directory = tmp_path / "inactive_highlights"
        directory.mkdir()
        save_interactions([_interaction(active=False)],
                          str(directory / "hydrogen_bond.h5"))
        client = _client()
        client.post("/select", data={"path": str(directory)})
        _open_structure(client, pdb_file)
        data = client.get("/api/structure/view",
                          params={"type": "hydrogen_bond", "idx": 0,
                                  "frame": 0}).json()
        assert data["active"] is False
        assert data["contact_pdb"] is None
        assert [len(_atom_lines(view["pdb"])) for view in data["highlights"]] == [
            2, 2]

    def test_highlight_pdb_is_parseable(self, client, pdb_file):
        _open_structure(client, pdb_file)
        data = client.get("/api/structure/view",
                          params={"type": "hydrogen_bond", "idx": 0,
                                  "frame": 0}).json()
        for view in data["highlights"]:
            parsed = mda.Universe(io.StringIO(view["pdb"]), format="PDB")
            assert parsed.atoms.n_atoms == len(_atom_lines(view["pdb"]))
        contact = mda.Universe(io.StringIO(data["contact_pdb"]), format="PDB")
        assert contact.atoms.n_atoms == 2


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
