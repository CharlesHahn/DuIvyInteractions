# -*- coding: utf-8 -*-
"""单对数据接口 GET /api/pair/{name}/{idx}：原子全局索引 + 页面同源 payload。

用真实数据：h5（测试数据目录）+ amber tpr/xtc（D927 结构）。
"""

import json
import re
from pathlib import Path

import numpy as np
import pytest

from DuIvyInteractions.core.datas import AtomData, Group, Interaction
from DuIvyInteractions.io.h5 import load_interactions
from DuIvyInteractions.visualizers import pair_summary
from DuIvyInteractions.visualizers.data import PT_TYPES
from DuIvyInteractions.visualizers.structure import (
    interaction_expressions, load_universe)

H5_DATA_DIR = Path(__file__).resolve().parents[1] / "interaction_h5data"
D927_DIR = Path(__file__).resolve().parents[1] / "test_MD_case_amber"
D927_TOPOLOGY = D927_DIR / "md.tpr"
D927_TRAJECTORY = D927_DIR / "md1ns.xtc"
PAIR_KEYS = {"pair_label", "idx", "n_pairs", "n_frames", "flags", "existence",
             "pt_types", "groups", "summary", "figure", "topology"}
GROUP_KEYS = {"label", "residue_selector", "is_water", "expression",
              "group_expression", "indices", "residue_indices", "atoms"}
ATOM_KEYS = {"name", "global_idx", "expression"}
FIGURE_PATTERN = r'<script id="pairFigure"[^>]*>(.*?)</script>'


@pytest.fixture(autouse=True)
def isolate_recent(monkeypatch, tmp_path):
    monkeypatch.setattr("DuIvyInteractions.visualizers.data.RECENT_PATH",
                        tmp_path / "recent.json")


def _client():
    from fastapi.testclient import TestClient
    from DuIvyInteractions.web.app import create_app
    return TestClient(create_app())


def _select(client, directory: Path) -> None:
    client.post("/select", data={"path": str(directory)})


@pytest.fixture
def bare_client():
    client = _client()
    _select(client, H5_DATA_DIR)
    return client


@pytest.fixture
def client():
    client = _client()
    _select(client, H5_DATA_DIR)
    response = client.post("/api/structure/select", data={
        "topology": str(D927_TOPOLOGY), "trajectory": str(D927_TRAJECTORY)})
    assert response.status_code == 200
    return client


@pytest.fixture(scope="module")
def amber_universe():
    return load_universe(str(D927_TOPOLOGY), str(D927_TRAJECTORY))


def _pair(client, type_name: str = "pi_stacking", idx: int = 0) -> dict:
    response = client.get(f"/api/pair/{type_name}/{idx}")
    assert response.status_code == 200
    return response.json()


def _h5_interaction(type_name: str):
    return load_interactions(str(H5_DATA_DIR / f"{type_name}.h5"))[0]


# ============================================================
# 载荷形状
# ============================================================


class TestPayloadShape:
    def test_key_set_is_exact(self, client):
        assert set(_pair(client)) == PAIR_KEYS

    def test_pair_metadata(self, client):
        data = _pair(client)
        interaction = _h5_interaction("pi_stacking")
        assert data["idx"] == 0
        assert data["n_pairs"] == interaction.n_pairs
        assert data["n_frames"] == interaction.n_frames
        assert data["pair_label"]
        assert data["flags"] == sorted(interaction.flags)

    def test_topology_block_with_structure(self, client, amber_universe):
        topology = _pair(client)["topology"]
        assert topology == {"has_structure": True,
                            "n_atoms": int(amber_universe.atoms.n_atoms),
                            "n_frames": len(amber_universe.trajectory)}

    def test_group_and_atom_key_sets(self, client):
        for group in _pair(client)["groups"]:
            assert set(group) == GROUP_KEYS
            for atom in group["atoms"]:
                assert set(atom) == ATOM_KEYS


# ============================================================
# 原子索引：indices / residue_indices / atoms[].global_idx
# ============================================================


class TestIndices:
    def test_group_indices_match_h5_order(self, client):
        data = _pair(client)
        interaction = _h5_interaction("pi_stacking")
        assert [group["indices"] for group in data["groups"]] == [
            [int(atom.atom_global_idx) for atom in group.atoms]
            for group in interaction.groups[0]]

    def test_atom_global_idx_matches_h5(self, client):
        data = _pair(client)
        interaction = _h5_interaction("pi_stacking")
        for group, source in zip(data["groups"], interaction.groups[0]):
            assert [atom["global_idx"] for atom in group["atoms"]] == [
                int(atom.atom_global_idx) for atom in source.atoms]

    def test_residue_indices_cover_one_residue(self, client, amber_universe):
        for group in _pair(client)["groups"]:
            indices = group["residue_indices"]
            assert indices is not None
            resindices = {int(amber_universe.atoms[index].resindex)
                          for index in indices}
            assert len(resindices) == 1
            assert set(group["indices"]).issubset(set(indices))

    def test_residue_indices_are_the_whole_residue(self, client,
                                                   amber_universe):
        for group in _pair(client)["groups"]:
            first = int(group["indices"][0])
            resindex = int(amber_universe.atoms[first].resindex)
            expected = sorted(
                int(index) for index in amber_universe.atoms[
                    amber_universe.atoms.resindices == resindex].indices)
            assert group["residue_indices"] == expected

    def test_water_group_has_three_atoms(self, client):
        data = _pair(client, "water_bridge")
        water = [group for group in data["groups"] if group["is_water"]]
        assert len(water) == 1
        assert len(water[0]["indices"]) == 3
        assert set(water[0]["residue_indices"]) == set(water[0]["indices"])


# ============================================================
# 帧数据：existence / pt_types
# ============================================================


class TestFrames:
    def test_existence_matches_h5(self, client):
        data = _pair(client)
        interaction = _h5_interaction("pi_stacking")
        assert data["existence"] == [
            int(value) for value in interaction.existence[0]]

    def test_pi_stacking_pt_types_per_frame(self, client):
        data = _pair(client)
        assert len(data["pt_types"]) == data["n_frames"]
        assert set(data["pt_types"]) <= set(PT_TYPES)

    def test_salt_bridge_has_no_pt_types(self, client):
        assert _pair(client, "salt_bridge")["pt_types"] == []


# ============================================================
# 与页面渲染同源：summary / figure
# ============================================================


class TestPageConsistency:
    def test_figure_matches_page_payload(self, client):
        data = _pair(client)
        html = client.get("/pair/pi_stacking/0").text
        match = re.search(FIGURE_PATTERN, html, re.S)
        assert match
        assert json.loads(match.group(1)) == data["figure"]

    def test_summary_matches_shared_builder(self, client):
        data = _pair(client)
        assert data["summary"] == pair_summary(_h5_interaction("pi_stacking"), 0)
        html = client.get("/pair/pi_stacking/0").text
        occupancy = "%.1f%%" % (data["summary"]["occupancy"] * 100)
        assert occupancy in html


# ============================================================
# 无结构退化
# ============================================================


class TestWithoutStructure:
    def test_topology_reports_no_structure(self, bare_client):
        topology = _pair(bare_client)["topology"]
        assert topology["has_structure"] is False
        assert topology["n_atoms"] == 0
        assert topology["n_frames"] == _h5_interaction("pi_stacking").n_frames

    def test_indices_still_available_residue_indices_none(self, bare_client):
        groups = _pair(bare_client)["groups"]
        assert groups
        for group in groups:
            assert group["indices"]
            assert group["residue_indices"] is None
            assert all(isinstance(atom["global_idx"], int)
                       for atom in group["atoms"])


# ============================================================
# 错误与既有接口兼容
# ============================================================


class TestErrorsAndCompat:
    def test_out_of_range_returns_404(self, client):
        interaction = _h5_interaction("pi_stacking")
        response = client.get(f"/api/pair/pi_stacking/{interaction.n_pairs}")
        assert response.status_code == 404

    def test_unknown_type_returns_404(self, client):
        assert client.get("/api/pair/nope/0").status_code == 404

    def test_structure_interactions_stays_compatible(self, client):
        response = client.get("/api/structure/interactions",
                              params={"type": "pi_stacking", "idx": 0})
        assert response.status_code == 200
        data = response.json()
        assert set(data) == {"pair_label", "n_frames", "existence",
                             "pt_types", "residues", "n_atoms"}
        group = data["residues"][0]
        assert GROUP_KEYS <= set(group)
        assert ATOM_KEYS <= set(group["atoms"][0])
        assert group["indices"] and group["residue_indices"]


# ============================================================
# group_expression：按 chain/resi + 原子名圈选（与 pair 无关）
# ============================================================


def _atom(global_idx: int, name: str) -> AtomData:
    return AtomData(atom_global_idx=global_idx, atom_idx_in_residue=0,
                    atom_name=name, atom_type="c3", atom_element="C",
                    atom_charge=0.0, atom_mass=0.0)


def _group(group_id: int, atoms) -> Group:
    return Group(group_id=group_id, group_type="hydrophobic", molecule="test",
                 residue_name="ALA", residue_id=1, atoms=list(atoms))


def _duplicate_name_pdb(path: Path) -> Path:
    """一个残基内 3 个同名原子（name 选择器无法区分）。"""
    lines = ["CRYST1   10.000   10.000   10.000  90.00  90.00  90.00"
             " P 1           1"]
    for serial, (name, x) in enumerate([("C1", 0.0), ("C1", 1.0),
                                        ("C1", 2.0)], start=1):
        lines.append(f"ATOM  {serial:5d} {name:>4s} ALA A   1    "
                     f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00          C")
    path.write_text("\n".join(lines + ["END"]) + "\n", encoding="utf-8")
    return path


def _duplicate_interaction() -> Interaction:
    groups = [(_group(0, [_atom(0, "C1"), _atom(1, "C1")]),
               _group(1, [_atom(2, "C1")]))]
    return Interaction(interaction_type="hydrophobic", groups=groups,
                       existence=np.ones((1, 1), dtype=bool),
                       metrics={"distance": np.full((1, 1), 3.0)},
                       times=np.zeros(1))


class TestGroupExpression:
    def test_starts_with_selector_and_lists_atom_names(self, client):
        for group in _pair(client)["groups"]:
            expression = group["group_expression"]
            assert expression.startswith(group["residue_selector"])
            assert "name " in expression
            assert "b >" not in expression

    def test_water_group_uses_full_number_selector(self, client):
        """mmCIF 不再截断：水组表达式即 chain W and resi 15689（无 b 区间）。"""
        water = [group for group in _pair(client, "water_bridge")["groups"]
                 if group["is_water"]]
        assert len(water) == 1
        assert water[0]["group_expression"] == water[0]["residue_selector"]
        assert water[0]["group_expression"] == "chain W and resi 15689"
        assert "b >" not in water[0]["group_expression"]

    def test_none_without_structure(self, bare_client):
        groups = _pair(bare_client)["groups"]
        assert groups
        assert all(group["group_expression"] is None for group in groups)

    def test_duplicate_atom_names_fall_back_to_b_bands(self, tmp_path):
        pdb = _duplicate_name_pdb(tmp_path / "dup.pdb")
        universe = load_universe(str(pdb), str(pdb))
        row = interaction_expressions(
            _duplicate_interaction(), 0, universe)["residues"][0]
        assert row["group_expression"] == row["expression"]
        assert "b >" in row["group_expression"]
