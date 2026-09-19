# -*- coding: utf-8 -*-
"""两个 Reader（GmxTprReader / GmxTprDumpReader）一致性单元测试。

同一体系（KRAS-RBD D927）分别经：
- GmxTprReader 读二进制 tpr（Tests/test_MD_case/md.tpr，MDAnalysis）
- GmxTprDumpReader 读 gmx dump 文本（Tests/original_draft/dump_md_D927.tpr.txt）

TestConsistentAcrossReaders：两 Reader 输出必须一致的维度（同源体系，
不一致说明某个 Reader 解析漂移）。
TestKnownDifferences：doc/future_optimizations.md 记录的已知差异，
用测试锁定现状，将来修复差异时需同步更新断言。
"""

import re
from collections import Counter
from pathlib import Path

import pytest

from DuIvyInteractions.system_readers import GmxTprReader, GmxTprDumpReader


TPR_FILE = Path(__file__).parent.parent / "test_MD_case" / "md.tpr"
DUMP_FILE = str(Path(__file__).parent.parent / "original_draft"
                / "dump_md_D927.tpr.txt").replace("\\", "/")

# 已知真实值（两个 Reader 已实测一致）
N_RESIDUES = 37551
N_ATOMS = 116383
N_INTER_BONDS = 307

# 已知差异的实测值
RBD_FIRST_PDB_NR_TPR = 1      # MDA 连续编号
RBD_FIRST_PDB_NR_DUMP = 157   # PDB 原始编号


def _strip_seg_prefix(name: str) -> str:
    """去掉 GmxTprReader 的 seg_N_ 前缀（如 seg_0_RBD_pro → RBD_pro）。"""
    return re.sub(r"^seg_\d+_", "", name)


@pytest.fixture(scope="module")
def tpr_sd():
    """GmxTprReader 读二进制 tpr。"""
    return GmxTprReader().read(str(TPR_FILE))


@pytest.fixture(scope="module")
def dump_sd():
    """GmxTprDumpReader 读 dump 文本。"""
    return GmxTprDumpReader().read(DUMP_FILE)


class TestConsistentAcrossReaders:
    """两 Reader 对同一体系必须一致的维度。"""

    def test_total_residues_match(self, tpr_sd, dump_sd):
        assert len(tpr_sd.residues) == len(dump_sd.residues) == N_RESIDUES

    def test_total_atoms_match(self, tpr_sd, dump_sd):
        n_tpr = sum(len(r.atoms) for r in tpr_sd.residues)
        n_dump = sum(len(r.atoms) for r in dump_sd.residues)
        assert n_tpr == n_dump == N_ATOMS

    def test_inter_residue_bond_count_match(self, tpr_sd, dump_sd):
        assert len(tpr_sd.inter_residue_bonds) == \
            len(dump_sd.inter_residue_bonds) == N_INTER_BONDS

    def test_molecule_composition_matches(self, tpr_sd, dump_sd):
        """去掉 seg_ 前缀后，分子类型与残基数分布一致。"""
        tpr_counts = Counter(_strip_seg_prefix(r.molecule_name)
                             for r in tpr_sd.residues)
        dump_counts = Counter(r.molecule_name for r in dump_sd.residues)
        assert dict(tpr_counts) == dict(dump_counts)

    def test_per_molecule_atom_counts_match(self, tpr_sd, dump_sd):
        """每个分子的原子数一致（去前缀匹配）。"""
        tpr_atoms = {}
        for r in tpr_sd.residues:
            name = _strip_seg_prefix(r.molecule_name)
            tpr_atoms[name] = tpr_atoms.get(name, 0) + len(r.atoms)
        dump_atoms = {}
        for r in dump_sd.residues:
            dump_atoms[r.molecule_name] = dump_atoms.get(r.molecule_name, 0) \
                + len(r.atoms)
        assert tpr_atoms == dump_atoms

    def test_residue_global_idx_range_match(self, tpr_sd, dump_sd):
        """全局残基索引均从 0 起连续。"""
        tpr_idx = sorted(r.residue_global_idx for r in tpr_sd.residues)
        dump_idx = sorted(r.residue_global_idx for r in dump_sd.residues)
        assert tpr_idx == dump_idx


class TestKnownDifferences:
    """doc/future_optimizations.md 记录的已知差异（锁定现状）。"""

    def test_molecule_name_prefix_differs(self, tpr_sd, dump_sd):
        """tpr 带 seg_N_ 前缀，dump 不带。"""
        tpr_names = {r.molecule_name for r in tpr_sd.residues}
        dump_names = {r.molecule_name for r in dump_sd.residues}
        assert all(n.startswith("seg_") for n in tpr_names)
        assert all(not n.startswith("seg_") for n in dump_names)

    def test_residue_numbering_differs(self, tpr_sd, dump_sd):
        """residue_idx_in_molecule：tpr 为 MDA 连续编号，dump 为 PDB 原始编号。"""
        tpr_rbd = [r for r in tpr_sd.residues
                   if _strip_seg_prefix(r.molecule_name) == "RBD_pro"]
        dump_rbd = [r for r in dump_sd.residues
                    if r.molecule_name == "RBD_pro"]
        assert tpr_rbd[0].residue_idx_in_molecule == RBD_FIRST_PDB_NR_TPR
        assert dump_rbd[0].residue_idx_in_molecule == RBD_FIRST_PDB_NR_DUMP

    def test_sol_bond_count_differs(self, tpr_sd, dump_sd):
        """SOL 键数：tpr 2 条（O-H1, O-H2），dump 3 条（SETTLE 展开含 H1-H2）。"""
        tpr_sol = next(r for r in tpr_sd.residues
                       if r.residue_name == "SOL" and "SOL" in r.molecule_name)
        dump_sol = next(r for r in dump_sd.residues
                        if r.residue_name == "SOL" and r.molecule_name == "SOL")
        assert len(tpr_sol.bonds) == 2
        assert len(dump_sol.bonds) == 3

    def test_bond_type_granularity_differs(self, tpr_sd, dump_sd):
        """键类型：tpr 全部 bond；dump 区分 bond/constrained/settle。"""
        tpr_types = {b.bond_type for r in tpr_sd.residues
                     for b in r.bonds}
        dump_types = {b.bond_type for r in dump_sd.residues
                      for b in r.bonds}
        assert tpr_types == {"bond"}
        assert dump_types >= {"bond", "constrained", "settle"}
