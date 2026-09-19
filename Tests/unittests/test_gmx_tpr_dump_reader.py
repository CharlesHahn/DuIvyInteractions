# -*- coding: utf-8 -*-
"""GmxTprDumpReader 单元测试：从真实 gmx dump 文本验证解析正确性。

测试数据：Tests/original_draft/dump_md_D927.tpr.txt
（RBD_D927_KRAS_GNP_Mg_MD 体系，gmx dump 输出，423092 行，已 git 追踪）

与 GmxTprReader（MDAnalysis）的已知差异（见 doc/future_optimizations.md）：
- 分子名不带 seg_ 前缀（如 "RBD_pro" vs "seg_0_RBD_pro"）
- 残基编号保留 PDB 原始编号（157~298），全局索引从 0 连续
- 能区分 bond / constrained / settle（GmxTprReader 全部标 "bond"）
"""

import pytest
from collections import Counter
from pathlib import Path

from DuIvyInteractions.system_readers import GmxTprDumpReader


# 测试数据路径（正斜杠：GmxTprDumpReader 用 split('/') 提取 system_name；
# 基于 __file__ 构造，保证任何 cwd 下均可运行）
DUMP_FILE = str(Path(__file__).parent.parent / "original_draft"
                / "dump_md_D927.tpr.txt").replace("\\", "/")
assert Path(DUMP_FILE).exists(), f"dump 文件缺失: {DUMP_FILE}"

# 真实数据已知值（由当前实现解析得到，数据文件已 git 追踪，稳定不变）
EXPECTED_MOLECULES = {
    "RBD_pro": 142, "D927": 1, "KRAS_pro": 167, "GNP_neg": 1,
    "Mg": 1, "SOL": 37021, "NA": 110, "CL": 108,
}
EXPECTED_N_RESIDUES = 37551
EXPECTED_N_ATOMS = 116383
EXPECTED_N_INTER_BONDS = 307
EXPECTED_BOND_TYPES = {"settle": 111063, "constrained": 2547, "bond": 2299}
EXPECTED_ATOM_COUNTS = {
    "RBD_pro": 2355, "D927": 53, "KRAS_pro": 2648, "GNP_neg": 45,
    "Mg": 1, "SOL": 111063, "NA": 110, "CL": 108,
}


@pytest.fixture(scope="module")
def system_data():
    """解析真实 dump 文本，返回 SystemData。"""
    return GmxTprDumpReader().read(DUMP_FILE)


class TestReaderMeta:
    """Reader 元信息。"""

    def test_name(self):
        assert GmxTprDumpReader().name == "gmx_tpr_dump"

    def test_returns_system_data(self, system_data):
        from DuIvyInteractions.core.datas import SystemData
        assert isinstance(system_data, SystemData)

    def test_system_name_from_filename(self, system_data):
        assert system_data.system_name == "dump_md_D927"


class TestSystemComposition:
    """体系组成（真实数据已知值）。"""

    def test_molecule_counts(self, system_data):
        counts = Counter(r.molecule_name for r in system_data.residues)
        assert dict(counts) == EXPECTED_MOLECULES

    def test_residue_count(self, system_data):
        assert len(system_data.residues) == EXPECTED_N_RESIDUES

    def test_atom_count(self, system_data):
        n = sum(len(r.atoms) for r in system_data.residues)
        assert n == EXPECTED_N_ATOMS

    def test_inter_residue_bond_count(self, system_data):
        assert len(system_data.inter_residue_bonds) == EXPECTED_N_INTER_BONDS

    def test_per_molecule_atom_counts(self, system_data):
        counts = {
            name: sum(len(r.atoms) for r in system_data.residues
                      if r.molecule_name == name)
            for name in EXPECTED_ATOM_COUNTS
        }
        assert counts == EXPECTED_ATOM_COUNTS


class TestBondTypes:
    """DumpReader 核心能力：区分 Bond / Constraint / Settle。"""

    def test_bond_type_distribution(self, system_data):
        counts = Counter(b.bond_type
                         for r in system_data.residues for b in r.bonds)
        assert dict(counts) == EXPECTED_BOND_TYPES

    def test_sol_has_three_settle_bonds(self, system_data):
        """SOL 水分子：SETTLE 展开为 3 条键（O-H1, O-H2, H1-H2）。"""
        sol = next(r for r in system_data.residues
                   if r.molecule_name == "SOL")
        assert len(sol.bonds) == 3
        assert all(b.bond_type == "settle" for b in sol.bonds)
        pairs = {tuple(sorted((b.atom1_idx_in_residue, b.atom2_idx_in_residue)))
                 for b in sol.bonds}
        assert pairs == {(0, 1), (0, 2), (1, 2)}

    def test_inter_residue_bonds_are_bond(self, system_data):
        """残基间键（肽键）为普通 bond 类型。"""
        assert all(ib.bond_type == "bond"
                   for ib in system_data.inter_residue_bonds)


class TestIndexInvariants:
    """全局索引唯一连续、局部索引与键引用有效。"""

    def test_atom_global_idx_unique_contiguous(self, system_data):
        idx = [a.atom_global_idx for r in system_data.residues
               for a in r.atoms]
        assert len(idx) == len(set(idx))
        assert min(idx) == 0
        assert max(idx) == EXPECTED_N_ATOMS - 1

    def test_residue_global_idx_unique_contiguous(self, system_data):
        idx = [r.residue_global_idx for r in system_data.residues]
        assert len(idx) == len(set(idx))
        assert sorted(idx) == list(range(EXPECTED_N_RESIDUES))

    def test_atom_local_idx_contiguous(self, system_data):
        for r in system_data.residues:
            local = sorted(a.atom_idx_in_residue for a in r.atoms)
            assert local == list(range(len(r.atoms)))

    def test_bond_refs_in_range(self, system_data):
        for r in system_data.residues:
            for b in r.bonds:
                assert b.atom1_idx_in_residue < len(r.atoms)
                assert b.atom2_idx_in_residue < len(r.atoms)

    def test_inter_bond_refs_in_range(self, system_data):
        n_res = len(system_data.residues)
        for ib in system_data.inter_residue_bonds:
            assert 0 <= ib.residue1_global_idx < n_res
            assert 0 <= ib.residue2_global_idx < n_res


class TestResidueDetails:
    """残基/原子字段抽查。"""

    def test_first_residue(self, system_data):
        r0 = system_data.residues[0]
        assert r0.residue_name == "ASN"
        assert r0.molecule_name == "RBD_pro"
        assert r0.residue_global_idx == 0
        assert len(r0.atoms) == 16

    def test_first_atom(self, system_data):
        a0 = system_data.residues[0].atoms[0]
        assert a0.atom_name == "N"
        assert a0.atom_type == "N3"
        assert a0.atom_element == "N"
        assert a0.atom_global_idx == 0

    def test_residue_idx_in_molecule_is_pdb_numbering(self, system_data):
        """dump 保留 PDB 原始编号（157 起），与 GmxTprReader 的 1~142 不同。"""
        rbd = [r for r in system_data.residues
               if r.molecule_name == "RBD_pro"]
        assert rbd[0].residue_idx_in_molecule == 157

    def test_molecule_names_no_seg_prefix(self, system_data):
        """与 GmxTprReader 的 seg_0_ 前缀不同（已知差异）。"""
        names = {r.molecule_name for r in system_data.residues}
        assert all(not n.startswith("seg_") for n in names)
