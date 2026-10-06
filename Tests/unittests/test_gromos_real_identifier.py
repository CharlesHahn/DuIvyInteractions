# -*- coding: utf-8 -*-
"""GROMOS 53A6 真实 tpr 基团识别集成测试。

使用 `Tests/test_MD_case_gromos/gromos53a6_md.tpr`（真实 GROMACS 拓扑：
130 残基蛋白 + 6 个配体 ZIN1-6），验证 GromosFFGroupIdentifier 对
真实力场数据的基团识别正确性（补充合成 SystemData 测试的覆盖）。

断言基线（已实测）：H_donor=200, H_acceptor=336, aromatic_ring=10,
charged_positive=10, charged_negative=15, hydrophobic=330, 总计 1584。
"""

import pytest
from pathlib import Path

from DuIvyInteractions.system_readers import GmxTprReader
from DuIvyInteractions.group_identifiers import GromosFFGroupIdentifier


TPR_FILE = Path(__file__).parent.parent / "test_MD_case_gromos" / "gromos53a6_md.tpr"
LIGAND_RESIDUES = {"1ZIN", "2ZIN", "3ZIN", "4ZIN", "5ZIN", "6ZIN"}


@pytest.fixture(scope="module")
def system_data():
    return GmxTprReader().read(str(TPR_FILE))


@pytest.fixture(scope="module")
def groups(system_data):
    return GromosFFGroupIdentifier().identify(system_data)


# ============================================================
# 体系读取
# ============================================================

class TestSystemRead:

    def test_residue_count(self, system_data):
        assert len(system_data.residues) == 136

    def test_atom_count(self, system_data):
        total = sum(len(r.atoms) for r in system_data.residues)
        assert total == 1311

    def test_ligand_residues_present(self, system_data):
        names = {r.residue_name for r in system_data.residues}
        assert LIGAND_RESIDUES <= names

    def test_no_water(self, system_data):
        """体系无水分（纯蛋白+配体）。"""
        names = {r.residue_name for r in system_data.residues}
        assert not names & {"SOL", "HOH", "WAT"}


# ============================================================
# 基团计数（真实 tpr 基线）
# ============================================================

class TestGroupCounts:

    def test_total_groups(self, groups):
        assert len(groups) == 1584

    def test_donors(self, groups):
        n = sum(1 for g in groups if g.group_type == "H_donor")
        assert n == 200

    def test_acceptors(self, groups):
        n = sum(1 for g in groups if g.group_type == "H_acceptor")
        assert n == 336

    def test_aromatic_rings(self, groups):
        n = sum(1 for g in groups if g.group_type == "aromatic_ring")
        assert n == 10

    def test_charged_positive(self, groups):
        n = sum(1 for g in groups if g.group_type == "charged_positive")
        assert n == 10

    def test_charged_negative(self, groups):
        n = sum(1 for g in groups if g.group_type == "charged_negative")
        assert n == 15

    def test_hydrophobic(self, groups):
        n = sum(1 for g in groups if g.group_type == "hydrophobic")
        assert n == 330


# ============================================================
# 蛋白/配体基团分布
# ============================================================

class TestProteinLigandSplit:

    def test_protein_ligand_group_split(self, groups):
        prot = [g for g in groups if g.residue_name not in LIGAND_RESIDUES]
        lig = [g for g in groups if g.residue_name in LIGAND_RESIDUES]
        assert len(prot) == 1410
        assert len(lig) == 174

    def test_ligand_has_no_aromatic_ring(self, groups):
        """配体芳香环不识别（GROMOS 芳香依赖蛋白残基白名单，配体 CR1 无强信号）。"""
        lig_aro = [g for g in groups
                   if g.residue_name in LIGAND_RESIDUES
                   and g.group_type == "aromatic_ring"]
        assert len(lig_aro) == 0

    def test_ligand_hydrophobic(self, groups):
        lig_hydro = [g for g in groups
                     if g.residue_name in LIGAND_RESIDUES
                     and g.group_type == "hydrophobic"]
        assert len(lig_hydro) == 30

    def test_ligand_acceptors(self, groups):
        lig_acc = [g for g in groups
                   if g.residue_name in LIGAND_RESIDUES
                   and g.group_type == "H_acceptor"]
        assert len(lig_acc) == 36


# ============================================================
# 带电基团语义抽查（方案A修复在真实 GROMOS 生效）
# ============================================================

class TestChargedSemantics:

    def test_nterm_leu_positive(self, groups):
        """N 端 LEU（NL 类型，+0.129）识别为正电（GROMOS UA：单原子即正电）。"""
        pos = [g for g in groups
               if g.group_type == "charged_positive" and g.residue_id == 0]
        assert len(pos) == 1
        assert pos[0].num_atoms == 1          # NL 单原子（UA 特性，3H 电荷分散不并入）

    def test_lys_nz_positive_not_acceptor(self, groups):
        """LYS 侧链铵（NL）为正电基团成员，不作受体（方案A修复）。"""
        lys_pos = [g for g in groups
                   if g.group_type == "charged_positive" and g.residue_name == "LYS"]
        assert len(lys_pos) == 5              # 5 段 × 1 LYS
        # 这些 NZ 不应出现在 H_acceptor 中
        lys_acc = [g for g in groups
                   if g.group_type == "H_acceptor" and g.residue_name == "LYS"
                   and g.atoms[0].atom_name == "NZ"]
        assert len(lys_acc) == 0

    def test_cterm_coo_negative(self, groups):
        """C 端 COO⁻ 识别为负电。"""
        neg = [g for g in groups
               if g.group_type == "charged_negative" and g.residue_id == 129]
        assert len(neg) == 1
        assert neg[0].num_atoms == 3          # C + 2O

    def test_glu_asp_negative(self, groups):
        neg = [g for g in groups
               if g.group_type == "charged_negative"
               and g.residue_name in ("GLU", "ASP")]
        assert len(neg) == 10                 # 5 段 × (1 GLU + 1 ASP)
