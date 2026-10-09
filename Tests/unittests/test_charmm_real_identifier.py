# -*- coding: utf-8 -*-
"""CHARMM36 真实 tpr 基团识别集成测试（SMO-BST）。

使用 Smoothened（SMO）+ β-谷甾醇（BST）体系，CHARMM36 + CGenFF + SPC/E 水
（GROMACS 2018.1，Omar et al. 2020, *Data in Brief*；Mendeley v94vzbwzf3）。

两套 tpr（注意差异）：
- ``SMO-BST_md_complex.tpr``：删水复合物——蛋白+配体合并为一个 molecule，**无水**；
- ``SMO-BST_md_fullatom.tpr``：全原子——蛋白/配体/水/离子分属独立 molecule，**含 63055 个 SOL 水**。

断言基线说明（2026-10-09）：
- 基团计数来自程序自动探测；其中**芳香环与蛋白残基组成精确对照**（PHE 32 + TYR 18
  + TRP 14 → 六元环 64；TRP 14 + HIS 8 → 五元环 22），可独立复核；
- 其余计数（供体/受体/疏水等）待人工核验，请谨慎引用。
"""

import pytest
from collections import Counter
from pathlib import Path

from DuIvyInteractions.system_readers import GmxTprReader
from DuIvyInteractions.group_identifiers.charmm_ff_identifier import (
    CharmmFFGroupIdentifier,
)


COMPLEX_TPR = Path(__file__).parent.parent / "test_MD_case_charmm36" / "SMO-BST_md_complex.tpr"
FULLATOM_TPR = Path(__file__).parent.parent / "test_MD_case_charmm36" / "SMO-BST_md_fullatom.tpr"

# 蛋白残基组成（与基团识别结果交叉对照）
PROTEIN_AROMATIC = {"PHE": 32, "TYR": 18, "TRP": 14, "HIS": 8}
CHARGED_RESIDUES = {"ARG": 28, "LYS": 14, "ASP": 16, "GLU": 23}


def _count(groups):
    return Counter(g.group_type for g in groups)


@pytest.fixture(scope="module")
def complex_data():
    sd = GmxTprReader().read(str(COMPLEX_TPR))
    return sd, CharmmFFGroupIdentifier().identify(sd)


@pytest.fixture(scope="module")
def fullatom_data():
    sd = GmxTprReader().read(str(FULLATOM_TPR))
    return sd, CharmmFFGroupIdentifier().identify(sd)


# ============================================================
# 删水复合物（SMO-BST_md_complex.tpr，无水）
# ============================================================

class TestComplexSystemRead:
    """复合物 tpr 解析：蛋白+配体，无水。"""

    def test_system_name(self, complex_data):
        sd, _ = complex_data
        assert "SMO-BST_md_complex" in sd.system_name

    def test_residue_count(self, complex_data):
        sd, _ = complex_data
        assert len(sd.residues) == 496          # 463 蛋白残基 + BST + 修饰残基

    def test_atom_count(self, complex_data):
        sd, _ = complex_data
        total = sum(len(r.atoms) for r in sd.residues)
        assert total == 7822

    def test_inter_residue_bonds(self, complex_data):
        sd, _ = complex_data
        assert len(sd.inter_residue_bonds) == 503

    def test_no_water(self, complex_data):
        """complex.tpr 为删水体系：无任何水残基、无水基团。"""
        sd, groups = complex_data
        resnames = {r.residue_name for r in sd.residues}
        assert not resnames & {"SOL", "HOH", "WAT", "TIP3"}
        counts = _count(groups)
        assert counts.get("water", 0) == 0

    def test_ligand_present(self, complex_data):
        """配体 BST（β-谷甾醇）存在。"""
        sd, _ = complex_data
        resnames = {r.residue_name for r in sd.residues}
        assert "BST" in resnames


class TestComplexGroupCounts:
    """复合物基团计数（对照蛋白残基组成）。"""

    def test_total_groups(self, complex_data):
        _, groups = complex_data
        assert len(groups) == 5677

    def test_donors(self, complex_data):
        _, groups = complex_data
        assert _count(groups)["H_donor"] == 824

    def test_acceptors(self, complex_data):
        _, groups = complex_data
        assert _count(groups)["H_acceptor"] == 752

    def test_aromatic_rings_match_protein_composition(self, complex_data):
        """芳香环 86 = 蛋白芳香残基精确对应：
        六元环 64 = PHE 32 + TYR 18 + TRP 14（苯环）；五元环 22 = TRP 14（吡咯）+ HIS 8（咪唑）。
        """
        _, groups = complex_data
        rings = [g for g in groups if g.group_type == "aromatic_ring"]
        assert len(rings) == 86
        sizes = Counter(len(g.atoms) for g in rings)
        assert sizes[6] == 64
        assert sizes[5] == 22

    def test_charged_positive(self, complex_data):
        _, groups = complex_data
        # 43 = ARG 28 + LYS 14 + PRO(0) N 端质子化氨基 1（已核实构成）
        assert _count(groups)["charged_positive"] == 43

    def test_charged_negative(self, complex_data):
        _, groups = complex_data
        # 40 = GLU 23 + ASP 16 + ARG(494) C 端羧基 1（已核实构成）
        assert _count(groups)["charged_negative"] == 40

    def test_hydrophobic(self, complex_data):
        _, groups = complex_data
        assert _count(groups)["hydrophobic"] == 1226

    def test_metal_binding(self, complex_data):
        _, groups = complex_data
        assert _count(groups)["metal_binding"] == 1374


class TestLigandBST:
    """配体 BST（β-谷甾醇，甾醇）识别。"""

    def test_ligand_atoms(self, complex_data):
        sd, _ = complex_data
        n = sum(1 for r in sd.residues if r.residue_name == "BST" for _ in r.atoms)
        assert n == 80

    def test_ligand_groups(self, complex_data):
        _, groups = complex_data
        lig = [g for g in groups if g.residue_name == "BST"]
        types = Counter(g.group_type for g in lig)
        # 甾醇：1 羟基供体、1 羟基受体、28 疏水碳、1 metal_binding(O)、1 halogen_acceptor(C1-O)
        assert types["H_donor"] == 1
        assert types["H_acceptor"] == 1
        assert types["hydrophobic"] == 28
        assert types["metal_binding"] == 1
        assert types["halogen_acceptor"] == 1

    def test_ligand_donor_structure(self, complex_data):
        """BST 唯一供体 = 甾醇羟基 O1-H38。"""
        _, groups = complex_data
        donors = [g for g in groups if g.residue_name == "BST"
                  and g.group_type == "H_donor"]
        assert len(donors) == 1
        assert donors[0].atoms[0].atom_name == "O1"      # D
        assert donors[0].atoms[1].atom_name == "H38"     # H

    def test_ligand_no_aromatic_ring(self, complex_data):
        """甾醇无芳香环——配体不应产生 aromatic_ring（正确）。"""
        _, groups = complex_data
        lig_aro = [g for g in groups if g.residue_name == "BST"
                   and g.group_type == "aromatic_ring"]
        assert len(lig_aro) == 0


# ============================================================
# 全原子体系（SMO-BST_md_fullatom.tpr，含水）
# ============================================================

class TestFullatomSystemRead:
    """全原子 tpr：蛋白/配体/水/离子独立分子，含水。"""

    def test_system_name(self, fullatom_data):
        sd, _ = fullatom_data
        assert "SMO-BST_md_fullatom" in sd.system_name

    def test_residue_count(self, fullatom_data):
        sd, _ = fullatom_data
        assert len(sd.residues) == 63555

    def test_atom_count(self, fullatom_data):
        sd, _ = fullatom_data
        total = sum(len(r.atoms) for r in sd.residues)
        assert total == 196991

    def test_molecule_split(self, fullatom_data):
        """全原子体系：蛋白/配体/水/离子分属独立 molecule。"""
        sd, _ = fullatom_data
        mols = sorted({r.molecule_name for r in sd.residues})
        assert any("BST" in m for m in mols)          # seg_1_BST
        assert any("SOL" in m for m in mols)          # seg_2_SOL
        assert any("CL" in m for m in mols)           # seg_3_CL
        assert any("Protein" in m for m in mols)      # seg_0_Protein_chain_A

    def test_water_residues(self, fullatom_data):
        sd, _ = fullatom_data
        waters = [r for r in sd.residues if r.residue_name in ("SOL", "HOH", "WAT", "TIP3")]
        assert len(waters) == 63055


class TestFullatomWater:
    """全原子体系的水识别（跨力场水排除/识别）。"""

    def test_water_groups(self, fullatom_data):
        """水被识别为 water 基团（63055 = SOL 残基数）。"""
        _, groups = fullatom_data
        assert _count(groups)["water"] == 63055

    def test_donors_include_water_h(self, fullatom_data):
        """H_donor = 蛋白/配体 824 + 水 O-H×2（63055×2）——精确吻合。"""
        _, groups = fullatom_data
        assert _count(groups)["H_donor"] == 824 + 63055 * 2

    def test_water_not_h_acceptor(self, fullatom_data):
        """CHARMM 水氧（OT 类型）不产生 H_acceptor——修复后的预期行为。"""
        _, groups = fullatom_data
        w_acc = [g for g in groups if g.group_type == "H_acceptor"
                 and g.residue_name in ("SOL", "HOH", "WAT", "TIP3")]
        assert len(w_acc) == 0

    def test_water_residues_in_identifier(self):
        """CHARMM 识别器 WATER_RESIDUES 含 TIP3（跨力场水排除基础）。"""
        assert CharmmFFGroupIdentifier.WATER_RESIDUES >= {"TIP3", "SOL", "HOH"}

    def test_protein_ligand_counts_match_complex(self, fullatom_data, complex_data):
        """全原子与删水复合物的蛋白/配体基团计数一致。

        注：水的 O-H 也会产生 H_donor 基团（残基名 SOL），故按"水残基"排除
        而非仅移除 water 类型。
        """
        _, g_full = fullatom_data
        _, g_cplx = complex_data
        water_res = CharmmFFGroupIdentifier.WATER_RESIDUES

        def non_water_counts(groups):
            return Counter(g.group_type for g in groups
                           if g.residue_name not in water_res)

        assert non_water_counts(g_full) == non_water_counts(g_cplx)