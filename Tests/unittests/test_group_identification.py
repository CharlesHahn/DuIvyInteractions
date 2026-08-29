# -*- coding: utf-8 -*-
"""基团识别集成测试：从真实 tpr 输入，验证完整的基团识别结果。

使用 D927 体系（RBD + D927 + KRAS + GNP + Mg + SOL + NA + CL）作为测试数据。
通过 GmxTprReader 读取 tpr 二进制文件（MDAnalysis），分子名带 seg_ 前缀。
"""

import pytest
from collections import Counter
from pathlib import Path

from DuIvyInteractions.input_readers import GmxTprReader
from DuIvyInteractions.group_identifiers import AmberFFGroupIdentifier


# 测试数据路径
TPR_FILE = Path(__file__).parent.parent / "test_MD_case" / "md.tpr"

# GmxTprReader 的分子名前缀
MOL_RBD = "seg_0_RBD_pro"
MOL_D927 = "seg_1_D927"
MOL_KRAS = "seg_2_KRAS_pro"
MOL_GNP = "seg_3_GNP_neg"
MOL_MG = "seg_4_Mg"
MOL_SOL = "seg_5_SOL"
MOL_NA = "seg_6_NA"
MOL_CL = "seg_7_CL"


@pytest.fixture(scope="module")
def system_data():
    """读取 tpr 文件，返回 SystemData。"""
    reader = GmxTprReader()
    return reader.read(str(TPR_FILE))


@pytest.fixture(scope="module")
def groups(system_data):
    """运行基团识别，返回 Group 列表。"""
    identifier = AmberFFGroupIdentifier()
    return identifier.identify(system_data)


@pytest.fixture(scope="module")
def group_counts(groups):
    """按类型统计基团数量。"""
    return Counter(g.group_type for g in groups)


@pytest.fixture(scope="module")
def mol_group_counts(groups):
    """按分子+类型统计基团数量。"""
    counts = {}
    for g in groups:
        counts.setdefault(g.molecule, Counter())[g.group_type] += 1
    return counts


# ============================================================
# SystemData 基础测试
# ============================================================

class TestSystemData:
    """验证 tpr 解析结果。"""

    def test_system_name(self, system_data):
        assert "md" in system_data.system_name

    def test_residue_count(self, system_data):
        assert system_data.n_residues == 37551

    def test_inter_residue_bonds(self, system_data):
        assert len(system_data.inter_residue_bonds) == 307

    def test_molecule_types(self, system_data):
        mol_names = {r.molecule_name for r in system_data.residues}
        expected = {
            MOL_RBD, MOL_D927, MOL_KRAS, MOL_GNP,
            MOL_MG, MOL_SOL, MOL_NA, MOL_CL,
        }
        assert mol_names == expected


# ============================================================
# 总体基团计数测试
# ============================================================

class TestGroupCounts:
    """验证各类型基团的总数。"""

    def test_total_groups(self, group_counts):
        total = sum(group_counts.values())
        assert total == 152385

    def test_H_donor(self, group_counts):
        assert group_counts["H_donor"] == 74623

    def test_H_donor_structure(self, groups):
        """H_donor atoms=[D, H]，第一个是 D，第二个是 H，metadata 为空。"""
        donors = [g for g in groups if g.group_type == "H_donor"]
        sample = donors[0]
        assert len(sample.atoms) == 2
        assert sample.atoms[0].atom_element in ("N", "O", "S", "F")  # D
        assert sample.atoms[1].atom_element == "H"                    # H
        assert sample.metadata == {}

    def test_H_acceptor(self, group_counts):
        assert group_counts["H_acceptor"] == 37954

    def test_water(self, group_counts):
        assert group_counts["water"] == 37021

    def test_aromatic_ring(self, group_counts):
        assert group_counts["aromatic_ring"] == 38

    def test_charged_positive(self, group_counts):
        assert group_counts["charged_positive"] == 42

    def test_charged_negative(self, group_counts):
        assert group_counts["charged_negative"] == 43

    def test_halogen_acceptor(self, group_counts):
        assert group_counts["halogen_acceptor"] == 881

    def test_halogen_donor(self, group_counts):
        assert group_counts["halogen_donor"] == 1

    def test_hydrophobic(self, group_counts):
        assert group_counts["hydrophobic"] == 734

    def test_metal(self, group_counts):
        assert group_counts["metal"] == 111


# ============================================================
# D927 分子测试
# ============================================================

class TestD927:
    """验证 D927 配体的基团识别。"""

    def test_aromatic_rings_count(self, mol_group_counts):
        assert mol_group_counts[MOL_D927]["aromatic_ring"] == 3

    def test_aromatic_ring_sizes(self, groups):
        rings = [g for g in groups
                 if g.molecule == MOL_D927 and g.group_type == "aromatic_ring"]
        sizes = sorted(len(g.atoms) for g in rings)
        assert sizes == [6, 6, 6]

    def test_aromatic_ring_atom_types(self, groups):
        """三个芳香环的原子类型应全部为芳香类型。"""
        from DuIvyInteractions.group_identifiers.amber_ff_identifier import STRONG_AROMATIC
        rings = [g for g in groups
                 if g.molecule == MOL_D927 and g.group_type == "aromatic_ring"]
        for ring in rings:
            types = [a.atom_type for a in ring.atoms]
            strong = sum(1 for t in types if t in STRONG_AROMATIC)
            assert strong >= len(types) - 1

    def test_halogen_donor(self, mol_group_counts):
        assert mol_group_counts[MOL_D927]["halogen_donor"] == 1

    def test_halogen_donor_structure(self, groups):
        """卤键供体 atoms=[C, X]，第一个是碳，第二个是卤素。"""
        donors = [g for g in groups
                  if g.molecule == MOL_D927 and g.group_type == "halogen_donor"]
        assert len(donors) == 1
        d = donors[0]
        assert len(d.atoms) == 2
        assert d.atoms[0].atom_element == "C"   # 碳
        assert d.atoms[1].atom_element == "F"   # 卤素
        assert d.metadata == {}

    def test_hydrophobic(self, mol_group_counts):
        assert mol_group_counts[MOL_D927]["hydrophobic"] == 12


# ============================================================
# 蛋白质芳香环测试
# ============================================================

class TestProteinAromaticRings:
    """验证蛋白质的芳香环识别。"""

    def test_RBD_ring_count(self, mol_group_counts):
        assert mol_group_counts[MOL_RBD]["aromatic_ring"] == 17

    def test_KRAS_ring_count(self, mol_group_counts):
        assert mol_group_counts[MOL_KRAS]["aromatic_ring"] == 18

    def test_RBD_ring_residues(self, groups):
        """RBD 芳香环应来自 HIS/TYR/TRP/PHE。"""
        rings = [g for g in groups
                 if g.molecule == MOL_RBD and g.group_type == "aromatic_ring"]
        residue_names = {g.residue_name for g in rings}
        assert residue_names == {"HIS", "TYR", "TRP", "PHE"}

    def test_RBD_HIS_count(self, groups):
        his_rings = [g for g in groups
                     if g.molecule == MOL_RBD
                     and g.group_type == "aromatic_ring"
                     and g.residue_name == "HIS"]
        assert len(his_rings) == 3

    def test_RBD_TYR_count(self, groups):
        tyr_rings = [g for g in groups
                     if g.molecule == MOL_RBD
                     and g.group_type == "aromatic_ring"
                     and g.residue_name == "TYR"]
        assert len(tyr_rings) == 11

    def test_RBD_TRP_has_two_rings(self, groups):
        trp_rings = [g for g in groups
                     if g.molecule == MOL_RBD
                     and g.group_type == "aromatic_ring"
                     and g.residue_name == "TRP"]
        assert len(trp_rings) == 2
        sizes = sorted(len(g.atoms) for g in trp_rings)
        assert sizes == [5, 6]

    def test_HIS_ring_size(self, groups):
        his_rings = [g for g in groups
                     if g.group_type == "aromatic_ring"
                     and g.residue_name == "HIS"]
        for ring in his_rings:
            assert len(ring.atoms) == 5

    def test_TYR_ring_size(self, groups):
        tyr_rings = [g for g in groups
                     if g.group_type == "aromatic_ring"
                     and g.residue_name == "TYR"]
        for ring in tyr_rings:
            assert len(ring.atoms) == 6


# ============================================================
# 蛋白质带电基团测试
# ============================================================

class TestProteinChargedGroups:
    """验证蛋白质的带电基团识别。"""

    def test_RBD_layer1_counts(self, groups):
        layer1 = [g for g in groups
                  if g.molecule == MOL_RBD
                  and g.group_type in ("charged_positive", "charged_negative")
                  and g.metadata.get("source") == "residue_name"]
        res_counts = Counter(g.residue_name for g in layer1)
        assert res_counts["ARG"] == 5
        assert res_counts["LYS"] == 14
        assert res_counts["ASP"] == 4
        assert res_counts["GLU"] == 9

    def test_KRAS_layer1_counts(self, groups):
        layer1 = [g for g in groups
                  if g.molecule == MOL_KRAS
                  and g.group_type in ("charged_positive", "charged_negative")
                  and g.metadata.get("source") == "residue_name"]
        res_counts = Counter(g.residue_name for g in layer1)
        assert res_counts["ARG"] == 10
        assert res_counts["LYS"] == 11
        assert res_counts["ASP"] == 15
        assert res_counts["GLU"] == 12

    def test_N_terminus_detected(self, groups):
        n_terms = [g for g in groups
                   if g.metadata.get("func_group") == "tertamine"]
        molecules = {g.molecule for g in n_terms}
        assert MOL_RBD in molecules
        assert MOL_KRAS in molecules

    def test_C_terminus_detected(self, groups):
        c_terms = [g for g in groups
                   if g.metadata.get("func_group") == "carboxylate"]
        molecules = {g.molecule for g in c_terms}
        assert MOL_RBD in molecules
        assert MOL_KRAS in molecules

    def test_KRAS_LYS309_dual_groups(self, groups):
        lys309 = [g for g in groups
                  if g.molecule == MOL_KRAS
                  and g.residue_name == "LYS"
                  and g.residue_id == 309]
        types = {g.group_type for g in lys309}
        assert "charged_positive" in types
        assert "charged_negative" in types

    def test_charged_group_net_charge_direction(self, groups):
        for g in groups:
            if g.group_type == "charged_positive":
                assert g.net_charge > 0.1, f"{g} net_charge={g.net_charge}"
            elif g.group_type == "charged_negative":
                assert g.net_charge < -0.1, f"{g} net_charge={g.net_charge}"


# ============================================================
# 金属离子测试
# ============================================================

class TestMetals:
    """验证金属离子检测。"""

    def test_Mg_detected(self, groups):
        metals = [g for g in groups
                  if g.group_type == "metal"
                  and g.atoms[0].atom_element == "Mg"]
        assert len(metals) == 1
        assert metals[0].atoms[0].atom_charge == 2.0

    def test_Na_count(self, groups):
        na_metals = [g for g in groups
                     if g.group_type == "metal"
                     and g.atoms[0].atom_element == "Na"]
        assert len(na_metals) == 110

    def test_total_metal(self, group_counts):
        assert group_counts["metal"] == 111


# ============================================================
# 水分子测试
# ============================================================

class TestWater:
    """验证水分子检测。"""

    def test_water_count(self, mol_group_counts):
        assert mol_group_counts[MOL_SOL]["water"] == 37021

    def test_water_H_donor(self, mol_group_counts):
        assert mol_group_counts[MOL_SOL]["H_donor"] == 74042

    def test_water_H_acceptor(self, mol_group_counts):
        assert mol_group_counts[MOL_SOL]["H_acceptor"] == 37021


# ============================================================
# 金属配位原子测试
# ============================================================

class TestMetalBinding:
    """验证金属配位原子检测。"""

    def test_metal_binding_exists(self, group_counts):
        """应检测到 metal_binding 基团。"""
        assert group_counts["metal_binding"] > 0

    def test_metal_binding_structure(self, groups):
        """每个 metal_binding 基团应有 1 个原子，元素为 O/N/S。"""
        bindings = [g for g in groups if g.group_type == "metal_binding"]
        for g in bindings:
            assert len(g.atoms) == 1
            assert g.atoms[0].atom_element in ("O", "N", "S")

    def test_metal_binding_source(self, groups):
        """metadata.source 应为 element。"""
        bindings = [g for g in groups if g.group_type == "metal_binding"]
        for g in bindings:
            assert g.metadata["source"] == "element"

    def test_water_excluded(self, groups):
        """水分子不应被识别为 metal_binding。"""
        water_binding = [g for g in groups
                         if g.group_type == "metal_binding"
                         and g.residue_name in ("SOL", "HOH", "WAT")]
        assert len(water_binding) == 0
