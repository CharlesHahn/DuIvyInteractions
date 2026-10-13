# -*- coding: utf-8 -*-
"""用户基团选择与过滤（--select / --pairing）单元 + 真实数据集成测试。

真实数据：
- Tests/test_MD_case_amber/md.tpr + md1ns.xtc（KRAS-RBD D927，8 segid）
- Tests/test_MD_case_charmm36/SMO-BST_md_fullatom.tpr（SMO-BST，4 segid）
"""

from pathlib import Path

import pytest

from DuIvyInteractions.core.datas import AtomData, Group
from DuIvyInteractions.selection import (
    WATER_SIDE, CONSTRAINT_POSITIONS,
    resolve_selection, validate_selections, assign_groups,
    filter_groups, build_tuple_filter,
)
from DuIvyInteractions.pipeline import Pipeline
from DuIvyInteractions.system_readers import GmxTprReader
from DuIvyInteractions.group_identifiers import (
    AmberFFGroupIdentifier, CharmmFFGroupIdentifier,
)
from DuIvyInteractions.io.h5 import load_interactions


AMBER_TPR = Path(__file__).parent.parent / "test_MD_case_amber" / "md.tpr"
AMBER_XTC = Path(__file__).parent.parent / "test_MD_case_amber" / "md1ns.xtc"
CHARMM_TPR = (Path(__file__).parent.parent
              / "test_MD_case_charmm36" / "SMO-BST_md_complex.tpr")
CHARMM_XTC = (Path(__file__).parent.parent
              / "test_MD_case_charmm36" / "SMO-BST_md_complex.xtc")


# ============================================================
# 合成数据辅助
# ============================================================

def _atom(gidx: int) -> AtomData:
    return AtomData(atom_global_idx=gidx, atom_idx_in_residue=0,
                    atom_name="X", atom_type="", atom_element="X",
                    atom_charge=0.0, atom_mass=0.0)


def _group(gid: int, gtype: str, res_name: str,
           idxs: list, molecule: str = "mol") -> Group:
    return Group(group_id=gid, group_type=gtype, molecule=molecule,
                 residue_name=res_name, residue_id=1,
                 atoms=[_atom(i) for i in idxs])


# ============================================================
# 单元测试：resolve / validate
# ============================================================

class TestResolveSelection:
    """选择串求值（基于真实 amber tpr）。"""

    @pytest.fixture(scope="class")
    def universe(self):
        import MDAnalysis as mda
        return mda.Universe(str(AMBER_TPR), str(AMBER_XTC))

    def test_resolve_syntax_error_raises(self, universe):
        with pytest.raises(ValueError, match="语法错误"):
            resolve_selection(universe, "no_such_keyword foo")

    def test_resolve_segid(self, universe):
        s = resolve_selection(universe, "segid seg_0_RBD_pro")
        assert isinstance(s, frozenset) and s
        assert all(idx < universe.atoms.n_atoms for idx in s)


class TestValidateSelections:
    """空/重叠校验。"""

    def test_empty_select_raises(self):
        with pytest.raises(ValueError, match="未选中任何原子"):
            validate_selections([frozenset(), frozenset({1})])

    def test_overlap_raises(self):
        with pytest.raises(ValueError, match="重叠"):
            validate_selections([frozenset({1, 2}), frozenset({2, 3})])

    def test_disjoint_ok(self):
        validate_selections([frozenset({1, 2}), frozenset({3, 4})])


# ============================================================
# 单元测试：assign / filter
# ============================================================

class TestAssignGroups:
    """基团归属：全部原子 ⊆ 某 select；水侧恒真；无归属不保留。"""

    def _groups(self):
        return [
            _group(0, "H_donor", "LYS", [0, 1]),
            _group(1, "H_acceptor", "ASP", [2, 3]),
            _group(2, "water", "SOL", [4, 5]),
            _group(3, "H_donor", "SER", [6, 7]),
        ]

    def test_assign_disjoint(self):
        sets = [frozenset({0, 1, 2, 3}), frozenset({6, 7})]
        a = assign_groups(self._groups(), sets, frozenset({"SOL"}))
        assert a == {0: 0, 1: 0, 2: WATER_SIDE, 3: 1}

    def test_water_always_water_side(self):
        # 即使 select 含 SOL 原子，water 基团仍恒属水侧
        sets = [frozenset({0, 1, 2, 3, 4}), frozenset({6, 7})]
        a = assign_groups(self._groups(), sets, frozenset({"SOL"}))
        assert a[2] == WATER_SIDE

    def test_partial_group_not_kept(self):
        # 基团只有部分原子在 select 内 → 不保留
        sets = [frozenset({0, 2})]  # g0 atoms={0,1} → 原子 1 不在 → 不全
        a = assign_groups(self._groups(), sets, frozenset({"SOL"}))
        assert 0 not in a

    def test_filter_groups_keeps_only_assigned(self):
        groups = self._groups()
        sets = [frozenset({0, 1, 2, 3}), frozenset({6, 7})]
        a = assign_groups(groups, sets, frozenset({"SOL"}))
        kept = filter_groups(groups, a)
        assert {g.group_id for g in kept} == {0, 1, 2, 3}


# ============================================================
# 单元测试：build_tuple_filter
# ============================================================

class TestBuildTupleFilter:
    """between 谓词：两端分属不同 select；水侧不参与。"""

    def _groups(self):
        return [
            _group(0, "H_donor", "LYS", [0, 1]),
            _group(1, "H_acceptor", "ASP", [2, 3]),
            _group(2, "water", "SOL", [4, 5]),
            _group(3, "H_donor", "SER", [6, 7]),
        ]

    def _assign(self):
        groups = self._groups()
        sets = [frozenset({0, 1, 2, 3}), frozenset({6, 7})]
        return assign_groups(groups, sets, frozenset({"SOL"}))

    def test_any_returns_none(self):
        assert build_tuple_filter("any", self._assign(),
                                  "hydrogen_bond") is None

    def test_between_cross_select(self):
        f = build_tuple_filter("between", self._assign(), "hydrogen_bond")
        g0, g3 = self._groups()[0], self._groups()[3]
        assert f((g0, g3)) is True      # LYS(sel0) × SER(sel1) → 跨侧

    def test_between_same_select(self):
        f = build_tuple_filter("between", self._assign(), "hydrogen_bond")
        g0, g1 = self._groups()[0], self._groups()[1]
        assert f((g0, g1)) is False     # 两端同属 sel0 → 非跨侧

    def test_between_water_side_not_participating(self):
        # pair 型检测器 water 已被 _filter_groups 排除，防御性验证
        f = build_tuple_filter("between", self._assign(), "hydrogen_bond")
        g2, g3 = self._groups()[2], self._groups()[3]
        assert f((g2, g3)) is False     # 水侧不参与配对

    def test_water_bridge_constraint_positions(self):
        assert CONSTRAINT_POSITIONS["water_bridge"] == (0, 2)
        f = build_tuple_filter("between", self._assign(), "water_bridge")
        g0, g2, g3 = self._groups()[0], self._groups()[2], self._groups()[3]
        # (donor=LYS sel0, water, acceptor=SER sel1) → 两端跨侧 → True
        assert f((g0, g2, g3)) is True

    def test_between_requires_interaction_type(self):
        with pytest.raises(ValueError):
            build_tuple_filter("between", self._assign(), None)


# ============================================================
# 真实数据集成测试（Amber KRAS-RBD）
# ============================================================

class TestRealAmber:
    """Amber 真实体系：RBD_pro × KRAS_pro 界面（--select ×2 + between）。"""

    @pytest.fixture(scope="class")
    def pipeline_ctx(self, tmp_path_factory):
        out = tmp_path_factory.mktemp("sel_amber") / "out"
        Pipeline("amber", "two_pass").run(
            str(AMBER_TPR), str(AMBER_XTC), str(out),
            interactions=["salt_bridge", "hydrogen_bond"],
            select=["segid seg_0_RBD_pro", "segid seg_2_KRAS_pro"],
            pairing="between")
        return out

    def test_salt_bridge_between_only_cross_select(self, pipeline_ctx):
        its = load_interactions(str(pipeline_ctx / "salt_bridge.h5"))
        assert len(its) == 1
        it = its[0]
        assert it.n_pairs > 0
        for pair in it.groups:
            mols = {g.molecule for g in pair}
            assert "seg_0_RBD_pro" in mols and "seg_2_KRAS_pro" in mols

    def test_hydrogen_bond_between_only_cross_select(self, pipeline_ctx):
        its = load_interactions(str(pipeline_ctx / "hydrogen_bond.h5"))
        it = its[0]
        for pair in it.groups:
            mols = {g.molecule for g in pair}
            assert len(mols) == 2
            assert "seg_0_RBD_pro" in mols and "seg_2_KRAS_pro" in mols


class TestRealAmberAny:
    """Amber：单 select 分子内分析（--select ×1 + any 默认）。"""

    def test_single_select_keeps_only_that_molecule(self, tmp_path):
        out = tmp_path / "out"
        Pipeline("amber", "two_pass").run(
            str(AMBER_TPR), str(AMBER_XTC), str(out),
            interactions=["salt_bridge"],
            select=["segid seg_0_RBD_pro"])
        it = load_interactions(str(out / "salt_bridge.h5"))[0]
        assert it.n_pairs > 0
        for pair in it.groups:
            assert all(g.molecule == "seg_0_RBD_pro" for g in pair)


class TestRealCharmm:
    """CHARMM 真实体系：蛋白 × BST 配体界面（complex tpr，单 segid）。"""

    def test_between_protein_ligand(self, tmp_path):
        out = tmp_path / "out"
        Pipeline("charmm", "two_pass").run(
            str(CHARMM_TPR), str(CHARMM_XTC), str(out),
            interactions=["salt_bridge", "hydrogen_bond"],
            select=["protein", "resname BST"],
            pairing="between")
        # 蛋白（protein 选择）与 BST 配体之间应存在相互作用
        it = load_interactions(str(out / "hydrogen_bond.h5"))[0]
        assert it.n_pairs > 0
        for pair in it.groups:
            resnames = {g.residue_name for g in pair}
            assert "BST" in resnames and resnames != {"BST"}


# ============================================================
# 报错路径（真实数据）
# ============================================================

class TestRealErrors:
    """选择歧义/空 → 报错。"""

    def test_overlapping_selects_raise(self, tmp_path):
        # 全 RBD_pro 与 RBD_pro 的 CA 子集必然重叠
        with pytest.raises(ValueError, match="重叠"):
            Pipeline("amber", "two_pass").run(
                str(AMBER_TPR), str(AMBER_XTC), str(tmp_path / "o"),
                interactions=["salt_bridge"],
                select=["segid seg_0_RBD_pro",
                        "segid seg_0_RBD_pro and name CA"])

    def test_between_with_one_select_raises(self, tmp_path):
        with pytest.raises(ValueError, match="至少 2 个"):
            Pipeline("amber", "two_pass").run(
                str(AMBER_TPR), str(AMBER_XTC), str(tmp_path / "o"),
                interactions=["salt_bridge"],
                select=["segid seg_0_RBD_pro"], pairing="between")

    def test_between_without_select_raises(self, tmp_path):
        with pytest.raises(ValueError, match="至少 2 个"):
            Pipeline("amber", "two_pass").run(
                str(AMBER_TPR), str(AMBER_XTC), str(tmp_path / "o"),
                interactions=["salt_bridge"], pairing="between")

    def test_syntax_error_raises(self, tmp_path):
        with pytest.raises(ValueError, match="语法错误"):
            Pipeline("amber", "two_pass").run(
                str(AMBER_TPR), str(AMBER_XTC), str(tmp_path / "o"),
                interactions=["salt_bridge"],
                select=["bogus keyword"])


# ============================================================
# 真实数据集成测试（GROMOS 53A6：同一蛋白两段残基区间之间）
# ============================================================

GROMOS_TPR = Path(__file__).parent.parent / "test_MD_case_gromos" / "gromos53a6_md.tpr"
GROMOS_XTC = Path(__file__).parent.parent / "test_MD_case_gromos" / "gromos53a6_md10ns.xtc"


class TestRealGromos:
    """GROMOS 真实体系：蛋白 1:60 × 61:120 两段残基区间（--select ×2 + between）。"""

    def test_between_residue_ranges(self, tmp_path):
        out = tmp_path / "out"
        Pipeline("gromos", "two_pass").run(
            str(GROMOS_TPR), str(GROMOS_XTC), str(out),
            interactions=["hydrophobic"],
            select=["protein and resid 1:60", "protein and resid 61:120"],
            pairing="between")
        it = load_interactions(str(out / "hydrophobic.h5"))[0]
        assert it.n_pairs > 0
        for a, b in it.groups:
            # 两端必须分属两个区间（resid 跨 60/61 分界）
            assert (a.residue_id <= 60 < b.residue_id) \
                or (b.residue_id <= 60 < a.residue_id)


# ============================================================
# CLI 层报错形态（ValueError → SystemExit 友好报错）
# ============================================================

class TestCliErrorUx:
    """CLI 把选择配置错误转成友好 SystemExit，而非裸 traceback。"""

    def _run_cli(self, monkeypatch, argv, tmp_path):
        import sys
        from DuIvyInteractions.DII import main
        monkeypatch.setattr(sys, "argv", argv)
        with pytest.raises(SystemExit) as ei:
            main()
        return ei.value

    def test_bad_select_exits_friendly(self, monkeypatch, tmp_path):
        e = self._run_cli(monkeypatch, [
            "dii", "run", "-t", str(AMBER_TPR), "-f", str(AMBER_XTC),
            "-o", str(tmp_path / "o"), "--ff", "amber",
            "--interactions", "salt_bridge", "--select", "bogus keyword"],
            tmp_path)
        assert "error:" in str(e) and "语法错误" in str(e)

    def test_overlap_exits_friendly(self, monkeypatch, tmp_path):
        e = self._run_cli(monkeypatch, [
            "dii", "run", "-t", str(AMBER_TPR), "-f", str(AMBER_XTC),
            "-o", str(tmp_path / "o"), "--ff", "amber",
            "--interactions", "salt_bridge",
            "--select", "segid seg_0_RBD_pro",
            "--select", "segid seg_0_RBD_pro and name CA"],
            tmp_path)
        assert "error:" in str(e) and "重叠" in str(e)

    def test_between_one_select_exits_friendly(self, monkeypatch, tmp_path):
        e = self._run_cli(monkeypatch, [
            "dii", "run", "-t", str(AMBER_TPR), "-f", str(AMBER_XTC),
            "-o", str(tmp_path / "o"), "--ff", "amber",
            "--interactions", "salt_bridge",
            "--select", "segid seg_0_RBD_pro", "--pairing", "between"],
            tmp_path)
        assert "error:" in str(e) and "至少 2 个" in str(e)


# ============================================================
# 不变式：CONSTRAINT_POSITIONS 覆盖全部相互作用类型
# ============================================================

class TestConstraintPositionsInvariant:
    """约束位置表必须覆盖全部 8 类相互作用（防未来类型漂移）。"""

    def test_covers_all_interactions(self):
        from DuIvyInteractions.pipeline import ALL_INTERACTIONS
        assert set(CONSTRAINT_POSITIONS) == set(ALL_INTERACTIONS)