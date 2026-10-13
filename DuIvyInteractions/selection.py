# -*- coding: utf-8 -*-
"""用户基团选择与过滤：--select / --pairing 参数的核心逻辑。

设计文档：doc/group_selection_design.md

规则：
1. 各 --select 串在 Universe 上求值为独立的原子集合（不合并）；
2. 基团保留 ⇔ 某个 select 集合完整覆盖该基团的全部原子；
3. water 基团恒属"水侧"，跳过 select 归属判定；
4. --pairing between 仅保留两端基团分属不同 select 的候选对（水侧不参与）；
5. 歧义/空状态抛 ValueError。
"""

from typing import TYPE_CHECKING, Callable, Dict, List, Optional, Union

from .core.datas import Group

if TYPE_CHECKING:
    import MDAnalysis as mda

# 水侧的恒定归属标记（隐式 select，用户不可关闭）
WATER_SIDE = "water"

# between 判定所需的不同侧最小数（跨侧 = 至少 2 个不同 select）
CROSS_MIN_SIDES = 2

# 基团归属值：select 索引（int）或 WATER_SIDE（str）
Side = Union[int, str]
Assignment = Dict[int, Side]

# 各检测器 tuple 中参与配对判定的成员位置（固定，不开放给用户）
# pair 型检测器均为两端 (0,1)；water_bridge 为 (donor, water, acceptor) → 判 (0,2)
CONSTRAINT_POSITIONS = {
    "hydrogen_bond": (0, 1),
    "pi_stacking": (0, 1),
    "salt_bridge": (0, 1),
    "hydrophobic": (0, 1),
    "halogen_bond": (0, 1),
    "metal_coordination": (0, 1),
    "water_bridge": (0, 2),
    "pi_cation": (0, 1),
}


def resolve_selection(universe: "mda.Universe", sel_str: str) -> frozenset:
    """将一条 --select 选择串求值为原子索引集合。

    Args:
        universe: 已加载轨迹、已 seek 到第一帧的 MDAnalysis Universe
        sel_str: MDAnalysis 选择语法字符串

    Returns:
        原子索引 frozenset。

    Raises:
        ValueError: 选择语法错误。
    """
    try:
        ag = universe.select_atoms(sel_str)
    except Exception as e:  # MDAnalysis SelectionError 等
        raise ValueError(f"select 语法错误: '{sel_str}' → {e}")
    return frozenset(ag.indices)


def validate_selections(sets: List[frozenset]) -> None:
    """校验 select 集合：空 / 两两重叠则报错。

    Raises:
        ValueError: 空 select 或 select 间原子重叠。
    """
    for i, s in enumerate(sets, 1):
        if not s:
            raise ValueError(f"select #{i} 未选中任何原子")
    for i in range(len(sets)):
        for j in range(i + 1, len(sets)):
            overlap = sets[i] & sets[j]
            if overlap:
                raise ValueError(
                    f"select #{i + 1} 与 select #{j + 1} 原子重叠 "
                    f"({len(overlap)} 个原子)，请用 and not 显式剔除后重试")


def _match_select(group: Group, sets: List[frozenset]) -> Optional[int]:
    """基团全部原子所在的 select 索引；无 select 完整覆盖则返回 None。"""
    indices = group.atom_indices
    for i, s in enumerate(sets):
        if all(idx in s for idx in indices):
            return i
    return None


def assign_groups(groups: List[Group], sets: List[frozenset],
                  water_residues: frozenset) -> Assignment:
    """确定每个基团的归属（select 索引或 WATER_SIDE）；无归属的基团不返回。"""
    assignment: Assignment = {}
    for g in groups:
        if g.residue_name in water_residues:
            assignment[g.group_id] = WATER_SIDE
            continue
        matched = _match_select(g, sets)
        if matched is not None:
            assignment[g.group_id] = matched
    return assignment


def filter_groups(groups: List[Group], assignment: Assignment) -> List[Group]:
    """按归属过滤基团池：仅保留有归属（含 WATER_SIDE）的基团。"""
    return [g for g in groups if g.group_id in assignment]


def _make_between(assignment: Assignment,
                  positions: tuple) -> Callable[[tuple], bool]:
    """组装 between 谓词：约束位置的基团须分属不同 select。"""
    def _between(t: tuple) -> bool:
        sides = set()
        for pos in positions:
            side = assignment.get(t[pos].group_id)
            if side is None or side == WATER_SIDE:
                return False
            sides.add(side)
        return len(sides) >= CROSS_MIN_SIDES
    return _between


def build_tuple_filter(pairing: str, assignment: Assignment,
                       interaction_type: Optional[str] = None
                       ) -> Optional[Callable[[tuple], bool]]:
    """组装 tuple 层过滤谓词（any → None，between → 谓词）。

    Args:
        pairing: 'any' 或 'between'
        assignment: assign_groups 的结果
        interaction_type: 相互作用类型名（between 时用于查 CONSTRAINT_POSITIONS）

    Returns:
        可调用谓词或 None。
    """
    if pairing != "between":
        return None
    if interaction_type is None:
        raise ValueError("build_tuple_filter: between 需要 interaction_type")
    return _make_between(assignment, CONSTRAINT_POSITIONS[interaction_type])