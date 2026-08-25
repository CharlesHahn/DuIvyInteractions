# -*- coding: utf-8 -*-
"""Group 基团数据结构定义。"""

from dataclasses import dataclass, field
from typing import List, Dict

from .constants import GROUP_TYPES


@dataclass
class Group:
    """一个可参与相互作用的基团。

    Attributes:
        group_id: 唯一标识
        group_type: 基团类型，必须是 GROUP_TYPES 中的值
        molecule: 所属分子名
        residue_name: 残基名
        residue_id: 全局残基号
        atom_indices: 全局原子索引列表
        atom_types: 力场原子类型，与 atom_indices 等长
        elements: 元素符号，与 atom_indices 等长
        charges: 原子电荷，与 atom_indices 等长
        metadata: 附加信息字典
    """

    group_id: int
    group_type: str
    molecule: str
    residue_name: str
    residue_id: int
    atom_indices: List[int]
    atom_types: List[str]
    elements: List[str]
    charges: List[float]
    metadata: Dict = field(default_factory=dict)

    def __post_init__(self):
        """验证数据完整性。"""
        if self.group_type not in GROUP_TYPES:
            raise ValueError(f"Invalid group_type: '{self.group_type}'")

        n = len(self.atom_indices)
        for name, lst in [("atom_types", self.atom_types),
                          ("elements", self.elements),
                          ("charges", self.charges)]:
            if len(lst) != n:
                raise ValueError(f"{name} length {len(lst)} != atom_indices length {n}")

    @property
    def num_atoms(self) -> int:
        """基团包含的原子数量。"""
        return len(self.atom_indices)

    def __repr__(self) -> str:
        return (f"Group(id={self.group_id}, type='{self.group_type}', "
                f"mol='{self.molecule}', res='{self.residue_name}{self.residue_id}', "
                f"atoms={self.num_atoms})")


def make_bond_key(idx1: int, idx2: int) -> str:
    """生成规范化的键字符串（较小索引在前）。"""
    return f"{min(idx1, idx2)}_{max(idx1, idx2)}"


def get_bonds_for_atom(group: Group, atom_idx: int) -> Dict[str, str]:
    """获取指定原子在 Group 中的所有键。"""
    bonds = group.metadata.get("bonds", {})
    return {k: v for k, v in bonds.items() if str(atom_idx) in k.split("_")}
