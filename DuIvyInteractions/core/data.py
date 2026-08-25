# -*- coding: utf-8 -*-
"""核心数据结构定义：Group（基团）和 Interaction（相互作用）。"""

from dataclasses import dataclass, field
from typing import List, Tuple, Dict
import numpy as np

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


# ============================================================
# Interaction 相互作用数据结构
# ============================================================

@dataclass
class Interaction:
    """一组基团之间的相互作用记录（全帧）。

    Attributes:
        interaction_type: 相互作用类型
        groups: 参与的基团元组（2个或更多）
        existence: bool 数组，每帧是否存在
        metrics: 几何指标字典，键为指标名，值为 numpy 数组
    """

    interaction_type: str
    groups: Tuple[Group, ...]
    existence: np.ndarray
    metrics: Dict[str, np.ndarray]

    def __post_init__(self):
        """验证数据完整性。"""
        n_frames = len(self.existence)
        for name, arr in self.metrics.items():
            if len(arr) != n_frames:
                raise ValueError(
                    f"{name} length {len(arr)} != existence length {n_frames}"
                )

    @property
    def n_frames(self) -> int:
        """帧数。"""
        return len(self.existence)

    @property
    def occupancy(self) -> float:
        """存在比例。"""
        return np.sum(self.existence) / self.n_frames

    @property
    def name(self) -> str:
        """相互作用名称，由基团信息组合而成。"""
        groups_str = "-".join(
            f"{g.molecule}_{g.residue_name}{g.residue_id}"
            for g in self.groups
        )
        return f"{self.interaction_type}_{groups_str}"

    def __repr__(self) -> str:
        groups_str = ", ".join(g.group_type for g in self.groups)
        return (f"Interaction(type='{self.interaction_type}', "
                f"groups=({groups_str}), frames={self.n_frames})")


@dataclass
class HydrogenBond(Interaction):
    """氢键相互作用。"""

    @classmethod
    def create(cls, donor: Group, acceptor: Group,
               existence: np.ndarray, distance: np.ndarray,
               angle: np.ndarray) -> "HydrogenBond":
        """创建氢键实例。"""
        return cls(
            interaction_type="hydrogen_bond",
            groups=(donor, acceptor),
            existence=existence,
            metrics={"distance": distance, "angle": angle}
        )

    @property
    def donor(self) -> Group:
        """供体基团。"""
        return self.groups[0]

    @property
    def acceptor(self) -> Group:
        """受体基团。"""
        return self.groups[1]

    @property
    def name(self) -> str:
        """氢键名称：供体→受体。"""
        return f"hbond_{self.donor.residue_name}{self.donor.residue_id}→{self.acceptor.residue_name}{self.acceptor.residue_id}"


@dataclass
class PiStacking(Interaction):
    """π-π 堆积相互作用。"""

    @classmethod
    def create(cls, ring1: Group, ring2: Group,
               existence: np.ndarray, distance: np.ndarray,
               angle: np.ndarray, offset: np.ndarray) -> "PiStacking":
        """创建 π-π 堆积实例。"""
        return cls(
            interaction_type="pi_stacking",
            groups=(ring1, ring2),
            existence=existence,
            metrics={"distance": distance, "angle": angle, "offset": offset}
        )

    @property
    def name(self) -> str:
        """π-π堆积名称：环1-环2。"""
        return f"pistack_{self.groups[0].residue_name}{self.groups[0].residue_id}-{self.groups[1].residue_name}{self.groups[1].residue_id}"


@dataclass
class WaterBridge(Interaction):
    """水桥相互作用。"""

    @classmethod
    def create(cls, water: Group, donor: Group, acceptor: Group,
               existence: np.ndarray, distance_donor_water: np.ndarray,
               distance_water_acceptor: np.ndarray,
               angle: np.ndarray) -> "WaterBridge":
        """创建水桥实例。"""
        return cls(
            interaction_type="water_bridge",
            groups=(water, donor, acceptor),
            existence=existence,
            metrics={
                "distance_donor_water": distance_donor_water,
                "distance_water_acceptor": distance_water_acceptor,
                "angle": angle
            }
        )

    @property
    def name(self) -> str:
        """水桥名称：供体-水-受体。"""
        return f"waterbridge_{self.groups[1].residue_name}{self.groups[1].residue_id}-{self.groups[0].residue_name}{self.groups[0].residue_id}-{self.groups[2].residue_name}{self.groups[2].residue_id}"
