# -*- coding: utf-8 -*-
"""核心数据结构定义：Group、Interaction、SystemData。"""

from dataclasses import dataclass, field
from typing import List, Tuple, Dict
import numpy as np

from .constants import GROUP_TYPES, BOND_TYPES


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


# ============================================================
# SystemData 体系分子数据结构
# ============================================================

@dataclass
class AtomData:
    """单个原子的信息。"""
    atom_global_idx: int        # 全局原子索引（整个体系唯一）
    atom_idx_in_residue: int    # 残基内索引（用于残基内键连接）
    atom_name: str              # 原子名（如 "CG"）
    atom_type: str              # 力场类型（如 "ca"）
    atom_element: str           # 元素符号（如 "C"）
    atom_charge: float          # 电荷
    atom_mass: float            # 原子质量（原子质量单位，-1.0 表示未设置）


@dataclass
class BondData:
    """残基内的一个键。"""
    atom1_idx_in_residue: int   # 残基内索引
    atom2_idx_in_residue: int   # 残基内索引
    bond_type: str              # 键类型，见 BOND_TYPES

    def __post_init__(self):
        """验证键类型。"""
        if self.bond_type not in BOND_TYPES:
            raise ValueError(f"Invalid bond_type: '{self.bond_type}'")


@dataclass
class ResidueData:
    """一个残基的数据。"""
    residue_name: str                   # 残基名（如 "TYR"）
    residue_global_idx: int             # 全局残基索引（整个体系唯一）
    residue_idx_in_molecule: int        # 分子内残基编号（PDB 编号）
    molecule_name: str                  # 所属分子名（如 "RBD_pro"）
    atoms: List[AtomData]               # 残基内的原子
    bonds: List[BondData]               # 残基内的键


@dataclass
class InterResidueBond:
    """残基间的共价键（如肽键、二硫键）。"""
    residue1_global_idx: int            # 残基 1 的全局索引
    atom_idx_in_residue1: int           # 原子在残基 1 内的索引
    residue2_global_idx: int            # 残基 2 的全局索引
    atom_idx_in_residue2: int           # 原子在残基 2 内的索引
    bond_type: str                      # 键类型

    def __post_init__(self):
        """验证键类型。"""
        if self.bond_type not in BOND_TYPES:
            raise ValueError(f"Invalid bond_type: '{self.bond_type}'")


@dataclass
class SystemData:
    """体系数据，连接 Reader 和 Identifier。"""
    system_name: str
    residues: List[ResidueData]
    inter_residue_bonds: List[InterResidueBond]

    def __post_init__(self):
        """验证数据完整性。"""
        # 验证 residue_global_idx 唯一性
        res_indices = [r.residue_global_idx for r in self.residues]
        if len(res_indices) != len(set(res_indices)):
            raise ValueError("Duplicate residue_global_idx found")

        # 验证 atom_global_idx 唯一性
        atom_indices = []
        for r in self.residues:
            for a in r.atoms:
                atom_indices.append(a.atom_global_idx)
        if len(atom_indices) != len(set(atom_indices)):
            raise ValueError("Duplicate atom_global_idx found")

    @property
    def n_residues(self) -> int:
        """残基数量。"""
        return len(self.residues)

    def get_residue_by_global_idx(self, idx: int) -> ResidueData:
        """根据全局索引获取残基。"""
        for r in self.residues:
            if r.residue_global_idx == idx:
                return r
        raise KeyError(f"Residue with global_idx={idx} not found")

    def get_atoms_by_molecule(self, molecule_name: str) -> List[AtomData]:
        """获取指定分子的所有原子。"""
        atoms = []
        for r in self.residues:
            if r.molecule_name == molecule_name:
                atoms.extend(r.atoms)
        return atoms

    def __repr__(self) -> str:
        return (f"SystemData(name='{self.system_name}', "
                f"residues={self.n_residues}, "
                f"inter_residue_bonds={len(self.inter_residue_bonds)})")
