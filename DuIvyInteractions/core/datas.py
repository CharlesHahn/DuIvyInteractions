# -*- coding: utf-8 -*-
"""核心数据结构定义：Group、Interaction、InteractionSparse、SystemData。"""

from dataclasses import dataclass, field
from typing import List, Tuple, Dict, Optional
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
        atoms: 基团内的原子列表
        metadata: 附加信息字典（键必须是字符串，值只能是 JSON 支持的类型：str, int, float, bool, None, list, dict）
    """

    group_id: int
    group_type: str
    molecule: str
    residue_name: str
    residue_id: int
    atoms: List["AtomData"]
    metadata: Dict = field(default_factory=dict)

    def __post_init__(self):
        """验证数据完整性。"""
        if self.group_type not in GROUP_TYPES:
            raise ValueError(f"Invalid group_type: '{self.group_type}'")
        if not self.atoms:
            raise ValueError("Group must contain at least one atom")

    @property
    def num_atoms(self) -> int:
        """基团包含的原子数量。"""
        return len(self.atoms)

    @property
    def atom_indices(self) -> List[int]:
        """全局原子索引列表。"""
        return [a.atom_global_idx for a in self.atoms]

    @property
    def net_charge(self) -> float:
        """基团的净电荷。"""
        return sum(a.atom_charge for a in self.atoms)

    def __repr__(self) -> str:
        return (f"Group(id={self.group_id}, type='{self.group_type}', "
                f"mol='{self.molecule}', res='{self.residue_name}{self.residue_id}', "
                f"atoms={self.num_atoms})")


# ============================================================
# Interaction 相互作用数据结构（per-type，矩阵存储）
# ============================================================

@dataclass
class Interaction:
    """一种相互作用类型的全部检测结果。

    按类型组织：一个 Interaction 对象包含该类型下所有基团对在全部帧上的结果。
    groups[i] 对应 existence[i] 和 metrics 中各数组的第 i 行。

    Attributes:
        interaction_type: 相互作用类型（如 "salt_bridge", "hydrogen_bond"）
        groups: 基团对列表，每对是一个 tuple
        existence: (n_pairs, n_frames) bool 数组
        metrics: 几何指标字典，值为 (n_pairs, n_frames) 数组
    """

    interaction_type: str
    groups: List[Tuple[Group, ...]]
    existence: np.ndarray
    metrics: Dict[str, np.ndarray]

    def __post_init__(self):
        """验证数据完整性。"""
        n_pairs = len(self.groups)
        n_frames = self.existence.shape[1] if self.existence.ndim == 2 else 0
        if self.existence.shape != (n_pairs, n_frames):
            raise ValueError(
                f"existence shape {self.existence.shape} != ({n_pairs}, {n_frames})"
            )
        for name, arr in self.metrics.items():
            if arr.shape != (n_pairs, n_frames):
                raise ValueError(
                    f"{name} shape {arr.shape} != ({n_pairs}, {n_frames})"
                )

    @property
    def n_pairs(self) -> int:
        """基团对数量。"""
        return len(self.groups)

    @property
    def n_frames(self) -> int:
        """帧数。"""
        return self.existence.shape[1]

    def occupancy(self) -> np.ndarray:
        """每对基团的存在比例，shape=(n_pairs,)。"""
        return np.sum(self.existence, axis=1) / self.n_frames

    def __repr__(self) -> str:
        return (f"Interaction(type='{self.interaction_type}', "
                f"pairs={self.n_pairs}, frames={self.n_frames})")


# ============================================================
# InteractionSparse 稀疏存储的相互作用检测结果（Pass1 输出）
# ============================================================

@dataclass
class InteractionSparse:
    """稀疏存储的相互作用检测结果（Pass1 输出）。

    以 (group_id1, group_id2, ...) 为键，存储每个基团组在哪些帧存在相互作用。
    不依赖内部索引（pair_idx），直接用 group_id 标识基团组。

    Attributes:
        interaction_type: 相互作用类型（如 "salt_bridge", "hydrogen_bond"）
        data: 稀疏数据，格式：
            {
                (group_id1, group_id2): {
                    "groups": (Group1, Group2),
                    "frames": [0, 1, 5, ...],
                    "metrics": {"distance": [...], "angle": [...]}
                }
            }
    """

    interaction_type: str
    data: Dict[Tuple[int, ...], dict]

    def __post_init__(self):
        """验证数据完整性。"""
        if not self.interaction_type:
            raise ValueError("interaction_type cannot be empty")

    @property
    def n_pairs(self) -> int:
        """基团组数量。"""
        return len(self.data)

    def __repr__(self) -> str:
        return (f"InteractionSparse(type='{self.interaction_type}', "
                f"pairs={self.n_pairs})")


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
