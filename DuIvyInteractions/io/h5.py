# -*- coding: utf-8 -*-
"""HDF5 格式的 Interaction 数据存储和加载。

支持将 Interaction 列表保存为单个 HDF5 文件，并支持无损加载。
"""

import h5py
import numpy as np
import json
from typing import List, Tuple, Dict
from pathlib import Path

from ..core.datas import Interaction, Group, AtomData


# HDF5 文件格式版本
FORMAT_VERSION = "1.0"


def save_interactions(interactions: List[Interaction], path: str, compress: bool = True) -> None:
    """保存 Interaction 列表到 HDF5 文件。

    Args:
        interactions: Interaction 列表
        path: 输出文件路径
        compress: 是否启用 gzip 压缩
    """
    path = str(path)
    
    with h5py.File(path, 'w') as f:
        # 写入格式版本
        f.attrs['format_version'] = FORMAT_VERSION
        f.attrs['n_interactions'] = len(interactions)
        
        # 为每个 Interaction 创建一个 group
        for i, interaction in enumerate(interactions):
            _write_interaction(f, i, interaction, compress)


def load_interactions(path: str) -> List[Interaction]:
    """从 HDF5 文件加载 Interaction 列表。

    Args:
        path: 输入文件路径

    Returns:
        Interaction 列表
    """
    path = str(path)
    interactions = []
    
    with h5py.File(path, 'r') as f:
        # 检查格式版本
        version = f.attrs.get('format_version', 'unknown')
        if version != FORMAT_VERSION:
            raise ValueError(f"Unsupported format version: {version}. Expected: {FORMAT_VERSION}")
        
        n_interactions = f.attrs['n_interactions']
        
        for i in range(n_interactions):
            interaction = _read_interaction(f, i)
            interactions.append(interaction)
    
    return interactions


def _write_interaction(f: h5py.File, idx: int, interaction: Interaction, compress: bool) -> None:
    """写入单个 Interaction。"""
    grp = f.create_group(f"interaction_{idx}")
    
    # 写入元数据
    grp.attrs['interaction_type'] = interaction.interaction_type
    grp.attrs['n_pairs'] = interaction.n_pairs
    grp.attrs['n_frames'] = interaction.n_frames
    
    # 写入 existence
    compression = 'gzip' if compress else None
    grp.create_dataset('existence', data=interaction.existence, compression=compression)
    
    # 写入 metrics
    metrics_grp = grp.create_group('metrics')
    for name, values in interaction.metrics.items():
        metrics_grp.create_dataset(name, data=values, compression=compression)
    
    # 写入 groups
    _write_groups(grp, interaction.groups, compress)


def _read_interaction(f: h5py.File, idx: int) -> Interaction:
    """读取单个 Interaction。"""
    grp = f[f"interaction_{idx}"]
    
    # 读取元数据
    interaction_type = grp.attrs['interaction_type']
    if isinstance(interaction_type, bytes):
        interaction_type = interaction_type.decode('utf-8')
    
    # 读取 existence
    existence = grp['existence'][:]
    
    # 读取 metrics
    metrics = {}
    metrics_grp = grp['metrics']
    for name in metrics_grp:
        metrics[name] = metrics_grp[name][:]
    
    # 读取 groups
    groups = _read_groups(grp)
    
    return Interaction(
        interaction_type=interaction_type,
        groups=groups,
        existence=existence,
        metrics=metrics
    )


def _write_groups(grp: h5py.Group, groups: List[Tuple[Group, ...]], compress: bool) -> None:
    """写入 groups 数据。"""
    groups_grp = grp.create_group('groups')
    
    # 统计总 group 数和总 atom 数
    n_pairs = len(groups)
    total_groups = sum(len(g_tuple) for g_tuple in groups)
    total_atoms = sum(
        len(g.atoms)
        for g_tuple in groups
        for g in g_tuple
    )
    
    # 创建数据集
    compression = 'gzip' if compress else None
    
    # group 元数据
    pair_indices = np.zeros(total_groups, dtype=np.int64)
    group_indices_in_pair = np.zeros(total_groups, dtype=np.int64)
    group_ids = np.zeros(total_groups, dtype=np.int64)
    group_types = []
    molecules = []
    residue_names = []
    residue_ids = np.zeros(total_groups, dtype=np.int64)
    metadata_jsons = []
    
    # atom 数据
    atom_pair_indices = np.zeros(total_atoms, dtype=np.int64)
    atom_group_indices_in_pair = np.zeros(total_atoms, dtype=np.int64)
    atom_global_indices = np.zeros(total_atoms, dtype=np.int64)
    atom_idx_in_residues = np.zeros(total_atoms, dtype=np.int64)
    atom_names = []
    atom_types = []
    atom_elements = []
    atom_charges = np.zeros(total_atoms, dtype=np.float64)
    atom_masses = np.zeros(total_atoms, dtype=np.float64)
    
    # 填充数据
    group_idx = 0
    atom_idx = 0
    
    for pair_i, g_tuple in enumerate(groups):
        for g_j, g in enumerate(g_tuple):
            # group 元数据
            pair_indices[group_idx] = pair_i
            group_indices_in_pair[group_idx] = g_j
            group_ids[group_idx] = g.group_id
            group_types.append(g.group_type)
            molecules.append(g.molecule)
            residue_names.append(g.residue_name)
            residue_ids[group_idx] = g.residue_id
            metadata_jsons.append(json.dumps(g.metadata, ensure_ascii=False))
            
            # atom 数据
            for atom in g.atoms:
                atom_pair_indices[atom_idx] = pair_i
                atom_group_indices_in_pair[atom_idx] = g_j
                atom_global_indices[atom_idx] = atom.atom_global_idx
                atom_idx_in_residues[atom_idx] = atom.atom_idx_in_residue
                atom_names.append(atom.atom_name)
                atom_types.append(atom.atom_type)
                atom_elements.append(atom.atom_element)
                atom_charges[atom_idx] = atom.atom_charge
                atom_masses[atom_idx] = atom.atom_mass
                atom_idx += 1
            
            group_idx += 1
    
    # UTF-8 字符串类型
    str_dtype = h5py.string_dtype(encoding='utf-8')
    
    # 写入 group 元数据
    groups_grp.create_dataset('pair_index', data=pair_indices, compression=compression)
    groups_grp.create_dataset('group_index_in_pair', data=group_indices_in_pair, compression=compression)
    groups_grp.create_dataset('group_id', data=group_ids, compression=compression)
    groups_grp.create_dataset('group_type', data=group_types, dtype=str_dtype, compression=compression)
    groups_grp.create_dataset('molecule', data=molecules, dtype=str_dtype, compression=compression)
    groups_grp.create_dataset('residue_name', data=residue_names, dtype=str_dtype, compression=compression)
    groups_grp.create_dataset('residue_id', data=residue_ids, compression=compression)
    groups_grp.create_dataset('metadata_json', data=metadata_jsons, dtype=str_dtype, compression=compression)
    
    # 写入 atom 数据
    atoms_grp = groups_grp.create_group('atoms')
    atoms_grp.create_dataset('pair_index', data=atom_pair_indices, compression=compression)
    atoms_grp.create_dataset('group_index_in_pair', data=atom_group_indices_in_pair, compression=compression)
    atoms_grp.create_dataset('atom_global_idx', data=atom_global_indices, compression=compression)
    atoms_grp.create_dataset('atom_idx_in_residue', data=atom_idx_in_residues, compression=compression)
    atoms_grp.create_dataset('atom_name', data=atom_names, dtype=str_dtype, compression=compression)
    atoms_grp.create_dataset('atom_type', data=atom_types, dtype=str_dtype, compression=compression)
    atoms_grp.create_dataset('atom_element', data=atom_elements, dtype=str_dtype, compression=compression)
    atoms_grp.create_dataset('atom_charge', data=atom_charges, compression=compression)
    atoms_grp.create_dataset('atom_mass', data=atom_masses, compression=compression)


def _read_groups(grp: h5py.Group) -> List[Tuple[Group, ...]]:
    """读取 groups 数据。"""
    groups_grp = grp['groups']
    
    # 读取 group 元数据
    pair_indices = groups_grp['pair_index'][:]
    group_indices_in_pair = groups_grp['group_index_in_pair'][:]
    group_ids = groups_grp['group_id'][:]
    group_types = _decode_strings(groups_grp['group_type'][:])
    molecules = _decode_strings(groups_grp['molecule'][:])
    residue_names = _decode_strings(groups_grp['residue_name'][:])
    residue_ids = groups_grp['residue_id'][:]
    metadata_jsons = _decode_strings(groups_grp['metadata_json'][:])
    
    # 读取 atom 数据
    atoms_grp = groups_grp['atoms']
    atom_pair_indices = atoms_grp['pair_index'][:]
    atom_group_indices_in_pair = atoms_grp['group_index_in_pair'][:]
    atom_global_indices = atoms_grp['atom_global_idx'][:]
    atom_idx_in_residues = atoms_grp['atom_idx_in_residue'][:]
    atom_names = _decode_strings(atoms_grp['atom_name'][:])
    atom_types = _decode_strings(atoms_grp['atom_type'][:])
    atom_elements = _decode_strings(atoms_grp['atom_element'][:])
    atom_charges = atoms_grp['atom_charge'][:]
    atom_masses = atoms_grp['atom_mass'][:]
    
    # 构建 atom 映射: (pair_i, group_j) -> List[AtomData]
    atom_map: Dict[Tuple[int, int], List[AtomData]] = {}
    for i in range(len(atom_pair_indices)):
        key = (int(atom_pair_indices[i]), int(atom_group_indices_in_pair[i]))
        atom = AtomData(
            atom_global_idx=int(atom_global_indices[i]),
            atom_idx_in_residue=int(atom_idx_in_residues[i]),
            atom_name=atom_names[i],
            atom_type=atom_types[i],
            atom_element=atom_elements[i],
            atom_charge=float(atom_charges[i]),
            atom_mass=float(atom_masses[i])
        )
        atom_map.setdefault(key, []).append(atom)
    
    # 构建 group 映射: (pair_i, group_j) -> Group
    group_map: Dict[Tuple[int, int], Group] = {}
    for i in range(len(pair_indices)):
        key = (int(pair_indices[i]), int(group_indices_in_pair[i]))
        group = Group(
            group_id=int(group_ids[i]),
            group_type=group_types[i],
            molecule=molecules[i],
            residue_name=residue_names[i],
            residue_id=int(residue_ids[i]),
            atoms=atom_map.get(key, []),
            metadata=json.loads(metadata_jsons[i])
        )
        group_map[key] = group
    
    # 构建 groups 列表
    n_pairs = int(pair_indices.max()) + 1 if len(pair_indices) > 0 else 0
    groups: List[Tuple[Group, ...]] = []
    
    for pair_i in range(n_pairs):
        # 找出这个 pair 的所有 group
        pair_groups = []
        g_j = 0
        while (pair_i, g_j) in group_map:
            pair_groups.append(group_map[(pair_i, g_j)])
            g_j += 1
        groups.append(tuple(pair_groups))
    
    return groups


def _decode_strings(arr: np.ndarray) -> List[str]:
    """将 h5py 字符串数组解码为 Python 字符串列表。"""
    result = []
    for s in arr:
        if isinstance(s, bytes):
            result.append(s.decode('utf-8'))
        elif isinstance(s, str):
            result.append(s)
        else:
            result.append(str(s))
    return result
