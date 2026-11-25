# -*- coding: utf-8 -*-
"""3D 结构可视化纯函数层：打开拓扑/轨迹、导出拓扑 PDB、选择表达式、接触对。"""

import threading
from collections import OrderedDict
from typing import Any, Dict, List, Optional, Sequence, Tuple

import MDAnalysis as mda
import numpy as np
from MDAnalysis.lib.distances import calc_bonds

STRUCTURE_CACHE_SIZE = 2
SAMPLE_ATOMS = 5
CONTACT_FRAME_LIMIT = 5000
PDB_RESID_DIGITS = 4
PDB_SERIAL_DIGITS = 5
PDB_END = "END"
TAG_BASE = 100.0  # 打标原子的 B-factor 起点（第 k 个打标原子 = TAG_BASE + k）
TAG_BAND = 0.25  # 选择表达式 B-factor 区间半宽
NO_TAG = 0.0  # 未打标原子的 B-factor
CHAIN_POOL = "ABCDEFGHJKLMNOPQRSTUVXYZ"  # 每个 segment 一个字母（不含 W/I）
CHAIN_WATER = "W"
CHAIN_ION = "I"
PDB_ATOM_LINE = (
    "ATOM  {serial:5d} {name}{altloc}{resname:>3s} {chain}{resseq:4d}"
    "{icode}   {x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{temp_factor:6.2f}"
    "          {element:>2s}")
WATER_RESIDUES = frozenset(
    {"SOL", "HOH", "WAT", "TIP3", "TIP4", "HO4", "HO5", "H2O"})
ION_ELEMENTS = frozenset({
    "NA", "K", "MG", "CA", "ZN", "FE", "MN", "CU", "CL", "BR", "I", "F",
    "LI", "CO", "NI", "CD", "HG"})

_UNIVERSES: "OrderedDict[Tuple[str, str], Any]" = OrderedDict()
_UNIVERSE_LOCK = threading.Lock()


def load_universe(topology: str, trajectory: str) -> "mda.Universe":
    """打开 (topology, trajectory) 对应的 Universe，最近 STRUCTURE_CACHE_SIZE 个复用。"""
    key = (topology, trajectory)
    with _UNIVERSE_LOCK:
        cached = _UNIVERSES.get(key)
        if cached is not None:
            _UNIVERSES.move_to_end(key)
            return cached
    try:
        universe = mda.Universe(topology, trajectory)
    except Exception as error:
        raise ValueError(f"cannot open structure: {error}") from None
    with _UNIVERSE_LOCK:
        _UNIVERSES[key] = universe
        while len(_UNIVERSES) > STRUCTURE_CACHE_SIZE:
            _UNIVERSES.popitem(last=False)
    return universe


def topology_pdb(universe: "mda.Universe", interaction: Any = None,
                 idx: Optional[int] = None) -> str:
    """第 0 帧全原子 PDB 文本；该对原子打 B-factor 标记，不写 CONECT。"""
    universe.trajectory[0]
    tags = _tagged_factors(interaction, idx)
    chains = residue_chains(universe)
    atom_chains = [chains[int(index)] for index in universe.atoms.resindices]
    positions = _wrapped_positions(universe, range(universe.atoms.n_atoms))
    lines = [_atom_line(universe, serial, index, atom_chains[index],
                        positions[index], tags.get(index, NO_TAG))
             for index, serial in enumerate(
                 range(1, universe.atoms.n_atoms + 1))]
    return "\n".join(lines + [PDB_END]) + "\n"


def residue_chains(universe: "mda.Universe") -> Dict[int, str]:
    """残基索引 → chainID：每个 segment 一个唯一字母；水 W、单原子离子 I。"""
    residues = universe.atoms.residues
    resnames = [str(name).strip().upper() for name in residues.resnames]
    counts = np.bincount(universe.atoms.resindices,
                         minlength=len(resnames)).tolist()
    letters = _segment_letters(residues)
    elements = _first_atom_elements(universe, len(resnames))
    return {index: _chain_letter(resnames[index], int(counts[index]),
                                elements[index],
                                letters[int(residues.segindices[index])])
            for index in range(len(resnames))}


def _segment_letters(residues: Any) -> List[str]:
    """segment 索引 → chain 字母（按 segment 顺序从 CHAIN_POOL 取）。"""
    n_segments = len(residues.segments)
    if n_segments > len(CHAIN_POOL):
        raise ValueError("too many segments for chain letters")
    return list(CHAIN_POOL[:n_segments])


def interaction_expressions(interaction: Any, idx: int,
                            universe: "mda.Universe"
                            ) -> Dict[str, List[Dict[str, Any]]]:
    """该对的 Mol* 选择表达式（每组一条残基级、每原子一条 B-factor 区间）。

    表达式只用 b 区间定位 topology_pdb 打标的原子，不依赖 chain/resi 编号。
    """
    tagged = _pair_atom_expressions(interaction, idx)
    chains = residue_chains(universe)
    residues: List[Dict[str, Any]] = []
    position = 0
    for group in interaction.groups[idx]:
        rows = tagged[position:position + len(group.atoms)]
        position += len(group.atoms)
        residues.append(_group_expression(interaction, group, rows, universe,
                                          chains))
    return {"residues": residues}


def frame_contacts(interaction: Any, idx: int, universe: "mda.Universe",
                   limit: int = CONTACT_FRAME_LIMIT) -> Optional[List]:
    """逐帧最近接触 [g1, a1, g2, a2]（无接触为 None）；超过 limit 帧返回 None。"""
    n_frames = len(universe.trajectory)
    if n_frames > limit:
        return None
    return [frame_contact(interaction, idx, universe, frame)
            for frame in range(n_frames)]


def frame_contact(interaction: Any, idx: int, universe: "mda.Universe",
                  frame: int) -> Optional[List[int]]:
    """单帧最近接触 [g1, a1, g2, a2]；该帧不活跃或无接触返回 None。"""
    if not _frame_active(interaction, idx, universe, frame):
        return None
    universe.trajectory[frame]
    groups = interaction.groups[idx]
    contact = closest_contact(universe, groups[0].atom_indices,
                              groups[-1].atom_indices)
    if contact is None:
        return None
    return [0, _atom_position(groups[0], contact[0]),
            len(groups) - 1, _atom_position(groups[-1], contact[1])]


def _first_atom_elements(universe: "mda.Universe",
                         n_residues: int) -> List[str]:
    """每个残基首原子的元素（拓扑未提供元素列时全部为空串）。"""
    try:
        elements = [str(element).strip().upper()
                    for element in universe.atoms.elements]
    except (AttributeError, ValueError):
        return [""] * n_residues
    starts = np.searchsorted(universe.atoms.resindices,
                             np.arange(n_residues))
    return [elements[int(start)] for start in starts]


def _chain_letter(resname: str, atom_count: int, element: str,
                  segment_letter: str) -> str:
    """残基 chainID：水 W / 单原子金属或卤素 I / 其余取所属 segment 字母。"""
    if resname in WATER_RESIDUES:
        return CHAIN_WATER
    if atom_count == 1 and element in ION_ELEMENTS:
        return CHAIN_ION
    return segment_letter


def _pair_atoms(interaction: Any, idx: int) -> List[Any]:
    """该对全部 group 原子（组序 + 组内序摊平，B-factor 打标序号即下标）。"""
    return [atom for group in interaction.groups[idx] for atom in group.atoms]


def _tagged_factors(interaction: Any,
                    idx: Optional[int]) -> Dict[int, float]:
    """该对原子的 global idx → B-factor 标记（无 interaction 时为空）。"""
    if interaction is None or idx is None:
        return {}
    return {int(atom.atom_global_idx): TAG_BASE + position
            for position, atom in enumerate(_pair_atoms(interaction, idx))}


def _b_factor_expression(position: int) -> str:
    """第 position 个打标原子的 B-factor 区间表达式。"""
    return (f"b > {TAG_BASE + position - TAG_BAND:.2f} and "
            f"b < {TAG_BASE + position + TAG_BAND:.2f}")


def _pair_atom_expressions(interaction: Any,
                           idx: int) -> List[Dict[str, str]]:
    """该对全部原子的原子名与 B-factor 表达式（顺序 = 打标序）。"""
    return [{"name": str(atom.atom_name),
             "expression": _b_factor_expression(position)}
            for position, atom in enumerate(_pair_atoms(interaction, idx))]


def _group_expression(interaction: Any, group: Any,
                      rows: List[Dict[str, str]],
                      universe: "mda.Universe",
                      chains: Dict[int, str]) -> Dict[str, Any]:
    """单组：标签 + 残基选择器 + 是否水 + 残基级/逐原子表达式。"""
    residue = _group_residue(universe, group)
    return {"label": group_label(interaction, group, universe),
            "residue_selector": _residue_selector(residue, chains),
            "is_water": _is_water_residue(residue),
            "expression": " or ".join(f"({row['expression']})"
                                      for row in rows),
            "atoms": rows}


def _group_residue(universe: "mda.Universe", group: Any) -> Any:
    """组首原子所在残基（组内原子视为同属一个残基）。"""
    return universe.atoms[int(group.atoms[0].atom_global_idx)].residue


def _residue_selector(residue: Any, chains: Dict[int, str]) -> str:
    """残基级选择器：chain 取 residue_chains + auth resi（截断与 topology 一致）。"""
    return (f"chain {chains[int(residue.resindex)]} and "
            f"resi {_pdb_resid(int(residue.resid))}")


def _is_water_residue(residue: Any) -> bool:
    """残基是否为水（resname 在 WATER_RESIDUES 内）。"""
    return str(residue.resname).strip().upper() in WATER_RESIDUES


def _frame_active(interaction: Any, idx: int, universe: "mda.Universe",
                  frame: int) -> bool:
    """该帧在结构/轨迹范围内且该对相互作用存在。"""
    if not 0 <= frame < len(universe.trajectory):
        return False
    if frame >= interaction.n_frames:
        return False
    return bool(interaction.existence[idx, frame])


def _atom_position(group: Any, global_index: int) -> int:
    """全局原子索引在 group.atoms 中的位置（缺失为 -1）。"""
    for position, atom in enumerate(group.atoms):
        if int(atom.atom_global_idx) == int(global_index):
            return position
    return -1


def involves_water(interaction: Any, idx: int) -> bool:
    """该对相互作用是否与水相关（基团类型为 water，或类型为 water_bridge）。"""
    if interaction.interaction_type == "water_bridge":
        return True
    return any(group.group_type == "water" for group in interaction.groups[idx])


def group_label(interaction: Any, group: Any,
                universe: Optional[Any] = None) -> str:
    """基团标签：h5 残基编号 + 基团类型，PDB 显示编号不同时补 (PDB N)。"""
    if universe is None:
        return f"{group.residue_name}{group.residue_id} · {group.group_type}"
    atom = universe.atoms[int(group.atoms[0].atom_global_idx)]
    auth_id = _pdb_resid(int(atom.resid))
    prefix = f"{group.residue_name}{group.residue_id}"
    if auth_id != group.residue_id:
        prefix += f" (PDB {auth_id})"
    return f"{prefix} · {group.group_type}"


def group_atoms(interaction: Any, idx: int) -> List[List[int]]:
    """该对每个基团的全局原子索引列表（按基团顺序）。"""
    return [list(group.atom_indices) for group in interaction.groups[idx]]


def contact_label(interaction: Any, idx: int,
                  contact: Optional[Tuple[int, int]],
                  universe: Optional[Any] = None) -> str:
    """接触对标签：首末基团的基团标签(+原子名)，编号取 universe（若给出）。"""
    groups = interaction.groups[idx]
    if contact is None:
        return (f"{group_label(interaction, groups[0], universe)} ┄ "
                f"{group_label(interaction, groups[-1], universe)}")
    return (f"{_atom_label(interaction, groups[0], contact[0], universe)} ┄ "
            f"{_atom_label(interaction, groups[-1], contact[1], universe)}")


def _atom_label(interaction: Any, group: Any, atom_index: int,
                universe: Optional[Any] = None) -> str:
    """基团标签 + 基团内第 atom_index 个原子的名字（越界时只到基团）。"""
    label = group_label(interaction, group, universe)
    if not 0 <= atom_index < len(group.atoms):
        return label
    return f"{label} {group.atoms[atom_index].atom_name}"


def _pdb_resid(resid: int) -> int:
    """PDB resSeq 四列显示值（超长取低 4 位，与 MDAnalysis 写器一致）。"""
    return _pdb_truncate(resid, PDB_RESID_DIGITS)


def _pdb_serial(serial: int) -> int:
    """PDB serial 五列显示值（超长取低 5 位，与 MDAnalysis 写器一致）。"""
    return _pdb_truncate(serial, PDB_SERIAL_DIGITS)


def _pdb_truncate(value: int, digits: int) -> int:
    """PDB 定宽整数列（取低 digits 位）。"""
    return int(str(int(value))[-digits:])


def _wrapped_positions(universe: Any,
                       indices: Sequence[int]) -> "np.ndarray":
    """当前帧所选原子的坐标（有盒时按最小镜像回绕，与主 PDB 同规则）。"""
    if not indices:
        return np.zeros((0, 3))
    atoms = universe.atoms[list(indices)]
    if universe.dimensions is None:
        return atoms.positions.copy()
    return atoms.wrap(inplace=False)


def _atom_line(universe: Any, serial: int, global_index: int, chain: str,
               xyz: Any, temp_factor: float = NO_TAG) -> str:
    """单条 ATOM 记录（原子名/残基名/resSeq 来自 Universe 的 auth 编号）。"""
    atom = universe.atoms[global_index]
    return PDB_ATOM_LINE.format(
        serial=_pdb_serial(serial), name=_pdb_atom_name(str(atom.name)),
        altloc=" ", resname=_pdb_resname(atom.resname), chain=chain,
        resseq=_pdb_resid(int(atom.resid)), icode=" ",
        x=float(xyz[0]), y=float(xyz[1]), z=float(xyz[2]),
        occupancy=1.0, temp_factor=temp_factor, element=_atom_element(atom))


def _pdb_atom_name(name: str) -> str:
    """PDB 原子名四列（短名右移一列，长名截断）。"""
    if len(name) >= 4:
        return name[:4]
    return (" " + name).ljust(4)


def _pdb_resname(resname: Any) -> str:
    """PDB 残基名三列（超长截断）。"""
    return str(resname)[:3]


def _atom_element(atom: Any) -> str:
    """PDB 元素列（Universe 未提供时留空，查看器按原子名推断）。"""
    return str(getattr(atom, "element", "") or "")[:2]


def alignment_errors(universe: "mda.Universe", interaction: Any) -> List[str]:
    """h5 原子索引/名称/力场类型与结构不一致的错误列表（空 = 一致）。"""
    atoms = _interaction_atoms(interaction)
    errors: List[str] = []
    if atoms:
        required = max(atom.atom_global_idx for atom in atoms) + 1
        if universe.atoms.n_atoms < required:
            errors.append(
                f"structure has {universe.atoms.n_atoms} atoms; h5 needs {required}")
    for atom in atoms[:SAMPLE_ATOMS]:
        message = _atom_error(universe, atom)
        if message:
            errors.append(message)
    return errors


def _interaction_atoms(interaction: Any) -> List[Any]:
    """interaction.groups 中按出现顺序摊平的全部原子。"""
    return [atom for pair in interaction.groups
            for group in pair for atom in group.atoms]


def _atom_error(universe: "mda.Universe", atom: Any) -> str:
    """单个抽样原子的名称/类型差异描述（一致时为空串）。"""
    index = atom.atom_global_idx
    if not 0 <= index < universe.atoms.n_atoms:
        return ""
    actual = universe.atoms[index]
    if actual.name != atom.atom_name:
        return (f"atom {index}: name {actual.name!r} != h5 "
                f"{atom.atom_name!r}")
    actual_type = _atom_type(actual)
    if atom.atom_type and actual_type and actual_type != atom.atom_type:
        return (f"atom {index}: type {actual_type!r} != h5 "
                f"{atom.atom_type!r}")
    return ""


def _atom_type(atom: Any) -> str:
    """原子的力场类型（拓扑未提供时为空串）。"""
    try:
        return str(atom.type)
    except (AttributeError, ValueError):
        return ""


def closest_contact(universe: "mda.Universe",
                    group1_atoms: Sequence[int],
                    group2_atoms: Sequence[int]) -> Optional[Tuple[int, int]]:
    """当前帧下两组原子的最近接触原子对（global idx），任一组为空返回 None。"""
    first, second = list(group1_atoms), list(group2_atoms)
    if not first or not second:
        return None
    coords1 = universe.atoms[first].positions
    coords2 = universe.atoms[second].positions
    best: Optional[Tuple[float, int, int]] = None
    for row, coord in enumerate(coords1):
        targets = np.broadcast_to(coord, coords2.shape)
        distances = calc_bonds(targets, coords2, box=universe.dimensions)
        column = int(distances.argmin())
        if best is None or distances[column] < best[0]:
            best = (float(distances[column]), row, column)
    return int(first[best[1]]), int(second[best[2]])
