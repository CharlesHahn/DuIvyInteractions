# -*- coding: utf-8 -*-
"""3D 结构可视化纯函数层：打开拓扑/轨迹、导出单帧 PDB、一致性校验、接触对。"""

import io
import threading
import warnings
from collections import OrderedDict
from dataclasses import dataclass
from typing import Any, FrozenSet, List, Optional, Sequence, Tuple

import MDAnalysis as mda
import numpy as np
from MDAnalysis.lib.distances import calc_bonds

STRUCTURE_CACHE_SIZE = 2
FRAME_CACHE_SIZE = 8
SAMPLE_ATOMS = 5
PDB_FORMAT = "PDB"
HIGHLIGHT_CHAINS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
PDB_RESID_DIGITS = 4
PDB_ATOM_LINE = (
    "ATOM  {serial:5d} {name}{altloc}{resname:>3s} {chain}{resseq:4d}"
    "{icode}   {x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{temp_factor:6.2f}"
    "          {element:>2s}")
WATER_RESIDUES = frozenset(
    {"SOL", "HOH", "WAT", "TIP3", "TIP4", "HO4", "HO5", "H2O"})

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


def frame_pdb(universe: "mda.Universe", frame: int) -> str:
    """第 frame 帧的 PDB 文本；有盒时回绕显示，原子顺序与 Universe 一致（索引映射前提）。"""
    _check_frame(universe, frame)
    universe.trajectory[frame]
    return _write_pdb(_display_copy(universe))


def frame_view(universe: "mda.Universe", frame: int, hide_water: bool = False,
               keep_atoms: FrozenSet[int] = frozenset()) -> Tuple[str, List[int]]:
    """第 frame 帧的 PDB 文本与逐行全局原子索引（hide_water 剔除水，keep_atoms 恒保留）。"""
    _check_frame(universe, frame)
    universe.trajectory[frame]
    shown = _shown_indices(universe, hide_water, keep_atoms)
    display = _display_copy(universe)
    return _write_pdb(display.atoms[shown] if hide_water else display), shown


def _check_frame(universe: "mda.Universe", frame: int) -> None:
    """帧号校验（越界报 ValueError）。"""
    if not 0 <= frame < len(universe.trajectory):
        raise ValueError(f"frame out of range: {frame}")


def _display_copy(universe: "mda.Universe") -> "mda.Universe":
    """回绕坐标后的 Universe 副本（原 Universe 不变）。"""
    display = universe.copy()
    if universe.dimensions is not None:
        display.atoms.positions = universe.atoms.wrap(inplace=False)
    return display


def _shown_indices(universe: "mda.Universe", hide_water: bool,
                   keep_atoms: FrozenSet[int]) -> List[int]:
    """保留原子的全局索引（升序；hide_water 时剔除水残基，keep_atoms 除外）。"""
    indices = np.arange(universe.atoms.n_atoms)
    if not hide_water:
        return [int(index) for index in indices]
    water = np.isin(universe.atoms.resnames, list(WATER_RESIDUES))
    keep = np.zeros(len(indices), dtype=bool)
    keep[[index for index in keep_atoms if 0 <= index < len(indices)]] = True
    return [int(index) for index in indices[~water | keep]]


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


def _pdb_resid(resid: int) -> int:
    """PDB resSeq 四列显示值（超长取低 4 位，与 MDAnalysis 写器一致）。"""
    return int(str(int(resid))[-PDB_RESID_DIGITS:])


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


@dataclass(frozen=True)
class Highlight:
    """单个基团的高亮小结构：组标签 + 仅含该组原子的 PDB 文本。"""

    label: str
    pdb: str


def highlights(interaction: Any, idx: int, universe: Any) -> List[Highlight]:
    """该对每个基团一个独立小结构（第 position 组用第 position 个 chain ID）。"""
    return [Highlight(label=group_label(interaction, group, universe),
                      pdb=_group_pdb(universe, group, _chain_id(position)))
            for position, group in enumerate(interaction.groups[idx])]


def contact_pdb(universe: Any,
                contact: Optional[Tuple[int, int]]) -> Optional[str]:
    """最近接触两个原子的 PDB 文本（两原子分别占 chain A/B）。"""
    if not contact:
        return None
    return _atoms_pdb(universe, [(int(contact[0]), _chain_id(0)),
                                 (int(contact[1]), _chain_id(1))])


def _chain_id(position: int) -> str:
    """第 position 个基团的 PDB chain ID（A/B/C…，超过 26 组循环）。"""
    return HIGHLIGHT_CHAINS[position % len(HIGHLIGHT_CHAINS)]


def _group_pdb(universe: Any, group: Any, chain: str) -> str:
    """该基团全部原子的 PDB 文本（组内顺序 = group.atoms 顺序）。"""
    return _atoms_pdb(universe, [(int(atom.atom_global_idx), chain)
                                 for atom in group.atoms])


def _atoms_pdb(universe: Any, rows: Sequence[Tuple[int, str]]) -> str:
    """(全局原子索引, chain) 序列 → PDB 文本（坐标取当前帧并回绕）。"""
    positions = _wrapped_positions(universe, [index for index, _ in rows])
    lines = [_atom_line(universe, serial, index, chain,
                        positions[serial - 1])
             for serial, (index, chain) in enumerate(rows, start=1)]
    return "\n".join(lines + ["END"]) + "\n"


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
               xyz: Any) -> str:
    """单条 ATOM 记录（原子名/残基名/resSeq 来自 Universe 的 auth 编号）。"""
    atom = universe.atoms[global_index]
    return PDB_ATOM_LINE.format(
        serial=serial, name=_pdb_atom_name(str(atom.name)), altloc=" ",
        resname=_pdb_resname(atom.resname), chain=chain,
        resseq=_pdb_resid(int(atom.resid)), icode=" ",
        x=float(xyz[0]), y=float(xyz[1]), z=float(xyz[2]),
        occupancy=1.0, temp_factor=0.0, element=_atom_element(atom))


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


def _write_pdb(universe: "mda.Universe") -> str:
    """把 Universe 当前帧写成 PDB 文本（不写 CONECT，序列号按 PDB 上限回绕）。"""
    buffer = io.StringIO()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        writer = mda.Writer(buffer, n_atoms=universe.atoms.n_atoms,
                            format=PDB_FORMAT, reindex=True, bonds=None)
        writer.write(universe)
        text = buffer.getvalue()
        writer.close()
    return text


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
