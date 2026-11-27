# -*- coding: utf-8 -*-
"""3D 结构可视化纯函数层：打开拓扑/轨迹、导出拓扑 mmCIF、选择表达式。"""

import re
import threading
from collections import OrderedDict
from typing import Any, Dict, List, Optional, Sequence, Tuple

import MDAnalysis as mda
import numpy as np

STRUCTURE_CACHE_SIZE = 2
SAMPLE_ATOMS = 5
CIF_BLOCK = "data_dii"
CIF_LOOP = "loop_"
CIF_END = "#"
CIF_INS_CODE = "?"
CIF_ALT_ID = "."
CIF_MODEL_NUM = "1"
CIF_OCCUPANCY = "1.00"
CIF_ATOM = "ATOM"
CIF_POLYMER = "polymer"
CIF_NON_POLYMER = "non-polymer"
CIF_WATER = "water"
CIF_SAFE_TOKEN = re.compile(r"[A-Za-z0-9_.+-]+")
CIF_ATOM_SITE_COLUMNS = (
    "_atom_site.group_PDB", "_atom_site.id", "_atom_site.type_symbol",
    "_atom_site.label_atom_id", "_atom_site.label_alt_id",
    "_atom_site.label_comp_id",
    "_atom_site.label_asym_id", "_atom_site.label_entity_id",
    "_atom_site.label_seq_id", "_atom_site.auth_seq_id",
    "_atom_site.auth_comp_id", "_atom_site.auth_asym_id",
    "_atom_site.auth_atom_id", "_atom_site.pdbx_PDB_ins_code",
    "_atom_site.pdbx_PDB_model_num",
    "_atom_site.Cartn_x", "_atom_site.Cartn_y", "_atom_site.Cartn_z",
    "_atom_site.occupancy", "_atom_site.B_iso_or_equiv")
TAG_BASE = 100.0  # 打标原子的 B-factor 起点（第 k 个打标原子 = TAG_BASE + k）
TAG_BAND = 0.25  # 选择表达式 B-factor 区间半宽
NO_TAG = 0.0  # 未打标原子的 B-factor
CHAIN_POOL = "ABCDEFGHJKLMNOPQRSTUVXYZ"  # 每个 segment 一个字母（不含 W/I）
CHAIN_WATER = "W"
CHAIN_ION = "I"
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


def topology_cif(interaction: Any, idx: Optional[int],
                 universe: "mda.Universe") -> str:
    """第 0 帧全原子拓扑 mmCIF；该对原子保留 B-factor 标记（旧表达式兜底）。

    auth_seq_id 原样使用残基编号（不再受 PDB 4 列截断），原子 id = 全局索引 + 1。
    """
    universe.trajectory[0]
    tags = _tagged_factors(interaction, idx)
    chains = residue_chains(universe)
    classes = _residue_classes(universe)
    entities = _entity_ids(chains, classes)
    resindices = universe.atoms.resindices
    atom_chains = [chains[int(index)] for index in resindices]
    atom_entities = [entities[(atom_chains[position], classes[int(index)])]
                     for position, index in enumerate(resindices)]
    positions = _wrapped_positions(universe, range(universe.atoms.n_atoms))
    rows = _cif_atom_rows(universe, atom_chains, atom_entities, positions,
                          tags)
    header = [CIF_BLOCK, CIF_END, CIF_LOOP]
    header += list(CIF_ATOM_SITE_COLUMNS)
    return "\n".join(header + rows + [CIF_END]) + "\n"


def _cif_atom_rows(universe: "mda.Universe", atom_chains: List[str],
                   atom_entities: List[int], positions: "np.ndarray",
                   tags: Dict[int, float]) -> List[str]:
    """全部原子的 _atom_site 行（serial = 全局索引 + 1，不截断）。"""
    rows = []
    for index, serial in enumerate(range(1, universe.atoms.n_atoms + 1)):
        rows.append(_cif_atom_row(
            universe, serial, index, atom_chains[index],
            atom_entities[index], positions[index],
            tags.get(index, NO_TAG)))
    return rows


def _cif_atom_row(universe: "mda.Universe", serial: int, global_index: int,
                  chain: str, entity: int, xyz: Any,
                  temp_factor: float = NO_TAG) -> str:
    """单行 _atom_site 记录（auth_seq_id 用完整残基编号）。"""
    atom = universe.atoms[global_index]
    name = _cif_token(str(atom.name))
    resname = _cif_token(str(atom.resname))
    return " ".join([
        CIF_ATOM, str(serial), _cif_element(atom), name, CIF_ALT_ID,
        resname, chain,
        str(entity), str(int(atom.resindex) + 1), str(_residue_number(atom)),
        resname, chain, name, CIF_INS_CODE, CIF_MODEL_NUM,
        f"{float(xyz[0]):.3f}", f"{float(xyz[1]):.3f}",
        f"{float(xyz[2]):.3f}", CIF_OCCUPANCY,
        f"{float(temp_factor):.2f}"])


def _residue_classes(universe: "mda.Universe") -> Dict[int, str]:
    """残基索引 → 实体类别：water / polymer（蛋白核酸）/ non-polymer。"""
    residues = universe.atoms.residues
    classes = {index: (CIF_WATER if str(name).strip().upper()
                       in WATER_RESIDUES else CIF_NON_POLYMER)
               for index, name in enumerate(residues.resnames)}
    for index in _polymer_resindices(universe):
        classes[index] = CIF_POLYMER
    return classes


def _polymer_resindices(universe: "mda.Universe") -> List[int]:
    """MDAnalysis 判定为蛋白/核酸的残基索引（用于 mmCIF 实体分类）。"""
    atoms = universe.select_atoms("protein or nucleic")
    return sorted({int(index) for index in atoms.resindices})


def _entity_ids(chains: Dict[int, str],
                classes: Dict[int, str]) -> Dict[Tuple[str, str], int]:
    """(chain, 类别) → 实体编号（同链同类共享一个实体，编号从 1 递增）。"""
    entities: Dict[Tuple[str, str], int] = {}
    for index in range(len(chains)):
        key = (chains[index], classes[index])
        if key not in entities:
            entities[key] = len(entities) + 1
    return entities


def _residue_number(item: Any) -> int:
    """残基 auth 编号原值（mmCIF 不截断）：水 resid 15689 保持 15689。"""
    return int(item.resid)


def _cif_element(atom: Any) -> str:
    """元素符号规范大小写（C / N / O / Cl / Na）；缺失时取原子名首字母。"""
    raw = str(getattr(atom, "element", "") or "").strip()
    if not raw:
        raw = str(atom.name).strip()[:1]
    return raw[:1].upper() + raw[1:].lower() if raw else CIF_INS_CODE


def _cif_token(text: str) -> str:
    """CIF 值：普通 token 原样输出，含空格/特殊字符时加引号。"""
    value = str(text)
    if value.startswith("_"):
        return "'" + value + "'"
    if CIF_SAFE_TOKEN.fullmatch(value):
        return value
    if "'" not in value:
        return "'" + value + "'"
    if '"' not in value:
        return '"' + value + '"'
    return "'" + value.replace("'", "") + "'"


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
                            universe: Optional["mda.Universe"] = None
                            ) -> Dict[str, List[Dict[str, Any]]]:
    """该对的 Mol* 选择表达式与原子全局索引（每组残基级 + 每原子 B-factor）。

    表达式只用 b 区间定位 topology_cif 打标的原子，不依赖 chain/resi 编号。
    universe 为 None（未选结构）时仍给出 indices，residue_indices 为 None。
    """
    tagged = _atom_rows(interaction, idx)
    chains = residue_chains(universe) if universe is not None else {}
    residues: List[Dict[str, Any]] = []
    position = 0
    for group in interaction.groups[idx]:
        rows = tagged[position:position + len(group.atoms)]
        position += len(group.atoms)
        residues.append(_group_expression(interaction, group, rows, universe,
                                          chains))
    return {"residues": residues}


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


def _atom_rows(interaction: Any, idx: int) -> List[Dict[str, Any]]:
    """该对全部原子的名称、全局索引与 B-factor 表达式（顺序 = 打标序）。"""
    return [{"name": str(atom.atom_name),
             "global_idx": int(atom.atom_global_idx),
             "expression": _b_factor_expression(position)}
            for position, atom in enumerate(_pair_atoms(interaction, idx))]


def _group_expression(interaction: Any, group: Any,
                      rows: List[Dict[str, Any]],
                      universe: Optional["mda.Universe"],
                      chains: Dict[int, str]) -> Dict[str, Any]:
    """单组：标签 + 残基选择器 + 是否水 + 全局索引 + 组/残基级表达式。"""
    residue = _group_residue(universe, group) if universe is not None else None
    selector = (_residue_selector(residue, chains)
                if residue is not None else None)
    expression = " or ".join(f"({row['expression']})" for row in rows)
    return {"label": group_label(interaction, group, universe),
            "residue_selector": selector,
            "is_water": (_is_water_residue(residue) if residue is not None
                         else group.group_type == "water"),
            "expression": expression,
            "group_expression": _group_expression_by_name(
                selector, rows, universe, group, expression),
            "indices": _group_indices(group),
            "residue_indices": _residue_indices(universe, group),
            "atoms": rows}


def _group_expression_by_name(residue_selector: Optional[str],
                              rows: List[Dict[str, Any]],
                              universe: Optional["mda.Universe"],
                              group: Any,
                              fallback: str) -> Optional[str]:
    """链/残基/原子名的精确组表达式；无结构或同名退化时回退 b 区间。"""
    if universe is None or residue_selector is None:
        return None
    if _covers_whole_residue(universe, group):
        return residue_selector
    names = _residue_atom_names(universe, group)
    if len(set(names)) != len(names):
        return fallback
    clause = " or ".join("name " + str(row["name"]) for row in rows)
    return residue_selector + " and (" + clause + ")"


def _covers_whole_residue(universe: "mda.Universe", group: Any) -> bool:
    """组原子集合是否等于该残基全部原子（按全局索引集合比较）。"""
    return (set(_group_indices(group))
            == set(_residue_indices(universe, group) or []))


def _residue_atom_names(universe: "mda.Universe", group: Any) -> List[str]:
    """组首原子所在残基的全部原子名（残基内顺序）。"""
    return [str(atom.name) for atom in _group_residue(universe, group).atoms]


def _group_indices(group: Any) -> List[int]:
    """组内全部原子的全局索引（顺序 = h5 group.atoms 顺序）。"""
    return [int(atom.atom_global_idx) for atom in group.atoms]


def _residue_indices(universe: Optional["mda.Universe"],
                     group: Any) -> Optional[List[int]]:
    """该组首原子所在残基的全部原子全局索引（升序）；无 universe → None。"""
    if universe is None:
        return None
    residue = _group_residue(universe, group)
    return sorted(int(index) for index in residue.atoms.indices)


def _group_residue(universe: "mda.Universe", group: Any) -> Any:
    """组首原子所在残基（组内原子视为同属一个残基）。"""
    return universe.atoms[int(group.atoms[0].atom_global_idx)].residue


def _residue_selector(residue: Any, chains: Dict[int, str]) -> str:
    """残基级选择器：chain 取 residue_chains + 完整 auth resi（mmCIF 不截断）。"""
    return (f"chain {chains[int(residue.resindex)]} and "
            f"resi {_residue_number(residue)}")


def _is_water_residue(residue: Any) -> bool:
    """残基是否为水（resname 在 WATER_RESIDUES 内）。"""
    return str(residue.resname).strip().upper() in WATER_RESIDUES


def involves_water(interaction: Any, idx: int) -> bool:
    """该对相互作用是否与水相关（基团类型为 water，或类型为 water_bridge）。"""
    if interaction.interaction_type == "water_bridge":
        return True
    return any(group.group_type == "water" for group in interaction.groups[idx])


def group_label(interaction: Any, group: Any,
                universe: Optional[Any] = None) -> str:
    """基团标签：h5 残基编号 + 基团类型，3D 编号不同时补 (PDB N)。"""
    if universe is None:
        return f"{group.residue_name}{group.residue_id} · {group.group_type}"
    atom = universe.atoms[int(group.atoms[0].atom_global_idx)]
    auth_id = _residue_number(atom)
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


def _wrapped_positions(universe: Any,
                       indices: Sequence[int]) -> "np.ndarray":
    """当前帧所选原子的坐标（有盒时按最小镜像回绕，与拓扑写出同规则）。"""
    if not indices:
        return np.zeros((0, 3))
    atoms = universe.atoms[list(indices)]
    if universe.dimensions is None:
        return atoms.positions.copy()
    return atoms.wrap(inplace=False)


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
