# -*- coding: utf-8 -*-
"""基团识别结果自动核验工具（TODO #3）。

以"残基名 + 原子名 + 键合图"推导的化学事实基准（独立于识别器代码），
逐残基 diff 识别器输出，报告缺失/多余的基团，供人工复核。

用法：
    python Tests/tools/verify_identification.py \
        --tpr Tests/test_MD_case_amber/md.tpr --ff amber \
        [--report /tmp/report.md]

设计：doc/group_identification_review.md
"""

import argparse
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Set, Tuple

# 脚本位于 Tests/tools/，仓库根为父目录的父目录——确保可 import 包
_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from DuIvyInteractions.system_readers import GmxTprReader
from DuIvyInteractions.group_identifiers import IDENTIFIER_CLASSES

# 标准氨基酸残基名（核验范围限定在此内；配体/水/离子跳过）
STANDARD_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "HID", "HIE", "HIP",
}

# 水/离子残基名（跳过）
NON_STANDARD = {"SOL", "HOH", "WAT", "TIP3", "CL", "NA", "MG",
                "NA+", "CL-", "K", "CA", "ZN", "MG2"}

# 受体判定的吡啶/亚胺型 N 类型（无 H 邻居、有孤对电子）
ACCEPTOR_N_TYPES = {"nb", "NB", "n", "n2", "na", "NA"}

# 正电/负电残基——期望"代表原子"（基团中心原子，与 group_atom_convention 对齐：
# charged 基团 atoms[0] 是中心原子，如 LYS=NZ、ARG=NE、ASP=CG、GLU=CD）
CHARGED_POS = {"LYS": {"NZ"}, "ARG": {"NE"}}
CHARGED_NEG = {"ASP": {"CG"}, "GLU": {"CD"}}

# 芳香环期望数（按残基名；HIS 三质子化态 HID/HIE/HIP 各 1 环）
AROMATIC_RINGS = {"PHE": 1, "TYR": 1, "TRP": 2, "HIS": 1,
                  "HID": 1, "HIE": 1, "HIP": 1}


def build_bond_graph(sd) -> Tuple[Dict[int, Set[int]], Dict[int, object]]:
    """从 SystemData 构建全局键合图（残基内 + 残基间）。

    Returns:
        atom_map: global_idx -> AtomData
        bond_graph: global_idx -> set(global_idx)
    """
    atom_map = {}
    bond_graph = defaultdict(set)
    for res in sd.residues:
        for a in res.atoms:
            atom_map[a.atom_global_idx] = a
        for b in res.bonds:
            a1 = res.atoms[b.atom1_idx_in_residue].atom_global_idx
            a2 = res.atoms[b.atom2_idx_in_residue].atom_global_idx
            bond_graph[a1].add(a2)
            bond_graph[a2].add(a1)
    for ub in sd.inter_residue_bonds:
        r1 = sd.residues[ub.residue1_global_idx]
        r2 = sd.residues[ub.residue2_global_idx]
        a1 = r1.atoms[ub.atom_idx_in_residue1].atom_global_idx
        a2 = r2.atoms[ub.atom_idx_in_residue2].atom_global_idx
        bond_graph[a1].add(a2)
        bond_graph[a2].add(a1)
    return atom_map, bond_graph


def _neighbor_elements(atom, atom_map, bond_graph) -> List[str]:
    """返回原子的邻居元素列表。"""
    return [atom_map[i].atom_element
            for i in bond_graph.get(atom.atom_global_idx, set())
            if i in atom_map]


def expected_groups(res, atom_map, bond_graph) -> Dict[str, Set[int]]:
    """按化学事实推导该残基的期望基团（重原子集合）。

    规则（独立于识别器代码，依据标准氨基酸化学）：
    - hydrophobic: C 且所有邻居 ∈ {C, H}
    - H_donor: N/O/S 与 ≥1 个 H 键连
    - H_acceptor: O 一律（含羟基 OH，有孤对电子）；N 无 H 邻居且为
      吡啶/亚胺型（ACCEPTOR_N_TYPES）——含 Pro 亚胺 N、His 无 H 吡啶 N
    - charged_positive/negative: 按残基名 + 侧链原子
    """
    exp = defaultdict(set)
    # 侧链原子名 -> global_idx 映射（供 charged/aromatic 用）
    name_to_idx = {a.atom_name: a.atom_global_idx for a in res.atoms}
    idx_to_name = {a.atom_global_idx: a.atom_name for a in res.atoms}

    for a in res.atoms:
        if a.atom_element == "H":
            continue
        elems = _neighbor_elements(a, atom_map, bond_graph)
        # 疏水
        if a.atom_element == "C" and elems and all(e in ("C", "H") for e in elems):
            exp["hydrophobic"].add(a.atom_global_idx)
        # 供体
        if a.atom_element in ("N", "O", "S") and "H" in elems:
            exp["H_donor"].add(a.atom_global_idx)
        # 受体
        if a.atom_element == "O":
            exp["H_acceptor"].add(a.atom_global_idx)
        elif a.atom_element == "N" and "H" not in elems \
                and a.atom_type in ACCEPTOR_N_TYPES:
            exp["H_acceptor"].add(a.atom_global_idx)

    # 带电基团（侧链原子）
    for p in CHARGED_POS.get(res.residue_name, set()):
        if p in name_to_idx:
            exp["charged_positive"].add(name_to_idx[p])
    for n in CHARGED_NEG.get(res.residue_name, set()):
        if n in name_to_idx:
            exp["charged_negative"].add(name_to_idx[n])

    # 蛋白 N 末端残基的主链 N：氨基 NH3⁺ 带正电（charged_positive）
    # 判定：主链 N 有 ≥3 个 H 邻居（N 末端三质子化）或残基是链首且 N 无肽键
    # 简化：主链 N 且 H 邻居数 ≥ 2（N 端质子化氨基）
    mainchain_n = next((a for a in res.atoms
                        if a.atom_name == "N" and a.atom_element == "N"), None)
    if mainchain_n is not None:
        h_neighbors = _neighbor_elements(mainchain_n, atom_map, bond_graph)
        if h_neighbors.count("H") >= 2:
            exp["charged_positive"].add(mainchain_n.atom_global_idx)

    # 蛋白 C 末端残基：末端羧基 COO⁻（C 与两个 O 键连、无肽键延续）
    # 简化判定：残基含 C 原子且其 O 邻居 ≥2、残基末尾无后续主链 C
    c_terminal = None
    for a in res.atoms:
        if a.atom_name == "C" and a.atom_element == "C":
            o_neighbors = sum(1 for i in bond_graph.get(a.atom_global_idx, set())
                              if i in atom_map and atom_map[i].atom_element == "O")
            if o_neighbors >= 2:
                c_terminal = a
                break
    if c_terminal is not None:
        exp["charged_negative"].add(c_terminal.atom_global_idx)

    return exp


def identified_repr(groups, group_type) -> Set[int]:
    """识别器输出中某基团类型的"代表原子"集合。

    按 group_atom_convention：H_donor/H_acceptor 的代表原子为 atoms[0]
    （D / A）；hydrophobic/charged 同理取 atoms[0]。
    """
    return {g.atoms[0].atom_global_idx
            for g in groups if g.group_type == group_type and g.atoms}


def run_verification(tpr_path: str, ff: str) -> Tuple[int, List[Tuple]]:
    """执行核验，返回 (核验残基数, 问题列表)。

    问题条目: (残基名, 残基号, 类型, 缺/多, 原子名, 原子类型)
    """
    reader = GmxTprReader()
    sd = reader.read(tpr_path)
    identifier = IDENTIFIER_CLASSES[ff]()
    groups = identifier.identify(sd)
    atom_map, bond_graph = build_bond_graph(sd)

    by_res: Dict[Tuple[str, int], List[object]] = defaultdict(list)
    for g in groups:
        by_res[(g.residue_name, g.residue_id)].append(g)

    issues: List[Tuple] = []
    checked = 0
    for res in sd.residues:
        if res.residue_name not in STANDARD_RESIDUES:
            continue
        checked += 1
        exp = expected_groups(res, atom_map, bond_graph)
        res_groups = by_res.get((res.residue_name, res.residue_global_idx), [])
        for gt in ("hydrophobic", "H_donor", "H_acceptor",
                   "charged_positive", "charged_negative"):
            e = exp.get(gt, set())
            a = identified_repr(res_groups, gt)
            for m in (e - a):
                am = atom_map[m]
                issues.append((res.residue_name, res.residue_global_idx, gt,
                               "缺", am.atom_name, am.atom_type))
            for x in (a - e):
                ax = atom_map[x]
                issues.append((res.residue_name, res.residue_global_idx, gt,
                               "多", ax.atom_name, ax.atom_type))
    return checked, issues


def write_report(checked: int, issues: List[Tuple], ff: str, tpr: str,
                 report_path: str) -> None:
    """写核验报告（markdown）。"""
    missing = [i for i in issues if i[3] == "缺"]
    extra = [i for i in issues if i[3] == "多"]
    lines = [
        "# 基团识别自动核验报告",
        "",
        f"- 体系: `{tpr}`",
        f"- 力场: `{ff}`",
        f"- 核验标准氨基酸残基数: {checked}",
        f"- 差异总数: {len(issues)}（缺失 {len(missing)} / 多余 {len(extra)}）",
        "",
        "## 缺失（识别器漏掉化学期望的基团）",
        "",
    ]
    if missing:
        for it in missing:
            lines.append(f"- {it[0]}({it[1]}) {it[2]}缺 {it[4]}({it[5]})")
    else:
        lines.append("无。")
    lines += ["", "## 多余（识别器识别了，需人工甄别）", ""]
    if extra:
        by_atom = Counter(f"{it[4]}({it[5]})" for it in extra)
        lines.append("按原子类型分布：" + str(dict(by_atom)))
        lines.append("")
        lines.append("明细（全部）：")
        for it in extra:
            lines.append(f"- {it[0]}({it[1]}) {it[2]}多 {it[4]}({it[5]})")
    else:
        lines.append("无。")
    lines += ["", "## 说明", "",
              "本报告为自动化初筛：\"多余\"条目需人工对照化学事实甄别"
              "（可能是基准规则偏窄，也可能是识别器误报）。",
              "判别依据见 `doc/group_identification_review.md`。", ""]
    Path(report_path).write_text("\n".join(lines), encoding="utf-8")
    print(f"报告已写入 {report_path}")


def main(argv=None) -> int:
    """CLI 入口。"""
    parser = argparse.ArgumentParser(description="基团识别自动核验")
    parser.add_argument("--tpr", required=True, help="GROMACS tpr 路径")
    parser.add_argument("--ff", required=True, choices=sorted(IDENTIFIER_CLASSES),
                        help="力场名称")
    parser.add_argument("--report", default="/tmp/identification_review.md",
                        help="报告输出路径")
    args = parser.parse_args(argv)

    checked, issues = run_verification(args.tpr, args.ff)
    write_report(checked, issues, args.ff, args.tpr, args.report)

    missing = sum(1 for i in issues if i[3] == "缺")
    extra = len(issues) - missing
    print(f"核验残基数: {checked}  缺失: {missing}  多余: {extra}")
    return 0 if missing == 0 else 1


if __name__ == "__main__":
    sys.exit(main())