#!/usr/bin/env python3
"""类型映射表自动对照验证。

方法：从 rtp 残基定义反推"每个原子类型出现在什么化学环境"，与 TYPE_FEATURES 对照。

核心思路：
1. 解析 rtp：每个残基的每个原子 → (原子名, 类型, 电荷, 连接伙伴类型)
2. 从连接伙伴类型推断化学环境：
   - 芳香性参考：若类型名出现在"已知芳香残基"（TYR/PHE/HIS/TRP + 嘌呤/嘧啶）的环中，则该类型应有芳香特征
   - 键长参考：从 rtp 的键长（若含）或 partner 判断
3. 与 TYPE_FEATURES 逐条对照，报告不一致

覆盖：用户力场（../amber14sb.ff）+ GROMACS 自带 7 个 amber 力场
"""
import os
import sys
import glob
from collections import defaultdict
from typing import Dict, List, Set, Tuple

# 已知芳香残基（化学事实）：其环内原子应为芳香
AROMATIC_RESIDUES = {"TYR", "PHE", "HIS", "HID", "HIE", "HIP", "TRP", "WAT",
                     "ADE", "GUA", "CYT", "URA", "THY", "ADE5", "GUA5", "CYT5", "URA5",
                     "A", "G", "C", "U", "T", "DA", "DG", "DC", "DT",
                     "RA", "RG", "RC", "RU"}

# 已知芳香残基的环原子名（化学事实）
# TYR: CG CD1 CD2 CE1 CE2 CZ
# PHE: CG CD1 CD2 CE1 CE2 CZ
# HIS: CG ND1 CD2 CE1 NE2 (咪唑 5 环)
# TRP: CG CD1 NE1 CE2 CD2 CE3 CZ3 CH2 CZ2 (双环)
# 嘌呤(A/G): C4 C5 C6 C8 + N7 N9 + (N1 C2 N3 吡啶环)
# 嘧啶(C/U/T): C2 C4 C5 C6 (+N1 N3)
RING_ATOMS = {
    "TYR": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "PHE": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "HIS": {"CG", "ND1", "CD2", "CE1", "NE2"},
    "HID": {"CG", "ND1", "CD2", "CE1", "NE2"},
    "HIE": {"CG", "ND1", "CD2", "CE1", "NE2"},
    "HIP": {"CG", "ND1", "CD2", "CE1", "NE2"},
    "TRP": {"CG", "CD1", "NE1", "CE2", "CD2", "CE3", "CZ3", "CH2", "CZ2"},
    # 嘌呤
    "ADE": {"C4", "C5", "C6", "C8", "N7", "N9", "N1", "C2", "N3"},
    "GUA": {"C4", "C5", "C6", "C8", "N7", "N9", "N1", "C2", "N3"},
    "A":   {"C4", "C5", "C6", "C8", "N7", "N9", "N1", "C2", "N3"},
    "G":   {"C4", "C5", "C6", "C8", "N7", "N9", "N1", "C2", "N3"},
    "DA":  {"C4", "C5", "C6", "C8", "N7", "N9", "N1", "C2", "N3"},
    "DG":  {"C4", "C5", "C6", "C8", "N7", "N9", "N1", "C2", "N3"},
    "RA":  {"C4", "C5", "C6", "C8", "N7", "N9", "N1", "C2", "N3"},
    "RG":  {"C4", "C5", "C6", "C8", "N7", "N9", "N1", "C2", "N3"},
    # 嘧啶
    "CYT": {"C2", "C4", "C5", "C6", "N1", "N3"},
    "URA": {"C2", "C4", "C5", "C6", "N1", "N3"},
    "THY": {"C2", "C4", "C5", "C6", "N1", "N3"},
    "C":   {"C2", "C4", "C5", "C6", "N1", "N3"},
    "U":   {"C2", "C4", "C5", "C6", "N1", "N3"},
    "T":   {"C2", "C4", "C5", "C6", "N1", "N3"},
    "DC":  {"C2", "C4", "C5", "C6", "N1", "N3"},
    "DT":  {"C2", "C4", "C5", "C6", "N1", "N3"},
    "DU":  {"C2", "C4", "C5", "C6", "N1", "N3"},
    "RC":  {"C2", "C4", "C5", "C6", "N1", "N3"},
    "RU":  {"C2", "C4", "C5", "C6", "N1", "N3"},
}


def parse_rtp(path: str) -> Dict[str, List[Dict]]:
    """解析 rtp 文件，返回 {残基名: [原子字典]}。原子字典含 name/type/charge/partners。"""
    residues = {}
    current_res = None
    in_atoms = False
    in_bonds = False
    atoms = []
    bond_pairs = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("[") and not line.startswith("[["):
                name = line.strip("[] ").strip()
                if in_atoms and current_res:
                    residues[current_res] = {"atoms": atoms, "bonds": bond_pairs}
                # 新段
                if name in ("atoms",):
                    in_atoms = True
                    in_bonds = False
                elif name in ("bonds",):
                    in_atoms = False
                    in_bonds = True
                elif name == "moleculetype":
                    # 开始新残基
                    if current_res and in_atoms:
                        residues[current_res] = {"atoms": atoms, "bonds": bond_pairs}
                    current_res = None
                    in_atoms = False
                    in_bonds = False
                else:
                    in_atoms = False
                    in_bonds = False
                continue

            if line.startswith(";") or line.startswith("#") or not line:
                continue

            # 残基名（[ XXX ] 段在 atoms 前出现）
            if line.isupper() and len(line) <= 5 and not in_atoms and not in_bonds:
                # 检查是否是残基标题（如 [ TYR ] 后的下一行）
                pass

            if in_atoms:
                parts = line.split()
                if len(parts) >= 4 and not parts[0].startswith("["):
                    atoms.append({
                        "name": parts[0],
                        "type": parts[1],
                        "charge": float(parts[2]),
                        "cg": int(parts[3]) if len(parts) > 3 else 0,
                    })
            elif in_bonds:
                parts = line.split()
                if len(parts) >= 2:
                    bond_pairs.append((parts[0], parts[1]))

    if current_res and in_atoms:
        residues[current_res] = {"atoms": atoms, "bonds": bond_pairs}
    return residues


def parse_rtp_fixed(path: str) -> Dict[str, List[Dict]]:
    """解析 rtp（修正版）：正确处理 [ 残基名 ] 段。"""
    residues = {}
    current_res = None
    section = None
    atoms = []
    bond_pairs = []

    with open(path) as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith(";") or stripped.startswith("#"):
                continue

            if stripped.startswith("["):
                # 段标题可能带行内注释：`[ URE ] ; urea added by EJS`
                name = stripped[1:].split("]")[0].strip()
                if name == "atoms":
                    section = "atoms"
                elif name in ("bonds", "angles", "dihedrals", "impropers", "pairs", "exclusions", "settles", "virtual_sites2", "virtual_sites3", "virtual_sites4", "position_restraints"):
                    section = None
                elif name == "bondedtypes":
                    section = None
                else:
                    # 新残基（如 TYR）或新段（如 moleculetype）
                    if name in ("moleculetype", "atomtypes", "bondtypes", "angletypes", "dihedraltypes", "constrainttypes", "defaults", "system", "molecules", "exclusions", "pairs", "settles", "virtual_sites2", "virtual_sites3", "virtual_sites4", "vsites2", "vsites3", "vsites4", "cmap", "position_restraints", "bondedtypes"):
                        section = None
                        # 保存上一个残基
                        if current_res:
                            residues[current_res] = {"atoms": atoms, "bonds": bond_pairs}
                            atoms = []
                            bond_pairs = []
                        current_res = None
                    else:
                        # 新残基名（大写）
                        if current_res:
                            residues[current_res] = {"atoms": atoms, "bonds": bond_pairs}
                            atoms = []
                            bond_pairs = []
                        current_res = name
                        section = None
                continue

            if section == "atoms" and current_res:
                parts = stripped.split()
                if len(parts) >= 4:
                    atoms.append({
                        "name": parts[0],
                        "type": parts[1],
                        "charge": float(parts[2]),
                        "cg": int(parts[3]) if len(parts) > 3 else 0,
                    })
            elif section == "bonds" and current_res:
                parts = stripped.split()
                if len(parts) >= 2:
                    bond_pairs.append((parts[0], parts[1]))

    if current_res:
        residues[current_res] = {"atoms": atoms, "bonds": bond_pairs}
    return residues


def get_type_environments(residues: Dict) -> Dict[str, Dict]:
    """统计每个类型的化学环境：出现在哪些残基、哪些原子名、partner 类型。"""
    type_env = defaultdict(lambda: {
        "residues": set(),
        "atom_names": set(),
        "in_aromatic_ring": False,
        "partners": defaultdict(int),
        "charge_sum": 0.0,
        "n": 0,
    })

    for resname, data in residues.items():
        if not data["atoms"]:
            continue
        atom_by_name = {a["name"]: a for a in data["atoms"]}
        ring_atoms = RING_ATOMS.get(resname, set())

        for a in data["atoms"]:
            t = a["type"]
            env = type_env[t]
            env["residues"].add(resname)
            env["atom_names"].add(a["name"])
            env["charge_sum"] += a["charge"]
            env["n"] += 1
            if a["name"] in ring_atoms:
                env["in_aromatic_ring"] = True

            # partners（通过 bonds）
            for b1, b2 in data["bonds"]:
                if b1 == a["name"] and b2 in atom_by_name:
                    env["partners"][atom_by_name[b2]["type"]] += 1
                elif b2 == a["name"] and b1 in atom_by_name:
                    env["partners"][atom_by_name[b1]["type"]] += 1

    return type_env


def main():
    # 力场列表：用户力场 + GROMACS 自带 7 个
    top_dir = "/home/hanyl/.micromamba/envs/DIP/share/gromacs/top"
    user_ff = "/mnt/work1/PMO/hanyl/JF_work/KRAS_MD/amber14sb.ff"
    force_fields = [("用户自定义 amber14sb", user_ff)]
    for ff in ["amber03.ff", "amber94.ff", "amber96.ff", "amber99.ff", "amber99sb.ff", "amber99sb-ildn.ff", "amberGS.ff"]:
        force_fields.append((ff, os.path.join(top_dir, ff)))

    # 汇总：类型 → 是否出现在芳香环
    type_in_arom: Dict[str, bool] = defaultdict(bool)
    type_partners: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
    type_residues: Dict[str, Set[str]] = defaultdict(set)
    type_charges: Dict[str, List[float]] = defaultdict(list)

    for ff_name, ff_path in force_fields:
        rtp_files = glob.glob(os.path.join(ff_path, "*.rtp"))
        for rtp in rtp_files:
            residues = parse_rtp_fixed(rtp)
            envs = get_type_environments(residues)
            for t, env in envs.items():
                if env["in_aromatic_ring"]:
                    type_in_arom[t] = True
                type_residues[t] |= env["residues"]
                for p, c in env["partners"].items():
                    type_partners[t][p] += c
                for _ in range(env["n"]):
                    type_charges[t].append(env["charge_sum"] / env["n"] if env["n"] else 0)

    # 导入我们的映射表
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from functional_groups import TYPE_FEATURES

    # 对照：类型在芳香环中出现，但映射表标非芳香
    print("=" * 70)
    print("【验证 1】芳香性对照：rtp 芳香环出现的类型 vs 映射表 aromatic 字段")
    print("=" * 70)
    conflicts_arom = []
    for t in sorted(type_in_arom):
        feat = TYPE_FEATURES.get(t)
        in_arom = type_in_arom[t]
        mapped_arom = bool(feat and feat[2])
        if in_arom and not mapped_arom:
            conflicts_arom.append(t)
            print(f"  ⚠️ {t}: rtp 中出现在芳香环，但映射表 aromatic={mapped_arom} (feat={feat})")
        elif not in_arom and mapped_arom and t not in ("ca", "cp", "cg", "ch", "cm", "cn", "cq", "c1", "na", "nb", "nh", "ni", "nj", "CA", "CB", "C5", "C6", "C7", "C*", "CW", "CR", "CN", "NA", "NB", "CV", "CQ", "N*"):
            # 映射为芳香但 rtp 未出现在芳香环 —— 可能是配体类型（不在蛋白 rtp 中）
            pass
    if not conflicts_arom:
        print("  ✅ 无冲突：所有在芳香环出现的类型都被映射为芳香")

    # 对照 2：类型在芳香环从未出现，但映射表标芳香（配体除外）
    print("\n" + "=" * 70)
    print("【验证 2】所有 amber 力场中类型总数 vs 映射表覆盖")
    print("=" * 70)
    all_ff_types = set(type_in_arom.keys())
    all_ff_types |= set(type_residues.keys())
    missing = sorted(t for t in all_ff_types if t not in TYPE_FEATURES)
    non_critical = {"C0", "Cs", "CS", "IB", "K", "Li", "LI", "MG", "Na", "Rb", "RB", "Zn", "ZN", "URE", "MW",
                    "OW_tip4p", "Cl", "OW", "HW", "OW_spc", "HW_spc", "VDW", "X", "Na+", "Cl-", "Ca"}
    critical_missing = [t for t in missing if t not in non_critical]
    print(f"  全部类型: {len(all_ff_types)}")
    print(f"  缺失: {missing}")
    print(f"  关键缺失: {critical_missing if critical_missing else '无 ✅'}")

    # 对照 3：电荷信息辅助（供体/受体判据可参考）
    print("\n" + "=" * 70)
    print("【验证 3】关键类型的电荷统计（供体/受体判据的参考）")
    print("=" * 70)
    for t in ["N", "N3", "N2", "NA", "NB", "O", "O2", "OH", "S", "SH", "C", "CA", "CX", "CT"]:
        if t in type_charges:
            cs = type_charges[t]
            avg = sum(cs) / len(cs)
            print(f"  {t}: 平均电荷 {avg:.3f} (n={len(cs)})")

    # 对照 4：partner 环境（诊断 C 类型双身份等）
    print("\n" + "=" * 70)
    print("【验证 4】C 类型的 partner 环境（验证双身份）")
    print("=" * 70)
    for t in ["C", "CA", "CX", "CT"]:
        if t in type_partners:
            partners = sorted(type_partners[t].items(), key=lambda x: -x[1])[:8]
            print(f"  {t}: partner 类型 {partners}")


if __name__ == "__main__":
    main()
