#!/usr/bin/env python3
"""官能团鉴定：从类型特征映射表 + 键连接 → 可参与相互作用的化学基团。

核心设计：
1. 类型特征映射表（type → 化学特征向量），覆盖 GAFF + amber14sb 蛋白类型
2. 环检测（SSSR 最小环）+ 芳香性验证 → 芳香环（π-π/π-阳/卤-π 判定基础）
3. H 键供体/受体鉴定（基于 D-H 键 + 类型特征 + 电荷）
4. 盐桥/带电基团、卤素（σ-hole）、金属中心

特征向量字段：
  element: 元素符号
  hyb: 杂化 (sp3/sp2/sp/d)
  aromatic: 是否芳香
  has_h: 是否带 H（通过连接确认）
  lone_pair: 是否有孤对（受体潜力）
  polarity: 极性类别 (nonpolar/polar/negative/positive)
"""
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Set, Optional
from parse_tpr_dump import parse_dump, fill_residues, MolType, AtomInfo

# ============================================================
# 1. 类型特征映射表（GAFF + amber14sb 蛋白类型）
# ============================================================
# 格式: 类型名 -> (element, hyb, aromatic, lone_pair, polarity)
#  polarity: n=nonpolar p=polar n-=负电 n+=正电 a=芳香
#  (has_h 由连接关系动态判定，不在此表)

TYPE_FEATURES = {
    # ------- GAFF 配体/非标残基类型 -------
    "c":    ("C", "sp2", False, False, "polar"),     # 羰基碳/酰胺碳
    "c3":   ("C", "sp3", False, False, "nonpolar"),  # sp3 饱和碳
    "c2":   ("C", "sp2", False, False, "nonpolar"),  # sp2 烯碳（非芳香）
    "ca":   ("C", "sp2", True,  False, "aromatic"),  # 芳香碳
    "cc":   ("C", "sp2", False, False, "polar"),     # 羧酸根碳
    "cd":   ("C", "sp",  False, False, "polar"),     # 腈碳
    "ce":   ("C", "sp2", False, False, "polar"),     # 酰胺碳
    "cf":   ("C", "sp2", False, False, "polar"),     # 氟代碳
    "cg":   ("C", "sp2", True,  False, "aromatic"),  # 芳香碳（吡咯 α）
    "ch":   ("C", "sp2", True,  False, "aromatic"),  # 芳香碳（咪唑）
    "cj":   ("C", "sp2", False, False, "polar"),     # 亚胺碳
    "ck":   ("C", "sp2", False, False, "polar"),     # 亚胺碳
    "cm":   ("C", "sp2", True,  False, "aromatic"),  # 吡啶 α 碳
    "cn":   ("C", "sp2", True,  False, "aromatic"),  # 吡啶 β 碳
    "cp":   ("C", "sp2", True,  False, "aromatic"),  # 吡咯/三唑碳
    "cq":   ("C", "sp2", True,  False, "aromatic"),  # 吡唑碳
    "cv":   ("C", "sp3", False, False, "polar"),     # 胍碳
    "c1":   ("C", "sp2", True,  False, "aromatic"),  # 杂芳环碳
    "c4":   ("C", "sp3", False, False, "nonpolar"),  # 季碳
    "n":    ("N", "sp2", False, False, "polar"),     # 酰胺氮/亚胺氮
    "n1":   ("N", "sp2", True,  False, "aromatic"),  # 吡咯氮
    "n2":   ("N", "sp2", True,  False, "aromatic"),  # 吡啶氮
    "n3":   ("N", "sp3", False, False, "polar"),     # 氨基氮（GAFF 注意：与蛋白 N3 同名不同义！）
    "n4":   ("N", "sp3", False, False, "positive"),  # 季铵氮
    "na":   ("N", "sp2", True,  True,  "aromatic"),  # 吡咯型芳香氮（带H为供体，无H为受体）
    "nb":   ("N", "sp2", True,  True,  "aromatic"),  # 吡啶型芳香氮（无H，受体）
    "nc":   ("N", "sp3", False, True,  "polar"),     # 胺氮
    "nd":   ("N", "sp2", False, False, "polar"),     # 酰胺氮
    "ne":   ("N", "sp3", False, True,  "polar"),     # 胺氮
    "nf":   ("N", "sp3", False, True,  "polar"),     # 胺氮
    "nh":   ("N", "sp2", True,  True,  "aromatic"),  # 吡咯型芳香氮（带H供体）
    "ni":   ("N", "sp2", True,  True,  "aromatic"),  # 咪唑氮
    "nj":   ("N", "sp2", True,  True,  "aromatic"),  # 咪唑氮
    "nk":   ("N", "sp2", False, True,  "polar"),     # 亚胺氮
    "n1":   ("N", "sp3", False, True,  "polar"),     # 氨氮
    "o":    ("O", "sp2", False, True,  "negative"),  # 羰基氧（受体）
    "o2":   ("O", "sp3", False, True,  "negative"),  # 羧酸根氧
    "oh":   ("O", "sp3", False, True,  "polar"),     # 羟基氧（供体+受体）
    "os":   ("O", "sp3", False, True,  "polar"),     # 醚氧（受体）
    "oe":   ("O", "sp2", False, True,  "negative"),  # 酯氧
    "o1":   ("O", "sp3", False, True,  "polar"),     # 水氧
    "ow":   ("O", "sp3", False, True,  "polar"),     # 水氧（TIP3P）
    "s":    ("S", "sp2", False, False, "nonpolar"),  # 硫（双键）
    "ss":   ("S", "sp3", False, True,  "nonpolar"),  # 硫醚硫
    "sh":   ("S", "sp3", False, True,  "polar"),     # 巯基硫（供体）
    "sx":   ("S", "sp3", False, True,  "polar"),     # 亚砜硫
    "s2":   ("S", "sp2", False, True,  "nonpolar"),  # 硫代羰基硫
    "p2":   ("P", "sp3", False, False, "polar"),     # 亚磷酸酯磷
    "p3":   ("P", "sp3", False, False, "polar"),     # 磷酸酯磷
    "p5":   ("P", "sp3", False, False, "polar"),     # 磷酸酯磷（磷酸基）
    "P":    ("P", "sp3", False, False, "polar"),     # 核酸骨架磷酸磷（amber RNA/DNA）
    "f":    ("F", "sp3", False, True,  "negative"),  # 氟（卤键 σ-hole）
    "cl":   ("Cl","sp3", False, True,  "negative"),  # 氯
    "br":   ("Br","sp3", False, True,  "negative"),  # 溴
    "i":    ("I", "sp3", False, True,  "negative"),  # 碘
    "hc":   ("H", "sp3", False, False, "nonpolar"),  # 烷基氢
    "h1":   ("H", "sp3", False, False, "nonpolar"),  # 烷基氢
    "h2":   ("H", "sp3", False, False, "nonpolar"),  # 烷基氢
    "h3":   ("H", "sp3", False, False, "nonpolar"),  # 烷基氢
    "h4":   ("H", "sp3", False, False, "nonpolar"),  # 烷基氢
    "h5":   ("H", "sp3", False, False, "nonpolar"),  # 烷基氢
    "ha":   ("H", "sp2", False, False, "nonpolar"),  # 芳香氢
    "hn":   ("H", "sp3", False, False, "polar"),     # 胺/酰胺氢
    "ho":   ("H", "sp3", False, False, "polar"),     # 羟基氢
    "hs":   ("H", "sp3", False, False, "polar"),     # 巯基氢
    "hw":   ("H", "sp3", False, False, "polar"),     # 水氢
    "hx":   ("H", "sp3", False, False, "polar"),     # 卤代氢
    # ------- amber14sb 蛋白类型（GROMACS 移植版方言，基于 rtp 实证）-------
    "C":    ("C", "sp2", None, False, "polar"),     # 双身份：主链羰基C / Tyr CZ 芳环碳（环境感知）
    "CA":   ("C", "sp2", True,  False, "aromatic"),  # 芳香碳（Phe/Tyr/Trp/His/Arg 侧链）
    "CB":   ("C", "sp2", True,  False, "aromatic"),  # Trp CD2（芳香）
    "CC":   ("C", "sp2", True,  False, "aromatic"),  # His CG（咪唑环碳，芳香）——rtp 实证修正
    "CD":   ("C", "sp3", False, False, "nonpolar"),  # 脂肪碳
    "CE":   ("C", "sp3", False, False, "nonpolar"),  # 脂肪碳
    "CF":   ("C", "sp3", False, False, "nonpolar"),  # 脂肪碳
    "CG":   ("C", "sp3", False, False, "nonpolar"),  # 脂肪碳
    "CH":   ("C", "sp3", False, False, "nonpolar"),  # 脂肪碳
    "CI":   ("C", "sp3", False, False, "nonpolar"),  # 脂肪碳
    "CJ":   ("C", "sp3", False, False, "nonpolar"),  # 脂肪碳
    "CK":   ("C", "sp2", True,  False, "aromatic"),  # 嘌呤 C8 芳香碳——rtp 实证修正
    "CL":   ("C", "sp3", False, False, "nonpolar"),  # 脂肪碳
    "CM":   ("C", "sp2", True,  False, "aromatic"),  # 嘧啶 C5/C6 芳香碳——rtp 实证修正
    "CN":   ("C", "sp2", True,  False, "aromatic"),  # Trp CE2（芳香）
    "CO":   ("C", "sp2", False, False, "negative"),  # Asp/Glu 羧基碳（连接 O2）
    "C4":   ("C", "sp2", True,  False, "aromatic"),  # amber03 核酸嘧啶 C5/C6 芳香碳
    "CP":   ("C", "sp2", True,  False, "aromatic"),  # amber03 核酸嘌呤 C8 芳香碳
    "CS":   ("C", "sp2", True,  False, "aromatic"),  # amber03 核酸嘧啶 C5/C6 芳香碳
    "CQ":   ("C", "sp2", True,  False, "aromatic"),  # 核酸嘧啶 C2 芳香碳（rna/dna）
    "CT":   ("C", "sp3", False, False, "nonpolar"),  # sp3 碳（amber99sb-ildn 等旧版 amber 的 α/侧链碳）
    "CV":   ("C", "sp2", True,  False, "aromatic"),  # 核酸碱基芳香碳（amber99sb-ildn rna/dna）
    "CX":   ("C", "sp3", False, False, "nonpolar"),  # α 碳（sp3）
    "CY":   ("C", "sp3", False, False, "nonpolar"),  # α 碳（sp3）
    "1C":   ("C", "sp3", False, False, "nonpolar"),  # 甲基碳
    "2C":   ("C", "sp3", False, False, "nonpolar"),  # 亚甲基碳
    "3C":   ("C", "sp3", False, False, "nonpolar"),  # 次甲基碳
    "4C":   ("C", "sp3", False, False, "nonpolar"),  # 季碳
    "C*":   ("C", "sp2", True,  False, "aromatic"),  # Trp CG（吡咯 3 位）
    "C5":   ("C", "sp2", True,  False, "aromatic"),  # 芳香碳（Trp）
    "C6":   ("C", "sp2", True,  False, "aromatic"),  # 芳香碳（His）
    "C7":   ("C", "sp2", True,  False, "aromatic"),  # 芳香碳（His）
    "C8":   ("C", "sp3", False, False, "nonpolar"),  # Arg/Lys 侧链碳（CB/CD/CG 等）
    "CW":   ("C", "sp2", True,  False, "aromatic"),  # His CD1/Trp CD1（芳香）
    "CR":   ("C", "sp2", True,  False, "aromatic"),  # His CE1（芳香）
    "N":    ("N", "sp2", False, True,  "polar"),     # 主链酰胺氮（供体）
    "N*":   ("N", "sp2", True,  True,  "aromatic"),  # 核酸嘌呤 N9 芳香氮（amber99sb-ildn）
    "NC":   ("N", "sp2", True,  True,  "aromatic"),  # 核酸嘧啶 N1/N3 芳香氮——rtp 实证修正
    "N3":   ("N", "sp3", False, True,  "positive"),  # Lys NZ/主链 N端（+。正电氨基）
    "N2":   ("N", "sp2", True,  False, "positive"),  # Arg 胍基氮（NE/NH1/NH2 正电）
    "NA":   ("N", "sp2", True,  True,  "aromatic"),  # His NE2/Trp NE1（吡咯型 供体）
    "NB":   ("N", "sp2", True,  True,  "aromatic"),  # His ND1（吡啶型 受体）
    "O":    ("O", "sp2", False, True,  "negative"),  # 羰基氧（受体）
    "O2":   ("O", "sp3", False, True,  "negative"),  # 羧酸根氧
    "OH":   ("O", "sp3", False, True,  "polar"),     # 羟基氧（Tyr/Ser/Thr）
    "OS":   ("O", "sp3", False, True,  "polar"),     # 醚氧
    "S":    ("S", "sp3", False, False, "nonpolar"),  # 硫（Met）
    "SH":   ("S", "sp3", False, True,  "polar"),     # 巯基硫（Cys）
    "H":    ("H", "sp3", False, False, "polar"),     # 酰胺氢
    "HC":   ("H", "sp3", False, False, "nonpolar"),  # 烷基氢
    "H1":   ("H", "sp3", False, False, "nonpolar"),  # 烷基氢
    "H2":   ("H", "sp3", False, False, "nonpolar"),  # 烷基氢
    "H3":   ("H", "sp3", False, False, "nonpolar"),  # 烷基氢
    "H0":   ("H", "sp3", False, False, "polar"),     # α 氢（amber03 的 HA1/HA2）
    "H4":   ("H", "sp2", False, False, "nonpolar"),  # 芳香氢（His/Trp）
    "H5":   ("H", "sp2", False, False, "nonpolar"),  # 芳香氢（His）
    "HP":   ("H", "sp3", False, False, "nonpolar"),  # 芳香氢（Phe/Tyr/Trp）
    "HA":   ("H", "sp2", False, False, "nonpolar"),  # 芳香氢（Phe/Tyr/Trp 环上）
    "HO":   ("H", "sp3", False, False, "polar"),     # 羟基氢
    "HS":   ("H", "sp3", False, False, "polar"),     # 巯基氢
    "HW":   ("H", "sp3", False, False, "polar"),     # 水氢
    # ------- 水/离子 -------
    "OW":   ("O", "sp3", False, True,  "polar"),     # TIP3P 水氧
    "HW":   ("H", "sp3", False, False, "polar"),     # TIP3P 水氢
    "Na":   ("Na","ion", False, False, "positive"),  # 钠离子
    "Cl":   ("Cl","ion", False, True,  "negative"),  # 氯离子
    "MG":   ("Mg","ion", False, False, "positive"),  # 镁离子
    # ------- 其他（未映射时降级为元素猜测）-------
}


def get_feature(type_name: str) -> Optional[Tuple]:
    """返回类型特征元组，未知类型返回 None"""
    return TYPE_FEATURES.get(type_name)


# ============================================================
# 2. 图论环检测（SSSR 最小环）
# ============================================================
def find_rings(bonds: List[Tuple[int, int]], max_ring_size: int = 12) -> List[List[int]]:
    """用深度优先搜索找最小环（SSSR 近似：对每个边找最短回路）。"""
    adj: Dict[int, Set[int]] = {}
    for i, j in bonds:
        adj.setdefault(i, set()).add(j)
        adj.setdefault(j, set()).add(i)

    rings: List[List[int]] = []
    seen_edges: Set[Tuple[int, int]] = set()

    for u, v in bonds:
        if (u, v) in seen_edges or (v, u) in seen_edges:
            continue
        # BFS 从 u 到 v 但不经过边(u,v)，找最短回路
        from collections import deque
        q = deque([(u, [u])])
        visited = {u}
        found = None
        while q and found is None:
            node, path = q.popleft()
            if len(path) > max_ring_size:
                continue
            for nb in adj[node]:
                if nb == u and len(path) > 1:
                    # 回到起点 → 环路
                    if nb == v:
                        pass
                    continue
                if nb not in visited:
                    visited.add(nb)
                    q.append((nb, path + [nb]))
        # 简化：改用找从 v 到 u 的路径
        # 用最短路 BFS（不带边(u,v)）
        q = deque([(v, [v])])
        visited = {v}
        while q:
            node, path = q.popleft()
            for nb in adj[node]:
                if nb == u and len(path) >= 2:
                    full = path + [u]
                    if sorted(full) not in [sorted(r) for r in rings]:
                        rings.append(full)
                    continue
                if nb not in visited and len(path) < max_ring_size:
                    visited.add(nb)
                    q.append((nb, path + [nb]))
        seen_edges.add((u, v))

    return rings


# 更简洁可靠的环检测：全部环枚举（环大小上限）
def find_all_rings(bonds: List[Tuple[int, int]], max_size: int = 10) -> List[List[int]]:
    """枚举所有长度 <= max_size 的简单环（去重，按顶点集）。

    算法：对每条边 (u,v)，删除该边后找 u→v 的最短路径；路径+边组成一个环。
    去重：按排序后的顶点元组。环大小限制 max_size。
    """
    from collections import deque
    adj: Dict[int, List[int]] = {}
    for i, j in bonds:
        adj.setdefault(i, []).append(j)
        adj.setdefault(j, []).append(i)

    rings: List[List[int]] = []
    seen: Set[Tuple[int, ...]] = set()

    for u, v in bonds:
        # BFS 从 u 到 v，不走边 (u,v) 本身
        q = deque([(u, [u])])
        visited = {u}
        while q:
            node, path = q.popleft()
            if len(path) > max_size:
                continue
            for nb in adj[node]:
                # 跳过边 (u,v) / (v,u)
                if (node == u and nb == v) or (node == v and nb == u):
                    continue
                if nb == v:  # 到达 v，成环
                    full = path + [v]
                    if len(full) - 1 <= max_size:  # 环边数
                        key = tuple(sorted(full))
                        if key not in seen:
                            seen.add(key)
                            rings.append(full)
                    continue
                if nb not in visited:
                    visited.add(nb)
                    q.append((nb, path + [nb]))
    # 冗余环去除：若环 R 的全部边都被更小的环覆盖，则 R 是稠合冗余环，删除
    def ring_edges(r):
        n = len(r)
        return {(min(r[i], r[(i + 1) % n]), max(r[i], r[(i + 1) % n])) for i in range(n)}
    rings_sorted = sorted(rings, key=len)
    final_rings = []
    coverage: Set[Tuple[int, int]] = set()  # 已被更小环覆盖的边
    for r in rings_sorted:
        e = ring_edges(r)
        if e.issubset(coverage):  # 全部边已被覆盖 → 冗余
            continue
        final_rings.append(r)
        coverage |= e
    return final_rings


# ============================================================
# 3. 官能团鉴定
# ============================================================
@dataclass
class AromaticRing:
    atoms: List[int]          # 环内原子索引
    aromatic: bool            # 是否全芳香
    elements: List[str]

# 强芳香类型（直接由类型名确定的芳香原子）
STRONG_AROMATIC = {"ca", "cg", "ch", "cm", "cn", "cp", "cq", "c1",
                   "n1", "n2", "na", "nb", "nh", "ni", "nj",
                   "CA", "CB", "CC", "CK", "CM", "C5", "C6", "C7", "C*", "CW", "CR", "CN", "CV", "CQ",
                   "NA", "NB", "NC", "N*"}


@dataclass
class Donor:
    atom: int
    h_atom: int
    type_name: str


@dataclass
class Acceptor:
    atom: int
    type_name: str


@dataclass
class ChargedGroup:
    atom: int
    charge: float
    sign: str
    type_name: str = ""


def identify_functional_groups(mt: MolType) -> Dict:
    """对一个 moltype 做官能团鉴定。返回结构化结果。"""
    # 连接图
    conn: Dict[int, Set[int]] = {}
    for i, j in mt.bonds:
        conn.setdefault(i, set()).add(j)
        conn.setdefault(j, set()).add(i)
    for i, j in mt.constraints:
        conn.setdefault(i, set()).add(j)
        conn.setdefault(j, set()).add(i)

    # 环检测（先用强芳香类型标记，再对环内歧义类型如 C 做环境感知升级）
    all_rings = find_all_rings(mt.bonds, max_size=8)
    # 去重（按原子集合）
    unique_rings = []
    seen_rings = set()
    for r in all_rings:
        k = tuple(sorted(r))
        if k not in seen_rings:
            seen_rings.add(k)
            unique_rings.append(r)

    rings_out = []
    for r in unique_rings:
        # 环内强芳香计数
        strong = sum(1 for idx in r if mt.atoms[idx].type_name in STRONG_AROMATIC)
        # 环内各原子最终芳香性
        atom_arom = {}
        for idx in r:
            a = mt.atoms[idx]
            if a.type_name in STRONG_AROMATIC:
                atom_arom[idx] = True
            else:
                # 歧义类型（如 C=CZ）：若环内 ≥4 个强芳香原子，则升级为芳香
                atom_arom[idx] = strong >= 4
        elems = [mt.atoms[idx].type_name for idx in r]
        aromatic = all(atom_arom.values())
        rings_out.append(AromaticRing(atoms=r, aromatic=aromatic, elements=elems))

    # 供体：D-H 键且 H 带正电（q(H) > 0）
    # 化学本质：供体的 H 是酸性氢，必须带 δ+；D 通常带 δ-（吸电子使 H 变正）
    donors = []
    for i, j in mt.bonds + mt.constraints:
        a1, a2 = mt.atoms[i], mt.atoms[j]
        if a1.z == 1 and a2.z != 1:
            d_atom, h_atom = a2, a1
        elif a2.z == 1 and a1.z != 1:
            d_atom, h_atom = a1, a2
        else:
            continue
        if d_atom.z in (7, 8, 16, 9) and h_atom.charge > 0:
            donors.append(Donor(atom=d_atom.idx, h_atom=h_atom.idx, type_name=d_atom.type_name))

    # 受体：N/O/F/S 原子本身带负电（或负电环境），有孤对可接受氢键
    # 化学本质：受体是富电子端（δ-），孤对电子与供体的 δ+ H 作用
    acceptors = []
    for a in mt.atoms:
        if a.z not in (7, 8, 9, 16):
            continue
        feat = get_feature(a.type_name)
        if feat is None:
            continue
        # 负电荷或电负性环境才有孤对可共享（排除正电 N 等）
        if feat[3] and a.charge < 0:  # 有孤对 + 负电
            acceptors.append(Acceptor(atom=a.idx, type_name=a.type_name))

    # 带电基团
    charged = []
    for a in mt.atoms:
        if abs(a.charge) > 0.3 and a.z != 1:  # 重原子强电荷
            sign = "+" if a.charge > 0 else "-"
            charged.append(ChargedGroup(atom=a.idx, charge=a.charge, sign=sign, type_name=a.type_name))

    # 卤素（σ-hole 潜在卤键）
    halogens = [a for a in mt.atoms if a.z in (9, 17, 35, 53)]

    # 金属
    metals = [a for a in mt.atoms if a.z in (3, 11, 12, 19, 20, 30, 26, 25, 29)]

    return {
        "rings": rings_out,
        "donors": donors,
        "acceptors": acceptors,
        "charged": charged,
        "halogens": halogens,
        "metals": metals,
        "n_atoms": len(mt.atoms),
    }


def format_report(mt: MolType, result: Dict, offset: int = 0) -> str:
    """格式化输出官能团报告。offset 为全局原子索引偏移（moltype 内 idx → 全局 idx）。"""
    lines = []
    lines.append(f"=== {mt.name} ({len(mt.atoms)} atoms) ===")
    lines.append(f"  芳香环 ({sum(1 for r in result['rings'] if r.aromatic)}):")
    for r in result["rings"]:
        if r.aromatic:
            names = [f"{mt.atoms[i].name}[{offset + i}]" for i in r.atoms]
            lines.append(f"    ring{len(r.atoms)}: {'-'.join(names)}")
    lines.append(f"  H键供体 ({len(result['donors'])}):")
    for d in result["donors"]:
        d_q = mt.atoms[d.atom].charge
        h_q = mt.atoms[d.h_atom].charge
        d_idx = offset + d.atom
        h_idx = offset + d.h_atom
        lines.append(f"    {mt.atoms[d.atom].name}({d.type_name})[{d_idx}] -H{mt.atoms[d.h_atom].name}[{h_idx}]  qD={d_q:.3f} qH={h_q:.3f}")
    lines.append(f"  H键受体 ({len(result['acceptors'])}):")
    for a in result["acceptors"]:
        a_idx = offset + a.atom
        lines.append(f"    {mt.atoms[a.atom].name}({a.type_name})[{a_idx}]  q={mt.atoms[a.atom].charge:.3f}")
    lines.append(f"  强电荷基团 ({len(result['charged'])}):")
    for c in result["charged"]:
        c_idx = offset + c.atom
        lines.append(f"    {mt.atoms[c.atom].name}({c.type_name})[{c_idx}] {c.sign}{c.charge:.2f}")
    if result["halogens"]:
        hal = [f"{h.name}({h.type_name})[{offset + h.idx}]" for h in result["halogens"]]
        lines.append(f"  卤素: {', '.join(hal)}")
    if result["metals"]:
        met = [f"{m.name}[{offset + m.idx}]" for m in result["metals"]]
        lines.append(f"  金属: {', '.join(met)}")
    return "\n".join(lines)


if __name__ == "__main__":
    import sys
    path = sys.argv[1] if len(sys.argv) > 1 else "dump_md_D927.tpr.txt"
    mts = parse_dump(path)
    # 计算全局原子偏移
    offsets = {}
    cumul = 0
    for mt in mts:
        offsets[mt.name] = cumul
        cumul += len(mt.atoms)
    for mt in mts:
        fill_residues(mt)
        if mt.name in ("RBD_pro", "D927", "KRAS_pro", "GNP_neg"):
            result = identify_functional_groups(mt)
            print(format_report(mt, result, offsets.get(mt.name, 0)))
            print()