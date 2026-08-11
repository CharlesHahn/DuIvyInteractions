#!/usr/bin/env python3
"""解析 gmx dump 文本，提取每分子类型（moltype）的原子、类型名、连接（Bond+Constraint）、残基。

输入：dump_md_D927.tpr.txt（gmx dump -s md.tpr 输出）
输出：结构字典（内存），供官能团鉴定使用。

gmx dump 固定格式（2024.x）：
  moltype (N):
     name="X"
     atoms:
        atom (n):
           atom[i]={type= t, typeB= t, ptype= Atom, m=..., q=..., mB=..., qB=..., resind= r, atomnumber= z}
        atom (n):            # 第二块：原子名
           atom[i]={name="XX"}
        type (n):            # 第三块：类型名
           type[i]={name="T",nameB="T"}
        residue (n):
           residue[i]={name="RN", nr=rr, ic=' '}
     excls: ...
     Bond: ... (type=X (BONDS) i j)
     Constraint: ... (type=X (CONSTR) i j)
"""
import re
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional

# 正则
RE_HEADER = re.compile(r"^topology:")
RE_MOLTYPE_HEAD = re.compile(r"^\s+moltype \((\d+)\):")
RE_MOLTYPE_NAME = re.compile(r'^\s+name="([^"]+)"')
RE_ATOM_COUNT = re.compile(r"^\s+atom \((\d+)\):")
RE_ATOM_PARAM = re.compile(
    r'^\s+atom\[ *(\d+)\]=\{type= *(\d+), typeB= *(\d+), ptype= *(\w+), m= *([\d.eE+-]+), q= *([\d.eE+-]+),'
    r" mB= *([\d.eE+-]+), qB= *([\d.eE+-]+), resind= *(\d+), atomnumber= *(\d+)\}"
)
RE_ATOM_NAME = re.compile(r'^\s+atom\[ *(\d+)\]=\{name="([^"]+)"\}')
RE_TYPE_COUNT = re.compile(r"^\s+type \((\d+)\):")
RE_TYPE_NAME = re.compile(r'^\s+type\[ *(\d+)\]=\{name="([^"]+)",nameB="([^"]+)"\}')
RE_RES_COUNT = re.compile(r"^\s+residue \((\d+)\):")
RE_RES_NAME = re.compile(r'^\s+residue\[ *(\d+)\]=\{name="([^"]+)", nr= *(\d+), ic=\' ?\'')
RE_ILIST_HEAD = re.compile(r"^\s+(Bond|G96Bond|Morse|Cubic Bonds|Connect Bonds|Harmonic Pot\.|FENE Bonds|Tab\. Bonds|Tab\. Bonds NC|Restraint Pot\.|Angle|G96Angle|Restr\. Angles|Lin\. Angle|Bond-Cross|BA-Cross|Urey-Bradley|Quartic Angles|Tab\. Angles|Proper Dih\.|Ryckaert-Bell\.|Restr\. Dih\.|CBT Dih\.|Fourier Dih\.|Improper Dih\.|Tab\. Dih\.|CMAP|Constraint|Constr\. No Co\.|Settle|Position Rest\.|Dis\. Rest\.|Ori\. Rest\.|Angle Rest\.|Dih\. Rest\.):\s*$")
RE_ILIST_ENTRY = re.compile(r"^\s+\d+ type=(\d+) \((\w+)\)\s+([\d\s]+)$")
RE_ILIST_NR = re.compile(r"^\s+nr: (\d+)$")


@dataclass
class AtomInfo:
    idx: int          # moltype 内原子索引 (0-based)
    type_par: int     # type 参数编号（ffparams 索引）
    ptype: str        # Atom/Shell/Dummy...
    mass: float
    charge: float
    resind: int       # 残基索引
    z: int            # 原子序数 (atomnumber)
    name: str = ""
    type_name: str = ""   # 类型名（如 c3 / ca / N3）
    resname: str = ""
    resid: int = 0     # 原始残基号 nr


@dataclass
class MolType:
    idx: int
    name: str
    atoms: List[AtomInfo] = field(default_factory=list)
    bonds: List[Tuple[int, int]] = field(default_factory=list)       # (i,j) moltype 内原子索引
    constraints: List[Tuple[int, int]] = field(default_factory=list) # 同上
    residues: List[Tuple[int, str, int]] = field(default_factory=list)  # (resind, resname, nr)
    # 其他段（Angle/Dihedral/4体）暂不解析，锚定 Bond/Constraint 段边界即可


def parse_dump(path: str) -> List[MolType]:
    moltypes: List[MolType] = []
    cur: Optional[MolType] = None
    section: str = ""  # "", "atoms_param", "atoms_name", "atoms_type", "res", "ilist"
    ilist_section: str = ""   # Bond/Constraint/...
    ilist_nr: int = 0
    ilist_seen: int = 0
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            # 顶层新 moltype
            m = RE_MOLTYPE_HEAD.match(line)
            if m:
                cur = MolType(idx=int(m.group(1)), name="")
                moltypes.append(cur)
                section = ""
                continue
            if cur is None:
                continue
            # moltype name
            m = RE_MOLTYPE_NAME.match(line)
            if m:
                cur.name = m.group(1)
                continue
            # atom 段计数
            m = RE_ATOM_COUNT.match(line)
            if m:
                n = int(m.group(1))
                # 判断是哪个子块：之前已看到 name= 则为参数块；没有则按顺序
                if section == "" or section == "ilist":
                    section = "atoms_param"
                elif section == "atoms_param":
                    section = "atoms_name"
                continue
            # 原子参数
            m = RE_ATOM_PARAM.match(line)
            if m:
                if section == "atoms_param":
                    idx = int(m.group(1))
                    a = AtomInfo(
                        idx=idx,
                        type_par=int(m.group(2)),
                        ptype=m.group(4),
                        mass=float(m.group(5)),
                        charge=float(m.group(6)),
                        resind=int(m.group(9)),
                        z=int(m.group(10)),
                    )
                    cur.atoms.append(a)
                continue
            # 原子名
            m = RE_ATOM_NAME.match(line)
            if m:
                if section == "atoms_name":
                    idx = int(m.group(1))
                    # 顺序对应 atoms_param
                    if idx < len(cur.atoms):
                        cur.atoms[idx].name = m.group(2)
                continue
            # type 段
            m = RE_TYPE_COUNT.match(line)
            if m:
                section = "atoms_type"
                continue
            m = RE_TYPE_NAME.match(line)
            if m:
                if section == "atoms_type":
                    idx = int(m.group(1))
                    if idx < len(cur.atoms):
                        cur.atoms[idx].type_name = m.group(2)
                continue
            # residue 段
            m = RE_RES_COUNT.match(line)
            if m:
                section = "res"
                continue
            m = RE_RES_NAME.match(line)
            if m:
                if section == "res":
                    ri = int(m.group(1))
                    cur.residues.append((ri, m.group(2), int(m.group(3))))
                continue
            # 交互作用列表段（Bond/Constraint 等）
            m = RE_ILIST_HEAD.match(line)
            if m:
                ilist_section = m.group(1)
                section = "ilist"
                ilist_nr = 0
                ilist_seen = 0
                continue
            if section == "ilist":
                m = RE_ILIST_NR.match(line)
                if m:
                    ilist_nr = int(m.group(1))
                    continue
                m = RE_ILIST_ENTRY.match(line)
                if m:
                    # m.group(1)=type 参数编号, m.group(2)=段名, m.group(3)=原子索引串
                    if ilist_section in ("Bond", "Constraint"):
                        atoms_str = m.group(3).split()
                        if len(atoms_str) >= 2:
                            a1, a2 = int(atoms_str[0]), int(atoms_str[1])
                            if ilist_section == "Bond":
                                cur.bonds.append((a1, a2))
                            else:
                                cur.constraints.append((a1, a2))
    return moltypes


def fill_residues(mt: MolType) -> None:
    """把 atom 的 resind 映射到 resname/resid"""
    res_by_idx = {ri: (rn, nr) for ri, rn, nr in mt.residues}
    for a in mt.atoms:
        if a.resind in res_by_idx:
            a.resname, a.resid = res_by_idx[a.resind]


if __name__ == "__main__":
    import sys
    path = sys.argv[1] if len(sys.argv) > 1 else "dump_md_D927.tpr.txt"
    mts = parse_dump(path)
    print(f"解析到 {len(mts)} 个 moltype:")
    for mt in mts:
        fill_residues(mt)
        print(f"  [{mt.idx}] {mt.name}: 原子 {len(mt.atoms)}, "
              f"键 {len(mt.bonds)}, 约束 {len(mt.constraints)}, "
              f"残基 {len(mt.residues)}")
        # 展示前 3 个原子
        for a in mt.atoms[:3]:
            print(f"      idx={a.idx} name={a.name} type={a.type_name} q={a.charge:.3f} z={a.z} res={a.resname}({a.resid})")