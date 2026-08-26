# 相互作用基团识别定义

> 创建日期：2026-08-26

---

## 1. H 键供体（H_donor）

**定义**：满足以下全部条件的 D—H 键对。

| 条件 | 判据 |
|------|------|
| 键存在 | 从 bond 列表（残基内键 + 跨残基键）中读取 |
| D 的元素 | N、O、S 或 F |
| H 的电荷 | q(H) > 0 |

**输出**：每个满足条件的 D—H 键生成一个 `Group(group_type="H_donor")`，`atom_indices` 存 D 原子索引，`metadata["h_atom"]` 存 H 原子索引。

---

## 2. H 键受体（H_acceptor）

**定义**：满足以下全部条件的单个原子。

| 条件 | 判据 |
|------|------|
| 原子类型 | 在 ACCEPTOR_TYPES 列表中 |
| 电荷 | q < 0 |

**ACCEPTOR_TYPES 列表**：

| 类别 | 类型 |
|------|------|
| GAFF 氧 | o, o2, oh, os, oe, o1, ow |
| GAFF 氮 | n, n2, n3, nb, ni, nj, nc, ne, nf, nk |
| GAFF 硫 | s, ss, sh, sx, s2 |
| 卤素 | f, cl, br, i |
| Amber 氧 | O, OH, O2, OS, OW |
| Amber 氮 | N, N2, N3, NA, NB, N*, NC |
| Amber 硫 | S, SH |

**排除的类型**：

| 类型 | 排除原因 |
|------|---------|
| na | 质子化吡啶/吡咯氮，孤对电子不可用 |
| nh | 带 H 的吡咯氮，孤对电子在芳香体系中 |

---

## 3. 卤键供体（halogen）

**定义**：满足以下全部条件的单个原子。

| 条件 | 判据 |
|------|------|
| 元素 | F、Cl、Br 或 I |
| 连接 | 必须连接到碳原子（C） |

---

## 4. 卤键受体（halogen_acceptor）

**定义**：满足以下全部条件的单个原子。

| 条件 | 判据 |
|------|------|
| 元素 | C、P 或 S |
| 连接 | 必须连接到 O、P、N 或 S |

---

## 5. 金属离子（metal）

**定义**：元素在金属离子列表中的单个原子。

**METAL_IONS 列表**（来自 PLIP config.py）：

| 类别 | 元素 |
|------|------|
| 碱金属 | Li, Na, K, Rb, Cs |
| 碱土金属 | Mg, Ca, Sr, Ba |
| 过渡金属 | Cr, Mn, Fe, Co, Ni, Cu, Zn, Ru, Rh, Pd, Ag, Cd, W, Os, Ir, Pt, Au, Hg |
| 镧系 | La, Ce, Pr, Sm, Eu, Gd, Tb, Yb, Lu |
| 其他 | Al, Ga, In, Sb, Tl, Pb |

---

## 6. 水分子（water）

**定义**：残基名在水分子残基名列表中的整个残基。

**WATER_RESIDUES 列表**：

| 残基名 | 来源 |
|--------|------|
| SOL | GROMACS 标准 |
| HOH | PDB 标准 |
| WAT | 部分力场 |

**输出**：每个水残基生成一个 `Group(group_type="water")`，`atom_indices` 包含该水分子的所有原子。

---

*文档结束*
