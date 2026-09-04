# 相互作用基团识别定义

> 创建日期：2026-08-26
> 更新日期：2026-08-28（补充缺失的基团类型，修正输出格式描述）

---

## 1. H 键供体（H_donor）

**定义**：满足以下全部条件的 D—H 键对。

| 条件 | 判据 |
|------|------|
| 键存在 | 从 bond 列表（残基内键 + 跨残基键）中读取 |
| D 的元素 | N、O、S 或 F |
| H 的电荷 | q(H) > 0 |

**输出**：每个满足条件的 D—H 键生成一个 `Group(group_type="H_donor")`，`atoms = [D, H]`（2 个原子，第一个是供体 D，第二个是氢 H）。`metadata` 为空。

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

**输出**：`atoms = [A]`（1 个原子）。`metadata` 为空。

---

## 3. 卤键供体（halogen_donor）

**定义**：卤素原子连接到碳原子，形成 C—X 对。

| 条件 | 判据 |
|------|------|
| 元素 | F、Cl、Br 或 I |
| 连接 | 必须连接到碳原子（C） |

**输出**：`atoms = [C, X]`（2 个原子，第一个是碳，第二个是卤素）。C→X 方向定义了 σ-hole 的方向。`metadata` 为空。

---

## 4. 卤键受体（halogen_acceptor）

**定义**：满足以下全部条件的原子及其邻居。

| 条件 | 判据 |
|------|------|
| 元素 | C、P 或 S |
| 连接 | 必须连接到 O、P、N 或 S |

**输出**：`atoms = [A, R1, R2, ...]`（1+n 个原子，第一个是受体原子 A，后续是 A 的 O/P/N/S 邻居 R）。检测时取使 X···A-R 角度最接近 120° 的 R。`metadata` 为空。

---

## 5. 芳香环（aromatic_ring）（2026-08-28 补充）

**定义**：满足以下两个条件的环。

| 条件 | 判据 |
|------|------|
| 条件 1：原子类型 | 环内至少 n-1 个原子的类型在 STRONG_AROMIC 中（n 为环大小） |
| 条件 2：兼容原子 | 环内不在 STRONG_AROMIC 中的原子，必须全部在 COMPATIBLE 中 |

**STRONG_AROMIC**：明确标记为芳香的原子类型（GAFF: ca, cg, ch, cm, cn, cp, cq, c1, na, nb, nh, ni, nj, n1, n2, pb; Amber: CA, CB, CC, CK, CM, C5, C6, C7, C*, CW, CR, CN, CV, CQ, NA, NB, NC, N*）。

**COMPATIBLE**：非芳香但在 n-1 个芳香原子"强制"下可参与共轭的类型（C, N, os, ss, cc, cd, pc, pd）。

**详细定义**：见 `aromatic_ring_definition.md`。

**输出**：`atoms = [环内所有原子]`（5-7 个原子），按环顺序排列（BFS 路径顺序）。`metadata` 为空。

---

## 6. 带电基团（charged_positive / charged_negative）（2026-08-28 补充）

**定义**：三层递进识别。

| 层级 | 适用范围 | 方法 |
|------|---------|------|
| 第一层 | 蛋白残基 | 残基名字典（ARG/LYS/HIP/ASP/GLU 等） |
| 第二层 | 非蛋白残基 | 官能团模式匹配（quartamine/tertamine/guanidine/sulfonium/phosphate/carboxylate 等） |
| 第三层 | 通用 | 部分电荷交叉验证（\|Σq\| > 0.1） |

**详细定义**：见 `charged_group_design.md`。

**输出**：`atoms = [电荷中心的所有原子]`（1-13 个原子不等）。

**metadata**：

| 键 | 值 | 说明 |
|----|----|------|
| source | `"residue_name"` 或 `"functional_group"` | 识别来源 |
| func_group | `"tertamine"`, `"carboxylate"`, ... 或 `None` | 官能团类型（第二层填写） |

---

## 7. 金属离子（metal）

**定义**：元素在金属离子列表中的单个原子。

**METAL_IONS 列表**（来自 PLIP config.py）：

| 类别 | 元素 |
|------|------|
| 碱金属 | Li, Na, K, Rb, Cs |
| 碱土金属 | Mg, Ca, Sr, Ba |
| 过渡金属 | Cr, Mn, Fe, Co, Ni, Cu, Zn, Ru, Rh, Pd, Ag, Cd, W, Os, Ir, Pt, Au, Hg |
| 镧系 | La, Ce, Pr, Sm, Eu, Gd, Tb, Yb, Lu |
| 其他 | Al, Ga, In, Sb, Tl, Pb |

**输出**：`atoms = [金属离子]`（1 个原子）。`metadata` 为空。

---

## 8. 水分子（water）

**定义**：残基名在水分子残基名列表中的整个残基。

**WATER_RESIDUES 列表**：

| 残基名 | 来源 |
|--------|------|
| SOL | GROMACS 标准 |
| HOH | PDB 标准 |
| WAT | 部分力场 |

**输出**：`atoms = res.atoms`（该残基的所有原子，通常为 O + H1 + H2，含虚拟位点时更多）。`metadata` 为空。

---

## 9. 疏水原子（hydrophobic）（2026-08-28 补充）

**定义**：满足以下两个条件的单个原子。

| 条件 | 判据 |
|------|------|
| 条件 1 | 原子是碳（element == 'C'） |
| 条件 2 | 所有邻居原子都是碳或氢（neighbor element ∈ {'C', 'H'}） |

**详细定义**：见 `hydrophobic_definition.md`。

**输出**：`atoms = [疏水碳原子]`（1 个原子）。`metadata` 为空。

---

## 10. 金属配位原子（metal_binding）（2026-08-28 补充）

**定义**：满足以下全部条件的单个原子。

| 条件 | 判据 |
|------|------|
| 元素 | O、N 或 S |
| 非水 | 残基名不在 WATER_RESIDUES 中 |

**输出**：`atoms = [配位原子]`（1 个原子）。`metadata = {"source": "element"}`。

---

*文档结束*
