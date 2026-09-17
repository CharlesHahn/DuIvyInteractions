# 力场类型映射

本项目直接从 GROMACS tpr 的力场原子类型判定化学基团。本文档列出 Amber 家族（Amber 蛋白 + GAFF 配体）类型到化学特征的映射规则，帮助科研用户理解判定依据。

## 映射原则

类型→化学特征的映射来自两条独立的、可交叉验证的证据链：

1. **GAFF 类型**：来自 antechamber 官方命名规则（Wang et al., J Comput Chem 2004）——类型名的后缀字母/数字编码化学特征（如 `ca` 的 `a` = aromatic）
2. **Amber 蛋白类型**：从 rtp 残基定义反推——已知残基的化学事实（如 TYR 苯环是芳香环）对照 rtp 文件中环原子使用的类型名（`CA`），推出该类型 = 芳香碳

两条链互相独立，结论一致，已验证对 Amber 全家族**零冲突**。

## 芳香类型（STRONG_AROMATIC）

这些类型的原子**由类型名直接确定**为芳香原子，是环检测的强证据。

| 类别 | 类型 | 说明 |
|:-----|:-----|:-----|
| GAFF 芳香碳 | `ca, cg, ch, cm, cn, cp, cq, c1` | `c`+芳香后缀 |
| GAFF 芳香氮 | `na, nb, nh, ni, nj, n1, n2` | `na`=吡咯型、`nb`=吡啶型 |
| GAFF 芳香磷 | `pb` | — |
| Amber 蛋白芳香碳 | `CA, CB, CC, CK, CM, C5, C6, C7, C*, CW, CR, CN, CV, CQ` | rtp 环原子类型 |
| Amber 蛋白芳香氮 | `NA, NB, NC, N*` | — |

## 兼容类型（COMPATIBLE_TYPES）

非芳香类型，但在环内 **n-1 个原子是芳香类型**的"强制"下可参与共轭。用于处理歧义类型与杂环。

| 类型 | 说明 |
|:-----|:-----|
| `C, N` | Amber 蛋白歧义类型（主链羰基碳/芳环碳 C；酰胺氮/芳香氮 N） |
| `os, ss` | GAFF 呋喃氧 / 噻吩硫 |
| `cc, cd` | GAFF 非纯芳香共轭环碳 |
| `pc, pd` | GAFF 共轭环内 sp2 磷 |

## 受体类型（ACCEPTOR_TYPES）

H 键受体（有孤对电子）的候选类型。

| 类别 | 类型 |
|:-----|:-----|
| GAFF 氧 | `o, o2, oh, os, oe, o1, ow` |
| GAFF 氮 | `n, n2, n3, nb, ni, nj, nc, ne, nf, nk`（排除 `na, nh`——吡咯型可作供体） |
| GAFF 硫 | `s, ss, sh, sx, s2` |
| 卤素 | `f, cl, br, i` |
| Amber 蛋白氧 | `O, OH, O2, OS, OW` |
| Amber 蛋白氮 | `N, N2, N3, NA, NB, N*, NC` |
| Amber 硫 | `S, SH` |

## 供体判定

供体原子 D（N/O/S/F）与其氢 H 之间的**键条目显式存在**（Bond + Constraint 合并，因 N-H 键在 Constraint 段），且 H 带正电荷。这是 100% 零推断的确定性判定。

## 金属离子（METAL_IONS）

来自 PLIP config.py：`Ca, Co, Mg, Mn, Fe, Cu, Zn, Li, Na, K, Rb, Sr, Cs, Ba, Cr, Ni, Ru, Rh, Pd, Ag, Cd, La, W, Os, Ir, Pt, Au, Hg, Ce, Pr, Sm, Eu, Gd, Tb, Yb, Lu, Al, Ga, In, Sb, Tl, Pb`。

## 水分子（WATER_RESIDUES）

残基名 `SOL, HOH, WAT`，氧原子 OW、氢原子 HW 由名字标识。

## 关键类型名坑

- **`C` 双身份**：既作主链羰基碳，也作 Tyr CZ 芳环碳。判定需"环内 ≥4 个强芳香邻居"升级。
- **`N3` vs `n3`**：大写 `N3` = 蛋白正电氨基，小写 `n3` = GAFF 中性氨基，含义不同。
- **`CA` 跨力场不同义**：amber = 芳香碳，GROMOS = α 碳。必须通过特征映射表区分。
- **N-H 键在 Constraint 段**：供体识别必须合并 Bond + Constraint 两类键。

## 兼容性范围

已验证的力场：amber03、amber94、amber96、amber99、amber99sb、amber99sb-ildn、amberGS、amber14sb + GAFF/GAFF2 配体。类型名跨版本差异（如 CX vs CT）已全部处理。若出现新版本新类型，只需在 `STRONG_AROMATIC`/`COMPATIBLE_TYPES`/`ACCEPTOR_TYPES` 等特征表补一条即可（特征空间映射设计目标）。