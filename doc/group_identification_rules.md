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

*文档结束*
