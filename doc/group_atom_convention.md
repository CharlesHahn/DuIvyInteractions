# 基团原子存储约定

> 状态：定义定稿

---

## 1. 总则

每个 Group 的 `atoms` 字段是一个 `List[AtomData]`，存储基团内的所有原子。不同基团类型的原子数量和排列顺序有不同的约定。

---

## 2. 各基团类型的原子存储

### 2.1 H_donor

`atoms = [D, H]`（2 个原子）

| 位置 | 语义 | 元素 |
|------|------|------|
| atoms[0] | D（供体原子） | N、O、S 或 F |
| atoms[1] | H（氢原子） | H |

**metadata**：`{}`

---

### 2.2 H_acceptor

`atoms = [A]`（1 个原子）

| 位置 | 语义 |
|------|------|
| atoms[0] | A（受体原子） |

**metadata**：`{}`

---

### 2.3 aromatic_ring

`atoms = [环内所有原子]`（5 或 6 个原子）

原子按残基内索引顺序排列，无特定语义顺序。

**metadata**：`{}`

---

### 2.4 charged_positive / charged_negative

`atoms = [电荷中心的所有原子]`（1~13 个原子不等）

原子按残基内索引顺序排列。不同残基的原子数量不同：

| 残基 | 原子数 | 代表原子 |
|------|-------:|:---------|
| ARG | 9 | CZ, NE, NH1, NH2, HE, HH11, HH12, HH21, HH22 |
| LYS | 4 | NZ, HZ1, HZ2, HZ3 |
| ASP | 3 | CG, OD1, OD2 |
| GLU | 3 | CD, OE1, OE2 |

**metadata**：

| 键 | 值 | 说明 |
|----|----|------|
| source | `"residue_name"` 或 `"functional_group"` | 识别来源 |
| func_group | `"tertamine"`, `"carboxylate"`, ... 或 `None` | 官能团类型（第二层填写） |

---

### 2.5 halogen_donor

`atoms = [C, X]`（2 个原子）

| 位置 | 语义 | 元素 |
|------|------|------|
| atoms[0] | C（连接卤素的碳） | C |
| atoms[1] | X（卤素原子） | F、Cl、Br 或 I |

C→X 方向定义了 σ-hole 的方向（卤键沿 C-X 键轴方向）。

**metadata**：`{}`

---

### 2.6 halogen_acceptor

`atoms = [受体原子]`（1 个原子）

| 位置 | 语义 |
|------|------|
| atoms[0] | 受体原子（C、P 或 S，连接到 O/P/N/S） |

**metadata**：`{}`

---

### 2.7 metal

`atoms = [金属离子]`（1 个原子）

| 位置 | 语义 |
|------|------|
| atoms[0] | 金属离子 |

**metadata**：`{}`

---

### 2.8 water

`atoms = [O, H1, H2]`（3 个原子）

| 位置 | 语义 | 元素 |
|------|------|------|
| atoms[0] | O（氧原子） | O |
| atoms[1] | HW1（氢 1） | H |
| atoms[2] | HW2（氢 2） | H |

**metadata**：`{}`

---

### 2.9 hydrophobic

`atoms = [疏水碳原子]`（1 个原子）

| 位置 | 语义 |
|------|------|
| atoms[0] | 疏水碳原子（C，所有邻居 ∈ {C, H}） |

**metadata**：`{}`

---

## 3. 汇总表

| 基团类型 | atoms 数量 | atoms[0] | atoms[1] | atoms[2:] | metadata |
|:---------|----------:|:---------|:---------|:----------|:---------|
| H_donor | 2 | D（供体） | H（氢） | — | 空 |
| H_acceptor | 1 | A（受体） | — | — | 空 |
| aromatic_ring | 5-6 | 环原子 | 环原子 | ... | 空 |
| charged_positive | 1-13 | 电荷中心原子 | ... | ... | source, func_group |
| charged_negative | 1-13 | 电荷中心原子 | ... | ... | source, func_group |
| halogen_donor | 2 | C（碳） | X（卤素） | — | 空 |
| halogen_acceptor | 1 | 受体原子 | — | — | 空 |
| metal | 1 | 金属离子 | — | — | 空 |
| water | 3 | O（氧） | H1 | H2 | 空 |
| hydrophobic | 1 | 疏水碳 | — | — | 空 |

---

*文档结束*
