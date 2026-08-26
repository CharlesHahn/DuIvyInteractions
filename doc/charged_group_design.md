# 带电基团（盐桥）识别设计

> 创建日期：2026-08-27
> 状态：设计定稿

---

## 1. 盐桥的文献定义

### 1.1 核心定义

> **"A salt bridge can be defined as an interaction between two groups of opposite charge in which at least one pair of heavy atoms is within hydrogen bonding distance."**
> — Donald et al. 2011, Proteins, PMC3069487

盐桥 = **相反电荷基团** + **重原子距离 < 阈值**（PLIP 取 5.5 Å，源自 Barlow & Thornton 1983 + 1.5 Å 扩展）。

盐桥不是纯粹的静电作用，而是**静电吸引 + 氢键的组合**：

> **"The salt bridge is a non-covalent interaction that combines an electrostatic attraction between oppositely charged chemical groups or atoms and a hydrogen bond."**
> — Spassov et al. 2023, Front Mol Biosci, PMC9871453

### 1.2 盐桥不限于氨基酸之间

三个层次的证据：

| 场景 | 证据 | 来源 |
|------|------|------|
| 蛋白内部 | ASP/GLU 与 LYS/ARG/HIS 之间 | Donald 2011, Kumar & Nussinov 2002 |
| 蛋白-配体 | PDB 中 1100+ 个配体-蛋白复合物含盐桥 | Kurczab et al. 2018, J Chem Inf Model |
| 超分子 | 离子对在超分子化学中广泛存在 | Wikipedia, Schneider 2009 |

**结论**：识别算法必须同时覆盖蛋白残基和非蛋白分子（配体、辅因子等）。

---

## 2. 设计策略：三层递进

```
第一层：残基名字典（蛋白残基）
   → 已知化学事实，零推断，100% 确定

第二层：官能团模式匹配（非蛋白残基）
   → 原子类型 + 键连接 → 识别带电官能团

第三层：部分电荷交叉验证（通用）
   → 识别出的官能团的净电荷是否符合预期
```

### 2.1 为什么不用部分电荷阈值？

| 方法 | 问题 |
|------|------|
| `\|q\| > 0.3` 逐原子阈值 | 魔法数字，无文献依据；逐原子丢失官能团整体语义 |
| 残基名 + 官能团模式 + 电荷验证 | 化学定义驱动；电荷仅作交叉验证，非主判据 |

### 2.2 信息论依据

tpr 文件中的部分电荷是力场参数化时的**化学判决的留存记录**。当 antechamber/sobtop 赋予某个 N 原子正的部分电荷时，已经隐含了"该 N 是质子化的"这一判断。我们直接利用这个已编码的信息，而非从头推断。

---

## 3. 第一层：蛋白残基名 → 带电基团

### 3.1 正电残基

| 残基 | 电荷 | 电荷中心原子 | 说明 |
|------|------|-------------|------|
| **ARG** | +1 | CZ, NE, NH1, NH2 | 胍基，永久正电 |
| **LYS** | +1 | NZ | 氨基，永久正电 |
| **HIP** | +1 | ND1, NE2 | 双质子化咪唑，正电 |

**HIS 质子化态说明**（Amber 力场）：

| 残基名 | 质子化态 | 电荷 | 能否形成盐桥 |
|--------|---------|------|-------------|
| **HIP** | 双质子化（ND1-H + NE2-H） | +1 | ✅ 可以 |
| **HID** | δ-质子化（ND1-H，NE2 中性） | 0 | ❌ 不能 |
| **HIE** | ε-质子化（NE2-H，ND1 中性） | 0 | ❌ 不能 |

tpr 文件通过残基名（HIP/HID/HIE）直接编码质子化态，无需推断。这比 PLIP（用 reduce 程序加 H 判断）更精确。

### 3.2 负电残基

| 残基 | 电荷 | 电荷中心原子 | 说明 |
|------|------|-------------|------|
| **ASP** | -1 | CG, OD1, OD2 | 羧基，永久负电 |
| **GLU** | -1 | CD, OE1, OE2 | 羧基，永久负电 |

### 3.3 不纳入的残基

| 残基 | 原因 |
|------|------|
| GLN / ASN | 酰胺基，中性（PLIP 列为 polar 但非 charged） |
| TYR | 酚羟基 pKa ~10，生理 pH 下中性 |
| HID / HIE | 单质子化 His，中性 |

### 3.4 N/C 末端

蛋白链的 N-末端（NH3+）和 C-末端（COO-）在 Amber 力场中：
- N-末端：残基名仍为标准名（如 ALA），但氨基已被质子化 → 会被**第二层官能团模式**（叔胺/季铵）自动捕获
- C-末端：同理，羧基去质子化 → 会被**第二层官能团模式**（羧酸盐）自动捕获

无需为末端做特殊处理。

---

## 4. 第二层：官能团模式匹配（非蛋白残基）

适用于配体、辅因子、非标准残基等无残基名字典的分子。

判据：**原子类型 + 键连接图（邻居原子的元素和数量）**。

### 4.1 正电官能团

| 官能团 | 模式 | 涉及原子 | 来源 |
|--------|------|---------|------|
| **季铵** `quartamine` | N 连 4 个非 H 邻居 | N + 4 邻居 | PLIP |
| **叔胺** `tertamine` | sp3 N 连 ≥3 个非 H 邻居 | N + 邻居 | PLIP |
| **胍基** `guanidine` | C 连 3 个 N，且至少一个 N 只连该 C | C + 3N | PLIP |
| **锍** `sulfonium` | S 连 3 个非 H 邻居 | S + 3 邻居 | PLIP |

### 4.2 负电官能团

| 官能团 | 模式 | 涉及原子 | 来源 |
|--------|------|---------|------|
| **羧酸盐** `carboxylate` | C 连 2 个 O + 1 个 C | C + 2O | PLIP |
| **磷酸盐** `phosphate` | P 的邻居全是 O | P + 所有 O | PLIP |
| **磺酸** `sulfonicacid` | S 连 3 个 O | S + 3O | PLIP |
| **硫酸盐** `sulfate` | S 连 4 个 O | S + 4O | PLIP |

### 4.3 判定函数参考（源自 PLIP `is_functional_group`）

```python
# 季铵：N 有 4 个非 H 邻居
def is_quartamine(atom, neighbors):
    return atom.element == 'N' and len(neighbors) == 4 and 'H' not in [n.element for n in neighbors]

# 叔胺：sp3 N 有 ≥3 邻居
def is_tertamine(atom, neighbors):
    return atom.element == 'N' and len(neighbors) >= 3

# 羧酸盐：C 连 2 个 O + 1 个 C
def is_carboxylate(atom, neighbors):
    if atom.element != 'C':
        return False
    o_count = sum(1 for n in neighbors if n.element == 'O')
    c_count = sum(1 for n in neighbors if n.element == 'C')
    return o_count == 2 and c_count >= 1

# 磷酸盐：P 的邻居全是 O
def is_phosphate(atom, neighbors):
    return atom.element == 'P' and all(n.element == 'O' for n in neighbors)

# 胍基：C 连 3 个 N，至少一个 N 只连 C
def is_guanidine(atom, neighbors):
    if atom.element != 'C' or len(neighbors) != 3:
        return False
    if not all(n.element == 'N' for n in neighbors):
        return False
    # 至少一个 N 只连该 C（可质子化）
    return any(len(n_neighbors) == 1 for n_neighbors in neighbor_neighbor_counts)

# 锍：S 连 3 个非 H 邻居
def is_sulfonium(atom, neighbors):
    return atom.element == 'S' and len(neighbors) == 3

# 磺酸：S 连 3 个 O
def is_sulfonicacid(atom, neighbors):
    return atom.element == 'S' and len(neighbors) == 3 and all(n.element == 'O' for n in neighbors)

# 硫酸盐：S 连 4 个 O
def is_sulfate(atom, neighbors):
    return atom.element == 'S' and len(neighbors) == 4 and all(n.element == 'O' for n in neighbors)
```

---

## 5. 第三层：部分电荷交叉验证

对第二层识别出的官能团，计算其**净部分电荷**进行验证：

| 官能团 | 预期净电荷 | 验证条件 |
|--------|-----------|---------|
| 季铵 / 叔胺 / 胍基 / 锍 | 正 | $\sum q_i > 0$ |
| 羧酸盐 / 磷酸盐 / 磺酸 / 硫酸盐 | 负 | $\sum q_i < 0$ |

**第三层不是主判据，而是交叉验证**。化学模式匹配是主判据，电荷验证是辅判据。

若化学模式匹配但电荷验证不通过（如一个羧酸盐的净电荷 ≈ 0），说明力场参数化时该基团可能未被去质子化，应**以力场为准**（不标记为带电基团）。

---

## 6. 电荷中心计算

### 6.1 公式

带电基团的电荷中心 = **带符号部分电荷加权的位置平均**：

$$\vec{r}_{center} = \frac{\sum_{i \in group} q_i \cdot \vec{r}_i}{\sum_{i \in group} q_i}$$

其中：
- $q_i$：原子 $i$ 的**带符号**部分电荷（不取绝对值）
- $\vec{r}_i$：原子 $i$ 的坐标
- $\sum q_i$：基团的净电荷（不为零，因为是带电基团）

### 6.2 物理意义

- 负电基团（$\sum q_i < 0$）：负电荷原子权重更大，中心偏向负电荷集中处
- 正电基团（$\sum q_i > 0$）：正电荷原子权重更大，中心偏向正电荷集中处
- 这比几何质心更准确地反映了静电作用的真实中心

### 6.3 为什么不用绝对值

以 ASP 羧基为例（OD1 和 OD2 对称分布于 CG 两侧）：

| 方法 | 计算 | 结果 |
|------|------|------|
| 带符号 q | $\frac{(+0.6)\vec{r}_{CG} + (-0.8)\vec{r}_{OD1} + (-0.8)\vec{r}_{OD2}}{-1.0}$ | 中心偏向 O 侧（负电荷主导） |
| 绝对值 \|q\| | $\frac{(0.6)\vec{r}_{CG} + (0.8)\vec{r}_{OD1} + (0.8)\vec{r}_{OD2}}{2.2}$ | 混合正负电荷，物理意义错误 |

绝对值把正电荷（CG 的 +0.6）和负电荷（OD 的 -0.8）混在一起加权，计算的是"电荷量的重心"而非"电荷的中心"。

---

## 7. 输出格式

### 7.1 Group 数据结构

```python
Group(
    group_id=gid,
    group_type="charged_positive",  # 或 "charged_negative"
    molecule=res.molecule_name,
    residue_name=res.residue_name,
    residue_id=res.residue_global_idx,
    atom_indices=[...],             # 官能团内所有原子的全局索引
    atom_types=[...],
    elements=[...],
    charges=[...],                  # 各原子的部分电荷
    center=(x, y, z),              # 电荷加权中心坐标
    metadata={
        "net_charge": -0.85,        # 净电荷
        "source": "residue_name",   # "residue_name" 或 "functional_group"
        "func_group": "carboxylate" # 官能团类型（第二层匹配时填写）
    }
)
```

### 7.2 每个带电残基/官能团生成一个 Group

- ASP → 一个 `charged_negative` Group，`atom_indices=[CG, OD1, OD2]`
- ARG → 一个 `charged_positive` Group，`atom_indices=[CZ, NE, NH1, NH2]`
- 配体羧酸盐 → 一个 `charged_negative` Group，`atom_indices=[C, O1, O2]`

---

## 8. 与 PLIP 的覆盖对比

### 8.1 蛋白侧

| PLIP 蛋白正电 | 我们的设计 | 差异 |
|:----|:----|:----|
| ARG | ✅ 第一层 | 无 |
| LYS | ✅ 第一层 | 无 |
| HIS（所有质子化态） | ✅ 仅 HIP | **我们更精确**（tpr 直接编码质子化态） |

| PLIP 蛋白负电 | 我们的设计 | 差异 |
|:----|:----|:----|
| ASP | ✅ 第一层 | 无 |
| GLU | ✅ 第一层 | 无 |

### 8.2 配体侧

| PLIP 配体正电 | 我们的设计 | 差异 |
|:----|:----|:----|
| `quartamine` | ✅ 第二层 | 无 |
| `tertamine` | ✅ 第二层 | 无 |
| `guanidine` | ✅ 第二层 | 无 |
| `sulfonium` | ✅ 第二层 | 无 |

| PLIP 配体负电 | 我们的设计 | 差异 |
|:----|:----|:----|
| `carboxylate` | ✅ 第二层 | 无 |
| `phosphate` | ✅ 第二层 | 无 |
| `sulfonicacid` | ✅ 第二层 | 无 |
| `sulfate` | ✅ 第二层 | 无 |

### 8.3 我们的优势

| 方面 | PLIP | 我们 |
|------|------|------|
| 输入 | PDB（无电荷、无力场） | tpr（有电荷、有力场类型） |
| His 质子化 | reduce 推断（可能不准） | 残基名直接编码（HIP/HID/HIE） |
| 配体电荷 | 无（靠化学模式推断） | 有（第三层交叉验证） |
| 电荷中心 | 质心或单原子 | 电荷加权中心 |

---

## 9. 实现路径

### 9.1 函数结构

```
_find_charged(res, bond_graph, start_id)
├── _identify_protein_charged(res)          # 第一层：残基名字典
├── _identify_functional_group_charged(res, bond_graph)  # 第二层：官能团模式
└── _verify_and_build_group(atoms, expected_sign, res, gid)  # 第三层：电荷验证 + 构建 Group

_charge_weighted_center(atoms, positions)   # 电荷加权中心
```

### 9.2 改动范围

| 函数 | 改动 |
|------|------|
| `_find_charged` | **重写**：删除 `\|q\|>0.3` 逻辑，改为三层递进 |
| 新增 `_identify_protein_charged` | 第一层残基名字典 |
| 新增 `_identify_functional_group_charged` | 第二层官能团模式匹配 |
| 新增 `_verify_and_build_group` | 第三层电荷验证 + Group 构建 |
| 新增 `_charge_weighted_center` | 电荷加权中心计算 |
| `_build_bond_graph` | **复用**（已有的全局连接图） |

### 9.3 常量定义

```python
# 蛋白带电残基（第一层）
POSITIVE_RESIDUES = {
    "ARG": {"center_atoms": ["CZ", "NE", "NH1", "NH2"], "charge": +1},
    "LYS": {"center_atoms": ["NZ"], "charge": +1},
    "HIP": {"center_atoms": ["ND1", "NE2"], "charge": +1},
}
NEGATIVE_RESIDUES = {
    "ASP": {"center_atoms": ["CG", "OD1", "OD2"], "charge": -1},
    "GLU": {"center_atoms": ["CD", "OE1", "OE2"], "charge": -1},
}
```

---

## 10. 参考文献

1. Donald JE, Kulp DW, DeGrado WF. Salt Bridges: Geometrically Specific, Designable Interactions. Proteins. 2011;79(3):898-915. PMC3069487.
2. Barlow DJ, Thornton JM. Ion-pairs in proteins. J Mol Biol. 1983;168(4):867-885.
3. Kumar S, Nussinov R. Close-range electrostatic interactions in proteins. ChemBioChem. 2002;3(7):604-617.
4. Kurczab R, et al. Salt Bridge in Ligand-Protein Complexes—Systematic Theoretical and Statistical Investigations. J Chem Inf Model. 2018;58(11):2224-2238.
5. Spassov DS, et al. A role of salt bridges in mediating drug potency. Front Mol Biosci. 2023;9:1066029. PMC9871453.
6. PLIP 源码: https://github.com/pharmai/plip — `structure/preparation.py` (`is_functional_group`, `append_if_charged_func_group`)
7. PLIP 文档: https://github.com/pharmai/plip/blob/master/DOCUMENTATION.md — "Charged Groups" 章节

---

*文档结束*
