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

| 层级 | 适用范围 | 方法 | 确定性 |
|------|---------|------|--------|
| 第一层 | 蛋白残基 | 残基名字典 | 100%（零推断） |
| 第二层 | 非蛋白残基 | 官能团模式匹配（元素+邻居） | 高 |
| 第三层 | 通用 | 部分电荷交叉验证（\|Σq\| > 0.1） | 辅助过滤 |

### 2.1 为什么不用部分电荷阈值？

| 方法 | 问题 |
|------|------|
| `\|q\| > 0.3` 逐原子阈值 | 魔法数字，无文献依据；逐原子丢失官能团整体语义 |
| 残基名 + 官能团模式 + 电荷验证 | 化学定义驱动；电荷仅作交叉验证，非主判据 |

### 2.2 信息论依据

tpr 文件中的部分电荷是力场参数化时的**化学判决的留存记录**。当 antechamber/sobtop 赋予某个 N 原子正的部分电荷时，已经隐含了"该 N 是质子化的"这一判断。我们直接利用这个已编码的信息，而非从头推断。

---

## 3. 第一层：蛋白残基名 → 带电基团

### 3.1 正电残基（+1）

| 残基名 | 电荷 | 电荷中心原子 | 说明 |
|--------|------|-------------|------|
| **ARG** | +1 | CZ, NE, NH1, NH2, HE, HH11, HH12, HH21, HH22 | 胍基，永久正电 |
| **LYS** | +1 | NZ, HZ1, HZ2, HZ3 | 氨基，永久正电 |
| **HIP** | +1 | ND1, NE2, HD1, HE2 | 双质子化咪唑，正电 |
| **ORN** | +1 | NE, HE1, HE2, HE3 | 鸟氨酸，伯胺 |
| **DAB** | +1 | ND, HD1, HD2, HD3 | 二氨基丁酸，伯胺 |
| **M3L** | +1 | NZ, CM1, CM2, CM3, HM11-33 | 三甲基化 Lys（季铵） |
| **MLY** | +1 | NZ, CH1, CH2, HH11-23 | 二甲基化 Lys（叔胺） |

### 3.2 负电残基（-1）

| 残基名 | 电荷 | 电荷中心原子 | 说明 |
|--------|------|-------------|------|
| **ASP** | -1 | CG, OD1, OD2 | 羧基，永久负电 |
| **GLU** | -1 | CD, OE1, OE2 | 羧基，永久负电 |
| **CYM** | -1 | SG | 去质子化 Cys（硫醇盐） |
| **KCX** | -1 | NZ, CX, OQ1, OQ2, HZ | 羧基化 Lys（氨基甲酸） |
| **PCA** | -1 | CA, C, O, N, CD, CG | 焦谷氨酸 |

### 3.3 磷酸化残基（-2 或 -1）

| 残基名 | 电荷 | 电荷中心原子 | 说明 |
|--------|------|-------------|------|
| **SEP** | -2/-1 | OG, P, O1P, O2P, O3P | 磷酸化 Ser |
| **TPO** | -2/-1 | OG1, P, O1P, O2P, O3P | 磷酸化 Thr |
| **PTR** | -2/-1 | OH, P, O1P, O2P, O3P | 磷酸化 Tyr |

磷酸化残基的电荷取决于质子化态：
- 完全去质子化：-2（两个 H 被移除）
- 部分质子化：-1（一个 H 保留）

### 3.4 不纳入的残基

| 残基 | 原因 |
|------|------|
| GLN / ASN | 酰胺基，中性 |
| TYR | 酚羟基 pKa ~10，生理 pH 下中性 |
| HID / HIE | 单质子化 His，中性 |
| LYN | 去质子化 Lys，中性 |
| ASH | 质子化 ASP，中性 |
| GLH | 质子化 GLU，中性 |
| CYS | 标准 Cys，中性 |

### 3.5 HIS 质子化态说明（Amber 力场）

| 残基名 | 质子化态 | 电荷 | 能否形成盐桥 |
|--------|---------|------|-------------|
| **HIP** | 双质子化（ND1-H + NE2-H） | +1 | ✅ 可以 |
| **HID** | δ-质子化（ND1-H，NE2 中性） | 0 | ❌ 不能 |
| **HIE** | ε-质子化（NE2-H，ND1 中性） | 0 | ❌ 不能 |

tpr 文件通过残基名（HIP/HID/HIE）直接编码质子化态，无需推断。

### 3.6 N/C 末端

蛋白链的 N-末端（NH3+）和 C-末端（COO-）在 Amber 力场中：
- N-末端：残基名仍为标准名（如 ALA），但氨基已被质子化 → 会被**第二层官能团模式**（叔胺/季铵）自动捕获
- C-末端：同理，羧基去质子化 → 会被**第二层官能团模式**（羧酸盐）自动捕获

无需为末端做特殊处理。

---

## 4. 第二层：官能团模式匹配（非蛋白残基）

适用于配配体、辅因子、非标准残基等无残基名字典的分子。

判据：**元素 + 邻居元素 + 邻居数量**（严格参照 PLIP `is_functional_group`）。

### 4.1 正电官能团

#### quartamine（季铵）

**定义**：N 有 4 个邻居，且**没有一个邻居是 H**。

| 条件 | PLIP | 我们 |
|------|------|------|
| 元素 | `atom.atomicnum == 7` | `atom.atom_element == 'N'` |
| 邻居数 | `len(n_atoms) == 4` | `len(neighbors) == 4` |
| 无 H 邻居 | `'1' not in n_atoms` | `all(n.atom_element != 'H' for n in neighbors)` |

**物理意义**：NR4+（季铵盐），永远带正电。

---

#### tertamine（叔胺 / 质子化胺）

**定义**：N 有 **≥3 个邻居（包括 H）**。

| 条件 | PLIP | 我们 |
|------|------|------|
| 元素 | `atom.atomicnum == 7` | `atom.atom_element == 'N'` |
| sp3 杂化 | `GetHyb() == 3` | 不检查（无此信息） |
| 邻居数 | `len(n_atoms) >= 3` | `len(neighbors) >= 3` |

**物理意义**：sp3 N 有 ≥3 个邻居 → NH3、NH4+、NR3、NR3H+ 等。包含 H 邻居，因为 NH4+（4 个 H）是正电。需要第三层电荷验证过滤中性胺。

**PLIP 设计意图**：在晶体结构中，质子化态不确定，所以标记所有 sp3 N ≥3 邻居为候选。在 MD 中，力场已确定质子化态，电荷验证会过滤中性胺。

---

#### guanidine（胍基）

**定义**：C 有 3 个 N 邻居，且**至少一个 N 只连了该 C**（可质子化）。

| 条件 | PLIP | 我们 |
|------|------|------|
| 元素 | `atom.atomicnum == 6` | `atom.atom_element == 'C'` |
| N 邻居数 | `n_atoms.count(7) == 3` | N 邻居数 = 3 |
| 总邻居数 | `len(n_atoms) == 3` | 总邻居数 = 3 |
| 可质子化 N | `min(nitro_partners) == 1` | 至少一个 N 只连该 C |

**物理意义**：胍基（如 ARG 侧链）的中心 C 连 3 个 N，其中一个 N 是终端 NH2，可质子化形成正电。

---

#### sulfonium（锍）

**定义**：S 有 3 个邻居，且**没有一个邻居是 H**。

| 条件 | PLIP | 我们 |
|------|------|------|
| 元素 | `atom.atomicnum == 16` | `atom.atom_element == 'S'` |
| 邻居数 | `len(n_atoms) == 3` | `len(neighbors) == 3` |
| 无 H 邻居 | `'1' not in n_atoms` | `all(n.atom_element != 'H' for n in neighbors)` |

**物理意义**：R3S+（锍盐），永远带正电。排除亚砜（R2S=O）和巯基（R-SH）。

---

### 4.2 负电官能团

#### phosphate（磷酸盐）

**定义**：P 的**邻居全是 O**。

| 条件 | PLIP | 我们 |
|------|------|------|
| 元素 | `atom.atomicnum == 15` | `atom.atom_element == 'P'` |
| 邻居 | `set(n_atoms) == {8}` | 所有邻居都是 O |

**物理意义**：PO4^3-、PO3^2- 等磷酸基，永远带负电。

---

#### sulfonicacid（磺酸）

**定义**：S 有 **3 个 O 邻居**。

| 条件 | PLIP | 我们 |
|------|------|------|
| 元素 | `atom.atomicnum == 16` | `atom.atom_element == 'S'` |
| O 邻居数 | `n_atoms.count(8) == 3` | O 邻居数 = 3 |

**物理意义**：R-SO3-（磺酸基），带负电。

---

#### sulfate（硫酸盐）

**定义**：S 有 **4 个 O 邻居**。

| 条件 | PLIP | 我们 |
|------|------|------|
| 元素 | `atom.atomicnum == 16` | `atom.atom_element == 'S'` |
| O 邻居数 | `n_atoms.count(8) == 4` | O 邻居数 = 4 |

**物理意义**：R-SO4^2-（硫酸基），带负电。

---

#### carboxylate（羧酸盐）

**定义**：C 有 **2 个 O 邻居**和**恰好 1 个 C 邻居**。

| 条件 | PLIP | 我们 |
|------|------|------|
| 元素 | `atom.atomicnum == 6` | `atom.atom_element == 'C'` |
| O 邻居数 | `n_atoms.count(8) == 2` | O 邻居数 = 2 |
| C 邻居数 | `n_atoms.count(6) == 1` | C 邻居数 = **恰好 1** |

**物理意义**：R-COO-（羧酸盐），带负电。

**关键**：PLIP 要求**恰好 1 个 C**，不是 ≥1 个。这排除了酯（R-COOR'，有 2 个 C 邻居）和酰胺（R-CONH2，有 2 个 C 邻居）。

---

---

## 5. 第三层：部分电荷交叉验证

对第二层识别出的官能团，计算其**净部分电荷**进行验证：

| 官能团 | 预期净电荷 | 验证条件 |
|--------|-----------|---------|
| quartamine | +1 | $\sum q_i > 0.1$ |
| tertamine | +1 或 0 | $\sum q_i > 0.1$（过滤中性胺） |
| guanidine | +1 | $\sum q_i > 0.1$ |
| sulfonium | +1 | $\sum q_i > 0.1$ |
| phosphate | -1/-2/-3 | $\sum q_i < -0.1$ |
| sulfonicacid | -1 | $\sum q_i < -0.1$ |
| sulfate | -2 | $\sum q_i < -0.1$ |
| carboxylate | -1 | $\sum q_i < -0.1$ |

**阈值**：$|\sum q_i| > 0.1$（排除数值误差）

**设计原则**：第二层是"候选生成"（可能有假阳性），第三层是"电荷验证"（过滤假阳性）。化学模式匹配是主判据，电荷验证是辅判据。

若化学模式匹配但电荷验证不通过，说明力场参数化时该基团可能未被去质子化，应**以力场为准**（不标记为带电基团）。

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

### 7.1 Group 字段

每个带电基团生成一个 Group，包含以下信息：

| 字段 | 说明 |
|------|------|
| group_type | `"charged_positive"` 或 `"charged_negative"` |
| atoms | `List[AtomData]`，基团内的所有原子（单一数据源，无冗余） |
| metadata.source | 来源：`"residue_name"`（第一层）或 `"functional_group"`（第二层） |
| metadata.func_group | 官能团类型（第二层匹配时填写，如 `"carboxylate"`） |

原子属性（索引、类型、元素、电荷）通过 `atoms` 列表中的 `AtomData` 对象访问，无需存储冗余副本。

### 7.2 每个带电残基/官能团生成一个 Group

- ASP → 一个 `charged_negative` Group，atoms 包含 CG, OD1, OD2
- ARG → 一个 `charged_positive` Group，atoms 包含 CZ, NE, NH1, NH2
- 配体羧酸盐 → 一个 `charged_negative` Group，atoms 包含 C, O1, O2

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

**我们额外支持的带电残基**（PLIP 不支持或需推断）：

| 残基 | 电荷 | PLIP 支持 | 我们 |
|------|------|----------|------|
| CYM | -1 | ❌ | ✅ 第一层 |
| SEP/TPO/PTR | -2/-1 | ❌ | ✅ 第一层 |
| M3L/MLY | +1 | ❌ | ✅ 第一层 |
| KCX | -1 | ❌ | ✅ 第一层 |
| PCA | -1 | ❌ | ✅ 第一层 |

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

### 9.1 架构概览

带电基团识别采用三层递进架构：

- **第一层**（残基名字典）：根据残基名直接查找已知带电残基，零推断
- **第二层**（官能团模式匹配）：遍历残基内原子，用元素+邻居模式识别带电官能团
- **第三层**（电荷验证）：对第二层结果计算净电荷，过滤假阳性

第一层和第二层结果合并后进行去重。

### 9.2 模块职责

| 模块 | 职责 |
|------|------|
| 带电基团主函数 | 三层递进调用 + 去重 |
| 第一层模块 | 残基名字典查找 |
| 第二层模块 | 官能团模式匹配 |
| 第三层模块（嵌入第二层） | 电荷验证 + Group 构建 |
| 键连接图 | 复用已有的全局连接图 |

### 9.3 去重策略

第一层和第二层可能对同一基团产生重复识别。例如：
- ARG 的胍基：第一层（残基名字典）和第二层（guanidine 模式）都会识别
- ASP 的羧基：第一层（残基名字典）和第二层（carboxylate 模式）都会识别
- N/C 末端：第二层会识别（tertamine/carboxylate），但第一层不会（残基名是标准名）

**去重规则**：

1. **去重 key**：`frozenset(atom_indices)`（原子全局索引集合）
2. **优先级**：`residue_name` 来源 > `functional_group` 来源
3. **逻辑**：
   - 若同一原子集合只被一层识别 → 直接保留
   - 若同一原子集合被两层识别 → 保留 `residue_name` 来源（第一层更确定）

**示例**：

- ARG 第一层识别：包含 CZ, NE, NH1, NH2 及所有 H 原子，来源为 residue_name
- ARG 第二层识别（guanidine 模式）：仅包含 CZ, NE, NH1, NH2（不含 H），来源为 functional_group
- 两个 Group 的原子集合不同，都保留

### 9.4 常量定义

常量定义见 §3.1（正电残基）、§3.2（负电残基）、§3.3（磷酸化残基）。每个残基名映射到其带电基团的原子名列表。

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
