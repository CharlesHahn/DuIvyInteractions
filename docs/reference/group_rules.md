# 基团识别规则

基团鉴定阶段（两段式架构的①，见[核心概念](concepts.md)）从 tpr 的原子类型、键合图、电荷和显式氢原子识别出可参与相互作用的化学基团，**只做一次、与帧无关**。本文档列出各基团的判定规则及 4 个力场识别器（Amber / GROMOS / CHARMM / OPLS）的实现差异，供理解工具行为与验证结果。各力场类型特征表的完整内容见[力场类型映射](force_field.md)。

## 基团类型全集（GROUP_TYPES）

工具可识别的基团类型（`core/constants.py` 的 `GROUP_TYPES`）：

| 基团类型 | 参与相互作用 |
|:---------|:------------|
| `H_donor` | 氢键、水桥 |
| `H_acceptor` | 氢键、水桥 |
| `aromatic_ring` | π-π 堆积、π-阳离子、卤键（π 受体） |
| `charged_positive` | 盐桥、π-阳离子 |
| `charged_negative` | 盐桥 |
| `halogen_donor` | 卤键 |
| `halogen_acceptor` | 卤键 |
| `metal` | 金属配位 |
| `metal_binding` | 金属配位 |
| `water` | 水桥 |
| `hydrophobic` | 疏水 |

## 判定依据总览

| 基团 | tpr 中的证据 | 确定性 |
|:-----|:-------------|:-------|
| 芳香环 | 类型 ∈ `STRONG_AROMATIC`/兼容类型 → 图论环检测（BFS）→ 两条件过滤（GROMOS 分级判定） | 高 |
| H 键供体 | D–H 键条目显式存在（Bond/Constraint/Settle 合并）＋ q(H)>0 | 100%（零推断） |
| H 键受体 | 类型表 ＋ 结构判据 ＋ q<0（排除正电成员/环内带 H 氮） | 高 |
| 带电基团 | 残基名字典/官能团模式 ＋ 净电荷验证（阈值 ±0.1） | 蛋白高 |
| 疏水 | C 且邻居 ∈ {C,H}（GROMOS：类型白名单 ＋ 无极性邻居） | 高 |
| 卤键 | 卤素元素/类型 ＋ 邻接碳（区分芳香/烷基卤） | 高 |
| 金属 | 元素 ∈ `METAL_IONS`（由读取器推断） | 100% |
| 水 | 残基名 ∈ `WATER_RESIDUES` ＋ OW/HW 原子名 | 100% |

## 识别流程

`identify()` 按固定顺序逐残基识别（基类 `AmberFFGroupIdentifier._identify_residue`；GROMOS/CHARMM/OPLS 为子类，覆写力场相关部分），跨残基键单独补检供体：

1. 芳香环（`_find_aromatic_rings`）→ 2. 供体（`_find_donors`）→ 3. 带电基团（`_find_charged`，提前识别以便为受体排除正电成员）→ 4. 受体（`_find_acceptors`，带排除集合）→ 5. 卤键供体 → 6. 卤键受体 → 7. 金属 → 8. 水 → 9. 疏水 → 10. 金属配位原子

**受体排除集合**（供受体判定用，`_identify_residue` 中构建）：① 正电基团成员（无孤对，不可作受体）；② 芳香环内带 H 的氮（吡咯型 NH，孤对参与芳香共轭，不可作受体）。

## 逐基团判定规则

### 芳香环（aromatic_ring）

**目的**：为 π-π / π-阳离子 / 卤键（π 受体）提供环几何。

- **环检测**：图论方法，遍历每条边做 BFS 找环（无大小限制，`_detect_rings`），随后按环大小排序去重——被已接受小环完全覆盖的大环剔除（`_deduplicate_aromatic_rings`）
- **芳香过滤（Amber / CHARMM / OPLS，两条件）**：环内 **≥ n−1 个原子**属 `STRONG_AROMATIC`，剩余原子必须全部属 `COMPATIBLE_TYPES`
- **GROMOS 分级判定**（类型粒度粗，`_filter_aromatic_rings` 覆写）：环成分必须 ⊆ `{C, CR1, NR}`；含 `NR`（芳香氮）直接判芳香；全 C 环回退残基名白名单（PHE/TYR/TRP/His 变体，见[力场类型映射](force_field.md)）
- 环内原子按环顺序存储（BFS 路径顺序），π-π 检测器据此用相邻原子叉积计算环法向量

**注意**：芳香判定**只用类型名与键合图，不依赖键序**（即使 tpr 不保留键序也可靠）；平面性不是基团鉴定的判据，而是 π-π 检测器的可选几何检查（[相互作用判据](criteria.md)）。

### H 键供体（H_donor）

**判据** = D-H 键条目存在（D ∈ {N, O, S, F}）+ q(H) > 0（`_classify_dh_pair`）。因 N-H 键位于 tpr 的 **Constraint 段**，读取器把 Bond / Constraint / Settle 合并为统一键列表，**必须依赖合并后的键图**（只扫 Bond 段会漏检）。跨残基键的 D-H 对单独补检（`_find_inter_residue_donors`）。判定不涉及坐标推断。

### H 键受体（H_acceptor）

**判据** = 类型 ∈ 受体类型表 + **部分电荷 q < 0** + 不在排除集合（正电成员、环内带 H 氮）。

**不使用"H 计数"**（从未实现）；带 H 环境的排除由各力场的类型表或结构判据完成：

| 力场 | 受体判定实现 |
|:-----|:-------------|
| Amber | 类型表 `ACCEPTOR_TYPES` ＋ 二义类型 `N` 的结构判据（带 H 邻居的普通酰胺 → 排除；无 H 的 Pro N → 保留；`N2`/`NB`/`NC` 由类型表直接保留）＋ 排除集合 |
| GROMOS | 类型表 `GROMOS_ACCEPTOR_TYPES` ＋ N 结构判据（`_find_acceptors` 覆写，对元素 N：**键数 ≥ 4（铵 RNH₃⁺）→ 排除**；**带 H 邻居（普通酰胺/侧链酰胺/带 H 吡咯/胍基的并集）→ 排除**；**无 H（Pro N / His 无 H 吡啶型 N）→ 保留**） |
| CHARMM | 类型表 `CHARMM_ACCEPTOR_TYPES` 已逐类型判定（蛋白剔 NH1/NH2/NH3/NC2/NY/NR1/NR3，留 N/NR2；CGenFF 剔 10 留 6），无二义 N 类型 |
| OPLS | 类型表 `OPLS_ACCEPTOR_TYPES` 已逐类型判定（氮只留 opls_239/511/900，剔 10） |

受体资格（9 类 N 环境 + 文献）详见 [力场类型映射](force_field.md) 的"受体资格判定"节与 `doc/acceptor_identification_evidence.md`。

### 带电基团（charged_positive / charged_negative）

三层递进（基类 `_find_charged`；GROMOS/CHARMM/OPLS 换字典或覆写验证）：

1. **第一层：残基名字典**——按残基名 + 原子名清单直接产出带电基团（`_identify_protein_charged`）

   | 力场 | 正电字典 | 负电字典 |
   |:-----|:---------|:---------|
   | Amber | ARG, LYS, HIP, ORN, DAB, M3L, MLY | ASP, GLU, CYM, KCX, PCA, SEP, TPO, PTR |
   | GROMOS | ARG, LYSH, HISH | ASP, GLU |
   | CHARMM | LYS, ARG, HSP（含咪唑碳 CG/CE1/CD2） | ASP, GLU, CYM |
   | OPLS | ARG, LYSH, HISH（含咪唑碳 CG/CE1/CD2） | ASP, GLU |

2. **第二层：官能团模式匹配**（参照 PLIP `is_functional_group`，力场无关）：正电——季铵（N 4 键无 H）、叔胺（N ≥3 键）、胍基（C 连 3 个 N 且至少一个 N 仅接该 C）、锍（S 3 键无 H）；负电——磷酸盐（P 邻居全 O）、磺酸（S 3 个 O 邻居）、硫酸盐（S 4 个 O 邻居）、羧酸盐（C 2 个 O + 恰 1 个 C 邻居）

3. **第三层：净电荷验证**——正电基团净电荷必须 > +0.1，负电基团必须 < −0.1（`CHARGE_THRESHOLD = 0.1`）

**CHARMM / OPLS 特殊**：铵氮（N 端 NH₃⁺、LYSH 侧链）部分电荷为负，单原子验证会失败；两力场覆写官能团层，仅对满足**末端铵结构判据**的氮（`_is_terminal_ammonium`：N 邻居 ⊆ {C,H}，且重原子邻居不连 ≥2 个 N），把电荷原子扩为 N + 键连 H 后验证。

**GROMOS 特殊**：其质子化残基的正电 N 本身带正电荷，官能团层会对每个 N 单原子额外产出正电基团，与字典层完整基团构成子集关系；`_deduplicate_charged` 覆写为删除同类型下被更完整基团严格包含的基团（Amber 无此问题）。

### 卤键供体 / 受体（halogen_donor / halogen_acceptor）

- **供体**：卤素（F/Cl/Br/I）键连碳原子 → 基团 atoms = [C, X]（第一个是碳，第二个是卤素）
- **受体**：中心原子 A ∈ {C, P, S}，其邻居 R ∈ {O, P, N, S} 非空 → 基团 atoms = [A, R₁, R₂, …]（几何判定用 X···A-R 角）

### 疏水原子（hydrophobic）

| 力场 | 判据 |
|:-----|:-----|
| Amber / CHARMM / OPLS | 碳原子，且**全部键连邻居 ∈ {C, H}**（无极性取代基） |
| GROMOS（反列举） | 类型 ∈ `GROMOS_HYDROPHOBIC_TYPES`（`C, CH0, CH1, CH2, CH3, CH4, CH2r`，排除 CH3p/CR1），且**邻居不含 O/N/S**——因联合原子力场脂肪族 C 无 H 邻居 |

### 金属配位（metal / metal_binding）

- **金属中心**：元素符号 ∈ `METAL_IONS`（42 种，由读取器推断元素，不依赖类型名）
- **配位原子**：元素 ∈ {O, N, S} 的原子（`_find_metal_binding`）；**水残基整体排除**（对残基名 ∈ `self.WATER_RESIDUES` 直接返回空）——因此"纯水分子配位"不会被报告，且按力场经 `WATER_RESIDUES` 正确排除（CHARMM TIP3 / OPLS HO4/HO5）

### 水（water）

残基名 ∈ `WATER_RESIDUES`（基类类属性，各力场集合见[力场类型映射](force_field.md)），基团 atoms = 整个水残基原子（OW/HW 由原子名标识）。

## 类型名歧义与处理

- **`C` 双身份**（Amber）：主链羰基碳 vs Tyr CZ 芳环碳。环内 ≥ n−1 个强芳香原子（六元环为 5）时 `C` 作为兼容类型参与共轭。
- **`N3` vs `n3`**：蛋白正电氨基（铵）vs GAFF 中性胺，大小写不同义。
- **`CA` 跨力场不同义**：Amber = 芳香碳，GROMOS = α 碳，必须经各力场独立特征表。
- **N-H 键在 Constraint 段**：供体识别必须依赖 Bond + Constraint（+ Settle）合并后的键图。

## 验证基准

在 Amber 真实测试体系（KRAS-RBD D927，`Tests/test_MD_case_amber`）上：

- **芳香环**：配体识别出 3 个芳香环 + 1 个非芳香含硫环（2,3-二氢噻吩式稠合杂环，C22=C23 双键非芳香噻吩）；RBD 的 27 个环（Pro×9 + Tyr×11 + His×3 + Trp×2 + Phe×1）判定结果与化学事实一致
- **供体**：RBD 的 263 个供体 H 电荷全部为 +0.19~+0.45，无 q(H) ≤ 0

真实集成测试见 `Tests/unittests/test_*_real_identifier.py` 与 `test_<ff>_interactions_{per_frame,two_pass}.py`（3 力场：Amber / GROMOS / CHARMM）。

## 支持边界

- **全原子力场**（Amber / CHARMM / OPLS）依赖显式氢原子做供体与疏水邻接判定，完整可用
- **GROMOS** 已支持（53A6/54A7），但为联合原子力场：极性 H（N-H、O-H）与芳香 C-H 显式，脂肪族非极性 H 并入 CH1/CH2/CH3——类型粗的部分用结构判据/分级判定/反列举补偿（受体、芳香环、疏水），依赖脂肪族 H 的判定（供体邻接细节、疏水取代基检查）受限；GROMOS 配体（ATB 参数化）不在第一版支持承诺内
- **键序缺失鲁棒**：tpr 即使不保留键序（全 func=1），类型名已编码芳香性，芳香判定可靠；键序/平面性不作为基团鉴定判据
