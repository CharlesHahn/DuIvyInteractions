# GROMOS 力场基团识别器设计方案

> 创建日期：2026-10-07
> 修订：2026-10-07（原子类型对抗性审查 + 53A6/54A7 家族兼容性审查后定稿）
> 分支：feature/gromos-identifier
> 状态：方案定稿（待实现）
> 用途：新力场入库——在现有 "新力场入库 = 填特征表" 架构承诺下，评估并设计 GROMOS 基团识别器（覆盖 53A6 / 54A7 家族）

> ⚠️ **真实数据测试状态（2026-10）**：GROMOS 53A6 真实数据集成测试
> （`Tests/test_MD_case_gromos/`，tpr+xtc）已添加，其数值断言（基团计数、
> 相互作用 pair 数）为**程序自动探测，尚未经人工核验**。测试结果未经
> 与 gmx dump 或人工分子分析的交叉确认，引用须谨慎。

---

## 1. 背景

DII 的基团鉴定依赖四类证据源（见 `amber_ff_identifier.py`）：

1. 力场类型名语义（STRONG_AROMATIC / ACCEPTOR_TYPES）
2. 键合图（环检测、官能团邻居、卤键连接、疏水邻居、D–H 对）
3. 元素 + 电荷（供体 D 元素、q(H)>0、q<0、净电荷、金属元素）
4. 残基/原子名（带电残基字典、水残基名）

当前映射表基于 Amber 家族（GAFF + ff14SB 等）。本设计评估 GROMOS 54A7 能否按"填特征表"接入，并给出逐项识别方案。

**前置调研**：`doc/force_field_compatibility_survey.md`（三力场兼容性调研）结论——GROMOS 是唯一"类型名不编码化学环境"的力场，对应需降级证据源；本设计为该结论的实现落地。

**实证依据**：本地 GROMACS 2019（`D:\CharlieAPP\gmx2019_05_GPU\share\gromacs\top\gromos54a7.ff\`）的 `atomtypes.atp`、`aminoacids.rtp`、`spc.itp` 实际文件内容，与 GitHub GROMACS master 逐字一致。

---

## 2. GROMOS 54A7 关键事实（本地 rtp/atp 实测）

### 2.1 原子类型表（atomtypes.atp）

| 类型 | 语义 | 备注 |
|---|---|---|
| O / OM / OA / OE / OW | 羰基 / 羧基 / 羟基·糖·酯 / 醚·酯 / 水 氧 | |
| N / NT / NL / NR / NZ / NE | 肽胺 / 终端 NH2 / 终端 NH3 / **芳香氮** / Arg 胍 NH2 / Arg NE | ⚠️ **类型名 ≠ 原子名**：LYS 侧链铵氮原子名 NZ 但类型是 **NT**；ARG 胍氮原子名 NH1/NH2 且类型是 **NZ** |
| C | **裸碳（bare carbon）** | 羰基碳、羧基碳、**标准蛋白芳香环碳全部同名** |
| CH0 / CH1 / CH2 / CH3 / CH4 / CH2r | sp3 脂肪碳（UA） | CH2r = 环内 CH2 |
| CH3p | **极性脂肪族 CH3（54A7 新增）** | 胆碱 N+(CH3)3 的甲基；53A6 无此类型 |
| CR1 | **芳香/烯 sp2 CH 组（united-atom aromatic CH）** | 用于芳香环碳（黄素核 FMNO 的 FC6/FC9，环内 gb_16 键）**也用于链状烯烃碳**（MEBMT 的 CE=CZ，非环）——**是否芳香必须靠键合图环检测判定，不能仅凭 CR1**；标准蛋白残基未采用（用 C） |
| HC / H | 碳上氢 / 非碳氢（极性 H） | 极性/芳香 H 显式 |
| S / P / F / CL / BR | 硫 / 磷 / 卤素 | |
| CU1+ / CU2+ / FE / ZN2+ / MG2+ / CA2+ / NA+ / CL- | 金属/离子 | 元素可解析 |

### 2.2 蛋白芳香残基实测（aminoacids.rtp）

| 残基 | 环碳类型 | 环氮类型 | 环氢类型 |
|---|---|---|---|
| PHE | **C ×6** | — | HC |
| TYR | **C ×6** | — | HC |
| TRP | **C ×8**（含五元环）| **NR ×1**（NE1）| HC + H(HE1 极性) |
| HISA/HISB/HISH/HIS1/HIS2 | **C ×3**（咪唑）| **NR ×2**（ND1/NE2）| HC + H(HD1 极性) |

**结论**：GROMOS **标准蛋白残基**的芳香碳类型退化为 `C`（与羰基/羧基碳同名）；但类型体系本身**存在**芳香标注——`CR1`（芳香/烯 sp2 CH）、`NR`（芳香氮）。注意：CR1 既用于真芳香环（黄素核 FMNO 的 FC6/FC9）**也用于链状烯烃**（MEBMT 的 CE=CZ 非环），**仅凭 CR1 不能判芳香，必须环检测先行**。因此"芳香性完全不在类型中编码"的表述**不准确**：应为"标准蛋白芳香碳 = C（类型信号缺失），类型体系另有 CR1/NR 可作信号但需键合图验证"。类型语义证据源在 GROMOS 下**部分保留**（NR 可靠；CR1 需环验证）。

### 2.3 其他命名

- 带电残基：**ARG / LYS / ASP / GLU**（标准名）+ **HISH**（质子化 His）+ 中性变体（ARGN/LYSH/ASPH/GLUH）
- 水：**SOL**（spc.itp 实测，OW/HW1/HW2；另有 spce/tip3p/tip4p itp）
- 离子残基：`NA+ / CL- / MG2+ / ZN2+ / CA2+` 等
- 显式 H：极性 H（N-H、O-H）与芳香 C-H 全部显式；仅脂肪族非极性 H 并入 CH1/CH2/CH3（UA）

### 2.4 53A6 / 54A7 家族兼容性（实测 + 权威论文）

**结论：两版本原子类型表唯一差异 = 54A7 新增 `CH3p`（极性脂肪族 CH3，胆碱部分）**；芳香识别相关的 `C`/`CR1`/`NR`/`HC` 及其在残基中的用法**完全一致**（含 MEBMT/FMNO 的 CR1 用法）。识别器**无需按版本分支**，注册表单键 `"gromos"` 覆盖家族即可。

| 证据源 | 内容 |
|---|---|
| 本地 `atomtypes.atp` 严格 diff | 53A6 共 56 个类型、54A7 共 57 个；仅 54A7 多 `CH3p 15.03500` |
| 本地 `aminoacids.rtp` 逐原子核对 | PHE/TRP/HISA 的 C/NR/HC 用法两版本逐字一致 |
| [Schmid 2011（54A7/54B7 定义论文）](https://pubmed.ncbi.nlm.nih.gov/21533652/) | "54A7 … introduces a charged -CH3 atom type"——印证唯一差异为带电 CH3 |
| [Oostenbrink 2005（53A6 定义论文）](https://www.research-collection.ethz.ch/bitstreams/f7b64d3e-73f2-4244-ac8f-0166afb2d675/download) | 53A6 参数集背景 |

> `CH3p` 对识别的影响：其连接原子是 N（胆碱 N+），疏水判据"邻居 ∉ {O,N,S}"天然将其排除，不会被误判为疏水；测试补一条胆碱断言即可。

---

## 3. 逐项基团识别方案

| # | 基团 | Amber 判据 | GROMOS 判据 | 做法 |
|---|---|---|---|---|
| 1 | H_donor | D–H 键 + q(H)>0，D∈{N,O,S,F} | **同 amber**（极性 H 显式） | 复用 `_classify_dh_pair`，零改动 |
| 2 | H_acceptor | 类型 ∈ ACCEPTOR_TYPES + q<0 | 类型 ∈ GROMOS 受体表 + q<0 | 新建 `GROMOS_ACCEPTOR_TYPES`（见 §4），逻辑同 amber |
| 3a | 蛋白带电 | 残基名字典 + 原子名 | 字典键改 GROMOS 名，原子列表按 GROMOS rtp 实测调整、只含 5 个残基（见下方） | 需扩展字典 |
| 3b | 配体带电 | 官能团模式（元素+邻居） | **同 amber**（力场无关） | 复用，零改动 |
| 4 | 卤键供/受体 | 卤素+邻居 | **同 amber**（力场无关） | 复用，零改动 |
| 5 | 金属 | 元素 ∈ METAL_IONS | **同 amber**（GNOMOS 离子类型元素可解析） | 复用，零改动 |
| 6 | 水 | 残基名 ∈ WATER_RESIDUES | **SOL**（已在集合内） | 零改动 |
| 7 | 芳香环 | 类型 ∈ STRONG_AROMATIC + 环检测 | **分级：环内含 CR1/NR 直接判芳香；环内全 C 回退残基名白名单**（见 §5） | 需新实现 |
| 8 | 疏水 | C + 邻居 ∈ {C,H} | **重写：C/CH* + 邻居 ∉ {O,N,S}**（反列举） | 需新实现 |
| 9 | 金属配位 | 元素 ∈ {O,N,S} 非水 | **同 amber**（力场无关） | 复用，零改动 |

GROMOS 带电残基原子列表（本地 `aminoacids.rtp` 逐原子实测 + `aminoacids.hdb` 加氢规则确认，净电荷 Python/PowerShell 双独立验证）：

| GROMOS 残基 | 原子列表 | 净电荷 | 说明 |
|---|---|---|---|
| ARG | CZ, NE, NH1, NH2, HE, HH11, HH12, HH21, HH22 | +1.000 | 胍基阳离子；与 Amber `ARG` 原子集一致 |
| **LYSH** | NZ, HZ1, HZ2, HZ3 | +1.000 | **质子化 LYS（ε-铵 NH3+），对应 Amber `LYS`（Amber 默认恒为 NH3+）** |
| **HISH** | ND1, NE2, HD1, HE2 | +1.000 | 质子化 His（咪唑双质子化）；原子集与 Amber `HIP` 兼容 |
| ASP | CG, OD1, OD2 | -1.000 | 羧酸根 |
| GLU | CD, OE1, OE2 | -1.000 | 羧酸根 |

> **质子化变体约定（关键）**：GROMOS 将可质子化残基拆为中性/质子化两个 rtp 条目——**LYS（NH2，净电荷 0）为中性，LYSH（NH3+，+1）才带正电**；同理 ASPH/GLUH/ARGN/HISA/B/1/2 均为中性变体，**不进带电字典**。CYS（-0.5 硫醇根）强度不足、PLIP 通行不含，**排除**。
>
> **字典裁剪声明**：Amber 字典中的修饰残基（ORN/DAB/M3L/MLY/CYM/KCX/PCA/SEP/TPO/PTR）**在 GROMOS 54A7 rtp 中全部不存在**（grep 零命中）——GROMOS 版本带电字典**只含 ARG/LYSH/HISH（正电）+ ASP/GLU（负电）**。

### 3.1 改动分类总览（核心：GROMOS 不是"只填类型表"）

对照 Amber 实现，GROMOS 版本按改动量分四类——**明确区别"填表"与"重写逻辑"**：

| 类别 | 基团 | 说明 |
|---|---|---|
| **A. 零改动复用** | H 键供体、配体带电（官能团）、卤键、金属、水、金属配位 | 判据基于元素/电荷/键合图/残基名，力场无关，直接复用 amber 方法 |
| **B. 仅换类型表** | H 键受体 | `ACCEPTOR_TYPES` → `GROMOS_ACCEPTOR_TYPES`（§4）；判定逻辑不变 |
| **C. 改字典+数据** | 蛋白带电 | 键名（IP→HISH）+ 原子列表（LYS 去 HZ3）+ 裁剪修饰残基；逻辑结构（残基名→原子列表）不变 |
| **D. 重写识别逻辑** | **芳香环、疏水** | C 类型一型多用（羰基/芳香/烯共用）→ 类型语义不足，必须组合"类型+残基名/键合图"判据（§5 分级判定、§6 反列举）——**不是填表能解决的** |

> **关键结论**：项目 CLAUDE.md 的"新力场入库 = 填特征表"承诺**对 GROMOS 不成立（完全填表）**——只有 A/B/C 类推此路径；**D 类（芳香环、疏水）因类型名不编码化学环境必须重写逻辑**。这也正是调研文档（`force_field_compatibility_survey.md` §5）"GROMOS 是唯一特例、需降级证据源"结论的代码层落实：**A/B/C 占 ~70% 代码量（直接复用），D 是新增的约 2 处逻辑**。

---

## 4. GROMOS 受体类型表（新建）

```python
GROMOS_ACCEPTOR_TYPES = frozenset({
    "O", "OM", "OA", "OE", "OW",       # 氧：羰基/羧基/羟基/醚/水
    "N", "NT", "NL", "NR", "NZ", "NE", # 氮：肽胺/终端/芳香/胍基
    "S",                                # 硫
    "F", "CL", "BR",                   # 卤素
})
```
判定逻辑与 amber `_find_acceptors` 相同（类型 ∈ 表 且 q < 0）。实测电荷支持：主链 O(-0.38)、ASP OD(-0.635)、ARG NH1(-0.26)、TRP NE1(-0.10) 均负。

---

## 5. 芳香环识别方案（核心难点）

**问题**：GROMOS 标准蛋白芳香碳类型 = `C`（与羰基/羧基碳同名），Amber 式 STRONG_AROMATIC 类型表对蛋白环失效；但类型体系另有 `CR1`（芳香 CH）、`NR`（芳香氮）可作强信号。

**方案（第一版）——分级判定**：
1. **环检测**：复用 amber `_detect_rings`（残基内 BFS）。
2. **环成分过滤**：环内原子类型全 ∈ {`C`, `CR1`, `NR`}（GROMOS 环碳/环氮全集；`CR1` 为修饰残基/黄素核的芳香碳，`NR` 为蛋白+核酸芳香氮）。
3. **强信号优先**：环内检测到 **`NR`** → **直接判定芳香环**（NR 语义专一，可靠）；环内含 **`CR1`** 需结合环成分确认（CR1 也用于链状烯烃，但**环检测先行**已保证此处是闭合环，环内 CR1 + 全 C/NR 成分即可判芳香）。此时无需限定残基——FMNO 等使用 CR1 的环自动覆盖。
4. **弱信号回退**：环内**全为 `C`**（无 CR1/NR）→ 回退残基名白名单：仅当残基 ∈ {PHE, TYR, TRP, HISA, HISB, HISH, HIS1, HIS2} 时判为芳香环（GROMOS 标准蛋白只有这 4 类残基含芳香环；其环碳恰好全 C）。

> 第一版不做非蛋白、非 CR1/NR 的纯 C 芳香环（GROMOS 配体参数化靠 ATB，第一版不支持配体；核酸碱基环含 NR，走强信号路径）。

**不采用**（第二版候选）：几何平面性兜底（PLIP 风格，环原子到拟合平面最大偏差 ≤ 阈值）——用于无任何类型/残基信号的极端场景。

**环序约定**：`Group.atoms` 仍按 BFS 路径序输出，检测器（PiStacking/π-cation 的 ring normal/center）与力场无关，零改动。

**判据分级与 Amber 分层思想对照**：强信号（CR1/NR）↔ STRONG_AROMATIC；弱信号（C+残基名）↔ COMPATIBLE_TYPES——识别器内部可抽象"强/弱"两层特征集合。

---

## 6. 疏水重写方案

- amber 判据："C 且所有邻居 ∈ {C, H}"——GROMOS UA 无碳上 H 邻居，失效。
- **GROMOS 判据**："C（含 CH0/CH1/CH2/CH3/CH4/CH2r）且所有邻居 ∉ {O, N, S}"（反列举排除极性取代）。
- **CH3p 处理**：54A7 的极性 CH3（胆碱 N+(CH3)3）邻居为 N → 被"邻居含 N"排除，天然不误判为疏水；显式纳入测试即可。
- **芳香环碳命中边界（实测）**：PHE/TYR 芳香环碳（类型 `C`）的邻居全为 C/HC → **会落入疏水判据**；TRP/HIS 环碳因邻居含 `NR`（芳香氮）被排除。此行为与 Amber 一致（Amber 的 PHE 环碳同样命中疏水判据）。**疏水-芳香去重属既有 TODO #2（未实现），GROMOS 版不引入新问题，但结果会含芳香环碳之间的疏水接触，需在测试中锁定现状**。
- 原子粒度按 UA 碳原子计，检测器距离判据不变。

---

## 7. 架构与注册

```
DuIvyInteractions/group_identifiers/gromos_identifier.py   # 新文件
DuIvyInteractions/group_identifiers/__init__.py             # 注册表加一项
```

```python
IDENTIFIER_CLASSES = {
    "amber": AmberFFGroupIdentifier,
    "gromos": GromosFFGroupIdentifier,   # 新增：覆盖 GROMOS 53A6 / 54A7 家族（见 §2.4 兼容性）
}
```

`dii run --ff gromos` 即可用；pipeline.py / DII.py / 全部检测器零改动。注册名 `gromos` 为家族键（两版本类型语义一致，无需按版本分支）。

**功能边界（明确声明）**：
- ✅ 可做：供体/受体/带电（蛋白）/卤键/金属/水/金属配位 + 降级芳香环 + 重写疏水
- ❌ 第一版不支持：GROMOS 配体（无 GAFF 等价物，需 ATB 工具链，属未来工作）

---

## 8. 可靠性与风险

| 项 | 置信度 | 说明 |
|---|---|---|
| 供体/受体/金属/水 | 高 | 显式 H + 电荷，与 amber 同证据强度 |
| 带电 | 高 | GROMOS 残基命名已实测（ARG/LYS/ASP/GLU/HISH）|
| 芳香环 | **中高** | 强信号（CR1/NR）类型语义有效；全 C 环靠残基名白名单（PHE/TYR/TRP/HIS）兜底，弱于 amber 但优于纯几何推断；非天然环不可识别 |
| 疏水 | 中 | 反列举语义，边缘情况（如 C 连卤素、CH3p）需测试 |
| 配体 | — | 第一版不支持 |
| 53A6/54A7 兼容 | 高 | 类型表唯一差异 CH3p，芳香识别相关类型零差异（§2.4）|

---

## 9. 里程碑（合并 dev 判据）

1. `gromos_identifier.py` + 注册表 → `dii run --ff gromos` 可跑通不报错
2. 非芳香类（供体/受体/盐桥/卤键/金属/水桥/金属配位）真实数据跑通
3. 芳香类（π-π / π-阳离子）分级方案跑通（强信号 CR1/NR 环 + 全 C 环 PHE/TYR/TRP/HIS 检出）
4. 单元测试绿（`Tests/unittests/test_gromos_identifier.py`：GROMOS 54A7 体系基团计数断言；**含 CH3p 胆碱疏水排除断言**；与 Amber 结果交叉验证）→ 合并 dev

---

## 10. 关联文档

- `doc/force_field_compatibility_survey.md` —— 三力场调研（§5 结论为本设计依据）
- `doc/paper_framework_v1.md` / `doc/why_force_field_group_identification.md` —— 论文框架与差异化论证
- `CLAUDE.md` —— "新力场入库 = 填特征表"架构承诺

## 11. 证据来源（原子类型对抗性审查记录）

本设计的原子类型判定均经本地力场文件 + anysearch 权威来源双重印证：

| 类型 | 本地证据 | 权威来源 |
|---|---|---|
| C（bare carbon） | `gromos54a7.ff/atomtypes.atp`；PHE/TYR/TRP/HISA rtp 环碳全 C | [GROMACS particle-type 官方文档](https://manual.gromacs.org/documentation/2025-rc/reference-manual/topologies/particle-type.html) |
| CR1（芳香/烯 sp2 CH） | `atomtypes.atp` 定义 `aromatic CH-group`；**环内使用**：FMNO 黄素核 FC6/FC9（gb_16 芳香键）；**非环使用**：MEBMT 链状烯烃 CE=CZ（gb_10）——证明不可单凭 CR1 判芳香 | [Horta 2016 JCTC（CR1 united-atom aromatic CH group）](https://pubs.acs.org/doi/abs/10.1021/acs.jctc.6b00187) |
| NR（aromatic nitrogen） | `atomtypes.atp`；TRP NE1 / HISA ND1,NE2 / ADE N 全 NR | [GROMACS particle-type 官方文档](https://manual.gromacs.org/documentation/2025-rc/reference-manual/topologies/particle-type.html) |
| HC（碳上氢）/ H（极性氢） | `atomtypes.atp`；PHE 环氢 HC、TRP HE1 H | [GROMACS particle-type 官方文档](https://manual.gromacs.org/documentation/2025-rc/reference-manual/topologies/particle-type.html) + [GROMOS vol5](https://www.gromos.net/gromos11_pdf_manuals/vol5.pdf) |
| CH3p（54A7 新增） | 53A6/54A7 `atomtypes.atp` 严格 diff（唯一差异） | [Schmid 2011（charged -CH3 atom type）](https://pubmed.ncbi.nlm.nih.gov/21533652/) |
| 53A6/54A7 家族一致性 | 两版本类型表 + 芳香残基 rtp 逐原子 diff | [Schmid 2011](https://pubmed.ncbi.nlm.nih.gov/21533652/) / [Oostenbrink 2005](https://www.research-collection.ethz.ch/bitstreams/f7b64d3e-73f2-4244-ac8f-0166afb2d675/download) |

本地力场路径：`D:\CharlieAPP\gmx2019_05_GPU\share\gromacs\top\gromos53a6.ff\` 与 `gromos54a7.ff\`。

---

*文档结束*