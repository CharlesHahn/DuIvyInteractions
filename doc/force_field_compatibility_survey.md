# 多力场兼容性调研：Amber 之外的 GROMOS / CHARMM / OPLS

> 创建日期：2026-09-28
> 调研方式：anysearch 检索官方手册/论文 + 领域常识交叉验证（见 §7 诚实声明）
> 用途：评估 DII 基于 Amber 家族的基团鉴定指标能否对应到 GROMOS / CHARMM / OPLS；精确化论文图 6"适用域"边界的表述
> 状态：**调研记录**（非 agent 开发指令；未修改任何代码/文档）

---

## 1. 调研背景

DII 的基团鉴定指标（见 `amber_ff_identifier.py`）依赖四类证据源：
① 力场类型名语义（STRONG_AROMATIC / ACCEPTOR_TYPES）
② 键合图（环检测、官能团邻居、卤键连接、疏水邻居、D–H 对）
③ 元素 + 电荷（供体 D 元素、q(H)>0、q<0、净电荷、金属元素）
④ 残基/原子名（带电残基字典、水残基名）

当前映射表基于 Amber 家族（GAFF + ff14SB 等）。本文调研 GROMOS、CHARMM36、OPLS-AA 三个主要力场，逐项评估能否对应，回答"新力场入库 = 填特征表"这一架构承诺的普适边界。

---

## 2. 三个力场的关键事实

### 2.1 GROMOS（54A7 / 54A8 等）
- **联合原子（UA）**：GROMACS 手册原话 "united atom force fields, i.e. without explicit aliphatic (non-polar) hydrogens"；
- **但极性/芳香 H 显式**：GROMOS 官方手册（vol5）"part of the hydrogen atoms **(polar, aromatic)** are explicitly treated"——即 N-H、O-H、S-H 与芳香 C-H 显式，仅脂肪族非极性 H 并入 CH1/CH2/CH3；
- **类型名不编码化学环境**：原子类型粒度粗（O/OM/OA/OE/OW、N/NT/NL/NR/NZ/NE、C、CH1/CH2/CH3…），54A7 电荷表中 Phe/Tyr/Trp 芳香碳的类型名就是 **C**——与羰基碳、羧基碳同名；
- **无芳香键序**：全部单键 func=1，芳香性靠特殊 improper 二面角维持，不在类型/键序中显式；
- **残基命名**：与 Amber 一致，HIS 用 HISA/HISB/HISH 三态；水配 SPC/SPC-E（残基名 SOL）；金属离子为单独类型（CA2+/MG2+/ZN2+…）。

### 2.2 CHARMM36（蛋白力场 + CGenFF 配体力场）
- **全原子、显式 H**（主流 CHARMM22/27/36；老版 CHARMM19/20 曾有联合原子变体，非现代主流）；
- **蛋白力场类型编码化学**：C（羰基）、**CA（sp2 芳香碳）**、CT1/CT2/CT3（sp3 脂肪碳，MATCH 论文确认）、N/O/OH/OS…——与 GAFF 风格类似；
- **CGenFF 配体类型强编码**：命名规则（检索确认）为"元素 + G(36 typing) + 杂化 + 环大小 + 连接数"，如 **CG2R61 = sp2 芳香六元环 C-H**；
- ⚠️ **CA 双义坑**：蛋白中**原子名 CA（α 碳）≠ 类型名 CA（芳香 sp2）**——同名不同义。DII 读 tpr `atom_type` 字段而非原子名，架构天然免疫；
- **残基命名**：HIS → HSD/HSE/HSP 三态；水 TIP3P（残基名 TIP3）。

### 2.3 OPLS-AA（/L）
- **全原子、显式 H**（公认；GROMACS 内类型名带 `opls_XXX` 数字前缀，底层基础类型名如下）；
- **基础类型编码化学，命名最接近 GAFF**：**CA（芳香碳）、CT（sp3）、C（羰基）、CM（sp2 烯）、HC/HA（烷/芳 H）、OS、N/O/OH、NA/NB（芳香 N）**（检索确认）；
- 芳香性靠 improper 二面角维持（同 GROMOS，无芳香键序），但**类型名 `CA` 直接标注芳香**——信号在类型里，不需要键序；
- **残基/水**：蛋白残基与标准三字母一致；水 TIP3P/TIP4P（残基名 TIP3/TIP4，待核对 oplsaa.ff）。

---

## 3. 逐项指标对照（10 项）

| # | 基团 | Amber 指标 | GROMOS | CHARMM36 | OPLS-AA |
|---|---|---|---|---|---|
| 1 | 芳香环 | 类型 ∈ STRONG_AROMATIC + 环检测 | ❌ 基本失效（芳香碳类型名=C，无专用类型；需 rtp 反推或几何平面性替代） | ✅（CA 蛋白 + CG2R61 系 CGenFF） | ✅（CA + NA/NB） |
| 2 | H 键供体 | D–H 键 + q(H)>0，D∈{N,O,S,F} | ✅ 实际可用（极性 H 显式） | ✅ 全原子 | ✅ 全原子 |
| 3 | H 键受体 | 类型 ∈ ACCEPTOR_TYPES + q<0 | ⚠️ 部分（无现成受体类型表，可改"元素+q<0+H 计数"） | ✅ 需建表（O/OH/OS/N/O2…） | ✅ 需建表 |
| 4a | 蛋白带电 | 残基名字典 + 原子名 | ✅ 需改名（HISA/B/C） | ✅ 需改名（HSD/E/P） | ✅（HIS 写法待核对） |
| 4b | 配体带电 | 官能团模式（元素+邻居） | ✅ 力场无关 | ✅ 力场无关 | ✅ 力场无关 |
| 5 | 卤键供体 | 卤素 + C 邻居 | ✅ 力场无关 | ✅ 力场无关 | ✅ 力场无关 |
| 6 | 卤键受体 | A∈{C,P,S} + O/P/N/S 邻居 | ✅ 力场无关 | ✅ 力场无关 | ✅ 力场无关 |
| 7 | 金属 | 元素 ∈ METAL_IONS | ✅（CA2+/MG2+/ZN2+…） | ✅ | ✅ |
| 8 | 水 | 残基名 ∈ {SOL,HOH,WAT} | ✅（SOL） | ⚠️ 需加 TIP3 | ⚠️ 需加 TIP3/TIP4 |
| 9 | 疏水 | C + 邻居∈{C,H} | ⚠️ 需重写（UA 无 C-H，改"C 且邻居无 O/N/S"） | ✅ 原判据 | ✅ 原判据 |
| 10 | 金属配位 | 元素 ∈ {O,N,S} 非水 | ✅ 力场无关 | ✅ 力场无关 | ✅ 力场无关 |

---

## 4. 四方总对照

| 维度 | Amber 家族 | GROMOS | CHARMM36 | OPLS-AA |
|---|---|---|---|---|
| 显式 H | ✅ 全 | ⚠️ 仅极性/芳香 | ✅ 全 | ✅ 全 |
| 类型名编码化学环境 | ✅ | ❌ | ✅（CGenFF 强编码） | ✅ |
| 芳香环判据可用性 | ✅ | ❌ 需替代 | ✅ | ✅ |
| H 键供体 | ✅ | ✅（极性 H 显式） | ✅ | ✅ |
| H 键受体 | ✅ | ⚠️ 改判据 | ✅ 建表 | ✅ 建表 |
| 疏水判据 | ✅ | ⚠️ 重写 | ✅ | ✅ |
| 水残基名 | SOL/HOH | SOL | TIP3 | TIP3/TIP4 |
| HIS 命名 | HIP | HISA/B/C | HSD/E/P | （待核对） |
| **对应难度** | 基准 | **高（需降级证据源）** | 中（填两张表） | **低（基本换表）** |

---

## 5. 核心结论

1. **CHARMM 与 OPLS-AA 都能对应上**，走的就是代码架构已设计好的"新力场入库 = 填特征表"路线（CLAUDE.md 原话）。工作量 = 各建一套 STRONG_AROMATIC + ACCEPTOR_TYPES 表 + 水残基名单（TIP3/TIP4）+ HIS 三态命名，识别器逻辑零改动；
2. **GROMOS 是唯一特例**：类型名不编码化学（`C` 通吃羰基/芳香），对应需降级证据源（① 类型语义 → 几何平面性 / rtp 反推）。这反证了"类型名语义"这条证据源的独立价值；
3. **CHARMM 的 CA 双义坑（原子名 α 碳 vs 类型名芳香碳）** 与项目 Amber `C` 双身份属同类问题；DII 读 tpr `atom_type` 的架构天然免疫；
4. **对论文图 6"适用域"表述的升级**：不是"只支持 Amber"，而是"**支持类型名编码化学环境的全原子力场**（Amber / CHARMM / OPLS 均可填表接入）；不支持的仅是 GROMOS 这类类型无语义或 UA 的力场"。同时修正此前"GROMOS 无显式 H → 供体失效"的说法——**GROMOS 极性 H 显式，供体鉴定实际可用**，真正障碍是类型无语义（芳香）与疏水 H 缺失。

---

## 6. 可行性与后续工作（如未来扩展）

| 步骤 | 内容 | 工作量 |
|---|---|---|
| 1 | CHARMM 特征表：蛋白 CA/CT* + CGenFF CG* 系 → STRONG_AROMATIC / ACCEPTOR_TYPES | 中 |
| 2 | OPLS 特征表：CA/CT/CM/NA/NB/OS… → 同上 | 低 |
| 3 | WATER_RESIDUES 扩展 TIP3/TIP4 | 一行 |
| 4 | HIS 三态命名（HSD/HSE/HSP、HISA/B/C）进带电字典 | 几行 |
| 5 | GROMOS 支持：降级证据源设计（几何平面性兜底 + 疏水判据重写 + rtp 反推确认芳香环） | 高（独立设计） |

---

## 7. 诚实声明（证据来源分级）

- **检索确认**：GROMOS UA 且极性/芳香 H 显式（GROMOS vol5）；CHARMM36 蛋白 CT* 脂肪碳（MATCH 论文）；CGenFF 命名规则含杂化/环/连接数（mattermodeling 回答 + Vanommeslaeghe 2012）；OPLS 基础类型名 CA/CT/HA/HC/OS 等（MCCCS Towhee 列表 + GROMACS 论坛）；OPLS-AA 芳香性靠 improper 二面角维持（Towhee）；
- **领域常识 + 检索交叉支撑**：CHARMM36 / OPLS-AA 全原子显式 H；CHARMM36 蛋白芳香碳类型 CA；OPLS 芳香碳 CA；TIP3/TIP4 残基名；HIS 三态命名；
- **未单独查证**：GROMACS `oplsaa.ff` 中蛋白 HIS 的确切残基名写法、CHARMM36 蛋白 rtf 中芳香碳类型名的逐残基细节——如需落地扩展，应直接提取 `atomtypes.atp` / rtp 文件核对。

---

## 8. 参考来源

- GROMACS 手册《Force fields》（GROMOS UA 定义）：<https://manual.gromacs.org/current/user-guide/force-fields.html>
- GROMOS 手册 vol5《Program Library Manual》（极性/芳香 H 显式）：<https://www.gromos.net/gromos11_pdf_manuals/vol5.pdf>
- GROMOS 手册 vol3《Force Fields and Topology Data Set》（54A7 类型与电荷表）：<https://www.gromos.net/gromos11_pdf_manuals/vol3.pdf>
- MATCH 原子类型工具集论文（CHARMM36 CT* 语义）：<https://pmc.ncbi.nlm.nih.gov/articles/PMC3228871/>
- CGenFF I 论文（类型命名与芳香判定）：<https://pmc.ncbi.nlm.nih.gov/articles/PMC3528824/>
- OPLS-aa 原子类型列表（MCCCS Towhee）：<https://towhee.sourceforge.net/forcefields/oplsaa.html>
- GROMACS oplsaa.ff / charmm36.ff 的 atomtypes.atp（GROMACS GitHub）：<https://github.com/gromacs/gromacs/tree/main/share/top>
- Mobley 2018《Escaping atom types》（原子类型 vs 直接化学感知，跨力场语境）：<https://pmc.ncbi.nlm.nih.gov/articles/PMC6245550/>

---

## 9. 关联文档

- `doc/paper_framework_v1.md` —— 论文框架 v1（图 6 适用域据此精确化）
- `doc/why_force_field_group_identification.md` —— 差异化论证与文献证据链
- `doc/competitor_landscape_survey.md` —— 竞品全景
- `CLAUDE.md` —— "新力场入库 = 填特征表"架构承诺

---

*文档结束*