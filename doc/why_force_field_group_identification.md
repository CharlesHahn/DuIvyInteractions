# 为什么做"基于力场的基团识别"——差异化论证与竞品空缺分析

> 创建日期：2026-09-28
> 调研方式：anysearch 检索文献 + 源码交叉验证（Mobley 2018 / GAFF 2004 / ProLIF GSoC / InterMap 基准）
> 用途：回答"我们的特异性竞争优势是什么"与"为什么竞品不做基于力场的基团识别"；对外论证的完整材料
> 状态：**论证记录**（非 agent 开发指令）

---

## 1. 结论先行

1. **KDTree 不是差异化**：InterMap（2025）同样用 k-d tree 加速距离查询，这是性能共性。
2. **特异性竞争优势 = "确定性基团鉴定"整包**，由三件事构成：
   - **确定性/零推断**：H 键供体（D–H 键条目 + q(H)>0）、水桥（SOL 残基名）、金属（元素+电荷）100% 确定；
   - **帧无关**：基团鉴定只做一次、全轨迹复用，与 KDTree 的几何加速叠加；
   - **与力场自洽**：结果与模拟参数同源，无"坐标推断化学"的系统偏差。
3. **本质**：我们承担了别人不愿付的成本（跨力场特征映射表、绑定 GROMACS + Amber 生态），换来了别人结构上做不到的零推断。

---

## 2. 为什么别人不做？——四层原因

### 2.1 历史路径依赖：从 PDB 起家的惯性

PLIP 诞生于晶体学/对接场景，**PDB 里没有力场原子类型字段**，OpenBabel 重建是唯一选择；后续 MD 工具（ProLIF、getContacts、MD-LR、InterMap）大多从静态工具或通用分析生态演化，延续"重建化学"惯性。**dock pose 也没有类型**——"读类型"只在"有 tpr 的 GROMACS MD 分析"场景可行。

### 2.2 主流思潮是"逃离原子类型"（最关键的原因）

Mobley et al. 2018《Escaping atom types in force fields using direct chemical perception》（JCTC）将类型法贬称为 **indirect chemical perception**，指出其缺点：
- 类型必须编码足够化学环境信息 → 类型爆炸、参数冗余；
- 类型定义基于化学直觉、无严谨基础，**跨力场不可移植**；
- 会引入错误（联苯桥头键被误判为芳香键）。

力场社区随之转向 **SMIRNOFF/SMIRKS direct chemical perception**（直接以 SMARTS 匹配化学子结构分配参数）。分析工具作者跟随此思潮——ProLIF 论文"SMARTS 比力场特定原子类型更普适"即此价值观的直接产物。**在"直接化学感知"被视为先进方向的语境下，没人愿意回头读被批判的"原子类型"。**

### 2.3 技术成本高

- 类型名跨力场不同义（GAFF `ca` vs CHARMM `CG2R61` vs OPLS `CA`）→ 需维护特征映射表；
- 联合原子力场（GROMOS UA）无显式 H → 供体鉴定失效；
- tpr 二进制格式随 GROMACS 版本变动，解析成本高；
- 蛋白 + 配体常来自不同力场（本项目：amber14sb 蛋白 + GAFF 配体，需双命名处理）。

做"普适工具"的人不愿付此成本，宁可维护一套力场无关的 SMARTS。

### 2.4 产品目标错位

ProLIF 卖点是"任何复杂（蛋白/核酸/DNA/RNA）、任何格式、力场无关"。"读 tpr 类型"天然绑定 GROMACS + 特定力场家族，与其产品哲学冲突——被他们设计为"不值得做"。

---

## 3. 文献证据链（间接证据充分，直接讨论缺失）

**直接、专门讨论"用力场原子类型做基团鉴定"的文章：没有**（这本身就是空位存在的证据）。间接文献构成完整证据链：

| 文献 | 支持论据 |
|:-----|:---------|
| GAFF 论文（Wang 2004）：*"similar chemical environments are encoded as atom types"* | **类型 = 化学判决的记录**——"直读类型 = 复用判决"的权威出处 |
| Mobley 2018（批判原子类型） | 反方视角：类型不可移植、会引错；但其中 OEAntechamber 用 SMARTS 表达 GAFF 类型、Foyer/OpenBabel 用 SMARTS 编码 UFF/GAFF 的证据，**恰好证明 GAFF 类型与化学子结构一一对应**——读类型本质是"读已经做好的 SMARTS 判决" |
| ProLIF 作者（cbouy）GSoC 报告 | 作者本人承认："MD 拓扑不保留键序/电荷，必须推断"，且**推断依赖原子读取顺序**（azathioprine 例：读序导致共轭碳带负电、硝基非标准，需事后标准化）——"推断易错"的第一手证据 |
| InterMap 基准（2025） | 指出 ProLIF 把**肽键 N 误判为 H 键受体**（酰胺 N 无孤对），SMARTS 受体模式过宽——零推断的 DII 结构上不会犯 |
| MDAnalysis RDKitConverter 文档 | *"键序/电荷必须从拓扑猜测"*——整个 ProLIF 技术栈建立在"猜测"上 |

---

## 4. 关键反方思辨：Mobley 的批判为什么打不到我们

Mobley 批判原子类型针对的是**参数分配**任务——该场景要求类型编码完整化学环境才能查参数表，故类型爆炸。**我们的任务不是参数分配，而是基团鉴定**：

| 任务 | 对原子类型信息的需求 | Mobley 批评是否成立 |
|:-----|:---------------------|:--------------------|
| 参数分配 | 需编码全部化学环境以查参数表 | 成立（类型爆炸） |
| 基团鉴定 | 只提取部分化学特征（杂化/芳香/极性/带H/孤对） | 不成立（读取已压缩的特征摘要） |

GAFF 类型与 SMARTS 化学子结构一一对应（OEAntechamber 的证明），读类型 = 读取一次免费的、与力场自洽的化学判决，避免"从坐标反推"的噪声（ProLIF 的读序依赖问题、肽键 N 受体误判）。

**一句话**：别人不做，是因为他们 (a) 被"逃离原子类型"思潮主导，(b) 追求力场无关的普适性，(c) 没有意识到"基团鉴定"与"参数分配"对原子类型信息的需求强度完全不同。本项目不是反潮流，而是在潮流覆盖不到的细分场景（GROMACS tpr + Amber 家族 + 全原子 MD）里，把"被批判的原子类型"从参数分配语境中解放出来，用作零推断的化学判决源。

---

## 5. 诚实的边界（论证必须包含的部分）

优势有明确适用域：**GROMACS tpr + Amber 家族 + 全原子显式 H**。超出此域（任意力场、任意格式、对接 PDB、核酸/脂质复杂体系），普适性是 ProLIF/InterMap 的地盘。

卖点叙事应为"**确定性 + 与模拟力场自洽**"，而非"更准/更快"——后者会被 InterMap 的性能与 ProLIF 的通用性夹击。

---

## 6. 参考来源

- Mobley et al. 2018《Escaping atom types in force fields using direct chemical perception》，JCTC 14:6076：<https://pmc.ncbi.nlm.nih.gov/articles/PMC6245550/>
- Wang et al. 2004《Development and testing of a general amber force field》，JCC 25:1157：<https://onlinelibrary.wiley.com/doi/10.1002/jcc.20035>
- Wang et al. 2006《Automatic atom type and bond type perception in molecular mechanical calculations》（antechamber），J Mol Graph Model 24:247：<https://www.sciencedirect.com/science/article/abs/pii/S1093326305001737>
- ProLIF 作者（cbouy）GSoC 报告（MDAnalysis ↔ RDKit 转换，推断键序的读序依赖）：<https://www.mdanalysis.org/2020/08/29/gsoc-report-cbouy/>
- MDAnalysis RDKitConverter 文档：<https://docs.mdanalysis.org/stable/documentation_pages/converters/RDKit.html>
- InterMap 预印本（2025，含 ProLIF 肽键 N 误判与四工具 benchmark）：<https://www.biorxiv.org/content/10.64898/2025.12.15.694195v1.full-text>
- ProLIF 论文（SMARTS 定位声明）：<https://pmc.ncbi.nlm.nih.gov/articles/PMC8466659/>
- OEAntechamber（SMARTS 编码 GAFF 类型，GitHub 镜像）：<https://github.com/choderalab/oeante>

---

## 7. 关联文档

- `doc/project_background.md` —— 动机与差异定位论证（本文档的上级归档）
- `doc/competitor_landscape_survey.md` —— 竞品全景（MD 轨迹/构象相互作用检测）
- `doc/plip_prolif_group_identification_survey.md` —— PLIP / ProLIF 基团鉴定机制细节

---

*文档结束*
