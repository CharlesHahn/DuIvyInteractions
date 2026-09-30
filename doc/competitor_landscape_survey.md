# DII 竞品全景调研（MD 轨迹/蛋白-配体构象次级相互作用检测）

> 创建日期：2026-09-28
> 调研方式：anysearch 检索官方文档 / 论文 / benchmark，交叉验证
> 用途：论证本项目"直接读取 tpr 力场原子类型做确定性基团鉴定"在竞争格局中的空位；与 `plip_prolif_group_identification_survey.md`（鉴定机制细节）衔接成系列
> 状态：**调研记录**（非 agent 开发指令）

---

## 1. 调研范围

围绕"**从 MD 轨迹或蛋白-配体构象中检测次级相互作用**"这一功能，梳理全部竞品。竞争格局分三圈：

| 圈层 | 定位 | 代表工具 | 与 DII 的关系 |
|:-----|:-----|:---------|:-------------|
| ① 直接竞品 | MD 轨迹 IFP / 相互作用分析 | ProLIF、getContacts、MD-IFP、MD-LR、InterMap、MD_PLIP | 正面竞争 |
| ② 相邻竞品 | 静态结构 / 对接构象相互作用 | PLIP、Arpeggio、BINANA 2、nAPOLI、PyPLIF HIPPOS、ODDT、IChem/SPLIF 系、NAPS | 功能重叠，但不支持 MD 轨迹 |
| ③ 替代方案 | 单种相互作用的 MD 分析 | gmx hbond、cpptraj hbond、MDAnalysis HydrogenBondAnalysis | 只覆盖氢键等单类 |

**核心判断：所有 ① 圈直接竞品都以"重建化学"为基团鉴定手段（RDKit SMARTS / VMD 命名约定 / OpenBabel），"直接用 tpr 力场原子类型做确定性鉴定"这一空位至今无人站。**

---

## 2. ① 直接竞品：MD 轨迹 IFP 类

| 工具 | 基团鉴定机制 | MD 支持 | 局限（据官方文档 / benchmark） |
|:-----|:-------------|:--------|:-------------------------------|
| **ProLIF**（2021） | MDAnalysis→RDKit（推断键级/电荷）→ **SMARTS** 逐残基匹配 | ✅ 原生 | 靠显式 H 推断化学；长轨迹/复杂体系**内存爆炸**（InterMap 基准：MPRO 上 35 分钟只跑完 1%、内存 >64GB） |
| **getContacts** | VMD 脚本 + **严格原子命名约定**选原子 | ✅ 原生 | **跨力场兼容性差**（命名约定随力场变）；部分相互作用未对配体定义；检测量级少（MPRO 只测出 InterMap 的 4.9%） |
| **MD-IFP**（Kokh/Wade） | RDKit（配体）+ 距离/角度 | ✅ 但专用 | 本质为 **RAMD 配体解离**设计；只支持蛋白-配体 / 蛋白A-蛋白B，**不能分子内、不能核酸** |
| **MD-LR**（2023） | **内部调用 PLIP**（OpenBabel 规则）并行化 | ✅ 但限 GROMACS | 只能 GROMACS 生态输入；配体必须是单个分子组（小分子限定）；**继承 PLIP 不能做分子内 IFP**；核酸不能作独立配体 |
| **InterMap**（2025） | MDAnalysis + **SMARTS**（RDKit 集成）+ k-d tree + Numba + bitarray | ✅ 原生 | 最新最快：比 ProLIF 快 2–3 个数量级、内存 91 倍节省；但**仍是 SMARTS 重建化学**、力场无关路线 |
| **MD_PLIP** | 采样帧喂 PLIP | ⚠️ 帧采样 | 抽帧分析，非全轨迹；继承 PLIP 的 PDB 重建弱点 |

**关键佐证——InterMap benchmark（2025-12，bioRxiv）**：
- ProLIF：MPRO 轨迹 64GB 内存打爆，线性外推需 ~3500 分钟；
- MD-LR：仅检测到 InterMap 的 0.6% 相互作用；
- getContacts：检测量为 InterMap 的 1.2%–4.9%，且内存更高。

→ 说明 **MD 轨迹级相互作用分析的性能与化学完整度仍是开放战场**。

---

## 3. ② 相邻竞品：静态结构 / 对接构象类

| 工具 | 基团鉴定机制 | 特点 |
|:-----|:-------------|:-----|
| **PLIP** | OpenBabel 感知 + 规则 + 残基字典 | 8 类相互作用；无 CONECT 的 PDB 判芳香失败（DII 的实证痛点） |
| **Arpeggio** | **SMARTS 原子类型**（继承 CREDO）+ OpenBabel + BioPython KDTree | 15 种相互作用子类型、SIFt 指纹；支持蛋白-蛋白/核酸；5 Å 截断 |
| **nAPOLI**（Arpeggio 团队） | 图论策略，基于 Arpeggio 接触 | 大尺度保守相互作用检测/可视化 |
| **BINANA 2** | 几何距离/角度规则（Python + JS） | docking 后分析、批量工作流 |
| **PyPLIF HIPPOS** | 原子对规则 | AutoDock Vina / PLANTS docking 后处理 |
| **ODDT** | OpenBabel / RDKit 规则 | StIF 与 PLEC 两种 IFP 方案 |
| **IChem / SPLIF / APIF / SIFt 系** | 药效团 / 结构 IFP（2004–2017 历史谱系） | docking 打分 / pose 排序的老一代 |
| **NAPS** | 网络分析 | 蛋白结构 / 核酸复合物网络统计，非接触级检测 |

---

## 4. ③ 替代方案：单类相互作用的 MD 分析

- `gmx hbond`：GROMACS 官方氢键分析；
- `cpptraj hbond`：AmberTools 氢键分析；
- MDAnalysis `HydrogenBondAnalysis`：Python 生态氢键分析。

共同局限：**只做氢键一类**，无 π-堆积/盐桥/卤键/水桥等精细分类，也不产出"基团-基团"级别的分类结果。

---

## 5. 对 DII 的启示（差异化校验）

1. **"力场类型直读 + MD 轨迹"仍旧是空位**：① 圈竞品全部走"重建化学"路线（ProLIF/InterMap 的 SMARTS-RDKit、getContacts 的命名约定、MD-LR 的 PLIP-OpenBabel），没有一个读 tpr 类型语义——`project_background.md` 的差异定位叙事再次得到外部证据支撑；
2. **性能军备竞赛已经开始**：InterMap（2025）证明 k-d tree 可将 IFP 检测提速 2–3 个数量级——DII 的 TwoPass + KDTree 优化路线与该趋势一致（水桥 65h→~5s），可引 InterMap 作为性能对标对象；
3. **可对标的可复制结论**：InterMap 基准中 ProLIF 的失败模式（长轨迹内存爆炸、肽键 N 被误判为受体）恰好是 DII"零推断 + 稀疏存储 + 供体零推断"能打的点；
4. **新威胁感知**：若未来有人把"SMARTS 匹配 + 力场拓扑缓存"组合（每帧只查几何、不重建化学），会逼近 DII 的定位——但当前无实现，属于前瞻性风险。

---

## 6. 参考来源

- ProLIF 论文（J Cheminform 2021）：<https://pmc.ncbi.nlm.nih.gov/articles/PMC8466659/>
- getContacts 官网（VMD 脚本、结构/动力学）：<https://getcontacts.github.io/>
- MD-IFP 仓库（HITS-MCM，RAMD 解离设计）：<https://github.com/HITS-MCM/MD-IFP>
- MD-LR 论文（MD–Ligand–Receptor，IJMS 2023，PLIP 集成 + HPC）：<https://www.mdpi.com/1422-0067/24/14/11671>
- InterMap 预印本（2025，k-d tree + Numba + SMARTS，含四工具 benchmark）：<https://www.biorxiv.org/content/10.64898/2025.12.15.694195v1.full-text>
- MD_PLIP 仓库（帧采样 + PLIP）：<https://github.com/cmwoodley/MD_PLIP>
- PLIP 官方算法文档：<https://github.com/pharmai/plip/blob/master/DOCUMENTATION.md>
- Arpeggio 论文（15 种相互作用子类型，JMB 2017）：<https://pmc.ncbi.nlm.nih.gov/articles/PMC5282402/>
- BINANA 2（JCIM 2022）：<https://pmc.ncbi.nlm.nih.gov/articles/PMC8889568/>
- PyPLIF HIPPOS（JCIM 2020）：<https://pubmed.ncbi.nlm.nih.gov/32687350/>
- ODDT（J Cheminform 2015）：<https://link.springer.com/article/10.1186/s13321-015-0078-2>
- NAPS（NAR 2016）：<https://pmc.ncbi.nlm.nih.gov/articles/PMC4987928/>
- gmx hbond 手册：<https://manual.gromacs.org/current/reference-manual/analysis/hydrogen-bonds.html>
- MDAnalysis HydrogenBondAnalysis：<https://docs.mdanalysis.org/stable/documentation_pages/analysis/hydrogenbonds.html>

---

## 7. 关联文档

- `doc/plip_prolif_group_identification_survey.md` —— PLIP / ProLIF 基团鉴定机制细节（本文档的机制级补充）
- `doc/project_background.md` —— 项目动机与差异定位论证（本文档的论据归档）

---

*文档结束*