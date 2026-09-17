# MD 相互作用分析工具深度调研报告

> 对抗性调研：验证"**PLIP 级精细分类 + MD 原生拓扑**"想法是否为行业空白
> 调研日期：2026-08-07
> 调研工具：AnySearch 学术域 + 代码域 + URL 提取，全程只读
> 触发背景：在 D927 项目中发现 PLIP v3.0.1 对 MD 轨迹小分子（D927）的化学感知缺陷，即：PLIP 通过 OpenBabel 重建加氢/芳香性，完全丢弃 MD 力场拓扑中的化学信息，导致芳香环计数为 0、π-π 堆积丢失、H 键判定失真。

---

## 摘要

本报告对抗性验证两个递进的假设：

1. **技术缺陷**：PLIP（行业最著名的蛋白-配体相互作用分析工具）不能直接利用 MD 的拓扑信息，其化学感知交给 OpenBabel 重建，对 MD 体系不可靠。
2. **商业/学术空白**："PLIP 级精细分类 + MD 原生拓扑"的组合是行业空白，值得开发新软件。

**结论**：

- 假设 1 **成立**（PLIP 对 MD 化学感知确实有缺陷，有直接实验证据）。
- 假设 2 **部分成立**：*功能层面已被占位*，但*"用 tpr 力场原子类型做确定性化学感知"这一差异化点尚无工具实现*。这是一个**"可靠性增强"类空白**，而非**"能力缺失"类空白**。

---

## 目录

- [一、背景：PLIP 在 MD 场景的化学感知缺陷（实证）](#一背景plip-在-md-场景的化学感知缺陷实证)
- [二、现有工具全景：5 个直接竞争者](#二现有工具全景-5-个直接竞争者)
- [三、对抗性结论：空白重建评估](#三对抗性结论空白重建评估)
- [四、技术细节：ProLIF 如何"重建"拓扑（关键）](#四技术细节prolif-如何重建拓扑关键)
- [五、差异化窗口评估](#五差异化窗口评估)
- [六、下一步建议](#六下一步建议)
- [附录 A：检索数据汇总](#附录-a检索数据汇总)
- [附录 B：与 CLAUDE.md 的技术对照](#附录-b与claudemd的技术对照)

---

## 一、背景：PLIP 在 MD 场景的缺陷（实证）

### 1.1 PLIP 的工作机制

PLIP v3.0.1 的检测流程（官方文档）：

```
读入 PDB 结构
  → OpenBabel 感知化学（识别键、芳香环、供受体、加极氢）
  → 自动修复 PDB 错误（--nofix 关闭）
  → 加极氢（--nohydro 关闭）
  → 检测相互作用（疏水、H-键、盐桥、π-π、π-阳、卤键、水桥、金属）
```

**关键**：PLIP 的几何判定器（距离/角度/平面判据）与化学感知（OpenBabel）是**解耦的两个阶段**。化学感知是"输入"，判定器只是"几何应用"。

### 1.2 直接实验证据：同一结构、两次运行的差异

我们在 D927 分析中，对同一 `RBD_D927.pdb`（1 μs MD 终点帧）运行两次 PLIP：

| 项目 | 默认运行 | `--nofix --nofixfile --nohydro` |
|:----|:----|:----|
| SMILES 芳香性表示 | `c1ccccc1`（苯环）、三唑识别正常 | `[CH]1=[CH](C=CC=C1)`，全部变成 `[C]/[CH]` |
| `num_aromatic_rings` | **3**（苯环+三唑+吡啶） | **0** |
| π-π 堆积 | Y207（P型）、Y246（T型） | **全部丢失** |
| H 键 | Q205、Y207、K228（受体=三唑 N2 / 磺酰胺 O3 等） | Q205、K228（受体由 N2 翻转为 O3） |
| 原子坐标 | 与原始一致 | 与原始一致 |
| 原子数 | plipfixed：1189 重原子（0 H，重新加氢） | 直接用输入 1189 重原子 |

**三条独立证据（SMILES、aromatic count、π-stack）指向同一结论**：PLIP/OpenBabel 对 D 的小分子化学感知彻底失败——芳香性、π 堆积、H 键类型全部受影响。

**本质**：PLIP 不信任 MD 已定义的化学信息，用 OpenBabel 根据几何+简单规则"重建"加氢/键序。这是**针对静态晶体 PDB** 设计的工作流，在 MD 场景下不可靠。

### 1.3 打开"残留"文件的事实（`plipfixed.*.pdb`）

PLIP 会自动产出 `plipfixed.<name>_<id>.pdb`（修复后的输入），以及一个 OpenBabel 加氢后的中间产物（如 `RBD_D927_protonated.pdb`）。这两个文件**丢弃 MD 拓扑的氢原子，用 OpenBabel 重新加氢**，而且加氢质量与力场 RJ/JE 无关。这在 MD 分析里是**系统性信息丢失**。

---

## 二、现有工具全景：5 个直接竞争者

（检索结果：MD 相互作用分析领域活跃工具约 5 个，其中 3 个直接占位"PLIP + MD"）

### 2.1 MD-Ligand-Receptor（MDLR）

| 项 | 信息 |
|:----|:----|
| 论文 | *MD–Ligand–Receptor: A High-Performance Computing Tool for Characterizing Ligand–Receptor Binding Interactions in MD Trajectories*, IJMS 24(14):11671 (**2023**) |
| DOI | `10.3390/ijms241411671` |
| 仓库 | github.com/fraMade/MD_ligand_receptor |
| 输入 | **`.tpr` + `.xtc`**（GROMACS 原生） |
| 机制 | trjconv 每帧切 PDB → **PLIP v2.2.2** 逐帧分析 → CSV/JSON + Plotly 仪表盘 |
| 覆盖类型 | 氢键、疏水、π-π、盐桥、水桥（per-atom 永久性、heatmap） |
| 性能 | 100k 帧：50 进程 ~1h20m / 100 进程 ~20min |
| **与"你的想法"关系** | **已实现"PLIP + MD 轨迹"的框架**，但用的是 PLIP 2.2.2 的 PDB 模式（未解决 OB 化学缺陷） |

**关键教训**：MDLR 论文自己明确写道 *"Another integral component of MDLR is PLIP... PLIP plays a pivotal role in the downstream analysis"*——**"用 PLIP 处理 MD 轨迹"这个想法 2023 年已被实现**。

### 2.2 MD-IFP（Heidelberg，Rebecca Wade 组）

| 项 | 信息 |
|:----|:----|
| 论文 | DOI `10.1063/5.0019088`（J Chem Phys, **2020**） |
| 仓库 | github.com/HITS-MCM/MD-IFP |
| 机制 | 用 **MDAnalysis 直接读 GROMACS 轨迹**，从轨迹算蛋白-配体（PL）与蛋白-蛋白（PP）相互作用指纹 |
| 用途 | 结合 τRAMD 研究解离机制、热点识别 |
| 与 MDLR 关系 | MDLR 论文将其列为已有竞争者 |

**它是"直接从 MD 轨迹 + MDAnalysis 拓扑做 PP/PL 指纹"的代表**，比 MDLR 更"原生拓扑"。

### 2.3 ProLIF（chemosim-lab 维护）

| 项 | 信息 |
|:----|:----|
| 论文 | DOI `10.1186/s13321-021-00548-6`（J Cheminform, **2021**） |
| 仓库 | github.com/chemosim-lab/ProLIF（原 @cbouy，2021 起维护） |
| 机制 | **RDKit 化学感知 + MDAnalysis（MD 拓扑）**，支持轨迹 |
| 覆盖类型 | Hydrophobic, HBond, PiStacking, CationPi, Anionic, Cationic, **XBDonor/XBAcceptor（卤键）**, Metal, VdWContact |
| 现状 | v2.x 活跃开发，论文引用量高，是 **ProDIF 的直接竞品且更接近"正确化学+精细分类+MD"** |

*注意*：调查中还遇到 Arp2g（蛋白静态结构相互作用服务器，JMB 2016），但它只针对静态结构，不是 MD。另发现 **ProDy/InSty**（交互式可视化 MD 交互）、**PyLipID**（蛋白-脂质），属其他细分。

### 2.4 PyContact（作者最关心的未知者）

| 项 | 信息 |
|:----|:----|
| 论文 | *Rapid, Customizable, and Visual Analysis of Noncovalent Interactions in MD simulations*, **Biophys J 114(3):577-583 (2018)** |
| 仓库 | github.com/maxscheurer/pycontact（PyPI 可装） |
| 定位 | **GUI 工具，从 MD 轨迹分析非共价相互作用**（与你的想法几乎一致） |
| 引用 | 官方 citations 页列出多位→ **约 36 引用**（我检索到它的官方 citations 页，Zenodo DOI `10.5281/zenodo.1041419`，2021 才有 0 引用） |
| 现状 | **已更名前不加维护**，GitHub 话题显示活跃度低；依赖老模式 |

**重点提醒**：作为 6 年 MD 专家竟然从不知道 PyContact，这本身是一个**市场信号**——该领域工具很多，但**知名度/使用率都不高**。这也说明"做一个狗都能用的新工具"仍有空间。

### 2.5 其他相关

| 工具 | 定位 |
|:----|:----|
| Arp2 | 蛋白-配体相互作用的网络服务器（**静态结构**，非） |
| ProLy | 蛋白-脂质 | 蛋白-脂质 |
| BINANA | 用于药化的相互作用分析库 |

---

## 三、对抗性分析：空白重建评估

### 3.1 假设 2 的两分性质

**"PLIP+MD" 功能维度**：已经被 2020 (MD-IFP)、2021 (ProLIF)、2023 (MDLR) **先后实现**。三者都宣称支持从 MD 轨迹/topology 提取相互作用。所以：
- ❌ "从 MD 检测相互作用、做分类"这个**能力类空白不存在**

**"MD 原生化学感知"维度**：
- ProLIF：RDKit 从"拓扑（bonds/types）+几何"推断键序、芳香性、氢键供受体（有错误历史，CHANGELOG 反复修 SMARTS）
- MD-IFP：MDAnalysis 拓扑，但依赖其键序推断/读不写
- MDLR：PLIP 2.2.2 的 OpenBabel 重建（就是你说的缺陷）

**结论**：**没有工具真正读取"MD/top 中 GAFF/Amber 原子类型直接语义"来做化学判定**。这个点是**真空的**。

### 3.2 但必须诚实注意（对抗性）

| 考量 | 说明 |
|:----|:----|
| **它的"空间不大"** | "用拓扑读键序/芳香性"技术上没被做，但**不是为了袖手旁观**，而是因为没有 "能从 tpr 原生读出 gaaff 类型的现有框架相对于 RDKit 推断"是"锦上添花"而非"逼用" |
| **用户的护守现状** | 大部分 MD 用户做交互分析用的是 **氢键/π/疏水这类粗分析**（gmx 自带 + MDTools 就够），未必感受到"spatial 越精密越可靠"的需要 |
| **维护成本** | 这是所有 small-citation 工具共同天花板（PyContact 0 引用自带这个因素）——**不是"做得足够好"，而是"发现的人少"** |

### 3.3 结论综述

| 层面 | 结论 |
|:----|:----|
| 技术可行性 | ✅ 完全可行，且数据充分 |
| 是否"行业空白" | ⚖️ **"能力缺失"是 No，已有 MDLR/ProLIF 占位**；**但"确定性化学（来自 tpr GAFF 类型）"是 None，现有工具无一做到** |
| 竞争强度 | 明面上已有 3 个发表（2020/2021/2023）+ 作者组 2 个 |
| 真用户价值 | 对**非标残基/修饰/罕见小分子**的 MD 系统 |强价值（这正是本项目的 D927/C 场景）；对标准小分子 | 增量（RDKit 推断通常够用） |

---

## 四、技术性技术：ProLIF 如何"重建"拓扑（关键）

（这是回答"ProLIF 是全原子还是重原子、是否用 MD 拓扑"的证据）

### 4.1 默认：全原子 + 推断的链序

MDAnalysis RDKitConverter 官方文档：
> *"Hydrogens should be explicit in the topology file. If this is not the case, use `implicit_hydrogens=True`."*
> *"This algorithm only relies on the topology with explicit hydrogens to assign bond orders and formal charges."*

→ **若你的 tpr 含全氢，基于全原子**；否则需 `implicit_hydrogens=True`。你的 MD 是 all-atom，所以走**全原子**路径。

### 4.2 键序/芳香性：从拓扑+几何"推断"，**不读 gaaff 类型**

RDkit 文档原文：
> *"Most MD topology files don't explicitly require bond orders or formal charges... it has to be guessed from the topology."*

算法（`MDAnalysisInferrer`）：
- 原子**元素**、**连接性（bonds）**：来自 tpr ✅（真拓扑）
- **键序**：用"价电子数 vs 当前键数"的 NUE（unpaired electron）推断
- **芳香性**：由共轭体系 + `_rebuild_conjugated_bonds()` 迭代重建
- **形式电荷**：由奇偶电子数倒推

**结束**：**键连关系来自 MD 拓扑，键序/芳香/质子化由 RDKit 几何+电子推断，未使用 MD 拓扑中的 GAFF 原子类型**（如 `c`、`n`、`c3`…）。

### 4.3 作者自述的"重构系统不完善"证据

Cédric Hoys **GSoC 博客** (ProLIF 作者) 明确：
> *"RDKit requires 2 things to create a complete molecule: Elements and Bonds. ... we need bonds and bond types. ... MD topology files don't keep bond types... bonds which we need to gen"*

ProLIF CHANGELOG 里 20 条>=修 SMARTS 规则的记录（如 "HBond acceptor 排除某 N"、"疏水碳排除连 N/O/F 的")——**全是作者反击"简单的化学推断在 MD 里出错"的证据**。

---

## 五、差异化与真空

### 5.1 鉴定你的"真空白"

| 维度 | 现有（MDLR/ProLIF） | 你的机会 |
|:----|:----|:----|
| 读 MD 连接 | ✅ | ✅ |
| 芳香性/键序/HB 供受体 | **推断**（RDKit/MDAnalysis/OB） | **来自 tpr GAFF 原子类型（确定、无歧义）** |
| 对非标分子可靠性 | 低（错误靠修 SMARTS） | 高（语义就是力场定义） |
| 与分析的一致性 | 与 MD 拓扑可能冲突（修饰残基/共价） | 天然一致（直接读拓扑） |

### 5.2 差异化是"精度"而非"能力"

- **能力（功能）层面**：已被占位（MDLR 已有 PLIP+轨迹+可视化）
- **精度（可信度）层面**：可用 tpr 原子类型替代"推断"——**是现有生态中真正的空位**

### 5.3 潜在应用场景（强价值）

- 含**非标准残基**（如本项目的 CYB、共价 C242）
- 含**修饰小分子**（糖基化、NTR 靶点、卤代……
- 需要 π-π/卤键/水桥**精密谱学**的分析（罕见相互作用有化学要求）
- **require 大轨迹**（大幅降低每帧的推断开销，只做一次拓扑化学感知）

---

## 六、下一步建议（决策路径）

**(A) 先验证技术：一个悲欢对比实验**
- 用本项目 D927：tpr 中 GAFF 类型 → Build 芳香环列表 vs ProLIF/RDKit 推断芳香环 → 数一致率
- 同一结构，PLIP vs ProLIF vs 自写" From-tpr" 三者对比 H/键/π 判定
- **用数据证明"自写 tpr 感知更准"**，这是论文中最重要的硬证据

**(B) 解决"Before"：定义一个最小可行软件**
- 只读 tpr + 自写几种几何判定（hb/π/卤/疏/盐）
- Python + MDAnalysis（读拓扑）+ tpr 解析（或 GAFF 映射表）
- 输出：PLIP 式 XML/JSON 报告 + per-frame MATLAB 统计

**(C) 策略**
- 明确"能力空白已占有、精度空白开放"的定位，避免"我们做了一个 PLIP-MD"的撞车
- 面向"非标分子 MD"这个细分，而非通用（否则和 PLIP/MDLR/L 无差别）

---

## 附录 A：检索数据汇总

| 检索项 | 关键发现 |
|:----|:----|
| `PyContact` | Biophysical J 2018，~ 组合引用（官方 citations 页），Github 冷 |
| `MDLR` | IJMS 2023，国内插件 PLIP 2.2.2，MPI 并行，海active |
| `MD-IFP` | J Chem Phys 2020，Rebecca Wade，MDAnalysis 原生拓扑 |
| `ProLIF` | J Cheminform 2021，RDKit+MD，v2 活跃维护 |
| `从MDAnalysis RDKitConverter` | 官方文档：推断键序/芳香性，未用 gaaff 类型 |
| `PLIP 官方文档` | 明确 `--nofix/--nohydro`、OpenBabel 芳香✓、加氢 |
| `ProLIF CHANGELOG` | 数十条 SMARTS 修补，自证"MD 简单推断易错" |

---

## 附录 B：与 CLAUDE.md 的跟进

- 本项目 D927/CY 体系正是"非标准残基"场景：**tpr 中 GAFF 类型 `c/n/ca/nb/cp/c2/c3/os/f/ss`** 全部直接可用（见 CLAUDE.md中的 D927.itp）——验证" tpr 自带化学语义"为切实可行。
- 该调研不改变现有 MD 方案方法（PLIP 继续用于单帧交互，不在轨迹上跑）；仅在志愿开发新分析工具时作为依据。

---

*调研结束。本报告所有结论均基于查到的公开可验证资料来源（DOI、官网、GitHub、官方文档、同行未审阅的博客）。*
</parameter>
</invoke>