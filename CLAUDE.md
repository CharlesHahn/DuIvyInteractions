# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## 编程哲学（项目规范）

### 核心三原则

1. **KISS（Keep It Simple, Stupid）**
   - 代码行数硬约束：函数 ≤ 20 行，类 ≤ 100 行
   - 不做"可能有用"的功能，只做当前需要的
   - 简单方案优先，复杂方案需论证

2. **显式优于隐式**
   - 数据流必须清晰可见，不隐藏复杂度
   - 不使用魔法数字，全部定义为常量
   - 类型注解必加，提高可读性

3. **可读性即文档**
   - 代码即文档，docstring 只写"是什么"，不解释"为什么"
   - "为什么"属于设计文档（doc/），不属于代码
   - 美丽优于丑陋，简单优于复杂

### Python 之禅（最相关 5 条）

```python
import this
# 1. Beautiful is better than ugly.
# 2. Simple is better than complex.
# 3. Complex is better than complicated.
# 4. Explicit is better than implicit.
# 5. Readability counts.
```

### 实践规矩

| 规矩 | 说明 |
|:----|:----|
| 函数 ≤ 20 行 | 强制拆分复杂逻辑 |
| 类 ≤ 100 行 | 单一职责 |
| 无魔法数字 | 全部定义为常量 |
| 类型注解必加 | `def func(x: int) -> str:` |
| docstring 简洁 | 只写"是什么"，不写"为什么" |
| 不重复设计文档 | 设计决策在 doc/，不在代码 |

---

## 项目：基于 MD 拓扑力场参数的相互作用判定工具

### 项目目标

构建一个**通用**的分子相互作用分析工具：
- **输入**：MD 模拟的拓扑（GROMACS tpr）＋ 轨迹（xtc）
- **输出**：PLIP 级精细分类的相互作用（π-π 堆积、H 键、盐桥、疏水、π-阳离子、卤键、水桥、金属配位等）
- **核心创新**：**直接读取 tpr 力场参数编码的化学语义**做基团鉴定，不用 OpenBabel（PLIP）/RDKit（ProLIF）重建化学

### 两段式架构

```
① 基团鉴定（本项目难点）
   tpr 力场参数（原子类型 + 键合图 + 显式 H + 电荷）
   → 确定性识别能参与相互作用的化学基团（芳香环、H 键供受体、带电基团…）
   —— 只做一次，与帧无关

② 几何判定（成熟件，可复用现有方法）
   基团带标签后，逐帧用 PLIP 式距离/角度/平面临近判据
   → 每帧相互作用列表 → 时间统计
```

### 动机与背景（调研结论，详见 plip_md_research_survey.md）

- **PLIP v3.0.1 缺陷（实证）**：通过 OpenBabel 重建加氢/芳香性，丢弃 MD 拓扑的化学信息。对 trjconv 导出的 `RBD_D927.pdb`（无 CONECT/键序）判 `num_aromatic_rings=0`、π-π 堆积全丢、H 键受体类型翻转。本质：为静态晶体 PDB 设计的工作流，对 MD 轨迹帧不可靠（早期调研记录的"默认模式 3 环"仅适用于含键序标注的 PDB）。
- **竞争者已占位**：MDLR（PLIP+轨迹，2023）、MD-IFP（MDAnalysis 原生拓扑，2020）、ProLIF（RDKit+MD，2021，最接近且活跃维护）。
- **差异定位**："能力空白"已被占位；"**用 tpr 力场原子类型做确定性化学感知**"是真正的空白——现有工具无一读取 GAFF/Amber 类型语义做化学判定，全部是推断（ProLIF 的 CHANGELOG 有 20+ 条 SMARTS 修补记录，自证"MD 简单推断易错"）。

## 第一性原理分析结论（已作批判性论证，成立）

1. **核心论据**：tpr 原子类型（GAFF 的 `ca/na/nb/os`…）是参数化时由 antechamber/sobtop 依据分子电子结构做出的**化学判决的留存记录**；RDKit/MDAnalysis 做的事是"从连接性+几何反推这个判决"。信息论上：重建 = 丢弃已编码信息 + 注入推断噪声；直接读取 = 零损失、零歧义。
2. **对"键序可能缺失"鲁棒**：即使 tpr `[bonds]` 不保留键序（全 func=1），芳香性判定仍可行——**类型名本身已编码芳香性**（`ca`=芳香碳、`na`=吡咯 N、`nb`=吡啶 N）。键序只是第二重交叉证据。
3. **全原子显式 H 是系统性优势**：H 键供体（D–H 键条目存在）、水桥（SOL 残基名）、金属（元素+电荷）全部**零推断、100% 确定**——这是 PLIP（重加氢）/ProLIF（推断键序）结构上做不到的。
4. **卖点叙事**："**确定性 + 与模拟力场自洽**"，而非"更准"——对标准药物分子 RDKit 推断往往也准，硬碰正确率是脆弱卖点；但"基团鉴定只做一次、结果与帧无关、与力场同源"是竞争工具结构上无法做到的。

## 类型映射表的来源（科学事实，两条独立证据链）

类型→化学特征的映射不是从命名"推导"的，也不是"我的知识"——而是来自两条独立的、可验证的来源：

### 来源 1：Amber 蛋白类型——从 rtp 残基定义反推

方法：**已知化学事实（TYR 的环是芳香环）→ 查 rtp 文件（TYR 环原子用 `CA` 类型）→ `CA` = 芳香**。

```
[ TYR ]                          ← 已知化学事实：酪氨酸的苯环是芳香环
 [ atoms ]                       
    CG    CA          ← CG 是环原子 → CA 是芳香碳
    CZ    C           ← CZ 是环原子 → C 也是芳香碳（但主链 C 也是 C，双重身份）
    OH    OH          ← OH 是酚羟基，不在环上 → OH 不是芳香
```

**这不是"我的知识"**，而是：
- 残基的化学（TYR 有芳香环）= 有机化学的**客观事实**
- 类型名（`CA`）= rtp 文件的**客观记录**
- 两者对照 = 可复现的**科学推导**

**自动验证器**（`verify_type_mapping.py`）用同样的逻辑检查了全部 7 个 GROMACS 自带 amber 力场 + 用户自定义力场的所有 rtp 文件，发现 16 个芳香类型 × 7 力场**零冲突**（详见"全面兼容性验证"）。

### 来源 2：GAFF 类型——从 GAFF 命名规则 + 键长交叉验证

GAFF 类型**不是**从 amber 蛋白映射推导的，而是来自两条独立的证据：

**证据 A：GAFF 命名规则（antechamber 官方文档，Wang et al. J Comput Chem 2004）**

GAFF 的类型名有系统命名规则：

| 类型 | 规则 | 含义 | 证据来源 |
|:----|:----|:----|:----|
| `c3` | c + 数字 3 | sp3 碳 | GAFF 论文/文档 |
| `c2` | c + 数字 2 | sp2 烯碳 | GAFF 论文/文档 |
| `ca` | c + a | **a = aromatic** | GAFF 论文/文档 |
| `c` | c 无后缀 | sp2 羰基碳 | GAFF 论文/文档 |
| `na` | n + a | 吡咯型芳香氮 | GAFF 论文/文档 |
| `nb` | n + b | 吡啶型芳香氮 | GAFF 论文/文档 |
| `nh` | n + h | 带 H 的吡咯氮 | GAFF 论文/文档 |
| `os` | o + s | s = single bond 醚氧 | GAFF 论文/文档 |
| `o` | o 无后缀 | sp2 羰基氧 | GAFF 论文/文档 |
| `oh` | o + h | 羟基氧 | GAFF 论文/文档 |
| `ss` | s + s | 硫醚硫 | GAFF 论文/文档 |
| `hc`/`ha` | h + c/a | 烷基氢/芳香氢 | GAFF 论文/文档 |

**证据 B：tpr 键长参数独立验证**

对 D927 的 tpr 中 `ca` 类型键长（力场参数 b0）的实测：

```
ca-ca: b0=0.13984 nm = 1.398 Å（芳香键典型值）
ca-cp: b0=0.14058 nm = 1.406 Å（芳香键）
ca-nb: b0=0.13390 nm = 1.339 Å（芳香杂环键）
ca-ss: b0=0.17806 nm = 1.781 Å（C-S 单键——不影响 ca 的芳香性）
```

1.398 Å 落在芳香键区间（1.34-1.41 Å），与单键（1.52 Å）完全不同。这是**独立于命名规则的物理证据**。

### 两条来源的关系

```
GAFF 映射：          GAFF 论文/命名规则 → 类型特征    ✅ 系统规则
                     + tpr 键长交叉验证              ✅ 独立物理证据

amber 蛋白映射：      rtp 文件 → 残基已知化学 → 类型特征  ✅ 力场定义反推
                     + 验证器自动核对（7 力场全验证）    ✅ 程序化验证

两者没有互相推导。GAFF 的 ca 标芳香，是因为命名规则说 a=aromatic 且键长支持；
amber 蛋白的 CA 标芳香，是因为 rtp 说 TYR 的环碳用 CA 且 TYR 环是芳香。
这是两条独立的证据链，结论一致，互相印证。
```

## 基团鉴定可行性矩阵（逐基团）

| 相互作用 | 所需基团特征 | tpr 中证据 | 确定性 |
|:----|:----|:----|:----|
| π-π / π-阳 / 卤-π | 芳香环 | 类型标记 `ca/na/nb/cp/cg` → 图论环检测（SSSR）→ 环内全为芳香类型；可选键序 func=4/5 与平面性交叉验证 | 高 |
| H 键供体 | D–H（D=N/O/S/F） | D–H 键条目显式存在（全原子）＋ H 部分电荷为正 | 100%（零推断） |
| H 键受体 | 孤对可用 | 类型规则表（`o/oh/os/nb/n/f`…）＋ H 计数（`nb` 无 H 必为受体，`na` 有 H 视供/受）＋ 电荷为负验证 | 高 |
| 盐桥 | 形式电荷对 | 蛋白残基名字典（LYS/ARG/ASP/GLU/HIP 编码质子化态）；配体靠类型＋H 计数 | 蛋白 100%，配体高 |
| 疏水 | 非极性 C/S/X | 类型集合（`c3/c2`…）＋无极性取代 | 高 |
| 卤键 | σ-hole 卤素 | 卤素类型（`f/cl/br/i`）＋邻接碳类型区分芳香卤/烷基卤 | 高 |
| 金属配位 | 金属中心 | 元素/类型＋电荷（本项目 Mg²⁺） | 100% |
| 水桥 | 水分子 | 残基名 SOL/HOH ＋ OW/HW 原子名 | 100% |

## 两个工程难点与对策

1. **力场类型语义映射（"通用"的代价）**：类型名无跨力场统一字典（GAFF `ca` vs CHARMM36 `CG2R61` vs OPLS `CA` 含义不同）。对策：不做"类型→官能团"扁平字典，建**特征空间映射**（类型→{杂化, 芳香性, 极性, 带H, 孤对} 特征向量），基团由特征组合而成；新力场入库 = 填特征表而非重写逻辑。首版范围：Amber 系（GAFF/GAFF2/ff14SB）共享命名体系，~50-80 个类型可控。**支持边界**：全原子力场（GROMOS 联合原子无显式 H → 供体鉴定失效，声明不支持或降级）。
2. **tpr 二进制读取**：GROMACS xdr 序列化，随版本演化。路线：MDAnalysis 垫底（原型）→ 自研解析（生产）→ `gmx dump` 文本兜底（验收基准）。原型阶段 MDAnalysis + gmx dump 交叉验证，不阻塞可行性判断。

## 已验证结论（2026-08-11 实证完成，详见 tpr_to_functional_groups.md）

✅ **从 tpr 到官能团全链路已打通**（D927 体系验证）：

1. **tpr 键序**：dump 实测无 func=4/5（conj/arom），全部 BONDS(func=1)+CONSTR(func=2)。键序靠 **b0 键长**（芳环 1.34-1.41 Å vs 单键 1.52 Å）+ 类型名联合判定。
2. **类型字典**：D927 的 GAFF 类型 `c3/c/ca/c2/cp/n/nh/nb/ss/os/f/hc/hn/ha/h4/h1` 全部在 tpr `type (N):` 段中（per-moltype）；蛋白为 amber14sb 方言 `N3/CX/2C/1C/C/O/N/CA/CW/CR/NA/NB/N2` 等。特征映射表已建（functional_groups.py）。
3. **环检测**：边删除+BFS + 稠合冗余去除。D927 找到 4 个环（3 个芳香环 + 1 个非芳香含硫环，见下），RBD 27 环（Pro×9+Tyr×11+His×3+Trp×2+Phe×1）全部化学正确。
   - **D927 芳香环 = 3 个**：6元苯环（C7-C8-C9-C10-C11-C12，全 ca）、6元含氮杂环（C15-C20-C19-C18-N17-N16，ca+cp+nb）、6元另一芳环（C24-C29-C28-C27-C26-C25）
   - **5 元含硫环（C19-S21-C22-C23-C20）判为非芳香**（正确）：S21=ss（sp3 硫醚）、C22/C23=c2（烯碳），仅 2 个 ca 强芳香 → 不满足芳香判据。b0 键长验证：C22=C23 双键（1.33 Å）、C-S 单键（1.74-1.78 Å）→ 是 2,3-二氢噻吩式稠合杂环，非芳香噻吩
4. **H 键供体**：判据 = D-H 键（**Bond+Constraint 合并**）+ **q(H)>0**。RBD 263 个供体 H 电荷全部 +0.19~+0.45，q(H)≤0 为 0。
5. **力场确认**：amber14sb 蛋白（GROMACS 移植版，fudgeQQ=0.8333 佐证）+ GAFF 配体，同属 Amber 家族。力场文件在 `../amber14sb.ff/`（不在 DIP 环境！）。

⚠️ **类型名坑**（务必记住）：
- `C` = 主链羰基碳 AND Tyr CZ 芳环碳（双身份，需"环内≥4强芳香邻居"升级）
- `N3`（蛋白正电氨基）vs `n3`（GAFF 中性氨基）——大小写不同义
- `CA` 跨力场不同义（amber=芳香碳，GROMOS=α碳）——必须特征映射表
- N-H 键在 `Constraint:` 段，不在 `Bond:` 段！供体识别必须合并两者

### 自定义力场 vs GROMACS 自带 amber 系列（兼容性结论，2026-08-11 全面实证）

**项目自定义力场**（`../amber14sb.ff/`）：所有修改都是**纯新增**（.bak 对比实证）：
- `aminoacids.rtp` 新增 `[ CYB ]` 残基（BBO 共价，GAFF 类型 ca/nb/os/ss/na/cp/n2/o/c2/n）
- `atomtypes.atp` / `ffbonded.itp` / `ffnonbonded.itp` 新增 GAFF 参数 + Sobtop GNP_neg（p5/oh）
- **未改动任何标准 amber 类型**（N/C/O/CA/CX/2C 等保持原生）——科学做法

**全面兼容性验证（7 个 GROMACS 自带 amber 力场全查）**：

| 力场 | 蛋白氨基酸类型 | 关键缺失 |
|:----|:----|:----|
| amber03.ff | 44 | 无 ✅ |
| amber94.ff | 43 | 无 ✅ |
| amber96.ff | 43 | 无 ✅ |
| amber99.ff | 43 | 无 ✅ |
| amber99sb.ff | 43 | 无 ✅ |
| amber99sb-ildn.ff | 43 | 无 ✅ |
| amberGS.ff | 43 | 无 ✅ |

**补充映射的类型**（本次为覆盖全部 amber 家族添加）：
- `CT` = sp3 碳（amber99sb-ildn 及旧版的 α/侧链碳，amber14sb 用 CX）
- `CV`/`CQ` = 核酸碱基芳香碳
- `N*` = 核酸嘌呤/嘧啶芳香氮、`NC` = 核酸氨基氮
- `P` = 核酸骨架磷酸磷
- `H0` = amber03 的 α 氢

**不影响基团鉴定的缺失类型**（全部为离子/水模型/虚拟位点）：`C0`(Ca²⁺)、`Cs`、`K`、`Li`、`Rb`、`Zn`、`Cl`、`MW`、`IB`、`OW_tip4p`、`URE`。**注意**：金属配位检测用原子序数 `a.z`（不依赖类型名），即使这些金属类型不在 TYPE_FEATURES 也能检测金属中心。

**结论**：方法对整个 Amber 家族（amber03/94/96/99/99sb/99sb-ildn/GS/14sb + GAFF 配体）**完全兼容**。类型名跨版本差异（如 CX vs CT）已全部处理；若出现新版本新类型，只需在 TYPE_FEATURES 补一条特征（设计目标即如此）。

## 常用命令（侦察/验证阶段）

```bash
# tpr dump（⚠️ 必须用 DIP 环境的 gmx，MDsoftware 的 gmx 跑不起来）
/home/hanyl/.micromamba/envs/DIP/bin/gmx dump -s <md.tpr> > dump.txt 2> dump.log

# 解析 dump + 官能团鉴定（本目录脚本）
/home/hanyl/.micromamba/envs/mamba/bin/python parse_tpr_dump.py dump_md_D927.tpr.txt
/home/hanyl/.micromamba/envs/mamba/bin/python functional_groups.py dump_md_D927.tpr.txt
```

## 相关体系速览（验证用，详细组分表见上级 KRAS_MD/CLAUDE.md）

三体系正式 MD 各 1 μs（均在 KRAS_MD/ 下）：
- `RBD_KRAS_GNP_Mg_MD/` — Control：RBD + KRAS + GNP + Mg²⁺
- `RBD_D927_KRAS_GNP_Mg_MD/` — D927：RBD + D927 + KRAS + GNP + Mg²⁺
- `RBD_BBO_KRAS_GNP_Mg_MD/` — BBO：RBD + CYB（C242 共价）+ KRAS + GNP + Mg²⁺

已知化学事实（用于基准对比的"标准答案"）：D927 含 **3 个芳香环**（苯环 + 含氮杂环 + 另一芳环）+ **1 个非芳香含硫环**（2,3-二氢噻吩式稠合杂环，C22=C23 双键，非芳香噻吩）；CYB 弹头 C=C 与 C242 SG 加成后已为 C—C 单键（sp²→sp³）。

### ⚠️ 与 PLIP 的实证对比（2026-08-11 核实）

PLIP v3.0.1 对同一 D927 结构（trjconv 导出的 PDB）判定 `num_aromatic_rings=0`、`pi_stacks` 为空：
- **根因**：PLIP 输入是 trjconv 纯坐标 PDB（无 CONECT/键序），OpenBabel 从几何重建化学失败
- **本方案**：从 tpr 直接识别 3 个芳香环 + 1 个非芳香含硫环，与化学事实一致
- **DIP 佐证**：你们自己的 PiStacking 分析（dip_PiStacking.yaml）必须**手动指定环原子索引**（Pi_rings_Index），因为 DIP 不能自动判定芳香环——这正是本项目工具的价值
- 注：早期调研文档（plip_md_research_survey.md §1.2）记录的"默认模式 3 环"结果是 PLIP 对含键序/芳香标注 PDB 的能力，对 trjconv 导出 PDB 不适用

---

## 项目架构设计（2026-08-12）

### 目录结构

```
DuIvyInteraction/
│
├── DuIvyInteractions/             # 主包
│   ├── core/                      # 领域模型 + 接口定义
│   │   ├── __init__.py
│   │   ├── datas.py               # Group, Interaction, SystemData 数据类
│   │   ├── interfaces.py          # Reader, GroupIdentifier, InteractionDetector ABC
│   │   └── constants.py           # 元素周期表、力场类型常量
│   │
│   ├── system_readers/            # 系统数据读取器（Reader 接口实现）
│   │   ├── __init__.py
│   │   ├── gmx_tpr_reader.py      # 从 tpr 二进制解析（MDAnalysis）
│   │   └── gmx_tpr_dump_reader.py # 从 gmx dump 文本解析
│   │
│   ├── group_identifiers/         # 基团识别器（策略模式，可插拔）
│   │   ├── __init__.py
│   │   └── amber_ff_identifier.py # Amber 力场识别器
│   │
│   ├── interaction_detectors/     # 相互作用判定器（策略模式，可插拔）
│   │   ├── __init__.py
│   │   ├── hydrogen_bond_detector_*.py
│   │   ├── pi_stacking_detector_*.py
│   │   ├── saltbridge_detector_*.py
│   │   ├── hydrophobic_detector_*.py
│   │   ├── halogen_bond_detector_*.py
│   │   ├── metal_coordination_detector_*.py
│   │   ├── water_bridge_detector_*.py
│   │   └── pi_cation_detector_*.py
│   │
│   ├── io/                        # 结果文件读写（Interaction 序列化/反序列化）
│   │   └── __init__.py
│   │
│   ├── utils/                     # 通用工具（无状态、可复用）
│   │   └── __init__.py
│   │
│   └── visualizers/               # 可视化
│       └── __init__.py
│
├── Tests/                         # 单元测试（镜像 src 结构）
├── doc/                           # 中文文档
├── docs_en/                       # 英文文档
├── test_MD_case/                  # 测试数据（已 gitignore）
├── .gitignore
├── CLAUDE.md
├── README.md
└── LICENSE
```

### 架构设计原则

1. **单一职责**：每个模块只做一件事
2. **依赖方向**：高层依赖低层，不反向（identifiers/detectors → core → utils）
3. **策略模式**：识别器和判定器都是可插拔的策略
4. **可测试性**：模块可独立测试

### core vs utils 的边界

| 目录 | 本质 | 特征 |
|:----|:----|:----|
| **core/** | 领域模型 + 接口定义 | "是什么" — 数据结构、ABC、常量 |
| **utils/** | 通用工具函数 | "怎么做" — 纯函数、无状态、可复用 |

**判断标准**：如果一个模块换了项目还能用 → utils；只在本项目有意义 → core

### 为什么没有 io/ 目录

- `tpr_parser.py` 是 Amber 特定的 → 合并进 `identifiers/amber.py`
- `output.py` 是通用工具 → 放在 `utils/output.py`
- `io/` 作为目录职责不清晰，不如按性质分散

### 核心接口设计

#### 数据类（core/data.py）

```python
@dataclass
class Group:
    """一个可参与相互作用的基团"""
    group_id: int                    # 唯一标识
    group_type: str                  # "aromatic_ring", "donor", "acceptor", ...
    molecule: str                    # 所属分子名（如 "D927", "RBD_pro"）
    residue_name: str                # 残基名
    residue_id: int                  # 残基号
    atom_indices: List[int]          # 全局原子索引列表
    center: Optional[Tuple] = None   # 质心坐标（动态计算）
    normal: Optional[Tuple] = None   # 法向量（芳香环用）
    properties: dict = None          # 额外属性

@dataclass
class Interaction:
    """一个检测到的相互作用"""
    interaction_type: str            # "hydrogen_bond", "pi_stacking", ...
    group1: Group
    group2: Group
    frame: int
    time_ps: float
    distance: float
    angle: Optional[float] = None
    is_active: bool = True
```

#### 基团识别器 ABC（core/interfaces.py）

```python
class GroupIdentifier(ABC):
    @property
    @abstractmethod
    def name(self) -> str: ...

    @abstractmethod
    def identify(self, topology) -> List[Group]: ...
```

#### 相互作用判定器 ABC（core/interfaces.py）

```python
class InteractionDetector(ABC):
    @property
    @abstractmethod
    def name(self) -> str: ...

    @property
    @abstractmethod
    def required_group_types(self) -> List[str]: ...

    @abstractmethod
    def detect_frame(self, groups, coordinates, frame, time_ps) -> List[Interaction]: ...
```

### 扩展性设计

| 扩展场景 | 实现方式 |
|:----|:----|
| **新力场识别器** | 继承 `GroupIdentifier`，实现 `identify()` |
| **新相互作用类型** | 继承 `InteractionDetector`，实现 `detect_frame()` |
| **新判据** | 同一类型可有多个 Detector（如 `HBondStrict`, `HBondLoose`） |
| **新输出格式** | 在 `utils/output.py` 添加新函数 |
| **新可视化** | 在 `visualize/plotter.py` 添加新方法 |

### 使用示例

```python
from DuIvyInteractions.identifiers.amber import AmberGroupIdentifier
from DuIvyInteractions.detectors.hydrogen_bond import HydrogenBondDetector
from DuIvyInteractions.detectors.pi_stacking import PiStackingDetector
from DuIvyInteractions.pipeline import Pipeline

# 配置
identifier = AmberGroupIdentifier()
detectors = [
    HydrogenBondDetector(distance_cutoff=3.5, angle_cutoff=150),
    PiStackingDetector(distance_cutoff=5.0, angle_cutoff=30),
]

# 运行
pipeline = Pipeline(identifier, detectors)
results = pipeline.run("md.tpr", "md.xtc")

# 输出
from DuIvyInteractions.utils.output import save_results_csv
save_results_csv(results, "interactions.csv")

# 可视化
from DuIvyInteractions.visualize.plotter import InteractionPlotter
plotter = InteractionPlotter(results)
plotter.plot_timeline("hydrogen_bond")
plotter.plot_heatmap("pi_stacking")
```

### 当前代码状态（2026-08-12）

已完成第一阶段（基团鉴定）的 D927 体系验证：
- ✅ `parse_tpr_dump.py` — tpr dump 解析器（200行）
- ✅ `functional_groups.py` — 特征映射 + 环检测 + 官能团鉴定（470行）
- ✅ `verify_type_mapping.py` — 映射表自动验证器（312行）
- ✅ 已使用全局原子索引（跨分子类型唯一）

待迁移：将现有代码重构到新架构中