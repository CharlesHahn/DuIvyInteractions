# PLIP 与 ProLIF 基团鉴定机制调研（竞品对照）

> 创建日期：2026-09-28
> 调研方式：anysearch 检索官方文档 / 论文 / 源码并交叉验证
> 用途：论证本项目（直接读取 tpr 力场原子类型做确定性基团鉴定）的差异化定位；补充 `project_background.md` 中"竞争者已占位"一节的具体机制细节
> 状态：**调研记录**（非 agent 开发指令）

---

## 1. 调研背景

本项目核心差异化：**直接读 tpr 力场原子类型做确定性基团鉴定**，不做化学重建。为验证"现有工具无一读取力场类型语义"，需要弄清两个最接近的竞品（PLIP / ProLIF）**具体怎么做基团鉴定**。本文档记录调研结论，含可复现的来源链接。

---

## 2. ProLIF：纯 SMARTS 模式匹配

**定位**：RDKit + MDAnalysis，力场无关，支持 MD 轨迹（最接近且活跃维护的竞品）。

### 2.1 数据转换路径

1. 输入：MDAnalysis Universe 或 RDKit 分子；
2. **MDAnalysis → RDKit 转换**：MDAnalysis 的 RDKitConverter 在"所有 H 原子显式存在"的前提下，从拓扑**推断键级与形式电荷**，重建出 RDKit 分子；
3. **按残基切分**：RDKit 母分子按残基名/残基号/链自动碎片化为逐残基的 RDKit 分子，便于逐残基编码相互作用。

### 2.2 基团鉴定 = 每个相互作用类内置 SMARTS 查询

论文 Table 3（当前源码 `prolif/interactions/interactions.py` 与之对应，注明部分 SMARTS 灵感来自 **Pharmit** 与 **RDKit BaseFeatures.fdef**）：

| 基团 | SMARTS |
|:-----|:-------|
| 阴离子 | `[−{1−}]`（形式电荷 -1） |
| 阳离子 | `[+{1−}]` |
| 芳香环 | `a1:a:a:a:a:a:1`（六元）、`a1:a:a:a:a:1`（五元） |
| H 键受体 | `[N,O,F,−{1−};! + {1−}]` |
| H 键供体 | `[#7,#8,#16][H]`（N/O/S–H） |
| 卤键受体 | `[#7,#8,P,S,Se,Te,a;! + {1−}][*]` |
| 卤键供体 | `[#6,#7,Si,F,Cl,Br,I]-[Cl,Br,I,At]` |
| 金属 | `[Ca,Cd,Co,Cu,Fe,Mg,Mn,Ni,Zn]` |
| 疏水 | `[#6,#16,F,Cl,Br,I,At; + 0]` |

### 2.3 关键点

- 每个相互作用是一个 Python 类，`detect()` 对两个残基 RDKit 分子做 SMARTS 子结构匹配 → 命中原子对 → 再套几何阈值（论文 Table 2：如 H 键 D–A ≤ 3.5 Å、角度 130°–180°；π-堆积环心 ≤ 6.0 Å 等）；
- 论文明确表述：SMARTS **"比依赖元素/原子量更精确，比依赖力场特定原子类型更普适"**——这是与本项目路线最直接的对立声明；
- 要求显式 H（v2.2+ 提供 `implicit_hydrogens` 模式）；
- **每帧都要跑一次 RDKit 转换 + SMARTS 匹配**（有并行优化，但化学信息仍是从坐标/拓扑重建）。

---

## 3. PLIP：规则系统 + OpenBabel 感知（混合路线）

**定位**：为静态晶体 PDB 设计、面向蛋白-配体复合物的规则系统（PLIP 2021 扩展到 DNA/RNA）。

### 3.1 数据转换路径

- 输入：PDB 文件；OpenBabel **加极性 H**（"非确定性"加氢可用 `--nohydro` 关闭）、去除 alt 构象/模型/位置；
- 过程分两步：**先在全结构里搜索可能参与相互作用的原子/原子团（基团鉴定），再套几何规则**。源码中鉴定在 `plip/structure/preparation.py`，几何检测在 `plip/structure/detection.py`。

### 3.2 各基团的鉴定方式（三种机制混合）

| 基团 | PLIP 的鉴定方式 |
|:-----|:---------------|
| 疏水原子 | **纯规则**：碳 + 全部邻居 ∈ {C, H}（与 DuIvy 完全一致） |
| 芳香环 | **OpenBabel**（SSSR 环感知 + 芳香性）；若 Babel 未报芳香 → 兜底**平面性检查**（环上每原子法向量间夹角 < `AROMATIC_PLANARITY`=5°） |
| H 键供体/受体 | **OpenBabel 直接识别**；卤素排除、单列 |
| 蛋白带电基团 | **残基名字典**：正电 = ARG/HIS/LYS 侧链 N；负电 = ASP/GLU 羧基（只对结合位点穷举） |
| 配体带电基团 | **官能团模式匹配**：正电 = 季铵/叔胺/锍/胍；负电 = 磷酸/磺酸/硫酸/羧酸 |
| 卤键供体 | 只在配体找：所有连在碳上的 F/Cl/Br/I |
| 卤键受体 | 只在蛋白找：连了 O/P/N/S 的 C/P/S |
| 水 | 氧在配体最大 extent + `BS_DIST_MAX` 内 |
| 金属 | >50 种离子清单；蛋白目标 = Cys(S)/His(N)/含氧侧链 + 主链 O；配体目标 = 醇/酚/羧酸/磷酸/硫醇/咪唑/吡咯/Fe-S 簇；再按配位几何（线性/三角/四面体/八面体等）拟合去多余配位原子 |

### 3.3 关键点

- **基团鉴定是 OpenBabel 化学感知（芳香、HBD/HBA）+ 规则/字典（疏水、电荷、卤键）的混合**；
- 每处理一个结构都要跑 OpenBabel；
- 对 `trjconv` 导出的 PDB（无 CONECT/键序），Babel 无法感知芳香 → 只能依赖平面性兜底——这正是 DuIvy 差异化痛点的实证来源（见 `project_background.md` §5 的 D927 对比：PLIP 判 `num_aromatic_rings=0`）。

---

## 4. 三方对照（PLIP / ProLIF / DuIvy）

| 维度 | PLIP | ProLIF | DuIvy（本项目） |
|:-----|:-----|:-------|:----------------|
| 化学信息来源 | PDB 坐标 + OpenBabel 重建 | 坐标/拓扑 + RDKit 重建（推断键级/电荷） | **tpr 力场原子类型直读** |
| 基团鉴定手段 | OpenBabel 感知 + 规则 + 残基字典 | **SMARTS 子结构匹配**（Pharmit/RDKit 启发） | 类型→特征映射 + 键合图 + 电荷/残基字典 |
| 力场耦合 | 无（力场无关） | 无（力场无关，自称"普适"） | **与力场自洽**（Amber 家族 + GAFF） |
| H 键供体 | OpenBabel 判别 | `[#7,#8,#16][H]` SMARTS | D–H 键显式存在 + q(H)>0（零推断） |
| 与帧的关系 | 每结构跑一次 OpenBabel | **每帧** RDKit 转换 + SMARTS 匹配 | 识别一次、帧无关 |
| 失败模式 | 无 CONECT 的 PDB 芳香误判（兜底平面性） | 需显式 H；无 H 拓扑需推断 | 非全原子力场（GROMOS UA）不支持 |

---

## 5. 对项目的启示（差异化机会点）

1. **对赌点**：ProLIF 论文宣称"SMARTS 比力场原子类型更普适"，本项目的反方论据是——**同一 tpr 下 100% 自洽、零推断、鉴定只做一次**；ProLIF 的 CHANGELOG 有 20+ 条 SMARTS 修补记录（自证"MD 简单推断易错"，见 `project_background.md` §1）。
2. **结构优势**：全原子显式 H 下，H 键供体（D–H 键条目）、水桥（SOL 残基名）、金属（元素+电荷）的零推断能力是 PLIP（重加氢）/ProLIF（推断键序）结构上做不到的。
3. **可借鉴的细节**：PLIP 的官能团模式（季铵/叔胺/胍/锍/磷酸/磺酸/硫酸/羧酸）与金属配位目标原子清单，本项目 `amber_ff_identifier.py` 已对照实现；PLIP 的"芳香环平面性兜底"可作为本项目 π-stacking 的可选交叉验证（`PiStackingDetector` 已有 `check_planarity` 开关）。

---

## 6. 参考来源

- PLIP 官方算法文档（Detection of possible interacting groups 一节）：<https://github.com/pharmai/plip/blob/master/DOCUMENTATION.md>
- PLIP 论文（NAR 2015）：<https://pmc.ncbi.nlm.nih.gov/articles/PMC4489249/>
- PLIP 源码——基团鉴定（preparation.py）：<https://github.com/pharmai/plip/blob/master/plip/structure/preparation.py>
- PLIP 源码——相互作用检测（detection.py）：<https://github.com/pharmai/plip/blob/master/plip/structure/detection.py>
- ProLIF 论文（J Cheminform 2021，含 Table 2/3 阈值与 SMARTS）：<https://pmc.ncbi.nlm.nih.gov/articles/PMC8466659/>
- ProLIF 源码——相互作用类与 SMARTS（interactions.py）：<https://github.com/chemosim-lab/ProLIF/blob/master/prolif/interactions/interactions.py>
- ProLIF 官方文档（Interaction fingerprint）：<https://prolif.readthedocs.io/en/latest/source/modules/interaction-fingerprint.html>

---

*文档结束*
