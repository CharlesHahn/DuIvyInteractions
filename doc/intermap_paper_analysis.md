# InterMap 论文拆解（实验设计 / 主要创新点 / 主要结果）

> 创建日期：2026-09-28
> 来源：bioRxiv 预印本《InterMap: Accelerated Detection of Interaction Fingerprints on Large-Scale Molecular Ensembles》
> DOI：<https://doi.org/10.64898/2025.12.15.694195>；全文：<https://www.biorxiv.org/content/10.64898/2025.12.15.694195v1.full-text>
> 作者：Fajardo-Díaz, Bignon, Dehez, Karami, González-Alemán（Université Paris Cité / CNRS LORIA / Université de Lorraine）
> 代码：<https://github.com/Delta-Research-Team/intermap>；文档：<https://delta-research-team.github.io/intermap/>
> 状态：**论文拆解记录**（非 agent 开发指令；直接竞品的性能对标参照系）

---

## 1. 一句话总结

InterMap 是一个面向**大规模 MD 轨迹**的 IFP（相互作用指纹）检测 Python 包，用 **k-d tree + Numba + bitarray 位压缩**把 IFP 计算提速 2–3 个数量级、内存降低至 1/91，同时检出更多相互作用；定位是"任意格式、SMARTS 灵活定义、力场无关"。

---

## 2. 研究动机与故事线（叙事）

### 2.1 时代背景
-  GPU 与算法进步 → MD 可跑更大体系、更长时间尺度 → 数据生产超越后处理能力；
-  FAIR 轨迹仓库（MDDB / MDRepo / DynaRepo）沉淀毫秒级、上百体系数据 → 需要规模化 IFP 后处理；
-  IFP 价值被下游放大：聚类亚态、别构分析，以及作为 ML/DL 输入特征（结合亲和预测、pose ranking、虚拟筛选）。

### 2.2 对手诊断（现有 MD 轨迹 IFP 工具的四项结构性局限）
| 工具 | 机制 | 被指出的局限 |
|:-----|:-----|:-------------|
| getContacts | VMD 脚本 + 严格原子命名约定 | **跨力场不兼容**；部分相互作用未定义在配体上 |
| MD-IFP | RDKit + 距离/角度 | 只能 ligand-protein / proteinA-proteinB；**不能分子内、不能核酸** |
| MD-LR（MD-Ligand-Receptor） | 并行执行 PLIP | **锁死 GROMACS 生态**；配体须单分子组（小分子限定）；继承 PLIP 不能做分子内 IFP |
| ProLIF | SMARTS + RDKit | 能解决上述功能局限，但**原子选择较多时性能崩溃** |

共同死穴：**长轨迹、复杂体系、分子内 IFP**（分子内因候选对爆炸被普遍施加"近似不可行"的性能惩罚）。

### 2.3 技术论点
-  IFP 检测瓶颈 = 两两距离计算 O(N²)；
-  k-d tree（1975, Bentley）可将邻居搜索降到 O(logN)，**"据我们所知，在 IFP 计算上是被低估的选择"**——全篇立论点。

---

## 3. 主要创新点（四个技术支柱）

1. **k-d tree 砍距离查询（核心）**：
   - 用 numba-kdtree 构建/查询（Numba 生态无原生 k-d tree）；
   - **每帧建两棵树**：一棵"芳香"、一棵"非芳香"——非芳香原子多但 cutoff 短，芳香原子少但 cutoff 长，分开建树以匹配 k-d tree 查询特性（cutoff 越大、邻居越多，查询性能越差）。
2. **Numba JIT + 多核**：关键循环编译为机器码，多核执行且无明显内存开销。
3. **bitarray 位压缩存储**：IFP 编码为位向量（对比 NumPy/SciPy 最小 byte 数组省 8 倍），稀疏指纹再压缩 → 内存友好、适合大规模。
4. **MDAnalysis + SMARTS 集成**：任意格式拓扑/轨迹、丰富选择语法、SMARTS 定义相互作用 → **与原子/残基命名约定无关**（力场无关卖点）。

### 设计亮点（易用性）
- Python module / CLI 双入口；多副本（replica）一次提交；
- **atom / residue 双分辨率**（细粒度 vs 高压缩）；
- 支持子集相互作用计算；拓扑自定义注释（对非标准链/片段的命名）；
- 本地 Shiny 可视化 **InterVis**（不依赖联网 Web 服务）。

---

## 4. 工作流（Materials & Methods）

1. **输入解析**：目录层级 + 输出文件（IFP 二进制 + 配置文件 + LOG）；
2. **读取拓扑/轨迹**：MDAnalysis（要求可识别的文件扩展名）；
3. **两个用户自定义选择**（selection）定义互作两侧，SMARTS 查询识别可参与各相互作用的原子（Table 2，默认几何 cutoff 可自定义）；
4. **分块处理**：轨迹按 chunk 载入 RAM；
5. **逐帧检测**：每帧构建芳香/非芳香两棵 k-d tree → 邻居搜索 + 几何判定 → bitarray 记帧；
6. **输出**：CLI 存 `.PICKLE`，Python 模式返回压缩位数组字典；InterVis 读取配置文件 + pickle 做交互可视化。

---

## 5. 实验设计（做了哪些验证、如何做）

### 5.1 一致性/正确性验证（IFP recovery，InterMap vs ProLIF）
-  **数据**：MPRO（SARS-CoV-2 主蛋白酶）轨迹，等间距取 **1000 帧**，全部**分子内**相互作用；
-  **方法**：双工具默认 residue 模式，逐相互作用类型做帧级比对；Fig 3 三组箱线分布：两者都检出（绿）/ 仅 InterMap（橙）/ 仅 ProLIF（红）；
-  **结论**：π-cation、cation-π、anionic、cationic 四类**检出完全一致**；其余差异可合理解释（见 §6.1）。

### 5.2 性能基准（Table 4，核心实验）
-  **对比软件及锁定版本**：getContacts（master `da14deb`）、MD-IFP v1.1、MD-LR v1.0、ProLIF v2.0.0；
-  **数据集**：9 条真实 MD 轨迹——MPRO、IgG3-M1（免疫球蛋白 G3 × 链球菌 M1 蛋白）、p53、NSP13（SARS-CoV-2）、Ec T4P（大肠杆菌 IV 型菌毛）、Spike-open（Spike 三聚体开放态）、ACE2-RBD、Topoiso I-CPT（DNA 拓扑异构酶 I × 喜树碱）、Nucleosome（147 bp DNA × 组蛋白八聚体）；覆盖 **protein-protein / protein-DNA / protein-ligand / 大溶剂化复合物**；
  -  其中 Topoiso I-CPT 与 Nucleosome **保留全水盒子**，专门评估水-溶质相互作用；
-  **硬件环境**：六核个人电脑、64 GB DDR4、AMD Ryzen 5 @3.2 GHz 超线程、Xubuntu 22.04 LTS；
-  **指标**：墙钟时间 + 峰值内存（`/usr/bin/time -v`）+ 检出相互作用数（#inters，绝对值和相对 InterMap 的百分比）；getContacts 内存用**自研子进程采集脚本**度量（`time -v` 抓不到其子进程峰值）；
-  **分辨率**：residue 与 atom 双模式分别测。

### 5.3 可扩展性消融实验（Figure 4）
-  **数据**：MPRO 单轨迹 **50,000 帧**（最苛刻场景：分子内 IFP 候选对极多），限用轨迹一部分，但趋势外推一致；
-  **变量**：chunk size（一次载入 RAM 的帧数）× 处理器数，双变量网格；
-  **结论**：chunk=1 性能差；residue 模式显著快于 atom 模式；chunk 增大 → 加速但内存涨；处理器数增加收益保守（存储是瓶颈，不是检测）；推荐 **chunk size ≈ 处理器数**。

### 5.4 可视化功能展示（InterVis）
-  仅需配置 + pickle 输出即可启动本地 Shiny 应用；
-  过滤维度：相互作用类型 / 拓扑注释 / prevalence 阈值 / MDAnalysis 原子选择；
-  5 种交互图：2D 热图、双侧 prevalence 柱状图、相互作用寿命箱线图、逐帧时间序列、残基-残基网络图（边宽 ∝ prevalence）；支持 PNG 导出与 CSV 数据导出。

---

## 6. 主要结果（定量汇总）

### 6.1 一致性验证结果
-  π-cation / cation-π / anionic / cationic：**双工具完全一致**；
-  仅 InterMap 检出的 vdW/疏水接触：距离落在 cutoff 内、化学合理（ProLIF 漏报）；
-  **ProLIF 偶发把肽键 N 判为 H 键受体**（酰胺 N 因共振无孤对，不应做受体）——对 ProLIF 化学定义缺陷的直接批评；
-  π-stacking 分类差异源于环法向量构建细节；持久相互作用大部分同帧检出；
-  综述："差异是 minor 的，不影响任一工具所得科学结论"。

### 6.2 性能基准结果（正文给出关键数字）
| 软件 | 时间 | 内存 | 检出量 |
|:-----|:-----|:-----|:-------|
| **InterMap** | 3 min（IgG3-M1）~ 91 min（Nucleosome），全轨迹全部完成 | 稳定 0.6–2.6 GB（多数体系） | 最多（基准基准值） |
| **ProLIF** | MPRO：35 min 仅完成 1%，外推 ~3500 min；IgG3-M1 3789 min；p53 8000 min；NSP13 72000 min；**Spike/Nucleosome 加载一天无法启动** | 反复超 64 GB 崩溃（MPRO 即 64 GB = InterMap 0.7 GB 的 91 倍） | 低于 InterMap |
| **getContacts** | MPRO 29 min → Nucleosome 882 min（比 InterMap 慢 4–10 倍） | 3.5–42.5 GB（Nucleosome 除外均高于 InterMap） | **MPRO 130.8M = InterMap 的 4.9%；Ec T4P 66.49M = 1.2%** |
| **MD-LR** | Topoiso I-CPT：43 min（作者自带轨迹）vs InterMap 7 min | 0.2 GB | **仅 InterMap 的 0.6%** |
| **MD-IFP** | 无法处理任何一条基准轨迹（定义所限），未入表 | — | — |

**总体结论**：InterMap = 通用性（双分辨率、异质体系兼容）+ 2–3 个数量级提速 + 低而可预测的内存占用 + 检出更多相互作用。

### 6.3 主要工程取舍（作者自述的局限）
-  三步流程（载帧 / 检测 / 存储）中**只有检测步并行化**，载帧与存储是瓶颈；
-  chunk size 与核数增大 → 内存同步增长（多份数据结构副本），内存受限时需降 chunk 与核数；
-  atom 模式比 residue 模式慢且费内存（更细粒度的存储代价）。

---

## 7. 对 DII 的启示

1. **性能路线正确性背书**：k-d tree 砍距离查询 + 稀疏/压缩存储 = 行业认可的 IO 优化方向。DII 的 TwoPass + KDTree（水桥 65h→~5s）与之一致，"性能"话语可直接引用本文数字。
2. **化学缺陷的外部佐证**：ProLIF 肽键 N 误判受体 = 推断式化学感知的实锤错误，支撑 DII"零推断"叙事（见 `doc/why_force_field_group_identification.md`）。
3. **定位互不重叠**：InterMap 快而通吃（任意格式、SMARTS、力场无关）；DII 确定而自洽（tpr 类型直读、与力场同源）。若做对比实验，InterMap 是**性能上限参照系**；它验证不了、也不打算做的"tpr 独占拓扑化学信息"正是 DII 的护城河。
4. **可借鉴点**：atom/residue 双分辨率、chunk 流水线、本地可视化——DII 的 `visualizers/`（未实现）可参考其 InterVis 交互设计。

---

## 8. 诚实声明（抓取限制）

-  Table 2（相互作用默认几何定义表）与 Table 4（完整数据表）在网页端为 AJAX 弹窗，本次无法提取全文；本拆解的定量数字均来自论文正文叙述，已尽量完整转述。
-  本文为预印本（bioRxiv "New Results"，正文注明补充材料将发表于 *Nucleic Acids Research*），结论以最终发表版为准。

---

## 9. 参考来源

- 全文：<https://www.biorxiv.org/content/10.64898/2025.12.15.694195v1.full-text>
- DOI：<https://doi.org/10.64898/2025.12.15.694195>
- 代码仓库：<https://github.com/Delta-Research-Team/intermap>
- 官方文档：<https://delta-research-team.github.io/intermap/>

---

## 10. 关联文档

- `doc/competitor_landscape_survey.md` —— 竞品全景（① 圈直接竞品含 InterMap）
- `doc/why_force_field_group_identification.md` —— 差异化论证（引用本论文的 ProLIF 缺陷证据）
- `doc/plip_prolif_group_identification_survey.md` —— PLIP/ProLIF 鉴定机制细节

---

*文档结束*