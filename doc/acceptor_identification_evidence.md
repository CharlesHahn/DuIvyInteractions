# H 键受体资格证据清单（N 类原子判定依据）

> 创建日期：2026-09-30
> 状态：**证据锚定记录**——本次修复（剔除带 H 的非受体 N 类型）的逐项判定依据，供追溯
> 原则：每个"剔/留"决定必须有力场定义（rtp/rtf 实测）+ 化学事实 + 文献引用；证据不足的项一律"不动 + 登记待查"

---

## 1. 判定规则（化学事实 + 文献）

受体资格取决于该环境下的 N **是否有可用孤对**：

| N 环境 | 孤对状态 | 受体？ | 证据 |
|---|---|---|---|
| 普通酰胺 N（带 H，连 C=O） | 参与酰胺共振（不可用） | **非受体** | Eildal 2013, JACS 135:12998（酰胺→酯突变，供体转弱受体）；InterMap 论文（ProLIF 肽键 N 误判批评）；化学共识 |
| 铵 N（4 键，RNH₃⁺） | 无孤对 | **非受体** | 教科书（铵四键无孤对） |
| 带 H 吡咯 N（芳香环内） | 参与芳香共轭（不可用） | **非受体** | 教科书：吡咯 N 孤对参与芳香（LibreTexts / MasterOrganicChemistry："lone pair involved in aromatic system, less available"） |
| 胍基 N（连 C(N)(N)） | 参与胍基共振 | **非受体** | 化学事实（胍基 Y-芳香共振）；检索到的胍基文献全部为"胍基作供体"（Ng 2025 等），**无任何"胍基 N 作受体"文献** |
| 吡啶型 N（无 H，芳香环内） | 孤对在 sp2 轨道（可用） | **受体** | 教科书（吡啶 N 可质子化/配位） |
| 中性胺 N（带 H，sp3） | 自由孤对 | **受体** | 教科书 + Luisi 1998（氨基接受 H 键） |
| **核酸氨基 N**（腺嘌呤 N6、鸟嘌呤 N2，带 2H） | 与芳环共轭但仍可接受 | **受体** | **Luisi et al. 1998**（J Mol Biol，被引 132："amino groups in guanine and adenine bases accept hydrogen bonds"）+ **Baik et al. 2003**（JACS，被引 398："C6 amino group of adenine has potential to act as acceptor"）——两篇独立 |
| **Pro N**（无 H，连羰基 C，非芳香环） | 碱性高于普通酰胺，可接受 N-H···N | **受体** | **Deepak et al. 2016**（Biophys J，被引 72：系统分析高分辨晶体+中子衍射结构，"clearly establishes the H-bonding capability of proline nitrogen"） |
| 嘌呤糖苷 N9（无 H，芳香环内，3 键） | 参与嘌呤芳香（推理） | **证据不足** | 无直接文献；杂环化学推理（3 键环内 N=吡咯型）。**按"证据不足不动"原则暂保留为受体候选，登记待查** |

## 2. 各力场类型映射（本次修改依据）

### Amber（amber14sb GROMACS 移植 rtp 实测）
| 类型 | 语义（rtp 证据） | 处理 |
|---|---|---|
| `N` | 主链/侧链酰胺（带 H，非受体）+ **Pro N（无 H，受体）——同一类型二义** | 保留；`_find_acceptors` 按 H 邻居数区分（带 H 排除、无 H 保留） |
| `N3` | LYS NZ 铵 | **剔除** |
| `NA` | TRP NE1 / His 带 H 态（吡咯） | **剔除** |
| `NB` | His 无 H 态 / 核酸 N7（吡啶） | 保留 |
| `NC` | 核酸环内无 H N（A N1/C N3） | 保留 |
| `N2` | 核酸氨基（A N6） | 保留（Luisi 1998; Baik 2003） |
| `N*` | 核酸糖苷 N（嘌呤 N9） | 保留（证据不足，待查） |
| `n`(GAFF) | 配体酰胺 | **剔除** |
| `n2`(GAFF) | 语义待查证 | 保留（证据不足，待查） |
| `n3`(GAFF) | 中性胺 | 保留（受体） |

### CHARMM（top_all36_prot.rtf MASS 段类型名自证）
剔 `NH1`(peptide)/`NH2`(amide)/`NH3`(ammonium)/`NC2`(guanidinium)/`NY`(pyrrole)/`NR1`,`NR3`(protonated his)；留 `N`(proline)/`NR2`(unprotonated his)+CGenFF `NG*`（待逐项查证）。

### OPLS（GROMACS oplsaa.ff aminoacids.rtp 实测）
剔 `opls_238`(主链N)/`opls_237`(酰胺N)/`opls_287`(LYSH铵)/`opls_300`,`opls_303`(Arg胍)/`opls_503`(带H吡咯)/`opls_512`(双质子化His)/`opls_749`,`opls_750`,`opls_751`(ARGN中性胍)；留 `opls_239`(Pro)/`opls_511`(无H吡啶)/`opls_900`(中性胺LYS)。

### GROMOS（类型粒度粗，结构判据）
`GROMOS_ACCEPTOR_TYPES` 无法纯类型区分；`_find_acceptors` 对 N 加结构判据：**键数≥4（铵）或带 H（酰胺/吡咯/胍）→ 排除；无 H（Pro/His 吡啶）→ 保留**。

---

## 3. 已知边界与待查项（登记）

| 项 | 状态 | 说明 |
|---|---|---|
| 嘌呤 N9（amber `N*`）受体资格 | **待查** | 无直接文献；当前按"证据不足不动"保留为受体候选 |
| GAFF `n2` 语义 | **待查** | 需查 gaff.dat 定义后定剔/留 |
| CGenFF `NG*` 逐项 H 数 | **待查** | 配体类型需逐个核对后再定 |
| GROMOS 核酸氨基 N | **已知不一致** | GROMOS 简化判据"带 H 即排除"会连带排除核酸氨基（Luisi 判为受体，amber 用 `N2` 保留）；GROMOS 版设计文档已声明配体/核酸场景受限，属已知边界 |
| 胍基 N 非受体 | **化学推理级** | 无"胍基作受体"文献，判定靠胍基共振化学事实；方向与全部胍基文献（供体角色）一致 |

---

## 4. 关联文档

- `doc/force_field_compatibility_survey.md` —— 多力场类型体系调研（本次修改的力场证据来源）
- `doc/intermap_paper_analysis.md` —— ProLIF 肽键 N 误判批评（酰胺非受体的工具层佐证）
- `doc/why_force_field_group_identification.md` —— "力场类型=化学判决记录"方法论

---

*文档结束*