# 力场类型映射

本项目直接从 GROMACS tpr 的力场原子类型判定化学基团（详见[核心概念](concepts.md)）。本文档列出**已入库的 4 个力场家族**——Amber 家族（Amber 蛋白 + GAFF 配体）、GROMOS 53A6/54A7、CHARMM36/C36m（蛋白 + CGenFF 配体）、OPLS-AA/L——的类型到化学特征的映射规则，帮助科研用户理解判定依据并评估适用范围。

所有特征表常量均定义在各识别器模块顶部：`DuIvyInteractions/group_identifiers/{amber,gromos,charmm,opls}_ff_identifier.py`。多力场兼容性调研的完整论证见 `doc/force_field_compatibility_survey.md`。

## 映射原则与特征空间设计

类型→化学特征的映射来自两条独立的、可交叉验证的证据链：

1. **GAFF 类型**：来自 antechamber 官方命名规则（Wang et al., *J Comput Chem* 2004）——类型名的后缀字母/数字编码化学特征（如 `ca` 的 `a` = aromatic）
2. **蛋白类型**：从 rtp/rtf 残基定义反推——已知残基的化学事实（如 TYR 苯环是芳香环）对照 rtp 文件中环原子使用的类型名（如 Amber 的 `CA`），推出该类型 = 芳香碳

设计上采用**特征空间映射**：每个类型映射到特征向量 {杂化状态, 芳香性, 极性, 是否带 H, 孤对可用性}，基团由特征组合判定（如"芳香环 = 环内 ≥ n−1 个强芳香原子"）。新力场入库 = 填一张特征表，无需改判定逻辑（见[扩展指南](extension.md)）。

**已入库与边界**：

| 力场家族 | 蛋白 | 配体 | 边界 |
|:---------|:-----|:-----|:-----|
| Amber（amber03/94/96/99/99sb/99sb-ildn/GS/14sb） | ✅ | GAFF/GAFF2，映射验证零冲突 | — |
| GROMOS 53A6/54A7 | ✅（类型粒度粗，结构判据/分级判定/反列举补偿） | ❌ 不支持（ATB 参数化，第一版已声明） | 类型粗（`N`/`C` 通吃多种化学环境）；联合原子力场，依赖脂肪族显式 H 的判定（供体、疏水邻接）受限 |
| CHARMM36/C36m | ✅ | CGenFF（`NG*` 受体资格逐类型判定） | — |
| OPLS-AA/L（2001） | ✅ | 类型表内 | 需 `opls_XXX` 编号映射 |

## Amber 家族特征表

常量定义于 `DuIvyInteractions/group_identifiers/amber_ff_identifier.py` 模块顶部。

### 芳香类型（STRONG_AROMATIC）

这些类型的原子**由类型名直接确定**为芳香原子，是环检测的强证据（环内 ≥ n−1 个原子属此集合即可判定芳香环，见[基团识别规则](group_rules.md)）：

| 类别 | 类型 | 说明 |
|:-----|:-----|:-----|
| GAFF 芳香碳 | `ca, cg, ch, cm, cn, cp, cq, c1` | `c` + 芳香后缀 |
| GAFF 芳香氮 | `na, nb, nh, ni, nj, n1, n2` | `na` = 吡咯型、`nb` = 吡啶型 |
| GAFF 芳香磷 | `pb` | — |
| Amber 蛋白芳香碳 | `CA, CB, CC, CK, CM, C5, C6, C7, C*, CW, CR, CN, CV, CQ` | rtp 环原子类型 |
| Amber 蛋白芳香氮 | `NA, NB, NC, N*` | — |

### 兼容类型（COMPATIBLE_TYPES）

非芳香类型，但在环内 **n−1 个原子是芳香类型**的"强制"下可参与共轭。用于处理歧义类型与杂环：

| 类型 | 说明 |
|:-----|:-----|
| `C, N` | Amber 蛋白歧义类型（主链羰基碳/芳环碳 `C`；酰胺氮/芳香氮 `N`） |
| `os, ss` | GAFF 呋喃氧 / 噻吩硫 |
| `cc, cd` | GAFF 非纯芳香共轭环碳 |
| `pc, pd` | GAFF 共轭环内 sp2 磷 |

### 受体类型（ACCEPTOR_TYPES）

H 键受体（有可用孤对）的候选类型。**注意**：此表经 2026-09-30 受体资格修复（A1）修订，剔除了带 H 的非受体氮类型（详见下文"[受体资格判定](#acceptor-qualification)"节）：

| 类别 | 类型 |
|:-----|:-----|
| GAFF 氧 | `o, o2, oh, os, oe, o1, ow` |
| GAFF 氮 | `n2, n3, nb, ni, nj, nc, ne, nf, nk`——`n2` 亚胺型（sp2 N 双取代，有孤对，受体，GAFF 官方类型表定义）、`n3` 中性胺；已剔 `n`（配体酰胺）、`na/nh`（吡咯型，孤对参与芳香） |
| GAFF 硫 | `s, ss, sh, sx, s2` |
| 卤素 | `f, cl, br, i` |
| Amber 蛋白氧 | `O, OH, O2, OS, OW` |
| Amber 蛋白氮 | `N, N2, NB, NC`——`N` 兼作普通酰胺（带 H，非受体）与 Pro N（无 H，受体），同一类型二义，在 `_find_acceptors` 中按 H 邻居数区分；`N2` 核酸氨基（腺嘌呤 N6，带 2H，受体）；`NB/NC` 无 H 吡啶/环内氮；已剔 `N3`（铵）、`NA`（带 H 吡咯）、`N*`（核酸糖苷 N，3 键吡咯型） |
| Amber 硫 | `S, SH` |

### 其他常量

- **METAL_IONS**：42 种金属元素集合（来自 PLIP config.py），见下文"金属离子"节
- **WATER_RESIDUES**：`{SOL, HOH, WAT}`（类属性，见下文"水分子"节）
- **CHARGE_THRESHOLD**：0.1（带电基团净电荷验证阈值）
- **正电残基字典**（残基名 → 原子名清单）：ARG、LYS、HIP、ORN、DAB、M3L、MLY
- **负电残基字典**：ASP、GLU、CYM、KCX、PCA、SEP、TPO、PTR

## GROMOS 特征表（53A6 / 54A7）

常量定义于 `DuIvyInteractions/group_identifiers/gromos_ff_identifier.py`（本地 gromos54a7.ff 实测，53A6/54A7 通用）。GROMOS 类型**粒度粗**（`N`/`C` 通吃多种化学环境），无法纯类型区分，因此受体、芳香环、疏水三类判定改用**结构判据/分级判定/反列举**补偿。

### 受体类型（GROMOS_ACCEPTOR_TYPES）与 N 结构判据

```text
GROMOS_ACCEPTOR_TYPES = {O, OM, OA, OE, OW,          # 氧：羰基/羧基/羟基/醚/水
                         N, NT, NL, NR, NZ, NE,     # 氮：肽胺/终端/芳香/胍基
                         S,                          # 硫
                         F, CL, BR}                  # 卤素
```

`_find_acceptors` 对元素 N 的原子加结构判据（化学事实 + 文献，详见 `doc/acceptor_identification_evidence.md`）：

- **键数 ≥ 4（铵 RNH₃⁺）→ 排除**（无孤对，非受体）
- **带 H 邻居 → 排除**（普通酰胺/侧链酰胺/带 H 吡咯/胍基，非受体）
- **无 H → 保留**（Pro N / His 无 H 吡啶型 N，受体）

### 芳香环分级判定

GROMOS 无逐类型芳香表，采用分级判定（`_filter_aromatic_rings` 覆写）：

| 集合 | 内容 | 作用 |
|:-----|:-----|:-----|
| `GROMOS_RING_TYPES` | `{C, CR1, NR}` | 环原子类型全集（环碳/环氮），环成分过滤 |
| `GROMOS_AROMATIC_STRONG` | `{NR}` | 强信号：芳香氮，环内含 NR 直接判芳香 |
| `GROMOS_AROMATIC_RESIDUES` | PHE, TYR, TRP, HISA, HISB, HISH, HIS1, HIS2 | 弱信号：全 C 环仅限标准蛋白芳香残基 |

判定顺序：环成分 ⊆ `{C, CR1, NR}` → 含 `NR` 直接判芳香 → 全 C 环回退残基名白名单。

### 疏水反列举

GROMOS 为联合原子力场（脂肪族 C 无 H 邻居），疏水判定不用 Amber 的"邻居全为 C/H"，改为**类型白名单 + 极性邻居排除**（`_find_hydrophobic` 覆写）：

- `GROMOS_HYDROPHOBIC_TYPES` = `{C, CH0, CH1, CH2, CH3, CH4, CH2r}`（**不含** CH3p 极性胆碱 N⁺、CR1 芳香/烯 sp2）
- `GROMOS_HYDROPHOBIC_EXCLUDED` = `{O, N, S}`——任一邻居属此集合即排除

### 水残基与带电残基

- **WATER_RESIDUES**：`{SOL}`（GROMOS 配 SPC/SPC-E）
- **正电残基字典**：ARG、LYSH（质子化 ε-铵 NH₃⁺）、HISH（双质子化 His）——GROMOS 把 LYS 拆为中性 LYS / 质子化 LYSH，仅后者带正电
- **负电残基字典**：ASP、GLU
- 带电基团去重覆写 `_deduplicate_charged`：删除同类型下被更完整基团严格包含的基团（GROMOS 质子化残基的正电 N 原子本身带正电荷，官能团层会对每个 N 单原子额外产出正电基团，见[基团识别规则](group_rules.md)）

## CHARMM 特征表（CHARMM36 / C36m + CGenFF）

常量定义于 `DuIvyInteractions/group_identifiers/charmm_ff_identifier.py`（与 CHARMM-GUI 官方 `top_all36_prot.rtf` / `top_all36_cgenff.rtf` 逐字一致）。

### 受体类型（CHARMM_ACCEPTOR_TYPES）

| 类别 | 保留（受体） | 剔除（非受体） |
|:-----|:-----|:-----|
| 蛋白氧 | `O, OB, OC, OH1, OS` | — |
| CGenFF 氧 | `OG2D1–5, OG2P1, OG2R50, OG301–304, OG311, OG312, OG3C51, OG3C61, OG3R60` | — |
| 蛋白氮 | `N`（Pro N，无 H，受体）、`NR2`（无 H 吡啶型 His N，受体） | `NH1`（肽键）、`NH2`（酰胺）、`NH3`（铵）、`NC2`（胍基）、`NY`（吡咯）、`NR1/NR3`（质子化 His） |
| CGenFF 氮 | `NG2D1`（中性亚胺/席夫碱）、`NG2R50`（嘌呤 N7）、`NG2R60/NG2R62`（6 元吡啶型）、`NG2S3`（环外胺/苯胺型）、`NG3N1`（肼 N，sp3 胺）——**留 6** | `NG2O1`（硝基苯 N）、`NG2P1`（质子化亚胺）、`NG2R51`（5 元单键 sp2 N=His/Trp 吡咯带 H）、`NG2R52`（质子化席夫碱/脒/胍）、`NG2R61`（6 元单键亚胺 N）、`NG2RC0`（桥头 N）、`NG2S0`（N,N-二取代酰胺）、`NG2S1`（肽键 N）、`NG2S2`（末端酰胺 N）、`NG2S4`（羟肟酸 N）——**剔 10** |
| 硫 | `S, SM, SS` | — |
| 卤素 | `FGA1–3, FGR1, CLGA1, CLGA3, CLGR1, BRGA1–3, BRGR1, IGR1`（CGenFF 命名） | — |

> CGenFF `NG*` 各类型按官方 MASS 段注释逐项判定（留 6 剔 10）；其中 `NG2S0`（N,N-二取代酰胺）虽含配体 Pro 类弱受体，但与一般叔酰胺同类型不可分，登记为已知例外（见 `doc/acceptor_identification_evidence.md` §3）。

### 芳香环类型

- **CHARMM_STRONG_AROMATIC**：蛋白芳香碳 `CA, CAI, CPH1, CPH2, CPT, CY`；蛋白芳香氮 `NR1, NR2, NR3, NY`；CGenFF 芳香碳 `CG2R51–53, CG2R57, CG2R61–64, CG2R66, CG2R67, CG2R71, CG2RC0`；CGenFF 芳香氮 `NG2R50–52, NG2R57, NG2R60–62, NG2R67, NG2RC0`；杂环 O/S `OG2R50, SG2R50`
- **CHARMM_COMPATIBLE_TYPES**：`C, CC, CD`（环内主链羰基/羧基/酰胺碳）、`CG2D1, CG2D2, CG2DC1–3`（共轭烯）、`CG2O1–6`（环内羰基）

### 水残基与带电残基

- **WATER_RESIDUES**：`{TIP3, HOH, SOL, WAT}`（CHARMM 默认 TIP3/HOH，兼容 GROMACS 常规 SOL）
- **正电残基字典**：LYS、ARG、HSP（双质子化 His——HSP 咪唑氮带负电，必须纳入咪唑碳原子 CG/CE1/CD2（类型 CPH1/CPH2）才使电荷中心为正）
- **负电残基字典**：ASP、GLU、CYM（硫醇根 Cys⁻）
- CHARMM 蛋白 N 原子带负电荷，N 端 NH₃⁺ 识别须用**末端铵结构判据**（`_is_terminal_ammonium`：N 邻居 ⊆ {C,H} 且重原子邻居不连 ≥2 个 N，电荷原子 = N + 键连 H），见[基团识别规则](group_rules.md)

## OPLS 特征表（OPLS-AA/L 2001）

常量定义于 `DuIvyInteractions/group_identifiers/opls_ff_identifier.py`（蛋白 rtp 实测 + ffnonbonded 符号映射，版本 = OPLS-AA/L 2001）。

### 受体类型（OPLS_ACCEPTOR_TYPES）

| 类别 | 保留（受体） | 剔除（非受体） |
|:-----|:-----|:-----|
| 氧 | `opls_236, opls_272, opls_154, opls_167, opls_268, opls_269`（O/O2/OH/O_3） | — |
| 氮 | `opls_239`（Pro N，无 H）、`opls_511`（无 H 吡啶型 His N）、`opls_900`（LYS 中性胺 NZ）——**留 3** | `opls_238`（主链 N）、`opls_237`（侧链酰胺 N）、`opls_287`（LYSH 铵）、`opls_300/opls_303`（Arg 胍基）、`opls_503`（带 H 吡咯：HISD ND1、TRP NE1）、`opls_512`（双质子化 His N）、`opls_749/750/751`（ARGN 中性胍）——**剔 10** |
| 硫 | `opls_202, opls_200`（S/SH） | — |
| 卤素 | 离子卤素 `opls_400–403`；有机卤素 `opls_123, 151, 164, 226, 264, 709, 719, 721, 722, 726, 728, 730, 732, 786, 956, 965` | — |

### 芳香环类型

- **OPLS_STRONG_AROMATIC**：CA 系 `opls_145, 166, 302, 752`；5 元环 `opls_500–502, 506–510, 514`；芳香氮 `opls_503, 511, 512`（NA/NB）
- **OPLS_COMPATIBLE_TYPES**：`opls_235, opls_267`（环内羰基碳）、`opls_271`（环内羧酸碳）

### 水残基与带电残基

- **WATER_RESIDUES**：`{HOH, HO4, HO5, SOL, WAT}`（HOH/SPC、HO4/TIP4P、HO5/TIP5P）
- **正电残基字典**：ARG、LYSH（质子化 ε-铵）、HISH（双质子化 His，含咪唑碳 CG/CE1/CD2）——OPLS 的 LYS 中性、LYSH 质子化
- **负电残基字典**：ASP、GLU
- OPLS 铵氮带负电荷，复用 CHARMM 的**末端铵结构判据**（`_is_terminal_ammonium`，逻辑相同、独立维护）识别 LYSH 侧链铵与 N 端 NH₃⁺

(acceptor-qualification)=
## 受体资格判定（A1 修复，2026-09-30）

H 键受体资格取决于该环境下的 N **是否有可用孤对**。修复剔除了带 H 的非受体 N 类型（普通酰胺/铵/带 H 吡咯/胍基），保留 Pro N、His 无 H 吡啶、中性胺、核酸氨基。九类 N 环境判定（完整证据清单见 `doc/acceptor_identification_evidence.md`）：

| N 环境 | 孤对状态 | 受体？ | 判定依据 |
|:-------|:---------|:-------|:---------|
| 普通酰胺 N（带 H，连 C=O） | 参与酰胺共振 | 非受体 | Eildal 2013, *JACS* 135:12998；InterMap 分析（ProLIF 肽键 N 误判批评，`doc/intermap_paper_analysis.md`） |
| 铵 N（4 键，RNH₃⁺） | 无孤对 | 非受体 | 教科书（铵四键） |
| 带 H 吡咯 N（芳香环内） | 参与芳香共轭 | 非受体 | 教科书（杂环化学） |
| 胍基 N（连 C(N)(N)） | 参与胍基共振 | 非受体 | 化学事实（检索文献仅见"胍基作供体"，无"胍基作受体"文献） |
| 吡啶型 N（无 H，芳香环内） | 孤对在 sp2 轨道 | 受体 | 教科书（吡啶可质子化/配位） |
| 中性胺 N（带 H，sp3） | 自由孤对 | 受体 | Luisi 1998, *J Mol Biol* |
| 核酸氨基 N（腺嘌呤 N6、鸟嘌呤 N2，带 2H） | 与芳环共轭仍可接受 | 受体 | Luisi 1998；Baik 2003, *JACS* |
| Pro N（无 H，连羰基 C） | 碱性高于普通酰胺 | 受体 | Deepak 2016, *Biophys J* |
| 嘌呤糖苷 N（N9/N1，3 键环内 sp2 N） | 参与嘌呤芳香（吡咯型） | 非受体 | 教科书杂环化学；碱基配对受体位为 N1/N3/N7（2026-09-30 对抗性审查修正，见证据清单 §2/§3） |

**各力场的实现方式**：

- **Amber**：`ACCEPTOR_TYPES` 表 + 类型 `N` 按 H 邻居数区分（带 H 的普通酰胺跳过、无 H 的 Pro N 保留）+ 排除集合（正电基团成员、芳香环内带 H 的 N）
- **GROMOS**：`GROMOS_ACCEPTOR_TYPES` 表 + N 结构判据（键数 ≥ 4 或带 H 排除、无 H 保留）
- **CHARMM / OPLS**：类型表已完成逐类型判定（无二义 N 类型），按表直接过滤

所有力场的受体判定还要求 **原子部分电荷 q < 0**。

**已知边界**（登记于证据清单 §3）：GROMOS 的简化结构判据"带 H 即排除"会连带排除核酸氨基（Luisi 判为受体、Amber 以 `N2` 保留），GROMOS 配体/核酸场景受限；CHARMM `NG2S0`（N,N-二取代酰胺）与配体 Pro 类弱受体不可分，剔表时丢失该类弱受体。

## 供体判定

供体原子 D（N/O/S/F）与其氢 H 之间的**键条目显式存在**，且 H 带正电荷（q(H) > 0）。tpr 中 N-H 键位于 **Constraint 段**（而非 Bond 段）；读取器将 Bond、Constraint（以及 Settle）段**合并为统一的残基键列表**（`gmx_tpr_dump_reader.py` 的 `_parse_molblock` 对三类 ilist 段统一写入 `bonds`；MDAnalysis 读取器经 `u.bonds` 同样合并）。供体识别必须依赖合并后的键图，不能只扫 Bond 段；跨残基键的 D-H 对由 `_find_inter_residue_donors` 单独补检。该判定不涉及坐标推断。

## 金属离子（METAL_IONS）

来自 PLIP config.py，共 42 种（首字母大写，与元素符号一致）：

`Ca, Co, Mg, Mn, Fe, Cu, Zn, Li, Na, K, Rb, Sr, Cs, Ba, Cr, Ni, Ru, Rh, Pd, Ag, Cd, La, W, Os, Ir, Pt, Au, Hg, Ce, Pr, Sm, Eu, Gd, Tb, Yb, Lu, Al, Ga, In, Sb, Tl, Pb`

金属中心按**元素符号**识别（由读取器经 `ATOMIC_NUMBER_TO_ELEMENT`（`core/constants.py`）推断元素），不依赖类型名。

## 水分子（WATER_RESIDUES）

`GroupIdentifier.WATER_RESIDUES` 类属性（基类默认为空 `frozenset()`，各识别器子类覆盖）是**跨力场统一入口**：pipeline 按力场取用该集合排除水，识别器内部 `_find_water`、`_find_metal_binding` 也统一经 `self.WATER_RESIDUES` 判定：

| 力场 | WATER_RESIDUES（类属性值） | 说明 |
|:-----|:---------------|:-----|
| Amber（`amber_ff_identifier.py`） | `{SOL, HOH, WAT}` | SOL 为主 |
| GROMOS（`gromos_ff_identifier.py`） | `{SOL}` | 配 SPC/SPC-E |
| CHARMM（`charmm_ff_identifier.py`） | `{TIP3, HOH, SOL, WAT}` | CHARMM 默认 TIP3/HOH |
| OPLS（`opls_ff_identifier.py`） | `{HOH, HO4, HO5, SOL, WAT}` | HOH/SPC、HO4/TIP4P、HO5/TIP5P |

因此 CHARMM TIP3、OPLS HO4/HO5 等非 `SOL` 水残基名也能被正确识别与排除（2026-10-08 修复）。

## 类型名歧义与处理

- **`C` 双身份**：Amber 中 `C` 既作主链羰基碳，也作 Tyr CZ 芳环碳。环内芳香原子数 ≥ n−1（六元环为 5）时，`C` 作为兼容类型被"升级"参与共轭。
- **`N3` vs `n3`**：大写 `N3` = Amber 蛋白正电氨基（铵，非受体），小写 `n3` = GAFF 中性胺（受体）——大小写不同义。
- **`CA` 跨力场不同义**：在 Amber 类型体系中 `CA` 是芳香碳；在 GROMOS 体系中 `CA` 指 α 碳。必须通过各力场独立的特征表区分，不能跨力场复用类型名。
- **N-H 键在 Constraint 段**：供体识别必须合并 Bond + Constraint（+ Settle）三类键。

## 兼容性范围

- **Amber 家族**：amber03、amber94、amber96、amber99、amber99sb、amber99sb-ildn、amberGS、amber14sb + GAFF/GAFF2 配体；类型名跨版本差异（如 CX vs CT）已全部处理，映射验证零冲突
- **GROMOS**：53A6 / 54A7 家族（本地 gromos54a7.ff 实测）
- **CHARMM**：CHARMM36 / C36m 蛋白 + CGenFF 配体（与 CHARMM-GUI 官方 rtf 逐字一致）
- **OPLS**：OPLS-AA/L 2001（GROMACS oplsaa.ff）

若新版本引入新类型，只需在 `STRONG_AROMATIC`/`COMPATIBLE_TYPES`/`ACCEPTOR_TYPES` 等特征表补一条即可（特征空间映射的设计目标；完整步骤见[扩展指南](extension.md)）。
