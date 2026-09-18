# 核心概念

本文档说明工具的工作原理：如何从 GROMACS tpr 拓扑读取力场原子类型来鉴定化学基团，并据此检测分子间相互作用。

## 方法概述

分子间相互作用分析需要先确定分子的化学结构——哪些原子构成芳香环、哪些原子是氢键供体/受体、哪些基团带电荷等。工具基于这些化学基团，再结合轨迹坐标做几何判定。

确定化学结构有两种路径：

1. **从坐标重建**：从原子坐标出发，推断键序、芳香性和氢原子位置。常用工具（如 PLIP 使用 OpenBabel、ProLIF 使用 RDKit）采用此路径。
2. **从力场参数读取**（本工具）：直接从 GROMACS 拓扑（tpr）读取力场原子类型。这些类型在力场参数化时已由 antechamber/sobtop 依据分子电子结构确定。

本工具采用第二种路径。

## 力场原子类型的化学语义

力场原子类型不是任意标签，其命名编码了化学特征。以 GAFF（General Amber Force Field）为例，类型名后缀有系统含义：

| 类型 | 化学含义 | 命名依据 |
|:-----|:---------|:---------|
| `c3` | sp3 碳 | c + 数字 3 |
| `c2` | sp2 烯碳 | c + 数字 2 |
| `ca` | 芳香碳 | c + a（aromatic） |
| `c` | sp2 羰基碳 | c 无后缀 |
| `na` | 吡咯型芳香氮 | n + a |
| `nb` | 吡啶型芳香氮 | n + b |
| `os` | 醚氧（单键） | o + s（single） |

Amber 蛋白类型类似：残基定义（rtp 文件）中环原子使用的类型名（如 `CA`）可对照已知残基化学反推其含义。

类型名编码了杂化状态、芳香性、极性、是否带氢等特征。因此读取类型即可确定化学基团，不需要从坐标推断。这保证了：基团鉴定只依赖拓扑、不依赖坐标，跨帧一致、可复现。

## 两段式架构

分析分为两个阶段：

```
① 基团鉴定
   tpr 力场参数（原子类型 + 键合图 + 显式 H + 电荷）
   → 识别化学基团（芳香环、H 键供受体、带电基团等）
   —— 只做一次，与帧无关

② 几何判定
   基团带标签后，逐帧用距离/角度/平面临近判据
   → 每帧相互作用列表 → 时间统计
```

- **基团鉴定**只依赖拓扑，不依赖坐标，因此结果与帧无关、可复现
- **几何判定**对每个基团对逐帧计算几何指标（距离、角度等），按判据判定相互作用是否存在

## 显式氢与化学信息

全原子力场（Amber 系）中，以下化学信息在拓扑中显式存在：

- **氢键供体**：供体原子 D 与氢 H 之间的键条目存在，且 H 带正电荷
- **水分子**：由残基名（SOL/HOH）标识，氧原子 OW、氢原子 HW 明确
- **金属中心**：由元素（如 Mg）与电荷标识

这些信息直接从拓扑读取，不涉及坐标推断。

## 支持范围

- **力场**：Amber 家族（amber03/94/96/99/99sb/99sb-ildn/GS/14sb）蛋白 + GAFF/GAFF2 配体，类型映射已验证零冲突
- **依赖显式氢**：工具依赖全原子力场的显式氢原子；对联合原子力场（如 GROMOS，无显式 H）供体鉴定不适用
- **键序缺失**：即使 tpr 不保留键序（全 func=1），芳香性判定仍基于类型名，键序仅作交叉验证

## 详细规则

- 类型→化学特征的完整映射：见[力场类型映射](force_field.md)
- 各基团的识别规则：见[基团识别规则](group_rules.md)
- 各相互作用的几何判据：见[相互作用判据](criteria.md)

## 参考文献

- GAFF 力场：Wang J, Wolf RM, Caldwell JW, Kollman PA, Case DA. Development and testing of a general amber force field. *J Comput Chem*. 2004;25(9):1157-1174.
- Amber 力场（ff14SB）：Maier JA, Martinez C, Kasavajhala K, Wickstrom L, Hauser KE, Simmerling C. ff14SB: Improving the accuracy of protein side chain and backbone parameters from ff99SB. *J Chem Theory Comput*. 2015;11(8):3696-3713.
- PLIP 相互作用定义：Salentin S, Schreiber S, Haupt VJ, Adasme MF, Schroeder M. PLIP: fully automated protein-ligand interaction profiler. *Nucleic Acids Res*. 2015;43(W1):W443-W447.
- MDAnalysis：Michaud-Agrawal N, Denning EJ, Woolf TB, Beckstein O. MDAnalysis: A toolkit for the analysis of molecular dynamics simulations. *J Comput Chem*. 2011;32(10):2319-2327.