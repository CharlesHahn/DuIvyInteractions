# 使用指南

本部分面向**最终用户**（计算化学与结构生物学研究者）：如何安装、快速上手、使用 `dii` 命令、解读结果。

## 概述

DuIvyInteractions 是一个基于 MD 拓扑力场参数的分子间相互作用判定工具。其核心思路是**直接从 GROMACS tpr 拓扑读取力场原子类型来识别化学基团**，而非像 PLIP/ProLIF 那样从坐标重建化学结构（键序、芳香性、加氢），因此不会丢失 MD 拓扑中已有的化学信息，结果**确定**且与所用**力场自洽**。

当前版本（包版本 v0.0.1，`pyproject.toml`）功能要点：

- **两段式架构**：① 基团鉴定——只做一次、与帧无关，从 tpr 的原子类型 + 键合图 + 显式 H + 电荷确定性地识别基团（供体/受体/芳香环/带电/疏水/卤素/金属/水）；② 几何判定——逐帧按距离/角度/平面临近判据判定相互作用，再做时间统计。
- **4 个力场识别器**（`DuIvyInteractions/group_identifiers/__init__.py` 的 `IDENTIFIER_CLASSES` 注册表）：Amber 家族（amber03/94/96/99/99sb/99sb-ildn/GS/14sb 蛋白 + GAFF/GAFF2 配体）、GROMOS 53A6/54A7、CHARMM36/C36m（含 CGenFF）、OPLS-AA/L。`dii run --ff` 的取值（`amber`/`gromos`/`charmm`/`opls`）与之一一对应。
- **8 类相互作用**：氢键、π-π 堆积、盐桥、疏水、卤键、金属配位、水桥、π-阳离子（`pipeline.py` 的 `ALL_INTERACTIONS`）。
- **检测策略**：`two_pass`（默认；两轮遍历 + KDTree 预筛 + 稀疏存储，性能最优）与 `per_frame`（逐帧向量化）为现行策略；`per_tuple`（逐候选组遍历）处于"可能被舍弃"状态，仅作对照参考。
- **结果存储与导出**：HDF5 序列化（格式版本 1.0，`io/h5.py`）+ xvg/xpm/CSV 导出（`dii export`，与 GROMACS/DuIvyTools 工具链兼容）。
- **真实测试案例（3 力场）**：仓库自带 3 套真实 GROMACS 数据——Amber KRAS-RBD D927（116,383 原子，含 GNP/Mg²⁺）、GROMOS 53A6（130 残基蛋白 + 6 配体）、CHARMM36 SMO-BST——配套 `Tests/unittests/` 集成测试覆盖全部 8 类相互作用；OPLS-AA 的基团识别由单测验证（暂未配套真实轨迹数据）。

可靠性说明（2026-10-09 代码状态）：

- **H 键受体判定（A1 修复）**：跨力场剔除带 H 的非受体 N 类型（普通酰胺/铵/带 H 吡咯/胍基），保留 Pro N / His 无 H 吡啶 / 中性胺 / 核酸氨基；逐项证据见 `doc/acceptor_identification_evidence.md`。
- **跨力场水排除**：水残基名由 `GroupIdentifier.WATER_RESIDUES` 类属性统一管理（Amber `{SOL,HOH,WAT}`、GROMOS `{SOL}`、CHARMM `{TIP3,HOH,SOL,WAT}`、OPLS `{HOH,HO4,HO5,SOL,WAT}`），水桥/金属配位检测据此正确排除水。
- 已知限制（PBC 未处理、GROMOS 配体不支持等）见{doc}`/reference/limitations`。

## 阅读路径

- **安装**：见[安装与环境配置](install)。
- **快速上手**：用真实测试数据跑通完整流程，见[快速上手](quickstart)。
- **命令参考**：`dii run` / `dii export` 的全部参数与行为，见[命令参考](command)。
- **结果解读**：xvg/xpm/csv 每类文件的格式与含义，见[结果解读](result)。

```{toctree}
:maxdepth: 1
:caption: 快速开始

install
quickstart
```

```{toctree}
:maxdepth: 1
:caption: 使用参考

command
result
```
