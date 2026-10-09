# DuIvyInteractions 文档

基于 MD 拓扑力场参数的分子间相互作用判定工具。

**核心思路**：直接从 GROMACS tpr 拓扑读取力场原子类型来识别化学基团——基团鉴定基于力场参数本身，而非从坐标重建化学结构。这一做法直接复用 MD 拓扑中的力场类型语义（如 GAFF `ca`、CHARMM `CG2R61`、OPLS `CA` 等），结果**确定**（同一输入必得同一基团集合）且与所用**力场自洽**（化学信息来自模拟力场本身，而非第三方推断）；与依赖 OpenBabel（PLIP）/ RDKit（ProLIF）从坐标重建键序、芳香性与加氢的路线相比，不丢弃 MD 拓扑化学信息，也不引入重建不确定性。

## 特性概览

**多力场支持（4 力场）**：基团识别覆盖 Amber 家族（amber03/94/96/99/99sb/99sb-ildn/GS/14sb + GAFF）、GROMOS 53A6/54A7、CHARMM36/C36m（含 CGenFF）、OPLS-AA/L。识别器由 `IDENTIFIER_CLASSES` 注册表统一管理（`DuIvyInteractions/group_identifiers/__init__.py`），`dii run --ff` 可选 4 力场；水残基名（`WATER_RESIDUES`）按力场区分（Amber SOL/HOH/WAT、GROMOS SOL、CHARMM TIP3、OPLS HO4/HO5），保证水桥/金属配位检测正确排除水。

**8 类相互作用**：氢键、π-π 堆积、盐桥、疏水相互作用、卤键、金属配位、水桥、π-阳离子。每类提供 per_frame 与 two_pass 两种现行检测策略（per_tuple 策略处于待定状态），结果统一为 `List[Interaction]` 并存储为 HDF5。

**两段式架构**：

1. **基团鉴定**（与帧无关、只做一次）：tpr 原子类型 + 键合图 + 显式 H + 电荷 → 确定性识别基团（供体/受体/芳香环/带电/疏水/卤素/金属/水）。H 键受体判定经 A1 修复（2026-09-30）：剔除带 H 的非受体 N 类型，证据清单见 `doc/acceptor_identification_evidence.md`。
2. **几何判定**（逐帧）：按距离/角度/平面临近判据 → 每帧相互作用列表 → 时间统计。

**真实测试案例（3 力场）**：Amber KRAS–RBD D927（含 GNP/Mg 金属中心）、GROMOS 53A6（130 残基蛋白 + 6 配体 ZIN1–6）、CHARMM36 SMO–BST（Smoothened–β-谷甾醇，Mendeley v94vzbwzf3，含删水复合物与全原子 63055 SOL 水两套 tpr，64 项测试）三套真实 MD 数据均有集成测试覆盖（`Tests/test_MD_case_*` 与 `Tests/unittests/`），8 类相互作用的检测均在真实数据上验证（SMO-BST 体系无卤素/金属中心，此两类以零结果断言覆盖）。

## 快速开始

```bash
dii run -t md.tpr -f md.xtc -o out/ --ff amber      # 检测 → h5（--ff 可选 amber/gromos/charmm/opls）
dii export -i out/hydrogen_bond.h5 -o out_export/  # 导出 xvg/xpm/csv + 概览
```

```{toctree}
:maxdepth: 1
:caption: 使用指南

guide/index
```

```{toctree}
:maxdepth: 1
:caption: 参考手册

reference/index
```

```{toctree}
:maxdepth: 1
:caption: 其他

changelog
```