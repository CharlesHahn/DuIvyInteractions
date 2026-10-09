# 术语表

本文档汇总工具使用中涉及的术语与符号，供快速查阅。基团类型与键类型的枚举与 `DuIvyInteractions/core/constants.py` 的 `GROUP_TYPES` / `BOND_TYPES` 保持一致。

## 力场与通用概念

| 术语 | 说明 |
|:-----|:-----|
| **tpr** | GROMACS 二进制拓扑文件，含原子类型、键合图、电荷、力场参数 |
| **xtc** | GROMACS 压缩轨迹文件，含每帧原子坐标 |
| **MD** | 分子动力学（Molecular Dynamics）模拟 |
| **力场** | 描述原子间相互作用的参数集。本项目支持 **4 个力场家族**：Amber 家族、GROMOS 53A6/54A7、CHARMM36/C36m、OPLS-AA/L，由 `IDENTIFIER_CLASSES` 注册表统一管理 |
| **Amber 家族** | amber03/94/96/99/99SB/99SB-ildn/GS/14SB 蛋白力场 + GAFF/GAFF2 配体力场；全原子力场，显式 H 完整 |
| **GAFF / GAFF2** | General Amber Force Field，Amber 家族的通用有机分子力场（配体常用） |
| **GROMOS 53A6/54A7** | GROMOS 联合原子（united atom）力场：脂肪族 C 的 H 并入重原子，但**极性 H（N/O/S 上的 H）显式**。配 SPC/SPC-E 水模型（残基名 `SOL`）。配体不支持（需 ATB 参数化）。类型粒度粗（`N`/`C` 通吃多种化学环境），靠结构判据补偿 |
| **CHARMM36 / C36m** | CHARMM 蛋白力场；配体用 CGenFF（CHARMM General Force Field）。默认 TIP3 水模型（残基名 `TIP3`/`HOH`） |
| **OPLS-AA/L** | OPLS 全原子力场（2001 版，GROMACS `oplsaa.ff`）。水模型 HOH/SPC、HO4/TIP4P、HO5/TIP5P（残基名 `HOH`/`HO4`/`HO5`） |
| **联合原子（united atom）** | 将脂肪族 C 上的 H 并入重原子的力场风格（如 GROMOS）。此类力场无脂肪族显式 H，依赖显式 H 的判定（供体、疏水邻接）受限 |
| **显式 H** | 力场中显式建模的氢原子。H 键供体判定依赖 D–H 键与 q(H)>0，故需显式 H |
| **水模型** | 描述水分子的力场参数：SPC/SPC-E（GROMOS）、TIP3（CHARMM/Amber 常用）、TIP4P/TIP5P（OPLS 的 HO4/HO5）。残基名：`SOL`/`HOH`/`WAT`/`TIP3`/`HO4`/`HO5`，各力场的 `WATER_RESIDUES` 类属性收录 |

## 数据结构

| 术语 | 说明 |
|:-----|:-----|
| **SystemData** | 从 tpr 解析出的体系数据（原子、残基、键、残基间键） |
| **Group** | 一个可参与相互作用的化学基团（如一个芳香环、一个 H 键供体） |
| **pair** | 一对（或多个，如水桥三元组）参与相互作用的基团 |
| **Interaction** | 一种相互作用类型的全部检测结果（矩阵式存储） |
| **InteractionSparse** | TwoPass 策略 Pass1 的稀疏中间结果（键为 `(group_id, ...)` 元组） |
| **existence** | 布尔矩阵 `(n_pairs, n_frames)`，某 pair 在某帧是否存在 |
| **metrics** | 几何指标字典，如 distance、angle，形状 `(n_pairs, n_frames)` |
| **occupancy** | 占位率 = 某 pair 存在的帧数 / 总帧数 |
| **h5** | HDF5 格式结果文件，`dii run` 的输出（格式版本 1.0） |

## 基团类型（GROUP_TYPES）

| 基团类型 | 含义 | 参与相互作用 |
|:---------|:-----|:------------|
| `H_donor` | 氢键供体（D-H） | 氢键、水桥 |
| `H_acceptor` | 氢键受体（有孤对） | 氢键、水桥 |
| `aromatic_ring` | 芳香环 | π-堆积、π-阳离子、卤键 |
| `charged_positive` | 正电基团 | 盐桥、π-阳离子 |
| `charged_negative` | 负电基团 | 盐桥 |
| `halogen_donor` | 卤键供体（C-X） | 卤键 |
| `halogen_acceptor` | 卤键受体 | 卤键 |
| `metal` | 金属中心 | 金属配位 |
| `metal_binding` | 金属配位原子 | 金属配位 |
| `water` | 水分子 | 水桥 |
| `hydrophobic` | 疏水原子 | 疏水 |

## 键类型（BOND_TYPES）

| 键类型 | 含义 |
|:-------|:-----|
| `bond` | 未知键级（默认值） |
| `single` / `double` / `triple` / `aromatic` | 化学键 |
| `constrained` / `settle` / `virtual` | 约束键（如 SHAKE/LINCS）、SETTLE 水分子键、虚拟原子键 |

> 注：N–H 供体键在 tpr 的 `Constraint:` 段而非 `Bond:` 段，供体识别需合并两段；`BOND_TYPES` 中的 `constrained` 对应此类。

## 相互作用类型

| 类型 | 英文 | 说明 |
|:-----|:-----|:-----|
| 氢键 | hydrogen_bond | D-H···A |
| π-π 堆积 | pi_stacking | 两芳香环堆积（T 型/平行），指标含 `pistacking_type`（P/T/N） |
| 盐桥 | salt_bridge | 正负电荷对 |
| 疏水 | hydrophobic | 非极性原子接触 |
| 卤键 | halogen_bond | C-X···A |
| 金属配位 | metal_coordination | 金属-配位原子 |
| 水桥 | water_bridge | D-H···Ow···A |
| π-阳离子 | pi_cation | 芳香环-阳离子 |

## 检测策略

| 策略 | 说明 |
|:-----|:-----|
| `two_pass` | 两轮遍历：Pass1 逐帧发现活跃对（KDTree 预筛，稀疏存储）→ Pass2 补全全帧指标。默认策略，性能最优 |
| `per_frame` | 逐帧向量化处理全部候选对。候选对极多（如水桥）时适用；长轨迹时预分配矩阵内存占用高 |
| `per_tuple` | 逐候选组加载全部帧并向量化。对照实现，当前处于"可能被舍弃"状态（测试不跑） |

## 符号

| 符号 | 含义 |
|:-----|:-----|
| Å | 埃，长度单位（1 Å = 0.1 nm） |
| ° | 度，角度单位 |
| D | 氢键供体原子（N/O/S/F） |
| A | 氢键受体原子 |
| H | 氢原子 |
| Ow | 水分子氧原子 |
| X | 卤素原子（F/Cl/Br/I） |
| q(H) | 氢原子部分电荷 |
| P 型 / T 型 | π 堆积平行型 / 边对面（T 形）型 |