# 术语表

本文档汇总工具使用中涉及的术语与符号，供快速查阅。

## 通用概念

| 术语 | 说明 |
|:-----|:-----|
| **tpr** | GROMACS 二进制拓扑文件，含原子类型、键合图、电荷、力场参数 |
| **xtc** | GROMACS 压缩轨迹文件，含每帧原子坐标 |
| **MD** | 分子动力学（Molecular Dynamics）模拟 |
| **力场** | 描述原子间相互作用的参数集（本项目支持 Amber 家族） |
| **GAFF** | General Amber Force Field，通用有机分子力场（配体常用） |
| **ff14SB** | Amber 蛋白力场（本项目测试体系用 amber14sb） |

## 数据结构

| 术语 | 说明 |
|:-----|:-----|
| **SystemData** | 从 tpr 解析出的体系数据（原子、残基、键） |
| **Group** | 一个可参与相互作用的化学基团（如一个芳香环、一个 H 键供体） |
| **pair** | 一对（或多个，如水桥三元组）参与相互作用的基团 |
| **Interaction** | 一种相互作用类型的全部检测结果（矩阵式存储） |
| **existence** | 布尔矩阵 `(n_pairs, n_frames)`，某 pair 在某帧是否存在 |
| **metrics** | 几何指标字典，如 distance、angle，形状 `(n_pairs, n_frames)` |
| **occupancy** | 占位率 = 某 pair 存在的帧数 / 总帧数 |
| **h5** | HDF5 格式结果文件，`dii run` 的输出 |

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

## 相互作用类型

| 类型 | 英文 | 说明 |
|:-----|:-----|:-----|
| 氢键 | hydrogen_bond | D-H···A |
| π-π 堆积 | pi_stacking | 两芳香环堆积（T 型/平行） |
| 盐桥 | salt_bridge | 正负电荷对 |
| 疏水 | hydrophobic | 非极性原子接触 |
| 卤键 | halogen_bond | C-X···A |
| 金属配位 | metal_coordination | 金属-配位原子 |
| 水桥 | water_bridge | D-H···Ow···A |
| π-阳离子 | pi_cation | 芳香环-阳离子 |

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