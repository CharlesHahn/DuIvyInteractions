# 相互作用判据

每种相互作用类型由对应检测器在**几何判定阶段**用距离/角度判据判定每帧是否存在（两段式架构的②，见[核心概念](concepts.md)）。本文档列出 8 类相互作用的判据与阈值，供理解结果与调参。

> 所有距离单位 Å（Ångström），角度单位度（°）。阈值常量定义在检测器文件顶部（`DuIvyInteractions/interaction_detectors/*_detector_{per_frame,two_pass}.py`）。
>
> **策略一致性说明**：`per_frame` 与 `two_pass` 的阈值一致（含水桥——两策略采用同一套显式-H 判据，见下文"水桥"节）。`per_tuple` 策略已不被维护（测试约定跳过），其水桥仍用旧 PLIP 判据。
>
> **周期边界（PBC）说明**：下列几何量在检测时按**最小镜像（MIC）**约定计算（`two_pass`/`per_frame` 策略，PBC 已启用）：距离用 `calc_bonds(box)`、角度用 `calc_angles(box)`（三原子角顶点见各类型定义），聚合坐标（电荷中心、环心）与投影偏移先对原子做 MIC 回绕再计算。判据阈值本身不受 PBC 影响——PBC 只改变距离/角度的**数值来源**（跨盒对取最近镜像），不改变判据。无盒体系（`ts.dimensions` 缺失或全零）自动退化为欧氏计算。详见 `doc/pbc_design.md`。

## 几何量定义

各判据涉及以下几何量（均为逐帧、逐基团对计算）：

| 几何量 | 定义 | 用于 |
|:-------|:-----|:-----|
| `distance` | 两基团参考点间距离。参考点随类型不同：氢键 = 供体原子 D 与受体原子 A；盐桥 = 正/负电基团的**电荷加权中心**（原子部分电荷加权的质心）；π-堆积 = 芳香环几何中心；π-阳离子 = 环心与阳离子电荷中心；卤键 = 卤素 X 与受体 A；金属配位 = 金属与配位原子 | 全部类型 |
| `angle`（氢键） | 供体 D、氢 H、受体 A 三点角（顶点 H），即 D-H···A 角 | 氢键 |
| `angle`（π-堆积） | 两芳香环法向量夹角（取最小角 `min(raw, 180°−raw)`，处理方向不确定性；环法向量由环内相邻原子叉积的平均得到） | π-堆积 |
| `offset`（π-堆积） | 一环中心相对另一环平面的横向偏移：将环中心沿另一环法向投影到其平面，取投影点与另一环中心的距离；两环互为参考，取较小值 | π-堆积、π-阳离子 |
| `planarity` | 环内各原子局部法向量（相邻两键叉积）两两夹角的最大值 | π-堆积（可选） |
| `don_angle` | 碳 C、卤素 X、受体 A 三点角（顶点 X），即 C-X···A 角 | 卤键 |
| `acc_angle` | 受体 A 邻接原子中，X···A-R 角最接近 120° 者 | 卤键 |
| `angle_dh_ow`（水桥） | 供体 D、氢 H、水氧 Ow 三点角（顶点 H），即 D-H···Ow | 水桥 |
| `angle_ow_hw_a`（水桥） | 水氧 Ow、水的氢 Hw、受体 A 三点角（顶点 Hw），即 Ow-Hw···A | 水桥 |
| `angle_h_ow_a`（水桥） | 供体 H、水氧 Ow、受体 A 三点角（顶点 Ow）——"水在中间"约束 | 水桥 |

## 氢键（hydrogen_bond）

**判据**：供体 D 与受体 A 距离 ≤ 4.1 Å，且 D-H···A 角 ≥ 100°。

| 常量 | 值 | 含义 |
|:-----|:--|:-----|
| `HBOND_DIST_MAX` | 4.1 | D-A 最大距离（含） |
| `HBOND_DON_ANGLE_MIN` | 100 | D-H···A 最小角度 |

## π-π 堆积（pi_stacking）

**判据**：两芳香环中心距离 ∈ (0.5, 5.5] Å，且堆积类型非 N（`pistacking_type != 'N'`，即分类为 T 型或 P 型）；平面性满足为可选条件（默认关闭，`check_planarity=False`）。

| 常量 | 值 | 含义 |
|:-----|:--|:-----|
| `PISTACK_MIN_DIST` | 0.5 | 最小距离（不含，排除自身） |
| `PISTACK_DIST_MAX` | 5.5 | 环心最大距离（含） |
| `PISTACK_OFFSET_MAX` | 2.0 | 环心偏移上限 |
| `PISTACK_PLANARITY` | 5.0 | 环平面性偏差上限（°），默认关闭 |
| `PISTACK_ANG_DEV` | 30 | T/P 分类的角度偏差 |

**堆积类型分类**（`pistacking_type` metric，两策略实现一致），基于两环法向量夹角 `angle` 与偏移 `offset`：

- **P 型（平行堆积）**：`angle ≤ 30°` 且 `offset < 2.0 Å`
- **T 型（边对面堆积）**：`angle ≥ 60°` 且 `offset < 2.0 Å`
- **N（无）**：其余情况

输出 XPM 中：0=无（白色 `#FFFFFF`）、1=T 型（粉 `#F67088`）、2=P 型（蓝 `#38A7D0`）。

## 盐桥（salt_bridge）

**判据**：正电基团与负电基团的**电荷加权中心**距离 ≤ 5.5 Å（电荷中心 = 基团原子按部分电荷加权的质心）。

| 常量 | 值 | 含义 |
|:-----|:--|:-----|
| `SALTBRIDGE_DIST_MAX` | 5.5 | 电荷中心最大距离（含） |

## 疏水（hydrophobic）

**判据**：两疏水原子距离 ∈ (0.5, 4.0) Å（两端均不含）。

| 常量 | 值 | 含义 |
|:-----|:--|:-----|
| `HYDROPH_MIN_DIST` | 0.5 | 最小距离 |
| `HYDROPH_DIST_MAX` | 4.0 | 最大距离 |

## 卤键（halogen_bond）

**判据**：卤素 X 与受体 A 距离 ≤ 4.0 Å，C-X···A 角 ∈ [135°, 195°]（165°±30°），X···A-R 角 ∈ [90°, 150°]（120°±30°，取 A 的 R 邻居中最接近 120° 者）。

| 常量 | 值 | 含义 |
|:-----|:--|:-----|
| `HALOGEN_DIST_MAX` | 4.0 | X···A 最大距离（含） |
| `HALOGEN_DON_ANGLE` | 165 | C-X···A 最优角度 |
| `HALOGEN_ACC_ANGLE` | 120 | X···A-R 最优角度 |
| `HALOGEN_ANGLE_DEV` | 30 | 角度偏差上限 |

## 金属配位（metal_coordination）

**判据**：金属中心（如 Mg²⁺）与配位原子距离 < 3.0 Å。配位原子来自 `metal_binding` 基团（元素 ∈ {O, N, S} 的原子）。

| 常量 | 值 | 含义 |
|:-----|:--|:-----|
| `METAL_DIST_MAX` | 3.0 | 金属-配位原子最大距离（不含） |

> 注：纯水分子配位不报告——`metal_binding` 基团在识别阶段就把水残基整体排除（按力场 `WATER_RESIDUES`）；第一版只做距离判据，未做配位几何构型匹配。

## 水桥（water_bridge）

**判据（per_frame / two_pass）**：基于**水的显式氢**，水桥的两侧各判一个氢键（Baker & Hubbard 1984）：

- 供体侧 `D–H···Ow`：`H···Ow ≤ 2.5 Å` 且 `D-H···Ow ≥ 120°`
- 受体侧 `Ow–Hw···A`：`Hw···A ≤ 2.5 Å` 且 `Ow-Hw···A ≥ 120°`
- 水在中间：`∠(供体H–Ow–受体A)` ∈ **[71°, 140°]**（保留 PLIP `omega` 语义，保证水确实桥接两端）

| 常量 | 值 | 含义 |
|:-----|:--|:-----|
| `WATER_BRIDGE_HB_DIST_MAX` | 2.5 | H···A 最大距离（两侧） |
| `WATER_BRIDGE_HB_ANGLE_MIN` | 120 | D-H···A 最小角度（两侧） |
| `WATER_BRIDGE_MIN_DIST` | 0.5 | 防原子重叠下界 |
| `WATER_BRIDGE_MAX_DONORS` | 2 | 每水每帧最多供体数（2 对孤对） |
| `WATER_BRIDGE_MAX_BRIDGES` | 2 | 每水每帧最多桥数（2 个 H） |
| `WATER_BRIDGE_OMEGA_MIN` / `MAX` | 71 / 140 | 供体H–Ow–受体A 角范围（水在中间） |

**三层配额（逐帧）**：① 每个水 H 最多对应 1 个受体（取角度最优）；② 每个水 O 最多对应 2 个供体（取 H···O 最近）；③ 每水最多 2 个桥。配额消除"同一水被组合爆炸重复计入"的物理不合理结果。

**策略说明**：

- **per_frame / two_pass**：采用上述显式-H 判据与三层配额，二者帧级结果一致。
- **per_tuple**（已不被维护）：保留 PLIP 旧判据（Ow 距离 2.5–4.1 Å + theta ≥ 100° + omega 71–140°），因该策略逐候选组遍历、架构上无法实现跨组配额。

**文献依据与设计推导**见 `doc/water_bridge_design.md`。

## π-阳离子（pi_cation）

**判据**：环心到阳离子（正电基团）**电荷中心**距离 ∈ (0.5, 6.0) Å（两端均不含），偏移 < 2.0 Å。

| 常量 | 值 | 含义 |
|:-----|:--|:-----|
| `PICATION_MIN_DIST` | 0.5 | 最小距离 |
| `PICATION_DIST_MAX` | 6.0 | 环-阳离子最大距离 |
| `PICATION_OFFSET_MAX` | 2.0 | 偏移上限 |

## 各类型的指标（metrics）

检测器输出的指标（存于 h5 的 metrics，导出为 xvg）：

| 类型 | 指标 | 单位 |
|:-----|:-----|:-----|
| 氢键 | distance, angle | Å, ° |
| π-堆积 | distance, angle, offset, pistacking_type（启用平面性时另有 planarity_ring1/2） | Å, °, Å, 类型 |
| 盐桥 | distance | Å |
| 疏水 | distance | Å |
| 卤键 | distance, don_angle, acc_angle | Å, °, ° |
| 金属配位 | distance | Å |
| 水桥 | dist_h_ow, angle_dh_ow, dist_d_ow, dist_hw_a, angle_ow_hw_a, dist_ow_a, angle_h_ow_a | Å, °, Å, Å, °, Å, ° |
| π-阳离子 | distance, offset | Å, Å |

## 判据来源

判据与阈值沿用 PLIP（Protein-Ligand Interaction Profiler，Salentin et al. 2015）的相互作用定义体系：距离/角度截断值（如氢键 D-A 4.1 Å、盐桥 5.5 Å、卤键 165°±30°）均与 PLIP 一致，便于与既有 PLIP 结果对照。π-堆积 T/P 型分类标准与 PLIP 的 `pi-stacking` 分类一致（方法学参考 McGaughey et al. 1998）；水桥：两侧氢键几何采用 Baker & Hubbard 1984 标准（H···A ≤ 2.5 Å、D-H···A ≥ 120°），"水在中间"的 omega 约束（71°–140°）沿用 PLIP/Jiang et al. 2005；卤键几何参考 Auffinger et al.（检测器内注释标注）。具体阈值以检测器文件顶部常量为准。

## 结果解读要点

- **existence**：某 pair 在每帧是否存在（经判据过滤后的布尔矩阵）
- **metrics**：每 pair 每帧的指标值（`Interaction` 内为 (n_pairs, n_frames) 数组；two_pass 的非活跃帧以 NaN 填充）
- **occupancy（占位率）**：某 pair 存在的帧数 / 总帧数（`Interaction.occupancy()`），衡量该相互作用的持久性
