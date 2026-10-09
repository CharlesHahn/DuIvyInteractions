# 相互作用判据

每种相互作用类型由对应检测器在**几何判定阶段**用距离/角度判据判定每帧是否存在（两段式架构的②，见[核心概念](concepts.md)）。本文档列出 8 类相互作用的判据与阈值，供理解结果与调参。

> 所有距离单位 Å（Ångström），角度单位度（°）。阈值常量定义在检测器文件顶部（`DuIvyInteractions/interaction_detectors/*_detector_{per_frame,two_pass}.py`）。
>
> **策略一致性说明**：`per_frame` 与 `two_pass` 的阈值一致，**唯一例外是水桥**——per_frame 额外施加 2.5 Å 的距离下界，two_pass 只有 4.1 Å 上界（见下文"水桥"节）。`per_tuple` 策略已不被维护（测试约定跳过），其水桥阈值与 per_frame 相同。

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
| `theta`（水桥） | 供体 D、氢 H、水氧 Ow 三点角（顶点 H），即 D-H···Ow | 水桥 |
| `omega`（水桥） | 受体 A、水氧 Ow、氢 H 三点角（顶点 Ow），即 A-Ow···H | 水桥 |

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

**判据（two_pass，默认策略）**：水氧 Ow 到极性原子（供体 D / 受体 A）距离 < 4.1 Å，且 theta ≥ 100°、71° ≤ omega ≤ 140°。

| 常量 | 值 | 含义 |
|:-----|:--|:-----|
| `WATER_BRIDGE_MAXDIST` | 4.1 | Ow 到极性原子最大距离 |
| `WATER_BRIDGE_THETA_MIN` | 100 | D-H···Ow 最小角度 |
| `WATER_BRIDGE_OMEGA_MIN` | 71 | A-Ow···H 最小角度 |
| `WATER_BRIDGE_OMEGA_MAX` | 140 | A-Ow···H 最大角度 |

**策略差异（如实说明）**：

- **per_frame**：`apply_threshold` 额外要求两个距离均有 **2.5 Å 下界**——`2.5 < dist_dw < 4.1` 且 `2.5 < dist_wa < 4.1`（`WATER_BRIDGE_MINDIST = 2.5`，排除重叠/过近原子对）；候选生成阶段还要求第一帧 D···A 距离 ≥ 2.5 Å。
- **two_pass**：**缺 Ow-D 与 Ow-A 距离下界**，只要求 `dist_dw < 4.1` 且 `dist_wa < 4.1`。角度阈值（theta ≥ 100、71 ≤ omega ≤ 140）两策略一致。
- **per_tuple**（已不被维护）：与 per_frame 相同（含 2.5 Å 下界）。

因此同一体系下，per_frame 与 two_pass 的水桥结果在"水分子与极性原子过近（< 2.5 Å）"的情形可能不同；默认策略 two_pass 不含此下界。

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
| 水桥 | dist_dw, dist_wa, theta, omega | Å, Å, °, ° |
| π-阳离子 | distance, offset | Å, Å |

## 判据来源

判据与阈值沿用 PLIP（Protein-Ligand Interaction Profiler，Salentin et al. 2015）的相互作用定义体系：距离/角度截断值（如氢键 D-A 4.1 Å、盐桥 5.5 Å、卤键 165°±30°）均与 PLIP 一致，便于与既有 PLIP 结果对照。π-堆积 T/P 型分类标准与 PLIP 的 `pi-stacking` 分类一致（方法学参考 McGaughey et al. 1998）；水桥角度定义参考 Jiang et al. 2005；卤键几何参考 Auffinger et al.（检测器内注释标注）。具体阈值以检测器文件顶部常量为准。

## 结果解读要点

- **existence**：某 pair 在每帧是否存在（经判据过滤后的布尔矩阵）
- **metrics**：每 pair 每帧的指标值（`Interaction` 内为 (n_pairs, n_frames) 数组；two_pass 的非活跃帧以 NaN 填充）
- **occupancy（占位率）**：某 pair 存在的帧数 / 总帧数（`Interaction.occupancy()`），衡量该相互作用的持久性
