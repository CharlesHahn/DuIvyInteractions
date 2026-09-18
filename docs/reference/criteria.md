# 相互作用判据

每种相互作用类型由对应检测器在**几何判定阶段**用距离/角度判据判定每帧是否存在。本文档列出 8 类相互作用的判据与阈值，供理解结果与调参。

> 所有距离单位 Å（Ångström），角度单位度（°）。阈值常量定义在检测器文件顶部，三种检测策略（two_pass / per_frame / per_tuple）使用一致的阈值。

## 几何量定义

各判据涉及以下几何量（均为逐帧、逐基团对计算）：

| 几何量 | 定义 | 用于 |
|:-------|:-----|:-----|
| `distance` | 两基团参考点间距离。参考点随类型不同：氢键 = 供体 D 与受体 A；盐桥 = 正/负电中心；π 相关 = 芳香环中心；金属配位 = 金属与配位原子 | 全部类型 |
| `angle`（氢键） | 供体 D、氢 H、受体 A 三点角（顶点 H），即 D-H···A 角 | 氢键 |
| `angle`（π-堆积） | 两芳香环法向量夹角 | π-堆积 |
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
| `HBOND_DIST_MAX` | 4.1 | D-A 最大距离 |
| `HBOND_DON_ANGLE_MIN` | 100 | D-H···A 最小角度 |

## π-π 堆积（pi_stacking）

**判据**：两芳香环中心距离 ∈ (0.5, 5.5] Å，且堆积类型非 N（`pistacking_type != 'N'`，即分类为 T 型或 P 型）；平面性满足为可选条件（默认关闭）。

| 常量 | 值 | 含义 |
|:-----|:--|:-----|
| `PISTACK_MIN_DIST` | 0.5 | 最小距离 |
| `PISTACK_DIST_MAX` | 5.5 | 环心最大距离 |
| `PISTACK_OFFSET_MAX` | 2.0 | 环心偏移上限 |
| `PISTACK_PLANARITY` | 5.0 | 环平面性偏差上限（°），默认关闭 |
| `PISTACK_ANG_DEV` | 30 | T/P 分类的角度偏差 |

**堆积类型分类**（`pistacking_type` metric），基于两环法向量夹角 `angle` 与偏移 `offset`：
- **P 型（平行堆积）**：`angle ≤ 30°` 且 `offset < 2.0 Å`
- **T 型（边对面堆积）**：`angle ≥ 60°` 且 `offset < 2.0 Å`
- **N（无）**：其余情况

输出 XPM 中：白=无、粉=T 型、蓝=P 型。

## 盐桥（salt_bridge）

**判据**：正电中心（LYS/ARG 等）与负电中心（ASP/GLU 等）距离 ≤ 5.5 Å。

| 常量 | 值 | 含义 |
|:-----|:--|:-----|
| `SALTBRIDGE_DIST_MAX` | 5.5 | 电荷中心最大距离 |

## 疏水（hydrophobic）

**判据**：两疏水原子距离 ∈ (0.5, 4.0) Å。

| 常量 | 值 | 含义 |
|:-----|:--|:-----|
| `HYDROPH_MIN_DIST` | 0.5 | 最小距离 |
| `HYDROPH_DIST_MAX` | 4.0 | 最大距离 |

## 卤键（halogen_bond）

**判据**：卤素 X 与受体 A 距离 ≤ 4.0 Å，C-X···A 角 = 165°±30°，X···A-R 角 = 120°±30°。

| 常量 | 值 | 含义 |
|:-----|:--|:-----|
| `HALOGEN_DIST_MAX` | 4.0 | X···A 最大距离 |
| `HALOGEN_DON_ANGLE` | 165 | C-X···A 最优角度 |
| `HALOGEN_ACC_ANGLE` | 120 | X···A-R 最优角度 |
| `HALOGEN_ANGLE_DEV` | 30 | 角度偏差上限 |

## 金属配位（metal_coordination）

**判据**：金属中心（如 Mg²⁺）与配位原子距离 < 3.0 Å。配位原子来自 `metal_binding` 基团（金属周围可配位的原子）。

| 常量 | 值 | 含义 |
|:-----|:--|:-----|
| `METAL_DIST_MAX` | 3.0 | 金属-配位原子最大距离 |

> 注：纯水分子配位不报告（所有配位原子均来自水时忽略），第一版只做距离判据，未做配位几何构型匹配。

## 水桥（water_bridge）

**判据**：水氧 Ow 到极性原子（供体 D/受体 A）距离 < 4.1 Å，且两个角度在范围内（theta ≥ 100°、71° ≤ omega ≤ 140°）。

| 常量 | 值 | 含义 |
|:-----|:--|:-----|
| `WATER_BRIDGE_MAXDIST` | 4.1 | Ow 到极性原子最大距离 |
| `WATER_BRIDGE_THETA_MIN` | 100 | 水O-供体H-供体D 最小角度 |
| `WATER_BRIDGE_OMEGA_MIN` | 71 | 受体-水O-供体H 最小角度 |
| `WATER_BRIDGE_OMEGA_MAX` | 140 | 受体-水O-供体H 最大角度 |

## π-阳离子（pi_cation）

**判据**：环心到阳离子（正电基团）距离 ∈ (0.5, 6.0) Å，偏移 < 2.0 Å。

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
| π-堆积 | distance, angle, offset, pistacking_type | Å, °, Å, 类型 |
| 盐桥 | distance | Å |
| 疏水 | distance | Å |
| 卤键 | distance, don_angle, acc_angle | Å, °, ° |
| 金属配位 | distance | Å |
| 水桥 | dist_dw, dist_wa, theta, omega | Å, Å, °, ° |
| π-阳离子 | distance, offset | Å, Å |

## 判据来源

判据与阈值沿用 PLIP（Protein-Ligand Interaction Profiler）的相互作用定义体系：距离/角度截断值（如氢键 D-A 4.1 Å、盐桥 5.5 Å、卤键 165°±30°）均与 PLIP 一致，便于与既有 PLIP 结果对照。π-堆积 T/P 型分类标准与 PLIP 的 `pi-stacking` 分类一致。

## 结果解读要点

- **existence**：某 pair 在每帧是否存在（经判据过滤后的布尔矩阵）
- **metrics**：每 pair 每帧的指标值（非活跃帧为 NaN）
- **occupancy（占位率）**：某 pair 存在的帧数 / 总帧数，衡量该相互作用的持久性