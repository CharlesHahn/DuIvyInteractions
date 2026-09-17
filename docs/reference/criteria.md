# 相互作用判据

每种相互作用类型由对应检测器在**几何判定阶段**用距离/角度判据判定每帧是否存在。本文档列出 8 类相互作用的判据与阈值，供理解结果与调参。

> 所有距离单位 Å（Ångström），角度单位度（°）。阈值常量定义在 `interaction_detectors/*_two_pass.py` 顶部。

## 氢键（hydrogen_bond）

**判据**：供体 D 与受体 A 距离 ≤ 4.1 Å，且 D-H···A 角 ≥ 100°。

| 常量 | 值 | 含义 |
|:-----|:--|:-----|
| `HBOND_DIST_MAX` | 4.1 | D-A 最大距离 |
| `HBOND_DON_ANGLE_MIN` | 100 | D-H···A 最小角度 |

## π-π 堆积（pi_stacking）

**判据**：两芳香环中心距离 ∈ (0.5, 5.5] Å，偏移 ≤ 2.0 Å，且平面性满足（可选）。

| 常量 | 值 | 含义 |
|:-----|:--|:-----|
| `PISTACK_MIN_DIST` | 0.5 | 最小距离 |
| `PISTACK_DIST_MAX` | 5.5 | 环心最大距离 |
| `PISTACK_OFFSET_MAX` | 2.0 | 环心偏移上限 |
| `PISTACK_PLANARITY` | 5.0 | 环平面性偏差上限（°），默认关闭 |

**堆积类型分类**（`pistacking_type` metric）：T 型（边对面，两环法向夹角 ≈90°±偏差）或平行（P 型）。输出 XPM 中：白=无、粉=T、蓝=P。

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

**判据**：金属中心（如 Mg²⁺）与配位原子距离 < 3.0 Å。

| 常量 | 值 | 含义 |
|:-----|:--|:-----|
| `METAL_DIST_MAX` | 3.0 | 金属-配位原子最大距离 |

## 水桥（water_bridge）

**判据**：水氧 Ow 到极性原子（供体 D/受体 A）距离 ≤ 4.1 Å，且两个角度在范围内。

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

## 结果解读要点

- **existence**：某 pair 在每帧是否存在（经判据过滤后的布尔矩阵）
- **metrics**：每 pair 每帧的指标值（非活跃帧为 NaN）
- **occupancy（占位率）**：某 pair 存在的帧数 / 总帧数，衡量该相互作用的持久性