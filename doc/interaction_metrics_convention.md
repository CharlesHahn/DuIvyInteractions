# 相互作用指标约定

> 状态：进行中
> 用途：记录每种相互作用在 Interaction 数据结构中保存的 metrics 键名和含义

---

## 通用约定

- 单位：Å（MDAnalysis 对 GROMACS 轨迹返回 Å）
- 所有 metrics 数组 shape = (n_pairs, n_frames)
- existence 数组 shape = (n_pairs, n_frames), dtype=bool

---

## 各相互作用的 metrics

### 1. salt_bridge（盐桥）

| 键名 | 类型 | 说明 |
|:-----|:-----|:-----|
| distance | float | 电荷中心之间的距离（Å） |

**阈值**：distance < 5.5 Å

---

### 2. hydrogen_bond（氢键）

| 键名 | 类型 | 说明 |
|:-----|:-----|:-----|
| distance | float | D-A 距离（Å） |
| angle | float | D-H···A 角度（°） |

**阈值**：distance < 4.1 Å 且 angle > 100°

---

### 3. pi_stacking（π-π 堆积）

| 键名 | 类型 | 说明 |
|:-----|:-----|:-----|
| distance | float | 两环环心距离（Å） |
| angle | float | 两环法向量夹角（°），已取 min(θ, 180-θ) |
| offset | float | 环心投影偏移量（Å），双向取较小值 |
| pistacking_type | str(U1) | 堆积类型：'N'（无）、'P'（平行）、'T'（T 型） |
| planarity_ring1 | float | 环 1 法向量两两最大夹角（°），仅 check_planarity=True |
| planarity_ring2 | float | 环 2 法向量两两最大夹角（°），仅 check_planarity=True |

**阈值**：
- 0.5 Å < distance < 5.5 Å
- P 型：angle ≤ 30° 且 offset < 2.0 Å
- T 型：angle ≥ 60° 且 offset < 2.0 Å
- 平面性（可选）：两个环的 planarity ≤ 5.0°

**参数**：
- `check_planarity: bool = False`：是否逐帧检验环平面性

---

### 4. pi_cation（π-阳离子 / 阳离子-π）

| 键名 | 类型 | 说明 |
|:-----|:-----|:-----|
| distance | float | 正电荷中心到环心距离（Å） |
| offset | float | 电荷中心投影偏移量（Å） |

**阈值**：distance < 6.0 Å，offset < 2.0 Å

---

### 5. hydrophobic（疏水）

| 键名 | 类型 | 说明 |
|:-----|:-----|:-----|
| distance | float | 两疏水原子之间距离（Å） |

**阈值**：distance < 4.0 Å

---

### 6. halogen_bond（卤键）

| 键名 | 类型 | 说明 |
|:-----|:-----|:-----|
| distance | float | 卤素 X 到受体 A 的距离（Å） |
| don_angle | float | C-X···A 角度（°） |
| acc_angle | float | X···A-R 角度（°） |

**阈值**：distance < 4.0 Å，don_angle 在 165°±30°，acc_angle 在 120°±30°

---

### 7. metal_coordination（金属配位）

| 键名 | 类型 | 说明 |
|:-----|:-----|:-----|
| distance | float | 金属离子到配位原子距离（Å） |

**阈值**：distance < 3.0 Å

---

### 8. water_bridge（水桥）

| 键名 | 类型 | 说明 |
|:-----|:-----|:-----|
| dist_donor_water | float | 供体到水 O 距离（Å） |
| dist_water_acceptor | float | 水 O 到受体距离（Å） |
| omega | float | 受体-水O-供体H 角度（°） |
| theta | float | 水O-供体H-供体D 角度（°） |

**阈值**：距离在 2.5~4.1 Å，omega 在 71°~140°，theta > 100°

---

*文档结束*
