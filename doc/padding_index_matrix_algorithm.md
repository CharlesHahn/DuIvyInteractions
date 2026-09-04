# Padding 索引矩阵算法设计文档

> 创建日期：2026-08-25
> 状态：定稿（对抗性审查问题已全部修复，见 §6）
> 用途：InteractionDetectorPerFrame 的向量化计算方案

---

## 1. 问题背景

### 1.1 核心矛盾

逐帧检测相互作用时，每帧需要对 N 个候选 tuple 计算指标。每个 tuple 由若干 group 组成，每个 group 包含不定数量的原子。

numpy 要求矩形数组（统一 shape），但 group 内原子数不统一：

| 基团类型 | 典型原子数 | 示例 |
|:---------|:----------|:-----|
| ARG 胍基（带电） | 9 | CZ, NE, NH1, NH2, HE, HH11, HH12, HH21, HH22 |
| LYS 氨基（带电） | 4 | NZ, HZ1, HZ2, HZ3 |
| ASP 羧酸（带电） | 3 | CG, OD1, OD2 |
| H 键供体 | 2 | D, H |
| H 键受体 | 1 | A |
| 芳香环 | 5~7 | 环原子 |
| 水分子 | 3 | Ow, Hw1, Hw2 |

直接放进 numpy 数组会因 shape 不一致而失败。

### 1.2 朴素方案的问题

**方案 A：逐 tuple Python 循环**

```
每帧:
  for each tuple (N 个):
      取坐标 → 算指标
```

问题：N 可能很大（水桥 24.9 万 tuple），逐个循环慢。

**方案 B：逐 group Python 循环 + numpy 组装**

```
每帧:
  for each unique group (M 个):    ← M << N
      算 group 级量（电荷中心等）
  numpy 索引组装 N 个 tuple 结果   ← 零循环
```

问题：group 级循环虽然 M << N，但仍有 Python 循环，且每个 group 内的计算也是循环。

---

## 2. Padding 索引矩阵方案

### 2.1 核心思想

**把不定长度的 group 原子列表，补齐（padding）到统一最大长度，变成矩形数组，从而用纯 numpy 广播消除所有 Python 循环。**

关键：padding 位置填入的值是"垃圾"，但通过 valid 掩码和零电荷双重保险，确保它对计算结果零贡献。

### 2.2 数据结构

对一组同类型的 group（如所有带电基团），构建三个矩阵：

| 矩阵 | shape | dtype | 说明 |
|:-----|:------|:------|:-----|
| `group_indices` | (n_groups, max_atoms) | int | 每 group 的全局原子索引，padding 位填 0 |
| `group_charges` | (n_groups, max_atoms) | float | 每原子电荷，padding 位填 0 |
| `group_valid` | (n_groups, max_atoms) | bool | 有效位，padding 位为 False |

其中 `max_atoms` = 该组 group 中原子数的最大值。

### 2.3 构建过程

构建在 `detect()` 入口处执行一次，与帧数无关。

输入：group 列表（已按 required_group_types 过滤）

过程：

1. 找出 max_atoms（该组 group 中最大的原子数）
2. 对每个 group i，将其原子的全局索引填入 `group_indices[i, 0:n_i]`，剩余位置填 0
3. 对应位置填入电荷到 `group_charges[i, 0:n_i]`，剩余位置填 0
4. 对应位置设 `group_valid[i, 0:n_i] = True`，剩余位置为 False

输出：三个矩阵，形状均为 (n_groups, max_atoms)

### 2.4 每帧计算过程

输入：`positions` (n_atoms_total, 3) 当前帧全部原子坐标

**第一步：一次性取所有 group 的原子坐标**

用 `group_indices` 矩阵对 `positions` 做高级索引，得到 (n_groups, max_atoms, 3) 的坐标矩阵。

原理：numpy 高级索引——用 (n_groups, max_atoms) 的整数矩阵去索引 (n_atoms_total, 3) 的数组，结果 shape 为 (n_groups, max_atoms, 3)。

padding 位置（index=0）取到 atom 0 的坐标，是垃圾值。

**第二步：加权 + 掩码清零**

对坐标矩阵乘以电荷矩阵和 valid 掩码：

- 电荷为 0 的位置：加权后自然为 0
- valid 为 False 的位置：乘以 0，强制清零

双重保险确保 padding 位置对后续求和零贡献。

**第三步：聚合**

沿原子维度（axis=1）求和，得到每个 group 的电荷中心 (n_groups, 3)。

分母：有效电荷之和 (n_groups,)
分子：有效加权坐标之和 (n_groups, 3)

### 2.5 安全性论证

padding 位置的安全性由三层保障确保：

| 保障层 | 机制 | 效果 |
|:-------|:-----|:-----|
| 第一层 | padding 位置 charge = 0 | 加权坐标 = 坐标 × 0 = 0 |
| 第二层 | padding 位置 valid = False | 乘以掩码后强制为 0 |
| 第三层 | 分母只含有效电荷 | 归一化不受 padding 影响 |

即使 padding 位置的 index 指向了真实的 atom 0（垃圾坐标），三层保障确保它对最终结果零贡献。

---

## 3. 两步向量化架构

### 3.1 整体流程

```
detect() 入口（一次性）：
  ├─ 构建 padding 索引矩阵（group 级）
  └─ 构建 tuple 索引数组（tuple 级）

每帧（零循环）：
  ├─ 第一步：group 级向量化
  │   positions + padding 矩阵 → group 特征（电荷中心、环心、法向量等）
  │
  └─ 第二步：tuple 级向量化
      group 特征 + tuple 索引 → tuple 指标（距离、角度等）
```

### 3.2 第一步：group 级向量化

目的：计算每个 group 的帧相关特征量。

输入：padding 索引矩阵（构建一次）+ 当前帧 positions（每帧更新）

输出：group 特征矩阵，shape 为 (n_groups, ...) 的矩形数组

操作：纯 numpy 广播 + 聚合，零 Python 循环。

### 3.3 第二步：tuple 级向量化

目的：用 group 特征组装每个 tuple 的指标。

输入：group 特征（第一步输出）+ tuple 索引数组

tuple 索引数组：(n_tuples,) 的整数数组，记录每个 tuple 引用了哪些 group。

操作：numpy 高级索引取出 group 特征，做向量运算，得到 (n_tuples,) 的指标。

示例（盐桥距离）：

| 步骤 | 操作 | shape |
|:-----|:-----|:------|
| 取 pos 中心 | `pos_centers = group_centers[pos_indices]` | (n_tuples, 3) |
| 取 neg 中心 | `neg_centers = group_centers[neg_indices]` | (n_tuples, 3) |
| 算距离 | `distances = norm(pos_centers - neg_centers, axis=1)` | (n_tuples,) |

整个过程零循环。

---

## 4. 各检测器适配方案

### 4.1 盐桥（Salt Bridge）

| 项 | 值 |
|:---|:---|
| group 类型 | charged_positive, charged_negative |
| max_atoms | 9（ARG 胍基） |
| padding 矩阵 | group_indices, group_charges, group_valid |
| group 级量 | 电荷中心 (n_groups, 3) |
| tuple 级量 | pos_center - neg_center → 距离 (n_tuples,) |

### 4.2 氢键（Hydrogen Bond）

| 项 | 值 |
|:---|:---|
| group 类型 | H_donor (atoms=[D, H]), H_acceptor (atoms=[A]) |
| padding 矩阵 | donor 矩阵 (n_donors, 2), acceptor 矩阵 (n_acceptors, 1) |
| group 级量 | D 坐标 (n_donors, 3), H 坐标 (n_donors, 3), A 坐标 (n_acceptors, 3) |
| tuple 级量 | D-A 距离, D-H-A 角度 (n_tuples,) |

注：donor 固定 2 原子，acceptor 固定 1 原子，padding 实际无浪费。

### 4.3 π-堆积（Pi Stacking）

| 项 | 值 |
|:---|:---|
| group 类型 | aromatic_ring |
| max_atoms | 7（Trp 双环取单环最大 6，但需考虑稠合环） |
| padding 矩阵 | group_indices, group_valid（无电荷，用 valid 掩码） |
| group 级量 | 环心 (n_rings, 3), 法向量 (n_rings, 3) |
| tuple 级量 | 环心距, 法向量夹角, 投影偏移 (n_tuples,) |

~~原问题：`np.roll` + padding 破坏法向量计算。~~ **已修复**：用显式 circular neighbor indices（`prev_idx`/`next_idx`）替代 `np.roll`。见 §6.1。

### 4.4 π-阳离子（Pi Cation）

| 项 | 值 |
|:---|:---|
| group 类型 | aromatic_ring, charged_positive |
| padding 矩阵 | 两套：ring 矩阵 + charged 矩阵 |
| group 级量 | 环心 + 法向量（ring 组）, 电荷中心（charged 组） |
| tuple 级量 | 环心-电荷中心距, 投影偏移 (n_tuples,) |

~~同 §4.3 的 np.roll 问题。~~ **已修复**：同 §4.3。

### 4.5 疏水（Hydrophobic）

| 项 | 值 |
|:---|:---|
| group 类型 | hydrophobic |
| max_atoms | 1（每 group 只有一个原子） |
| padding 矩阵 | 不需要，直接用 group 原子索引数组 |
| group 级量 | 直接取坐标 (n_groups, 3) |
| tuple 级量 | 距离 (n_tuples,) |

~~原问题：全排列组合数量巨大，距离矩阵内存爆炸。~~ **已修复**：用 `KDTree.query_pairs()` 做空间预筛选。见 §6.2。

### 4.6 卤键（Halogen Bond）

| 项 | 值 |
|:---|:---|
| group 类型 | halogen_donor (atoms=[C, X]), halogen_acceptor (atoms=[A, R1, R2, ...]) |
| padding 矩阵 | donor 矩阵 (n_donors, 2), acceptor 矩阵 (n_acceptors, max_r+1) |
| group 级量 | C/X 坐标 (n_donors, 2, 3), A/R 坐标 (n_acceptors, max_r+1, 3) |
| tuple 级量 | X-A 距离, C-X-A 角度, X-A-R 角度 (n_tuples,) |

~~原问题：acc_angle 极值选择被 padding 污染。~~ **已修复**：acceptor R 部分用环形填充（circular fill），padding 位坐标是真实 R 原子的重复，`argmin` 结果正确。见 §6.3。

### 4.7 金属配位（Metal Coordination）

| 项 | 值 |
|:---|:---|
| group 类型 | metal, metal_binding |
| max_atoms | 1（每 group 只有一个原子） |
| padding 矩阵 | 不需要 |
| group 级量 | 直接取坐标 |
| tuple 级量 | 距离 (n_tuples,) |

### 4.8 水桥（Water Bridge）

| 项 | 值 |
|:---|:---|
| group 类型 | H_donor (2 原子), water (3 原子), H_acceptor (不定) |
| padding 矩阵 | 三套：donor (n_donors, 2), water (n_waters, 3), acceptor (n_acceptors, max_a) |
| group 级量 | D/H 坐标, Ow 坐标, A 坐标 |
| tuple 级量 | dist_dw, dist_wa, theta, omega (n_tuples,) |

~~原问题：候选三元组生成需要 KDTree 空间预筛选。~~ **已修复**：`detect()` 入口用 `KDTree.query_ball_point()` 对水分子做空间预筛选（水到供/受体 < 8.2Å）。见 §6.4。

---

## 5. 与 PerTuple 策略的对比

| 维度 | PerTuple | PerFrame (padding) |
|:-----|:---------|:-------------------|
| 轨迹读取 | n_tuples 次（串行） | 1 次 |
| group 坐标获取 | 每 tuple 独立取 | 每帧一次性取所有 group |
| 计算方式 | 每 tuple 独立循环 | 零循环 numpy 广播 |
| 内存 | 每 tuple 存 (F, n_atoms, 3) | 每帧存 (n_groups, max_atoms, 3) |
| 适用场景 | tuple 少、帧多 | tuple 多（如水桥） |

---

## 6. 对抗性审查：已知问题（全部已修复）

### 6.1 ✅ 已修复：`np.roll` + padding 破坏芳香环法向量

**影响**：§4.3 π-堆积、§4.4 π-阳离子

**问题**：`np.roll` 是环形移位，不知道哪些位置是 padding。对一个 5 元环 padding 到 7：

```
有效原子:  [a0, a1, a2, a3, a4,  -,  -]
np.roll(..., -1) 得到后邻居:
           [a1, a2, a3, a4,  -,  -,  a0]
                                      ↑ a4 的后邻居是 padding 位 ✗
                                     ↑ padding 位的后邻居是 a0 ✗
```

- a4 的法向量因错误邻居而错误
- valid=False 能消除 padding 位自身的贡献，但不能消除 a4 的错误

**修复方案**：用显式 circular neighbor indices 替代 `np.roll`。

**实现**（`pi_stacking_detector_per_frame.py` 的 `_build_circular_padding`）：

```python
# 对每个位置 j，显式计算它在环内的前/后邻居（只在有效原子数 n 内取模）
for j in range(max_atoms):
    prev_j = (j - 1) % n          # 只在有效原子数 n 内取模
    next_j = (j + 1) % n
    prev_idx[i, j] = atom_indices[prev_j]   # 存全局原子索引
    next_idx[i, j] = atom_indices[next_j]
```

`_ring_normals` 用预计算的 `self._ring_prev` / `self._ring_next` 索引取坐标，不依赖 `np.roll`。padding 位 j=5 的 `prev_j = (5-1) % 5 = 4`，仍指向有效原子 a4。padding 位自身有 `valid=False` 掩码清零，不影响环法向量。

### 6.2 ✅ 已修复：疏水距离矩阵内存爆炸

**影响**：§4.5 疏水

**问题**：D927 体系 ~10,000 个疏水原子。全排列 C(10000, 2) ≈ 5000 万对，距离矩阵 ~1.2 GB。

**修复方案**：用 `KDTree.query_pairs()` 做空间预筛选。

**实现**（`hydrophobic_detector_per_frame.py`）：

```python
tree = KDTree(h_pos)
pairs = tree.query_pairs(r=PREFILTER_CUTOFF, output_type='ndarray')
```

`PREFILTER_CUTOFF = 8.0 Å`（阈值 4.0 Å × 2），只返回距离 ≤ 8Å 的原子对，将 5000 万对降至可控数量。

### 6.3 ✅ 已修复：卤键 acc_angle 被 padding 污染

**影响**：§4.6 卤键

**问题**：acc_angle 对每个 R 计算角度后取最接近 120° 的。零填充方案中 padding 位坐标指向 atom 0（垃圾），计算出的角度可能恰好接近 120°，被 `argmin` 选中。

**修复方案**：acceptor R 部分用**环形填充**（circular fill）而非零填充。

**实现**（`halogen_bond_detector_per_frame.py` 的 `_build_acceptor_padding`）：

```python
# position 0 = A，position 1+ = R 部分环形填充
indices[i, 0] = atom_indices[0]
r_indices = atom_indices[1:]
for j in range(1, max_atoms):
    indices[i, j] = r_indices[(j - 1) % n_r]   # 循环重复真实 R 原子
```

padding 位的坐标是真实 R 原子的循环重复（如 `[R1, R2, R1, R2]`），不是 atom 0 的垃圾坐标。`argmin` 选到 padding 位的结果与选到对应的原始 R 位完全相同，角度值正确。

### 6.4 ✅ 已修复：水桥丢失 KDTree 空间预筛选

**影响**：§4.8 水桥

**问题**：水桥的性能瓶颈是候选三元组生成（24.9 万个），不是每帧计算。纯 padding 方案不解决候选生成问题。

**修复方案**：`detect()` 入口用 KDTree 做空间预筛选。

**实现**（`water_bridge_detector_per_frame.py`）：

```python
donor_tree = KDTree(d_pos)
acceptor_tree = KDTree(a_pos)

for wi in range(len(waters)):
    nearby_d = donor_tree.query_ball_point(wp, PREFILTER_CUTOFF)   # 8.2 Å
    nearby_a = acceptor_tree.query_ball_point(wp, PREFILTER_CUTOFF)
    for di in nearby_d:
        for ai in nearby_a:
            # 过滤 D-A 距离过近的情况
            if np.linalg.norm(d_pos[di] - a_pos[ai]) < WATER_BRIDGE_MINDIST:
                continue
            triple_d.append(di); triple_w.append(wi); triple_a.append(ai)
```

先用 KDTree 找每个水分子附近的供体和受体（< 8.2Å），再组合成三元组。将 24.9 万全排列降至实际候选数量。

---

## 7. 基类设计策略调整

### 7.1 问题：基类抽象程度过高

原方案设想基类封装 padding 矩阵构建，子类只用 tuple 级向量化。但不同检测器的差异太大：

| 差异点 | 盐桥 | π-堆积 | 卤键 | 水桥 |
|:-------|:-----|:-------|:-----|:-----|
| group 级量 | 电荷中心 | 环心+法向量 | C/X/A/R 坐标 | D/H/Ow/A 坐标 |
| 邻居关系 | 无 | 有（显式索引） | 有（环形填充） | 无 |
| 极值选择 | 无 | 无 | 有（argmin） | 无 |
| 空间预筛选 | 不需要 | 不需要 | 不需要 | 必须（KDTree） |

强行统一会导致基类过度复杂或子类无法表达特有逻辑。

### 7.2 调整方案：子类驱动

基类只提供帧迭代骨架，子类控制全部计算逻辑：

```
detect() 基类骨架：
  ├─ get_candidate_tuples(groups)        # 子类实现
  ├─ tuple_filter                        # 用户过滤
  ├─ filter_candidate_tuples(...)        # 子类可选覆盖
  │
  ├─ for each frame:
  │    └─ compute_metrics_for_frame(tuples, positions, frame)  # 子类全权控制
  │
  ├─ _post_process(results)              # 子类可选覆盖
  └─ _build_interaction(results)         # 基类固化
```

子类自行决定：
- 是否用 padding 矩阵
- 如何计算 group 级量
- 如何组装 tuple 级量
- 是否需要空间预筛选（KDTree 等）

padding 索引矩阵是子类可选的**实现技巧**，不是基类强制的接口。

### 7.3 基类 vs 子类职责

| 职责 | 归属 | 说明 |
|:-----|:-----|:-----|
| 帧迭代骨架 | 基类 | `detect()` 的 for 循环 + 结果累积 |
| 候选 tuple 生成 | 子类 | `get_candidate_tuples()` |
| 坐标预过滤 | 子类 | `filter_candidate_tuples()` |
| 每帧指标计算 | 子类 | `compute_metrics_for_frame()` |
| 阈值判定 | 子类 | `apply_threshold()` |
| 后处理 | 子类 | `_post_process()` |
| Interaction 构建 | 基类 | `_build_interaction()` |
| padding 矩阵 | 子类自行决定 | 可用可不用，基类不感知 |

---

## 8. 构建开销

| 项 | 复杂度 | 说明 |
|:---|:-------|:-----|
| group_indices | O(n_groups × max_atoms) | 一次构建，每帧复用 |
| group_charges | O(n_groups × max_atoms) | 同上 |
| group_valid | O(n_groups × max_atoms) | 同上 |
| tuple 索引数组 | O(n_tuples) | 一次构建 |

内存：盐桥 ~20 groups × 9 atoms × 3 cols × 8 bytes ≈ 4 KB，忽略不计。

---

## 9. 设计决策

| 决策 | 选择 | 理由 |
|:-----|:-----|:-----|
| padding 位置的 index 填什么 | 0（指向 atom 0） | 任意值均可，valid 掩码确保无影响 |
| padding 位置的 charge 填什么 | 0 | 配合 valid 掩码双重保障 |
| 矩阵在何时构建 | 子类自行决定 | 可在 `__init__` 或 `detect` 入口 |
| 矩阵存放在哪 | 子类实例属性 | 子类全权控制 |
| padding 是否是基类接口的一部分 | 否 | 是子类可选的实现技巧 |

---

*文档结束*
