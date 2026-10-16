# PBC（周期边界）处理设计文档（pbc_design）

> 创建日期：2026-10-14
> 分支：feature/pbc-handling
> 状态：方案定稿（待实现）——仅覆盖"纯距离类"相互作用（疏水、金属配位），
>       两种策略（two_pass / per_frame）；per_tuple 不在范围内（将被舍弃）。
> 参考：MDAnalysis lib.distances 官方文档（本地 2.10.0 实测验证）
> 范围演进：本版只含纯距离类；氢键/水桥/卤键/盐桥/π 类另开分支或追加章节。

---

## 1. 背景与问题

### 1.1 现状

现有检测器的距离/角度计算全部使用 **scipy KDTree / 手写 np.linalg.norm / 手写向量夹角**，
未传盒子尺寸——即**纯欧氏距离，无 PBC 处理**。当两个原子位于盒子边界两侧时
（如 x=0.1 与 x=L-0.1，盒长 L），欧氏距离 ≈ L-0.2Å 而被判定"无相互作用"，
但真实最小镜像距离 ≈ 0.2Å 应判定"有相互作用"——**跨盒相互作用被系统性漏检**。

### 1.2 实证（本地真实数据 KRAS-RBD，盒 105.78Å）

- 20000 原子子集：欧氏 ≤4.1Å 对 264,942 个；PBC(MIC) ≤4.1Å 对 271,986 个
- **PBC 新增 3,522 个跨盒近邻对（欧氏漏检）**——真实体系存在大量跨盒近邻
  （水分子在盒边缘聚集），漏检是实质性的，非浮点误差。

### 1.3 正确做法

最小镜像约定（MIC）：对每对原子取其在周期性重复中距离最近的那个像。
MDAnalysis `MDAnalysis.lib.distances` 的 `capped_distance` / `self_capped_distance` /
`calc_bonds` / `calc_angles` 传 `box=ts.dimensions` 即自动应用 MIC，
官方支持正交与三斜盒（本地 2.10.0 实测可用）。

## 2. 核心 API 与语义（本地实测验证）

| API | 语义 | 实测确认 |
|---|---|---|
| `capped_distance(ref, conf, max_cutoff, min_cutoff=None, box=None, method=None, return_distances=True)` | 返回 `(pairs[N,2] int64, distances[N] float64)`；区间 **(min, max]** 闭区间含 max；`pairs[:,0]=ref 索引, pairs[:,1]=conf 索引` | 4.1 精确被捕获（含边界） |
| `self_capped_distance(ref, max_cutoff, min_cutoff=None, box=None, method=None, return_distances=True)` | 自查询；**自动去重 i<j、无自对**；返回 `(pairs, distances)` | pairs=[[0,1],[0,2],[1,2]] |
| `calc_bonds(coords1, coords2, box=None)` | 逐对距离（**一对一，两数组长度必须相等**），MIC | 长度不等报错（防误用） |
| `calc_angles(c1, c2, c3, box=None)` | 三原子夹角，顶点是 c2，返回**弧度**；MIC 应用于构成角的向量 | 180° 与手写一致 |

**关键决策**：
- **显式 `method='pkdtree'`**——不依赖自动猜（nsgrid 有历史可靠性问题，
  MDAnalysis 1.0.1 曾因不可靠临时移除；pkdtree 与 scipy cKDTree 同源、可靠）。
- **box 防护**（无盒体系退化欧氏，见 §4 `_get_box`）。

## 3. 范围（本版）

本版只覆盖**纯距离类**相互作用：

| 检测器 | two_pass | per_frame | 距离计算点 |
|---|---|---|---|
| 疏水 hydrophobic | ✅ | ✅ | 4 处（two_pass Pass1/Pass2 + per_frame 预筛/逐帧） |
| 金属配位 metal_coordination | ✅ | ✅ | 4 处（two_pass Pass1/Pass2 + per_frame 预筛/逐帧） |

- **per_tuple：不动**（将被舍弃；`--strategy per_tuple` 保持非 PBC，结果与
  per_frame/two_pass 在跨盒对上有差异，limitations 注明）。
- 氢键/水桥/卤键/盐桥/π-堆积/π-阳离子：**不在本版**（角度类、聚合坐标类、
  派生几何量类另开分支/追加章节）。

## 4. 统一 box 获取函数（新增）

```python
def _get_box(ts):
    """取帧盒尺寸；无盒/退化盒返回 None（退化为欧氏距离，PBC 不生效）。"""
    dims = ts.dimensions
    if dims is None or dims[0] <= 0:
        return None
    return dims
```

- 有盒（amber/charmm/gromos 真实数据均实测有盒）→ MIC PBC 生效
- 无盒 → 退化为欧氏（与现状一致），不报错

## 5. 疏水相互作用 PBC 方案

### 5.1 判据（两种策略一致，apply_threshold 零改动）

```
0.5 Å < distance < 4.0 Å    （严格开区间，两端都不含）
HYDROPH_DIST_MAX = 4.0, HYDROPH_MIN_DIST = 0.5
```

### 5.2 数据流差异（决定改动面）

| 环节 | two_pass | per_frame |
|---|---|---|
| 候选发现 | 逐帧 `_kdtree_discover`（每帧建树、每帧发现近邻对） | 第一帧一次 `query_pairs(8.0)` 预筛固定候选对 |
| 指标计算 | Pass1 记「活帧+距离」→ Pass2 逐帧重算全帧距离 | 固定候选对逐帧向量化算全帧距离 |
| 去重 | 未实现（TODO 记录，同残基多原子只留最近；**不在本任务范围**） | `_deduplicate`（已实现） |

### 5.3 two_pass 替换

**Pass1 `_kdtree_discover(all_positions, box)`**：

```python
# 现实现（建树 → sparse 自查询 → 去重 i<j → 最小距离过滤）
from scipy.spatial import cKDTree
coords = all_positions[self._atom_indices]
tree = cKDTree(coords)
sparse_dist = tree.sparse_distance_matrix(tree, HYDROPH_DIST_MAX)
coo = sparse_dist.tocoo()
mask = coo.row < coo.col
pairs = (coo.row[mask], coo.col[mask])
distances = coo.data[mask]
dist_mask = distances > HYDROPH_MIN_DIST
return (pairs[0][dist_mask], pairs[1][dist_mask]), distances[dist_mask]

# 替换（一次调用：自查询 + 距离区间 + 去重 i<j + 无自对，全部自带）
from MDAnalysis.lib.distances import self_capped_distance
pairs, distances = self_capped_distance(
    coords, HYDROPH_DIST_MAX, min_cutoff=HYDROPH_MIN_DIST,
    box=box, method='pkdtree')
return (pairs[:, 0], pairs[:, 1]), distances
```

**逐点核对**：
- ✅ 去重：`self_capped_distance` 自动 i<j、无自对（实测）——与原 `mask = row < col` 一致
- ✅ 距离区间：`(0.5, 4.0]`（闭区间含 4.0）。原 sparse 也是 `<=4.0` 含边界；
  Pass2 `apply_threshold` 用 `< 4.0` 严格否决边界 → 新老最终判定一致
- ⚠️ `_kdtree_discover` 需新增 `box` 参数（由 `run_pass1` 逐帧传 `_get_box(ts)`）

**Pass2 `compute_pair_metrics`**：

```python
# 现实现
pos1 = all_positions[self._pair_idx1]
pos2 = all_positions[self._pair_idx2]
return {"distance": np.linalg.norm(pos1 - pos2, axis=1)}

# 替换
return {"distance": calc_bonds(pos1, pos2, box=box)}
```

### 5.4 per_frame 替换（自写 detect，无需改基类）

**候选预筛（第一帧一次）**：

```python
# 现实现
first_pos = trajectory[0].positions
h_pos = first_pos[self._hydro_idx]
tree = KDTree(h_pos)
pairs = tree.query_pairs(r=PREFILTER_CUTOFF, output_type='ndarray')  # 8.0 Å

# 替换
# 替换（return_distances=False 只返回 pairs 单个数组，不解包二元组）
pairs = self_capped_distance(
    h_pos, PREFILTER_CUTOFF,
    box=_get_box(trajectory[0]), method='pkdtree', return_distances=False)
pair_a, pair_b = pairs[:, 0], pairs[:, 1]
```

**逐帧距离计算**：

```python
# 现实现
dist = np.linalg.norm(pos_a - pos_b, axis=1)
# 替换
dist = calc_bonds(pos_a, pos_b, box=_get_box(ts))
```

**边界说明**：`query_pairs(r=8.0)` 与 `self_capped_distance(max=8.0)` 在恰 8.0 边界
可能差 1 对（预筛 cutoff 8.0 是安全裕度，实际判距 <4.0），Pass2 逐帧判距还会过滤，
无科学影响；文档注明即可。

## 6. 金属配位 PBC 方案

### 6.1 判据（两种策略一致）

```
distance < 3.0 Å    （严格小于）
METAL_DIST_MAX = 3.0
```

### 6.2 数据流差异（重要：per_frame 不用 KDTree）

| 环节 | two_pass | per_frame |
|---|---|---|
| 候选发现 | 逐帧 `_kdtree_discover`（cKDTree 交叉查询） | **笛卡尔积全候选**（无 KDTree）+ 第一帧 `norm` 预筛 |
| 指标计算 | Pass1 记「活帧+距离」→ Pass2 逐帧重算 | 固定候选逐帧向量化 |

金属 per_frame 的 PBC 化**无需候选发现替换**（笛卡尔积本身跨盒正确，候选对是全组合），
只需**两处距离计算 PBC 化**（第一帧预筛 norm + 逐帧 norm）。

### 6.3 two_pass 替换

**Pass1 `_kdtree_discover(all_positions, box)`**：

```python
# 现实现
m_coords = all_positions[self._metal_indices]
b_coords = all_positions[self._binding_indices]
tree_m = cKDTree(m_coords)
tree_b = cKDTree(b_coords)
sparse_dist = tree_m.sparse_distance_matrix(tree_b, METAL_DIST_MAX)
if sparse_dist.nnz == 0:
    return np.array([]), np.array([]), np.array([])
coo = sparse_dist.tocoo()
return coo.row, coo.col, coo.data

# 替换
pairs, distances = capped_distance(
    m_coords, b_coords, METAL_DIST_MAX,
    box=box, method='pkdtree')
# 空对判断：capped_distance 空对返回 pairs.shape=(0,2)（实测），用 len(pairs)==0
if len(pairs) == 0:
    return np.array([]), np.array([]), np.array([])
return pairs[:, 0], pairs[:, 1], distances
```

**逐点核对**：
- ✅ 顺序：`capped_distance(m_coords, b_coords)` → `pairs[:,0]=m 索引`（与
  原 `sparse_distance_matrix(tree_m, tree_b)` 的 `coo.row=m` 一致）
- ✅ 边界：`<=3.0` 含 3.0；Pass2 `apply_threshold` 用 `<3.0` 严格否决 → 新老一致

**Pass2 `compute_pair_metrics`**：

```python
# 现实现
m_pos = all_positions[self._pass2_metal_idx]
b_pos = all_positions[self._pass2_binding_idx]
return {"distance": np.linalg.norm(m_pos - b_pos, axis=1)}

# 替换
return {"distance": calc_bonds(m_pos, b_pos, box=box)}
```

### 6.4 per_frame 替换（自写 detect，无需改基类）

**第一帧预筛**：

```python
# 现实现
dist = np.linalg.norm(m_pos - b_pos, axis=1)
mask = dist < PREFILTER_CUTOFF

# 替换
dist = calc_bonds(m_pos, b_pos, box=_get_box(trajectory[0]))
mask = dist < PREFILTER_CUTOFF
```

**逐帧计算**：

```python
# 现实现
dist = np.linalg.norm(m - b, axis=1)
# 替换
dist = calc_bonds(m, b, box=_get_box(ts))
```

## 6A. 氢键相互作用 PBC 方案（2026-10-14 追加）

### 6A.1 判据（两种策略一致，apply_threshold 零改动）

```
D-A 距离 ≤ 4.1 Å 且 D-H···A 角度 ≥ 100°
HBOND_DIST_MAX = 4.1, HBOND_DON_ANGLE_MIN = 100.0
```

注意与疏水/金属不同：氢键距离是 **≤ 4.1 闭区间**（非 strict <），
`capped_distance`/`calc_bonds` 边界含 4.1 → 判定一致。

### 6A.2 数据流（PBC 化计算点共 6 处）

| 环节 | two_pass | per_frame |
|---|---|---|
| 候选发现 | 逐帧 `_kdtree_discover`（cKDTree 交叉查询 sparse 4.1） | 第一帧 `query_ball_point(a_pos, 12.3)` 预筛固定候选 |
| 指标计算 | Pass1 记「活帧+距离/角度」→ Pass2 逐帧重算（norm + 手写角度） | 固定候选逐帧向量化（norm + 手写角度） |

### 6A.3 角度 PBC 化（calc_angles 顶点核对）

```
现有手写（两策略一致）：hd_vec = d - h; ha_vec = a - h; cos = ... → degrees(arccos(clip))
替换：ang = np.degrees(calc_angles(d, h, a, box=box))
```
- `calc_angles(c1, c2, c3)` 顶点是 **c2** → `calc_angles(d, h, a)` 计算 H 处 D-H···A 角，
  与现有定义完全一致；PBC 下自动对构成角的两向量做 MIC。
- ⚠️ 退化构型（D==H 或 A==H）：现有手写产生 NaN，`calc_angles` 返回 0（实测）——
  两者 `apply_threshold` 均判 False（NaN≥100 与 0≥100 均假），**判定一致**；
  h5 metrics 数值从 NaN 变 0。氢键 D-H 是键连的（供体 N-H/O-H），原子重合物理
  不可能，实际不可达；实现后用既有测试确认无回归。

### 6A.4 two_pass 替换

**Pass1 `_kdtree_discover(all_positions, box)`（逐帧）**：

```python
# 现：tree_d/tree_a cKDTree + sparse_distance_matrix(tree_d, tree_a, 4.1)
#     coo.row=d 索引, coo.col=a 索引；角度手写
# 改
pairs, distances = capped_distance(
    d_coords, a_coords, HBOND_DIST_MAX, box=box, method='pkdtree')
active_donors, active_acceptors = pairs[:, 0], pairs[:, 1]
angles = np.degrees(calc_angles(
    d_coords[active_donors],
    all_positions[self._donor_h_idx[active_donors]],
    a_coords[active_acceptors], box=box))
```

- 顺序核对：`capped_distance(d_coords, a_coords)` → `pairs[:,0]=d`（与
  `sparse_distance_matrix(tree_d, tree_a)` 的 `row=d` 一致）✓

**Pass2 `compute_pair_metrics(..., box)`（逐帧，run_pass2 自写循环传 box）**：

```python
# 现：dist = norm(d-a)；手写角度
# 改
dist = calc_bonds(d, a, box=box)
ang = np.degrees(calc_angles(d, h, a, box=box))
```

`run_pass2` 自写循环同疏水/金属模式（基类零改动）：逐帧 `box=get_box(ts)` 传入，
结果构建与基类一致（不过滤、`_post_process` + `_build_interaction`）。

### 6A.5 per_frame 替换

**候选预筛（第一帧一次）**：

```python
# 现：tree = KDTree(d_pos); pairs = tree.query_ball_point(a_pos, PREFILTER)
#     pair_donor.append(di); pair_acceptor.append(ai)  # 按 acceptor 组织
# 改
pairs, _ = capped_distance(
    a_pos, d_pos, PREFILTER_CUTOFF,
    box=get_box(trajectory[0]), method='pkdtree', return_distances=False)
pair_acceptor = pairs[:, 0]   # 注意方向：ref=a → pairs[:,0]=acceptor
pair_donor = pairs[:, 1]      #               conf=d → pairs[:,1]=donor
```

⚠️ **方向核对**：`query_ball_point(a_pos, ...)` 语义是"每 acceptor 的邻近 donor"；
替换为 `capped_distance(a_pos, d_pos)` → `pairs[:,0]=acceptor, pairs[:,1]=donor`，
后续变量赋值需对应（现在 `pair_donor.append(di)` 对应 pairs[:,1]）。

**逐帧 `_compute_metrics(positions, box)`**：

```python
# 现：dist = norm(d-a)；手写角度
# 改
dist = calc_bonds(d_pos, a_pos, box=box)
ang = np.degrees(calc_angles(d_pos, h_pos, a_pos, box=box))
```

调用处逐帧 `_compute_metrics(positions, get_box(ts))`。

### 6A.6 依赖与兼容

- 基类接口 **零改动**（沿用疏水/金属模式：two_pass 自写 run_pass2）
- `apply_threshold` / 阈值常量 / `metric_names` **零改动**（`≤4.1 & ≥100` 原样）
- 其他检测器零影响；per_tuple 不动

### 6A.7 测试策略

1. **PBC 单测（test_pbc.py 追加）**：
   - 跨盒氢键：D 与 A 跨盒（MIC ≤4.1）→ 存在 True；欧氏 >4.1 → 旧漏检
   - 跨盒角度：H 在盒边、D/A 分居两侧 → `calc_angles(box)` 与回绕后手算一致
   - 边界：恰 4.1 → `≤4.1` 保留（闭区间语义）
2. **既有测试**：`test_hydrogen_bond_per_frame.py` / `test_hydrogen_bond_two_pass.py`
   跑通；真实数据（KRAS-RBD 86 氢键对）新旧逐 pair 逐帧 existence 对比应一致
3. **性能**：donor/acceptor 数量级小（KRAS-RBD 数百），pkdtree 每帧开销可接受

## 6B. 水桥相互作用 PBC 方案（2026-10-14 追加）

### 6B.1 判据（两种策略一致，apply_threshold 零改动）

```
dist_dw < 4.1 且 dist_wa < 4.1      （水 O 到供体 D / 受体 A 距离）
theta ≥ 100°                        （水O-供体H-供体D）
71° ≤ omega ≤ 140°                  （受体-水O-供体H）
WATER_BRIDGE_MAXDIST=4.1, THETA_MIN=100, OMEGA_MIN=71, OMEGA_MAX=140
```

### 6B.2 数据流（PBC 化计算点）

| 环节 | two_pass | per_frame |
|---|---|---|
| 候选发现 | 逐帧三棵 cKDTree（d/a/w）→ 两两 sparse（w↔d、w↔a）→ **按水分子分组**组合三元组 → 角度初筛 | 第一帧逐水 `query_ball_point(8.2)` 收集近邻 d/a + `norm(d-a)<2.5` 排除 → 固定候选三元组 |
| 指标存储 | Pass1 存 dist_dw/dist_wa/theta/omega 入 sparse；Pass2 重算全帧 4 指标 | 逐帧 `_compute_metrics` 向量化 4 指标 |
| 角度定义 | theta: `_vec_angle(d-h, ow-h)` H 顶点；omega: `_vec_angle(a-ow, h-ow)` Ow 顶点 | 同左 |

### 6B.3 角度顶点核对（calc_angles）

| 指标 | 现有手写 | calc_angles 替换 | 顶点 |
|---|---|---|---|
| theta | `_vec_angle(d-h, ow-h)` | `np.degrees(calc_angles(ow, h, d, box))` | **h**（c2）✓ 实测一致 |
| omega | `_vec_angle(a-ow, h-ow)` | `np.degrees(calc_angles(a, ow, h, box))` | **ow**（c2）✓ |

### 6B.4 two_pass 替换

**Pass1 `_kdtree_discover(positions, box)`（逐帧）**：

```python
# 现：三棵 cKDTree + dw_sparse(tree_w, tree_d) + wa_sparse(tree_w, tree_a)
# 改
dw_pairs, dw_dists = capped_distance(
    w_coords, d_coords, WATER_BRIDGE_MAXDIST, box=box, method='pkdtree')
wa_pairs, wa_dists = capped_distance(
    w_coords, a_coords, WATER_BRIDGE_MAXDIST, box=box, method='pkdtree')
if len(dw_pairs) == 0 or len(wa_pairs) == 0:
    return tuple(np.array([]) for _ in range(7))

# 分组逻辑保留（dw_pairs[:,0]=w 索引，与 dw_coo.row 语义一致，实测确认）
dw_by_w = defaultdict(list)
for i in range(len(dw_pairs)):
    dw_by_w[dw_pairs[i, 0]].append((dw_pairs[i, 1], dw_dists[i]))
# wa_by_w 同；组合三元组逻辑不变
# 角度: theta = degrees(calc_angles(ow, h, d, box)); omega = degrees(calc_angles(a, ow, h, box))
# 阈值筛选不变（mask 保留）
```

⚠️ **方向核对**：现 `tree_w.sparse_distance_matrix(tree_d)` → row=w、col=d；
改 `capped_distance(w, d)` → pairs[:,0]=w、[:,1]=d ✓（实测确认唯一对、无重复像）

⚠️ **Pass1 角度初筛从欧氏变 MIC**：跨盒帧的角度可能不同 → 初筛存活集可能变化
（正是 PBC 修复目的）；Pass2 重算全帧 MIC 指标并 apply_threshold → 最终判定自洽。

**Pass2 `compute_pair_metrics(..., box)`（逐帧，run_pass2 自写循环）**：

```python
dist_dw = calc_bonds(d, ow, box=box)
dist_wa = calc_bonds(ow, a, box=box)
theta = np.degrees(calc_angles(ow, h, d, box=box))
omega = np.degrees(calc_angles(a, ow, h, box=box))
```

`run_pass2` 自写循环同前（逐帧 get_box(ts)），结果构建与基类一致。

### 6B.5 per_frame 替换

**候选预筛（第一帧一次）——结构调整**：

```python
# 现：for wi: donor_tree.query_ball_point(wp, 8.2) → 逐水收集近邻 d/a
# 改（求全体 w↔d、w↔a 对，按 w 分组）
dw_pairs = capped_distance(
    w_pos, d_pos, PREFILTER_CUTOFF, box=first_box, method='pkdtree', return_distances=False)
wa_pairs = capped_distance(
    w_pos, a_pos, PREFILTER_CUTOFF, box=first_box, method='pkdtree', return_distances=False)
# 按 w 分组 → 组合三元组（与 two_pass 分组逻辑同构）
# D-A 最短距排除: calc_bonds(d, a, box) < WATER_BRIDGE_MINDIST → 排除
```

⚠️ **结构调整**：现实现是"逐水 query_ball_point 嵌入循环"，替换为"全部对 + 按 w 分组"
——与 two_pass 分组逻辑同构，需重构候选生成循环（triple_d/w/a 列表构建方式变化，
但后续逐帧 `_compute_metrics` 只依赖 `_triple_*` 索引，无需改）。

**逐帧 `_compute_metrics(positions, box)`**：

```python
dist_dw = calc_bonds(d, ow, box=box)
dist_wa = calc_bonds(ow, a, box=box)
theta = np.degrees(calc_angles(ow, h, d, box=box))
omega = np.degrees(calc_angles(a, ow, h, box=box))
```

调用处逐帧 `_compute_metrics(positions, get_box(ts))`。

### 6B.6 性能实测（关键：37021 水）

KRAS-RBD 水分子 **37021**、供体 581 —— 水桥是 pkdtree 慢基准的主要场景，实测：

| 方法 | 耗时（w↔d 一帧） | 结论 |
|---|---|---|
| 旧三树 cKDTree sparse | 0.0093s | 基准 |
| **capped_distance(pkdtree, box)** | **0.0118s** | ✅ **仅慢 1.3×，可接受** |
| capped_distance 自动 method | 0.0398s | ❌ 慢 4.3×——**必须显式 pkdtree** |
| apply_PBC+scipy（备选） | 0.0884s | 不必要 |

**结论**：显式 `method='pkdtree'` 是水桥正确选择；候选对 1425 与旧完全一致；
备选方案（apply_PBC 回绕+scipy）实测更慢，不需要。

### 6B.7 测试策略

1. **PBC 单测（test_pbc.py 追加）**：
   - 跨盒水桥：水 O 与供体 D 跨盒（MIC dist_dw ≤4.1）→ 三元组存在
   - 跨盒 theta/omega：H/Ow 在盒边 → `calc_angles(box)` MIC 正确
2. **既有测试**：`test_water_bridge_per_frame.py` / `test_water_bridge_two_pass.py`
3. **真实数据对比**：KRAS-RBD（水桥候选多）新旧逐三元组 existence 对比——判定应一致
4. **性能实测**：水桥两策略全量耗时（PBC vs 旧）——重点项

## 6C. 卤键相互作用 PBC 方案（2026-10-14 追加）

### 6C.1 判据（两种策略一致，apply_threshold 零改动）

```
distance ≤ 4.0 Å                     （X···A 距离，闭区间）
don_angle ∈ [135°, 195°]             （C-X···A 角，最优 165°±30°）
acc_angle ∈ [90°, 150°]              （X···A-R 角，最优 120°±30°；R 从 acceptor
                                      多个邻居原子中取最接近 120° 的——argmin 择优）
HALOGEN_DIST_MAX=4.0, DON_ANGLE=165, ACC_ANGLE=120, ANGLE_DEV=30
```

### 6C.2 数据流（关键：无 KDTree，候选发现无需替换）

| 环节 | two_pass | per_frame |
|---|---|---|
| 候选发现 | `initialize_candidates` 笛卡尔积全候选（**无 KDTree**） | 自写 detect：笛卡尔积 + 第一帧 `norm < 12.0` 预筛 |
| 指标计算 | Pass2 逐帧 `compute_pair_metrics`（norm + 手写角度 + argmin） | 逐帧 `_compute_metrics`（同左） |
| acceptor R | `_build_acceptor_padding` 环形填充（A, R1, R2...）→ 每对多 R 候选 | 同左 |

**卤键两策略候选都是笛卡尔积全组合（无 KDTree）**：候选本身跨盒正确（不依赖
距离预筛），**无需候选发现替换**——只需三处指标 PBC 化（distance、don_angle、
acc_angle）+ per_frame 预筛距离。

### 6C.3 距离与 don_angle 替换（同前类型模式）

```python
# 距离
distance = calc_bonds(x, a, box=box)          # 替代 norm(x-a)

# don_angle：C-X···A（顶点 X）
don_angle = np.degrees(calc_angles(c, x, a, box=box))   # 替代 _vec_angle(c-x, a-x)
```

### 6C.4 ⚠️ acc_angle 替换（本类型独特：多 R 候选 argmin 展开）

现实现（两策略同构）：acceptor 有 `max_r` 个候选 R 邻居（padding，第 0 列是 A
自身需跳过），对每个 R 算 X···A-R 角，取最接近 120° 的作为 acc_angle。

```python
# 现（two_pass compute_pair_metrics）：
r_only = r[:, 1:, :]                      # (n_pairs, max_r, 3) 跳过 A 自身
ar = r_only - a[:, None, :]               # A→R 向量
ax = x - a                                # A→X 向量
cos = np.sum(ax[:, None, :] * ar, axis=2) / (norm_ax * norm_ar)
all_acc_angles = np.degrees(np.arccos(np.clip(cos, -1, 1)))   # (n_pairs, max_r)
best_r = np.argmin(np.abs(all_acc_angles - HALOGEN_ACC_ANGLE), axis=1)
acc_angle = all_acc_angles[np.arange(len(best_r)), best_r]

# 改：展开 padding → calc_angles 批量 MIC 角 → reshape → argmin 还原
x_rep = np.repeat(x, max_r, axis=0)                   # 每 pair 的 X 复制 max_r 次
a_rep = np.repeat(a, max_r, axis=0)
r_flat = r_only.reshape(-1, 3)                          # (n_pairs*max_r, 3)
all_acc = np.degrees(calc_angles(x_rep, a_rep, r_flat, box=box))  # 顶点 A
all_acc_angles = all_acc.reshape(n_pairs, max_r)        # 还原 (n_pairs, max_r)
best_r = np.argmin(np.abs(all_acc_angles - HALOGEN_ACC_ANGLE), axis=1)
acc_angle = all_acc_angles[np.arange(len(best_r)), best_r]
```

**关键细节（对抗性审查实测确认）**：
- **顶点核对**：X···A-R 角顶点是 **A**（c2）→ `calc_angles(x, a, r)` ✓（与手写
  `ar = r - a`（A→R）、`ax = x - a`（A→X）夹角一致；实测 maxdiff 1e-5）
- **展开-还原等价**：`np.repeat` + `r_only.reshape(-1,3)` + `calc_angles`
  （MIC 向量化）→ reshape 回 `(n_pairs, max_r)` → **argmin 结果与逐 R 手算
  完全一致**（实测 best_r 一致、最终 acc_angle 一致）
- **padding 第 0 列**：展开前同样 `r_only = r[:, 1:, :]`（跳过 A 自身），逻辑保留
- **argmin 择优语义**：保留（择优不变，仅角度来源 PBC 化）
- **跨盒 acc 角**：`calc_angles(x,a,r,box)` MIC 正确（实测 180° vs 欧氏错误 0°）

### 6C.5 two_pass 替换

- `run_pass2` 自写循环（同前模式）：逐帧 `box=get_box(ts)` 传入
  `compute_pair_metrics(..., box)`，结果构建与基类一致
- `compute_pair_metrics(group_tuples, all_positions, box=None)`：三处按 §6C.3/6C.4 替换

### 6C.6 per_frame 替换

- **第一帧预筛**（现 `norm < 12.0`）：
  ```python
  dist = calc_bonds(x_pos, a_pos, box=get_box(trajectory[0]))
  mask = dist < PREFILTER_CUTOFF
  ```
- **逐帧 `_compute_metrics(positions, box)`**：按 §6C.3/6C.4 三处替换，
  调用处逐帧 `box=get_box(ts)` 传入

### 6C.7 测试策略

1. **PBC 单测（test_pbc.py 追加）**：
   - 跨盒卤键：X 与 A 跨盒（MIC ≤4.0）→ 判定存在；欧氏 >4.0 旧漏检
   - 跨盒 acc_angle：A 在盒边、R 跨盒 → `calc_angles(x,a,r,box)` MIC 正确 + argmin 择优
   - 展开-还原正确性：`calc_angles` 展开结果 vs 逐 R 手算，argmin 一致
2. **既有测试**：`test_halogen_bond_per_frame.py` / `test_halogen_bond_two_pass.py`
3. **真实数据**：KRAS-RBD 无卤素（结果 0 对，断言不变）；跨盒需合成含卤素小体系验证
4. **合成 MIC 用例**：构造含卤素跨盒体系验证 distance/don_angle/acc_angle

## 6D. 盐桥相互作用 PBC 方案（2026-10-14 追加）

### 6D.1 判据（两种策略一致，apply_threshold 零改动）

```
正电基团电荷中心 ↔ 负电基团电荷中心 距离 ≤ 5.5 Å
SALTBRIDGE_DIST_MAX = 5.5（闭区间）
```

### 6D.2 数据流（关键：无 KDTree，候选发现无需替换）

| 环节 | two_pass | per_frame |
|---|---|---|
| 候选发现 | `initialize_candidates` 笛卡尔积全候选（**无 KDTree**） | 自写 detect：笛卡尔积 + 第一帧 `all_dist < PREFILTER_CUTOFF` 预筛 |
| 指标计算 | Pass2 逐帧 `compute_pair_metrics`：`_charge_centers`（**电荷中心加权平均**）→ norm | 逐帧 `_charge_centers` → norm |
| 电荷中心 | `_charge_centers(positions, idx, q, valid)`：加权平均位点 `Σ(q·r)/Σq` | 同左 |

**盐桥两策略候选都是笛卡尔积全组合（无 KDTree）**：候选本身跨盒正确，**无需候选
发现替换**——只需电荷中心聚合 + 中心距离 PBC 化。

### 6D.3 核心难点：电荷中心 PBC 化（对抗性审查实测确认）

**科学定义**：电荷中心是**同一残基内多原子**加权平均位点。同一残基原子共价键连
（间距 ~1.5Å），**正常残基不跨盒**；真正的跨盒风险是**两电荷中心之间的距离**。
两层处理：
1. **电荷中心本身**：先以**基团第一个原子为参考**做 MIC 回绕（保基团连续），再聚合
   ——防御性保证：残基横跨盒边界时原子被拉回同一像，中心正确；正常残基回绕无影响
   （实测：回绕中心 = 欧氏中心）
2. **两中心间距**：跨盒时 MIC → `calc_bonds(pos_c, neg_c, box)`

**⚠️ 实测发现：`minimize_vectors` 不支持三维**——只接受 `(3,)` 或 `(n, 3)`，
不接受 `(n, max_atoms, 3)`（实测 ValueError）。必须**展平-回绕-还原**：

```python
# _charge_centers PBC 化（核心）
def _charge_centers(positions, group_idx, group_q, group_valid, box=None):
    coords = positions[group_idx]                    # (n, max_atoms, 3)
    if box is not None:
        n, ma = coords.shape[:2]
        ref = coords[:, 0:1, :]                      # 每基团第一个原子作参考
        rel = coords - ref                           # 相对位矢 (n, max_atoms, 3)
        rel_mic = minimize_vectors(rel.reshape(-1, 3), box).reshape(n, ma, 3)
        coords = ref + rel_mic                       # 回绕后坐标（保基团连续）
    weighted = coords * group_q[:, :, None] * group_valid[:, :, None]
    q_sum = np.sum(group_q * group_valid, axis=1)
    return np.sum(weighted, axis=1) / q_sum[:, None]
```

**对抗性审查实测确认**：
- 跨盒残基（原子 x=1.0 与 9.9、L=10）：回绕后中心 x=0.4（物理正确），
  不回绕欧氏中心 x=2.77（错误）✓
- 正常残基（间距 <半盒）：回绕与欧氏中心**完全一致**（MIC 不触发）✓
- `minimize_vectors` 展平 vs 逐组手算 maxdiff 0.0（完全等价）✓
- 负坐标中心（回绕可能产出负坐标）：`calc_bonds` MIC 正确处理
  （-0.133 vs 9.8 → 0.067Å）✓
- **参考原子选择无关性**（关键）：同一跨盒残基对，参考原子取 0（0.2）或 1（9.9）
  回绕 → 两中心 MIC 距离**完全一致**（4.3667Å 两种参考均同）——判定不依赖参考选择
  （`calc_bonds` 距离对中心像的位置不变，像无关性）✓

### 6D.4 two_pass 替换

- `run_pass2` 自写循环（同前模式）：逐帧 `box=get_box(ts)` 传入
  `compute_pair_metrics(..., box)`
- `compute_pair_metrics(group_tuples, all_positions, box=None)`：
  `_charge_centers(..., box)` 两处 + `calc_bonds(pos_c, neg_c, box)` 距离

### 6D.5 per_frame 替换

- **第一帧预筛**（现 `all_dist < PREFILTER_CUTOFF`）：
  ```python
  all_dist = calc_bonds(pos_centers[pos_flat], neg_centers[neg_flat],
                        box=get_box(trajectory[0]))
  mask = all_dist < PREFILTER_CUTOFF
  ```
- **逐帧**：`_charge_centers(positions, ..., box)` + `calc_bonds` 距离

### 6D.6 测试策略

1. **PBC 单测（test_pbc.py 追加）**：
   - 跨盒残基回绕：残基原子跨盒（参考 1.0、另一原子 9.9、L=10）→ 回绕后中心正确；
     不回绕中心错误
   - 跨盒盐桥距离：正/负电荷中心跨盒（MIC ≤5.5）→ 判定存在；欧氏 >5.5 旧漏检
   - 正常残基回绕无影响：原子间距 <半盒 → 回绕中心 = 欧氏中心
   - 参考原子无关性：不同参考原子回绕 → 中心间距离一致
   - `minimize_vectors` 展平：`(n,ma,3)` 展平-回绕-还原与逐组手算一致
2. **既有测试**：`test_saltbridge_per_frame.py` / `test_saltbridge_two_pass.py`
3. **真实数据对比**：KRAS-RBD 新旧逐 pair 逐帧 existence 对比应一致
   （盐桥中心距离 ≤5.5Å ≪ 半盒 52.9Å，MIC 不触发，判定不变）

## 7. box 传递方案（对抗性审查修正：基类零改动）

> **2026-10-14 对抗性审查修正（P0）**：原方案"基类 `compute_pair_metrics` 加
> 可选 box 参数、基类 `run_pass2` 传 box"是**错误的**——8 个 two_pass 检测器全部
> 覆盖 `compute_pair_metrics(self, group_tuples, all_positions)`（无 box 参数），
> 基类调用处传 `box=` 会对未修改子类抛 **TypeError**（实测确认）。Python 动态分派
> 要求基类调用参数被子类签名兼容，改基类调用 = 改全部子类。

**修正方案：基类 `core/interfaces.py` 零改动；由 PBC 化的子类内部缓存当前帧 box。**

### 7.1 two_pass 子类模式（疏水/金属）

子类覆盖 `run_pass2`，在逐帧循环前/中缓存 box，再调基类（基类签名不变）：

```python
# 子类（疏水/金属 two_pass）
def run_pass2(self, sparse, trajectory):
    if not sparse.data:
        return []
    self._build_indices_from_sparse(sparse)
    # 逐帧预取全部 box（NPT 系综每帧不同）→ 基类 run_pass2 循环里按帧取
    self._frame_boxes = [(_get_box(ts) if hasattr(ts, 'dimensions') else None)
                         for ts in trajectory]
    return super().run_pass2(sparse, trajectory)
```

基类 `run_pass2` **保持不动**：`frame_metrics = self.compute_pair_metrics(group_tuples, ts.positions)`。

子类 `compute_pair_metrics` 内用哪一帧的 box？——**关键设计**：基类调
`compute_pair_metrics(group_tuples, all_positions)` 时**不暴露帧号**，因此
`compute_pair_metrics` 无法知道当前帧。两个可靠方案：

| 方案 | 做法 | 评价 |
|---|---|---|
| **A. 逐帧循环放子类（推荐）**：不依赖基类 `run_pass2`，子类自写 Pass2 循环（`for f, ts in enumerate(trajectory): metrics = self.compute_pair_metrics(gt, ts.positions, _get_box(ts))`），完全掌控帧与 box | 需重写基类 run_pass2 逻辑（约 30 行），但**显式、无隐式状态** | ✅ 符合"显式优于隐式"；代码量小 |
| B. 子类缓存 `self._current_box`，但基类不传帧号 → 无法同步帧 | ❌ 不可行（基类不暴露帧） |

**选用方案 A**：疏水/金属 two_pass 的 `run_pass2` 自写循环（复制基类模板约 30 行），
逐帧 `box=_get_box(ts)` 传入 `compute_pair_metrics(..., box)`；synergy：不破坏基类、
不引入隐式状态。

```python
# 子类 run_pass2（方案 A 骨架）
def run_pass2(self, sparse, trajectory):
    if not sparse.data:
        return []
    self._build_indices_from_sparse(sparse)
    group_tuples = [entry["groups"] for entry in sparse.data.values()]
    n_groups, n_frames = len(group_tuples), trajectory.n_frames
    existence = np.zeros((n_groups, n_frames), dtype=bool)
    metrics = {name: np.full((n_groups, n_frames), np.nan)
               for name in self.metric_names}
    times = np.empty(n_frames)
    for f, ts in enumerate(trajectory):
        times[f] = ts.time
        fm = self.compute_pair_metrics(group_tuples, ts.positions,
                                       box=_get_box(ts))
        existence[:, f] = self.apply_threshold(fm)
        for name in self.metric_names:
            metrics[name][:, f] = fm[name]
    ...  # 与基类 run_pass2 相同的 results 构建
```

> 替代：把"逐帧循环"提为基类新 protected 方法（`_run_pass2_loop(sparse, trajectory,
> get_box_fn)`），但改动基类风险同上；本版采用子类自写，保持基类完全不动。

### 7.2 per_frame 子类（自写 detect，天然可控）

per_frame 疏水/金属自写 `detect()`，逐帧循环里直接 `box=_get_box(ts)` 传入
`calc_bonds(..., box=box)`——**无需任何基类改动**。

### 7.3 受影响面（修正后）

- `core/interfaces.py`：**零改动** ✓
- 其他 6 个 two_pass（氢键/水桥/卤键/盐桥/π-堆积/π-阳离子）：**零改动** ✓
- 其他 per_frame（氢键/水桥等）：**零改动** ✓

## 8. 实现落点

```
新增 DuIvyInteractions/core/pbc.py（或 utils/）：
  _get_box(ts) — 统一 box 获取与防护
改动 DuIvyInteractions/core/interfaces.py：
  compute_pair_metrics 加可选 box=None；run_pass2 传 box
改动 DuIvyInteractions/interaction_detectors/：
  hydrophobic_detector_two_pass.py    — _kdtree_discover(box) + compute_pair_metrics(box)
  hydrophobic_detector_per_frame.py   — 预筛 self_capped_distance + 逐帧 calc_bonds
  metal_coordination_detector_two_pass.py — _kdtree_discover(box) + compute_pair_metrics(box)
  metal_coordination_detector_per_frame.py — 预筛 calc_bonds + 逐帧 calc_bonds
  hydrogen_bond_detector_two_pass.py  — _kdtree_discover(box) + compute_pair_metrics(box)
  hydrogen_bond_detector_per_frame.py — 预筛 capped_distance + 逐帧 calc_bonds/calc_angles
新增 Tests/unittests/test_pbc.py — MIC 正确性单测（跨盒距离/角度/下界/无盒退化）
新增 doc/pbc_design.md — 本文档
```

检测器阈值（`apply_threshold`）与判定逻辑**零改动**——只换距离计算来源。

> **实现状态（2026-10-14）**：
> - ✅ 疏水 two_pass + per_frame 已实现（commit 12832b8）
> - ✅ 金属配位 two_pass + per_frame 已实现（commit 12832b8）
> - ✅ 氢键 two_pass + per_frame 已实现（commit d336c73，方案 §6A）
> - ✅ 水桥 two_pass + per_frame 已实现（commit 4dd6c70，方案 §6B）
> - ✅ 卤键 two_pass + per_frame 已实现（commit 8979827，方案 §6C）
> - ✅ 盐桥 two_pass + per_frame 方案定稿（§6D，本次实现）
> - ⏳ π-堆积 / π-阳离子：待后续分支/章节

## 9. 测试策略

### 9.1 MIC 正确性单测（test_pbc.py 新增）

1. **跨盒近邻**：构造疏水两原子 x=0.1 与 x=L-0.1（L=10）→ PBC 距离 ≈ 0.2Å
   判为疏水（存在性 True）；欧氏 ≈ 9.8Å 判无 → 断言 PBC 化后存在
2. **下界**：构造 dist < 0.5 的原子对 → 断言 `> 0.5` 仍排除（下界不破）
3. **无盒退化**：`_get_box` 返回 None → 距离欧氏、行为与旧一致、不报错
4. **金属跨盒**：金属与配位原子跨盒（MIC < 3.0）→ 存在性 True

### 9.2 真实数据回归

- 既有 `test_hydrophobic_per_frame.py` / `test_hydrophobic_two_pass.py` /
  `test_metal_coordination_per_frame.py` / `test_metal_coordination_two_pass.py`
  跑通
- **PBC 化后跨盒新增对会改变真实数据基线**（KRAS-RBD 疏水/金属的新增对需实测），
  断言更新为**重测后的 PBC 基线**，并记录新旧差异供人工核验
  （呼应 doc/TODO.md #3 基团识别人工核验）

### 9.3 性能

- 疏水/金属规模小（KRAS-RBD hydrophobic 734、metal 111/binding 937），
  实测（734 点 × 101 帧）：`self_capped_distance(pkdtree)` 0.037s vs
  scipy cKDTree 0.057s——**本规模下反而更快**（小规模 pkdtree 有优势）；
  实现后实测对比（PBC vs 旧），记录耗时。

## 10. 兼容性与已知边界

| 项 | 结论 |
|---|---|
| 依赖 | `MDAnalysis>=2.0`（pyproject）；`self_capped_distance(method='pkdtree')` 从 1.0.1 起存在，实现时在最低支持版本验证 |
| API 兼容 | **基类零改动**（§7 修正方案：子类 run_pass2 自写 / per_frame 自写 detect）→ 所有未 PBC 化检测器行为零变化 ✓；CLI/Pipeline 无签名变化 ✓；PBC 化的两检测器 `compute_pair_metrics(..., box=None)` 加可选参数 |
| 结果兼容 | **不兼容（预期内）**：PBC 化改变跨盒结果——这是修复目的；升级路径在 limitations 注明（旧结果与 v0.0.3 不可直接对比） |
| 性能 | 小规模（疏水 734 点）实测 pkdtree 快于 scipy；5000 点带 box 慢 2.3× 的基准仅适用于大体系（水桥等，本版范围外） |
| 精度 | `calc_bonds`/`capped_distance` 内部转 float32：`4.0000001` → `4.0`（实测）。对 `dist < 4.0` 开区间判定**边界恰好安全**（4.0 本身判 False）；Pass1/Pass2 同源 float32 保持一致；补精度锁定测试 |
| 三斜盒 | MDAnalysis 官方支持正交/三斜；现有测试数据全正交（90°），三斜仅依赖官方保证，本版不补三斜用例 |
| 无盒体系 | `_get_box` 防护：None/全 0 → box=None 退化欧氏 ✓ |
| NPT 系综 | 逐帧取 `ts.dimensions`（two_pass Pass1/Pass2、per_frame 逐帧均用当前帧 box）→ 天然支持 |

## 11. 不在本版范围（后续分支/章节）

- **氢键/水桥**：含角度（`calc_angles(box)` 可 PBC 化），水桥还有三元组分组逻辑
- **卤键**：`acc_angle` 为"多邻居 argmin 择优"结构，需展开 padding 批量算角再 argmin
- **盐桥**：电荷中心为聚合坐标，需先 `apply_PBC` 回绕残基原子群再聚合
- **π-堆积/π-阳离子**：环心/法向量/投影偏移为聚合与派生几何量，三步链 PBC 化
- **per_tuple**：全部不动（将被舍弃）

---

*文档结束*
