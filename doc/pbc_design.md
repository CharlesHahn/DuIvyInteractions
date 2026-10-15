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
新增 Tests/unittests/test_pbc.py — MIC 正确性单测（跨盒距离/下界/无盒退化）
新增 doc/pbc_design.md — 本文档
```

检测器阈值（`apply_threshold`）与判定逻辑**零改动**——只换距离计算来源。

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
