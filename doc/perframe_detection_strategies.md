# PerFrame 检测策略设计文档

> 创建日期：2026-08-25
> 状态：草案
> 用途：InteractionDetectorPerFrame 的三种实现策略对比与设计

---

## 1. 背景

### 1.1 核心问题

逐帧检测相互作用时，需要确定"哪些基团组（pair/triple）形成了相互作用"，并计算全部帧的指标。

关键矛盾：
- **候选基团组数量**可能巨大（氢键 52 万全排列）
- **构象变化**可能导致某些 pair 只在部分帧出现
- **数据完整性**要求每个 pair 有全帧的 metric 数据（绘图、分析）

### 1.2 三种策略

| 策略 | 候选对来源 | 全帧数据 | 遗漏风险 |
|:-----|:-----------|:---------|:---------|
| 策略一（当前） | 第一帧 prefilter 一次 | ✅ 预分配 | ⚠️ 有 |
| 策略二 | 第一帧 prefilter 一次 | ✅ 预分配 | ⚠️ 有 |
| 策略三 | 每帧动态发现 | ✅ Pass2 补全 | ✅ 无 |

---

## 2. 策略一：PerTuple（已实现，参考）

### 2.1 流程

```
for each tuple:
    for each frame:
        取坐标 → 算指标
```

### 2.2 特点

- 每个 tuple 独立遍历轨迹
- 轨迹被读 n_tuples 次
- 适用场景：tuple 少

### 2.3 问题

- 轨迹 I/O 开销大
- tuple 多时不可接受（水桥 ~65h）

---

## 3. 策略二：PerFrame 单轮预筛选（当前实现）

### 3.1 流程

```
第一帧: prefilter（KDTree / 全排列距离筛选）→ 确定候选 pair 集合
预分配: existence (n_pairs, n_frames), metrics (n_pairs, n_frames)
逐帧:   对候选 pair 算指标 → 填入数组
末尾:   过滤从未存在的 pair → Interaction
```

### 3.2 数据结构

```
detect() 入口:
  pair_set = prefilter(frame0)        # 确定候选对
  existence = zeros((n_pairs, n_frames))  # 预分配
  metrics = {k: zeros((n_pairs, n_frames)) for k in metric_names}

每帧:
  existence[:, f] = ...
  metrics[k][:, f] = ...
```

### 3.3 特点

- 轨迹只读 1 次
- 预分配矩阵，零动态开销
- 候选对在第一帧确定后固定

### 3.4 问题

- **遗漏风险**：prefilter 基于第一帧的空间距离，构象变化大的体系可能遗漏后续帧出现的 pair
- PREFILTER_CUTOFF = 3× 阈值，对蛋白热涨落已很保守，但无法保证 100% 覆盖

### 3.5 性能

| 检测器 | 耗时 | 加速比（vs PerTuple） |
|:-------|-----:|----------------------:|
| 盐桥 | 0.19s | 218× |
| 氢键 | 0.23s | >260× |
| π-堆积 | 0.19s | 130× |
| 水桥 | 4.9s | ~48,000× |

---

## 4. 策略三：PerFrame 两轮遍历 + 稀疏存储

### 4.1 设计目标

- **零遗漏**：不依赖空间预筛选，每帧独立检测
- **全帧数据完整**：每个 pair 有全部帧的 metric
- **按需计算**：Pass2 可选，统计分析只需 Pass1
- **内存低**：Pass1 稀疏存储（只存 active 的帧），长轨迹场景下比策略二的预分配稠密矩阵节省 99%+ 内存（典型：1,000 pairs × 500,000 帧，策略二需 4GB，策略三 Pass1 仅需 ~4MB）

### 4.2 流程

```
Pass 1（发现 + 记录）:
  results = {}
  for each frame:
      级联筛选（距离 → 角度 → ...）→ 发现 existence=True 的 pair
      对发现的 pair 记录 metrics
      不存在的帧不记录

  输出: results 字典（稀疏存储）

[可选] Pass 2（补全）:
  用户选择是否执行
  基于 Pass1 的 pair 集合，遍历全帧算 metrics
  填充 NaN → 完整 (n_pairs, n_frames) 数组

输出: List[Interaction]
```

### 4.3 Pass1 数据结构

#### 逐帧处理

```python
results = {}  # {pair_key: data}

for f, ts in enumerate(trajectory):
    frame_pairs = detect_frame(ts.positions)

    for pair_key, metrics in frame_pairs:
        if pair_key not in results:
            results[pair_key] = {
                "groups": (g1, g2),
                "frames": [],
                "metrics": {k: [] for k in metric_names},
            }
        results[pair_key]["frames"].append(f)
        for k, v in metrics.items():
            results[pair_key]["metrics"][k].append(v)
```

#### Pass1 结束后的数据

```python
results = {
    (101, 205): {
        "groups": (ARG73_group, ASP175_group),
        "frames": [0, 1, 2, 3, 5, 6, 7, ...],
        "metrics": {
            "distance": [3.2, 3.5, 3.1, 3.8, 3.3, 3.0, ...],
            "angle":    [150, 140, 155, 130, 148, 160, ...],
        }
    },
    (101, 308): {
        "groups": (ARG73_group, GLU200_group),
        "frames": [0, 1, 2, ...],
        "metrics": {
            "distance": [4.5, 4.2, 4.8, ...],
            "angle":    [120, 125, 118, ...],
        }
    },
}
```

**只有 existence=True 的帧被记录**。不存在的帧不占空间。

### 4.4 统计分析（直接用 Pass1 数据）

```python
for pair_key, data in results.items():
    occupancy = len(data["frames"]) / n_frames
    avg_dist = np.mean(data["metrics"]["distance"])
    # 只有存在的帧纳入统计
```

### 4.5 Pass2 补全（可选，绘图用）

```python
for pair_key, data in results.items():
    existence = np.zeros(n_frames, dtype=bool)
    existence[data["frames"]] = True

    for k in metric_names:
        full = np.full(n_frames, np.nan)
        full[data["frames"]] = data["metrics"][k]
        data["metrics"][k] = full

    data["existence"] = existence
    del data["frames"]
```

Pass2 后的数据：

```python
results = {
    (101, 205): {
        "groups": (ARG73_group, ASP175_group),
        "existence": [True, True, True, True, False, True, ...],  # (101,)
        "metrics": {
            "distance": [3.2, 3.5, 3.1, 3.8, NaN, 3.3, ...],    # (101,)
            "angle":    [150, 140, 155, 130, NaN, 148, ...],      # (101,)
        }
    },
}
```

- NaN 表示该帧不存在相互作用（matplotlib 原生支持，断点绘制）
- np.nanmean / np.nanstd 直接过滤 NaN

### 4.6 Pass1 级联筛选示例（氢键）

```
全排列 52 万对
  → 算 D-A 距离（向量化）→ 筛 distance < 4.1Å → 剩 ~638 对（99.9% 淘汰）
  → 算 D-H-A 角度（仅对存活对）→ 筛 angle ≥ 100° → 剩 ~119 对
  → 记录到 results
```

Pass1 内部可使用 KDTree 等空间索引加速级联筛选，框架不限定具体算法。

### 4.7 内存估算

#### Pass1 稀疏存储

| 场景 | 活跃帧次 | 每帧次存储 | 总计 |
|:-----|:---------|:-----------|-----:|
| 盐桥 | 2,473 | ~16 B | 39 KB |
| 氢键 | ~2,000 | ~16 B | 32 KB |
| 水桥 | 10,138 | ~32 B | 317 KB |

#### Pass2 稠密存储

| 场景 | n_pairs × n_frames | 总计 |
|:-----|:-------------------|-----:|
| 盐桥 101 帧 | 47 × 101 | 0.4 MB |
| 氢键 101 帧 | 119 × 101 | 1.0 MB |
| 氢键 1μs | 119 × 500,000 | 476 MB |

Pass2 的内存开销与策略二相同。Pass1 的内存开销远小于 Pass2。

### 4.8 计算量估算（氢键）

| | 计算内容 | 耗时 |
|:--|:---------|-----:|
| Pass1 每帧 | 全排列距离(52万) → 筛选(638) → 角度(638) | ~56ms |
| Pass1 合计 | 56ms × 101 帧 | ~5.7s |
| Pass2 每帧 | 119 对完整 metric | ~2.7ms |
| Pass2 合计 | 2.7ms × 101 帧 | ~0.3s |
| **Pass1 only** | | **~5.7s** |
| **Pass1 + Pass2** | | **~6.0s** |

对比策略二：0.23s。策略三慢 26×，但零遗漏。

---

## 5. 三种策略对比

| 维度 | 策略一（PerTuple） | 策略二（单轮预筛选） | 策略三（两轮遍历） |
|:-----|:-------------------|:--------------------|:-------------------|
| 轨迹 I/O | n_tuples 次 | 1 次 | 2 次 |
| 遗漏风险 | 有 | 有 | **无** |
| 全帧数据 | ✅ | ✅ | ✅（Pass2 后） |
| 计算效率 | 低 | **最高** | 中 |
| 内存 | 低 | 高（预分配） | Pass1 低 / Pass2 高 |
| 统计分析 | 需全量数据 | 需全量数据 | **Pass1 即可** |
| 接口 | `List[Interaction]` | `List[Interaction]` | `List[Interaction]` |

---

## 6. 策略三用户交互

```python
# Pass1：总是执行
detector = HydrogenBondDetectorPerFrame(strategy="two_pass")
results = detector.detect(groups, trajectory)

# 此时 results 已可用于统计分析（稀疏存储）

# 用户选择是否补全（绘图需要）
if need_plotting:
    results = detector.fill_missing(trajectory)
    # results 变为完整 (n_pairs, n_frames) 数组
```

---

## 7. 实现优先级

| 策略 | 状态 | 优先级 |
|:-----|:-----|:-------|
| 策略一（PerTuple） | ✅ 已实现 | 保留为参考 |
| 策略二（单轮预筛选） | ✅ 已实现 | 当前默认 |
| 策略三（两轮遍历） | ❌ 未实现 | 按需实现 |

策略三是策略二的"无遗漏"升级版。当用户需要 100% 覆盖保证时使用。

---

*文档结束*
