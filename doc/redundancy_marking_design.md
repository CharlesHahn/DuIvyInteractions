# 疏水-芳香重叠标记设计（TODO #2）

> 状态：设计定稿（待实现，本轮仅文档）
> 日期：2026-11-06
> 分支：dev

## 1. 背景与问题

### 1.1 现状

基团识别阶段**不排除芳香环原子**（`doc/hydrophobic_definition.md` §4）：芳香环碳在所有邻居均为
C/H 时，同时被识别为 `hydrophobic` 基团。实测（KRAS-RBD D927）：

| 项 | 数值 |
|---|---|
| `hydrophobic` 基团总数 | 734 |
| 其中属于 `aromatic_ring` 原子的 | **153（20.8%）** |

因此当两个芳香环贴近时，同一物理接触（环间 van der Waals 力）会被**报告两次**：

```
PHE-A 的 CZ  ──(3.8 Å)──  PHE-B 的 CZ
   ↑ 两端均为 hydrophobic 基团   → 输出"疏水相互作用"
   ↑ 两端均在 aromatic_ring 内   → 同时输出"π-π 堆积"
```

### 1.2 为什么不能靠改基团识别解决

若在识别阶段排除芳香环原子，则：
- 芳香环与脂肪族碳之间的疏水接触（占多数，实测 1144 对）会丢失
- 芳香环-芳香环之间**未形成 π-π 堆积**（距离/角度不符）但仍为疏水接触的情形也会丢失

→ 疏水基团必须保留环原子；重复计数只能在**检测之后**处理。

## 2. 设计原则

| # | 原则 | 说明 |
|---|---|---|
| **P1** | **不干涉两个检测器** | `pi_stacking` 与 `hydrophobic` 检测器代码零改动，各自只依赖自己的基团类型 |
| **P2** | **标记而非删除** | 重叠信息以**标记**（flag）形式附加，原始记录完整保留，可追溯、可复核 |
| **P3** | **默认不标记** | 仅当用户显式请求时执行（`--mark-redundant`），保证默认行为与既有结果完全一致 |
| **P4** | **逐帧标记** | 与 `existence`/`metrics` 同构的 `(n_pairs, n_frames)` bool 矩阵；对级汇总由导出层派生 |
| **P5** | **不改检测器契约** | 后处理在编排层（pipeline）完成；`Interaction` 仅新增**可选**字段（默认空） |

## 3. 判定语义

### 3.1 原子级匹配（不做残基对级）

对某一帧 f，疏水对 `(h1, h2)` 被标记 ⇔ **同帧**存在某个 π-π 堆积对 `(R1, R2)`，使得

```
(h1 ∈ R1 且 h2 ∈ R2)  或  (h1 ∈ R2 且 h2 ∈ R1)
```

其中 `R1`/`R2` 为两个芳香环的原子集合（`aromatic_ring` 基团的 `atom_indices`），
`h1`/`h2` 为疏水基团的代表原子。

**为什么不做残基对级**：残基对级（"若残基 A、B 的环已形成 π-π，则移除 A-B 间全部疏水对"）
会误标环外碳（如 PHE 的 CB、TRP 的 CB）与对方环的接触。实测差异：

| 口径 | 被标记的疏水对 | 占比（/2366） |
|---|---|---|
| 残基对级 | 107 | 4.5% |
| **原子级（本设计）** | **67** | **2.8%** |

原子级更精确，且天然处理 TRP 双环、HIS 咪唑环等多环残基（无需特判）。

### 3.2 逐帧而非对级

| 维度 | 逐帧矩阵 | 对级单一 bool |
|---|---|---|
| 精度 | 精确到帧 | 任一帧重叠即整对标记 |
| 与现有结构一致性 | 与 `existence` 相同 | 异质结构 |
| 可统计性 | 可算"冗余帧占比" | 二值 |
| 下游过滤 | 可按帧剔除 | 只能整对剔除（过度） |

**实测说明**（D927，two_pass）：至少一帧被标记的疏水对为 67，但逐帧合计仅 749 个
「对-帧」（占疏水对-帧 0.6%）。若按对级标记，这 67 对会被**整对**视为冗余，
而其中大部分帧并不冗余。

**结论**：存储逐帧矩阵；对级信息由导出层派生为"冗余帧占比"。

## 4. 数据结构

### 4.1 `Interaction` 新增可选字段

```python
# DuIvyInteractions/core/datas.py
@dataclass
class Interaction:
    interaction_type: str
    groups: List[Tuple[Group, ...]]
    existence: np.ndarray
    metrics: Dict[str, np.ndarray]
    times: np.ndarray
    flags: Dict[str, np.ndarray] = field(default_factory=dict)   # 新增，默认 {}
```

- 键名约定：`"redundant_pi_stacking"`；值：`(n_pairs, n_frames)` bool
- `__post_init__` 追加校验：每个 flag 的形状必须等于 `existence.shape`
- **39 处 `Interaction(...)` 构造点全部使用关键字参数** → 带默认值的字段不破坏任何现有调用

### 4.2 为什么不放进 `metrics`

`io/h5.py` 读取时**遍历 `metrics` 组下全部键**，`io/interaction_exporter.py` 的
`summary.csv` 对**每个 metric** 生成 `avg_<name>`/`std_<name>` 并输出独立 xvg。
若把标记放进 `metrics`：
- 会生成 `avg_redundant` 等聚合列与 xvg，污染几何指标语义
- 后续若增加更多标记类型，`metrics` 将被分类标签塞满

→ 标记必须走独立通道（`flags`）。

### 4.3 为什么不用 side-car 文件 / Group.metadata

| 备选 | 否决理由 |
|---|---|
| 独立 side-car h5（如 `hydrophobic_redundant.h5`） | 数据分裂为两文件，导出与消费需 join，用户易忽略 |
| `Group.metadata` | 粒度错误——冗余是「基团对 × 帧」信息，非基团属性 |

## 5. CLI 与执行流程

### 5.1 命令行

```bash
# 默认：不标记（行为与现状完全一致）
dii run -t md.tpr -f md.xtc -o out/ --ff amber --interactions hydrophobic

# 显式请求标记：若未请求 pi_stacking，则自动补跑（结果保留）
dii run -t md.tpr -f md.xtc -o out/ --ff amber \
        --interactions hydrophobic --mark-redundant
```

| 参数 | 类型 | 默认 | 语义 |
|---|---|---|---|
| `--mark-redundant` | flag | `False` | 标记与 π-π 堆积重叠的疏水对-帧 |

### 5.2 参数行为矩阵

| 请求的 `--interactions` | `--mark-redundant` | 行为 |
|---|---|---|
| 含 `hydrophobic` | `False` | 仅跑请求的类型，**不标记**（现状） |
| 含 `hydrophobic` | `True` 且**不含** `pi_stacking` | 自动补跑 `pi_stacking`（结果保留于输出目录），完成标记 |
| 含 `hydrophobic` | `True` 且含 `pi_stacking` | 直接使用已有 π-π 结果完成标记 |
| 不含 `hydrophobic` | `True` | 无疏水可标 → 提示 `[WARN] --mark-redundant requires hydrophobic` 并忽略 |
| — | — | `pi_stacking` 单独请求时不涉及标记 |

### 5.3 自动补跑的一致性要求（关键）

补跑 `pi_stacking` 必须使用与主流程**完全相同**的：
- `--ff`（力场识别器）
- `--strategy`（策略：`two_pass` / `per_frame`）
- `--select` / `--pairing`（参与池与配对规则 → `assignment` / `tuple_filter`）

否则会标记出"用户限定范围之外"的冗余对，或与用户实际看到的 π-π 结果不一致。

### 5.4 pipeline 流程

```
1. 组装待运行类型列表 names（现状逻辑）
2. 若 mark_redundant 且 "hydrophobic" ∈ names 且 "pi_stacking" ∉ names：
       追加 "pi_stacking"（标记为自动补跑，输出目录保留其 h5）
3. 逐类型检测并保存（现状逻辑不变）
4. 若 mark_redundant 且 hydrophobic 已产出：
       读取 hydrophobic.h5 与 pi_stacking.h5 → 计算逐帧标记
       → 写回 hydrophobic.h5（保留其余字段不变）
5. 控制台提示：
       "[INFO] marked N hydrophobic pairs (M pair-frames) overlapping pi-stacking (kept, not removed)"
```

**实现取舍**：采用"先存后读回标记"（步骤 4）而非"内存中标记后保存"——
后者需重构 pipeline 的"检测一个存一个"流程，前者改动更小、风险更低（多一次 I/O）。

## 6. 导出与呈现

| 输出 | 内容 |
|---|---|
| `hydrophobic_summary.csv` | 新增列 `avg_redundant_pi_stacking`：该对**活跃帧**中冗余帧的比例（0–1） |
| `hydrophobic_redundant_pi_stacking.xvg`（可选） | 逐帧 0/1 时序，便于观察冗余随时间的分布 |
| 控制台概览 | 标记的对数与「对-帧」数，并明确"未删除" |
| 文档 | `limitations` 由"疏水-芳香未去重"更新为"以标记形式提供，用户可自行过滤" |

**导出器改动**：在 `InteractionExporter` 基类增加一个可选的 flags 处理钩子
（默认不输出），仅疏水导出器覆写。其余 7 类导出器行为不变。

## 7. h5 格式兼容性

### 7.1 方案：新增可选组 `flags/`，格式版本保持 `1.0`

```
interaction_0/
├── existence      (n_pairs, n_frames) bool
├── times          (n_frames,)
├── metrics/       ...
├── groups/        ...
└── flags/         ← 新增（可选）
    └── redundant_pi_stacking  (n_pairs, n_frames) bool
```

- **旧代码读新文件**：`_read_interaction` 只读取 `existence`/`times`/`metrics`/`groups`，
  **忽略未知的 `flags/`** → 不报错（标记不可见）
- **新代码读旧文件**：`flags/` 不存在 → `flags = {}`

### 7.2 已知代价（如实记录）

- 旧版本代码读取并**重存**新文件时，`flags/` 会丢失
- 因此**不升版本号**（`1.0` 保持）意味着"新增可选组"属**非破坏性扩展**；
  若未来出现其他破坏性变更，再统一升版

> 备选：升版 `1.1` 且 `load_interactions` 同时接受 `{"1.0","1.1"}`。更规范，但需额外兼容测试；
> 本设计按"非破坏性扩展"处理，保留 1.0。

## 8. 影响面

| 项 | 影响 |
|---|---|
| 8 个检测器 | **零改动**（`flags` 有默认值） |
| `hydrophobic` 结果 | 仅在 `--mark-redundant` 时新增 `flags`，`existence`/`metrics` 完全不变 |
| `pi_stacking` 结果 | 不受影响（仅被读取） |
| 默认行为 | 完全不变（默认不标记） |
| h5 格式 | 新增可选组，前后兼容 |
| 导出 | 疏水 summary 新增 1 列；其余不变 |
| 性能 | 后处理成本 O(活跃对-帧)（预建 原子→环 索引，避免嵌套循环） |

### 实测数据（D927，101 帧，two_pass）

| 指标 | 数值 |
|---|---|
| 疏水对 | 2366 |
| π-π 堆积对 | 10 |
| 原子级标记：至少一帧被标记的疏水对 | **67（2.8%）** |
| 原子级标记：标记的「对-帧」 | **749（占疏水对-帧 0.6%）** |

## 9. 备选方案对比（均否决）

| 方案 | 做法 | 否决理由 |
|---|---|---|
| A. 删除式去重 | 直接从疏水结果删除重叠对 | 丢失原始信息，不可追溯；与 P2 冲突 |
| B. 检测器内部集成 | 疏水检测器自行加载 aromatic_ring 并剔除 | 违反单一职责，改变检测器契约（只依赖自身基团类型） |
| C. 残基对级标记 | 按残基对匹配 | 过度标记（4.5% vs 2.8%），误标环外碳 |
| D. 标记入 `metrics` | `metrics["redundant"]` | 污染几何指标语义与导出 |
| E. 识别阶段排除环原子 | 疏水基团排除芳香碳 | 丢失环-脂肪碳、未成堆积的环-环疏水接触 |

## 10. 验证计划

| 项 | 方法 |
|---|---|
| 标记正确性（合成） | 构造两芳香环 + 其环碳疏水对，断言对应「对-帧」被标记，其余不标记 |
| 反向验证（非重叠） | 构造环-脂肪碳疏水对，断言**不**被标记 |
| 逐帧语义 | 构造仅部分帧形成 π-π 的体系，断言标记矩阵随帧变化 |
| 默认零影响 | `--mark-redundant` 未给出时，hydrophobic 的 `existence`/`metrics` 与改动前逐位一致 |
| 自动补跑 | `--interactions hydrophobic --mark-redundant` 时输出目录中存在 `pi_stacking.h5` |
| 一致性 | 自动补跑的 π-π 与显式请求 `pi_stacking` 的结果一致（同 ff/strategy/select/pairing） |
| h5 兼容 | 旧文件（无 `flags/`）可读；新文件（有 `flags/`）round-trip 后标记一致 |
| 真实数据 | D927 / CHARMM 全原子体系跑通，记录标记对数与对-帧数 |

## 11. 待确认项（实现前锁定）

| # | 事项 | 当前建议 |
|---|---|---|
| 1 | CLI 参数名 | `--mark-redundant` |
| 2 | flag 键名 | `redundant_pi_stacking` |
| 3 | 是否导出逐帧 xvg | 可选（默认导出，便于核查） |
| 4 | 自动补跑时打印级别 | `[INFO]`（非 WARN，因属用户显式请求的派生行为） |

## 12. 关联文档

- `doc/hydrophobic_definition.md` §4（疏水基团不排除环原子 + 去重规则原型）
- `doc/criteria.md`（检测判据）、`docs/reference/limitations.md`（已知限制）
- `DuIvyInteractions/core/datas.py`（`Interaction` 定义）
- `DuIvyInteractions/io/h5.py`（序列化格式）、`io/interaction_exporter.py`（导出契约）
- `DuIvyInteractions/pipeline.py`（编排层后处理位置）
- `doc/TODO.md` #2（本设计的来源）
