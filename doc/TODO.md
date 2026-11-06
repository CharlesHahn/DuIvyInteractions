# TODO 清单

> 记录推迟到后续阶段实现的任务，避免遗忘。
>
> 状态标记：✅ 已完成 | ⏳ 进行中 | ❌ 未开始

---

## 相互作用检测阶段需要实现的任务

以下任务在基团识别阶段有意推迟，需要在实现 `interaction_detectors/` 时一并完成。

### 1. 芳香环平面性检验 ✅

**推迟原因**：初始帧构象可能未弛豫，平面性可能很差；热涨落导致瞬时平面性波动；力场原子类型已编码芳香性，不需要几何验证。

**实现位置**：`interaction_detectors/pi_stacking_detector_*.py`

**判据**：
- 对环内每个原子，计算相邻两键的法向量
- 所有法向量两两夹角 ≤ 5.0°（PLIP `AROMATIC_PLANARITY`）
- 不满足 → 该帧此环不参与 π-π 堆积

**参考**：`doc/aromatic_ring_definition.md` §3.3

**实现状态**（2026-08-28 更新）：已在所有 pi_stacking 检测器中通过 `check_planarity: bool = False` 参数实现，默认关闭。

### 2. 疏水-芳香去重 ❌

**推迟原因**：π-π 堆积已包含疏水接触（芳香环之间的 van der Waals 力），如果同时报告 π-π 堆积和疏水相互作用，会重复计数同一物理接触。

**实现位置**：`interaction_detectors/` 的主流程编排

**规则**：
- 先检测 π-π 堆积
- 再检测疏水相互作用
- 移除与 π-π 堆积重叠的疏水接触（两个芳香环的原子之间的疏水接触）

**参考**：`doc/hydrophobic_definition.md` §4


### 3. 基团识别结果需要人为审查！ ✅（2026-11-03 标记完成）

目前基于真实数据识别得到的基团的数据，还需要人去手工一个一个核对一下确认一下。单元测试的真实数据测试是基于现有代码输出结果做的，不一定对。以人工审查过的结果为准！

**完成状态（2026-11-03）**：以"自动核查 + 文献核验 + 对抗性审查"替代人工逐个核对完成：
- 自动 diff 工具 `Tests/tools/verify_identification.py`（独立化学基准，4 力场）
- 三力场真实体系 + OPLS 合成体系核查：缺失=0（`doc/group_identification_review.md` §4）
- 8 个子 agent 文献核验"多余"项：PRO/MET/CYS/HIS 受体全部确认正确
- 对抗性审查：F1 卤素、F2 嘌呤环均确认"忠实力场语义非缺陷"
- 实际修复：3c GROMOS 中性胺受体判据（`d1dc7fa`）

**优先级**：已关闭（2026-11-03）

### 4. 可选的候选对预过滤优化 ✅

**检测算法**：遍历基团对，每对向量化计算全部帧的距离。预加载基团原子坐标到内存（F×G×3），每对直接得到一行 existence 和 distance 向量。

**问题**：候选对数量可能很大（如氢键 ~74000×38000），虽然每对计算很快，但总对数多时仍耗时。

**优化方案**：提供可选的预过滤参数（如 `--prefilter`），用户在确认体系构象变化不大的情况下，可在第一帧用粗略距离 cutoff 去除明显不可能的候选对。

**注意**：
- 此优化**默认关闭**，因为构象变化大的体系会遗漏真实相互作用（科学错误）
- 是否启用由用户通过命令行参数自行判定
- cutoff 应远大于相互作用距离阈值（如氢键 4.1 Å → 预过滤用 15 Å）

**实现状态**（2026-08-28 更新）：已在 PerTuple 和 PerFrame 检测器中通过 `PREFILTER_CUTOFF` 常量和 `filter_candidate_tuples()` 方法实现。PerTuple 用第一帧距离预过滤，PerFrame 用 KDTree 预筛选。

### 5. 基团组用户自定义过滤（tuple_filter） ✅

**需求**：用户需要指定相互作用的基团来源，如"只检测 RBD 和 KRAS 之间的盐桥"。

**实现**：`detect` 方法的 `tuple_filter` 参数，接受 `(Tuple[Group,...]) -> bool` 函数。

**执行顺序**：
```
get_candidate_tuples → tuple_filter → filter_candidate_tuples(距离) → 遍历计算
```

**使用示例**：
```python
# 只检测不同蛋白之间的盐桥
detector.detect(groups, trajectory,
    tuple_filter=lambda gt: gt[0].molecule != gt[1].molecule)

# 只检测 RBD 和 KRAS 之间的盐桥
detector.detect(groups, trajectory,
    tuple_filter=lambda gt: {gt[0].molecule, gt[1].molecule} == {"RBD_pro", "KRAS_pro"})
```

**状态**：已实现。

### 6. 官能团分类粗糙问题 ✅（2026-11-03 对抗性审查：设计如此，不修改）

**原问题**：PLIP 的 `is_functional_group` 对氮基团分类粗糙。"tertamine" 实际包含所有
sp3 氮（NR₃、NR₃H⁺、NH₃⁺ 等），不是化学意义上的"叔胺"。

**对抗性审查（2026-11-03）→ 关闭，不修改代码**：
1. **标签粗糙但功能正确**：N 末端 NH₃⁺ 被标 `func_group='tertamine'`（化学是伯胺），但
   判成 `charged_positive` 是**正确且必需**的（NH₃⁺ 本就带正电，参与盐桥/π-cation）。
2. **修改方案有破坏性副作用**：拟议的"重原子（非 H）≥3"判据会丢失蛋白 N 末端
   NH₃⁺ 的正电基团——实测 KRAS-RBD 42→40（ASN0 RBD、MET143 KRAS 两个 N 端丢失），
   因为 N 末端残基不在 `POSITIVE_RESIDUES` 字典，唯一依赖 tertamine 路径判正电。
   改判据 → 盐桥/π-cation 漏检。
3. **标签无下游消费**：`func_group` 在检测器/pipeline/selection 层**零消费者**
   （grep 实测）；pi_cation 的 tertamine 角度检查（#7）尚未实现。当前标签不影响结果。
4. **忠实 PLIP**：PLIP 用 `tertamine` 粗粒度覆盖所有 sp3 正电胺，属其既定行为；
   项目"忠实复现 PLIP 定义"是设计声明（AGENTS.md），非缺陷。

**结论**：#6 是"标签命名与化学命名不一致"，**非功能缺陷**；修改判据会破坏功能而被否决。
若日后实现 #7（tertamine π-cation 角度检查）时需区分 N 末端，再单独设计（须保留
N 末端正电资格）。

**优先级**：已关闭（2026-11-03，设计如此）

### 7. PLIP 对 tertamine 的 pi-cation 特殊处理 ✅（决策：不实现，设计如此）

**PLIP 的实际操作**（源码 `plip/structure/detection.py` → `pication()`，2026-11 核对）：
计算胺的三个取代基平面法向量 `amine_normal = cross(c0→c1, c2→c0)`，与环法向量
`ring.normal` 取夹角 `a = min(b, 180-b)`；**`a ≤ 30°` 时保留该 pi-cation 对，`a > 30°` 时排除**
（源码写作 `if not a > 30.0: append(...)`）。

> ⚠️ 本条目原始描述（"若夹角 ≤ 30° 排除"）**方向写反**，已按源码更正。

**物理含义**：`amine_normal` ≈ N 的孤对轴方向。夹角小（≤30°）表示孤对轴指向环心
→ 阳离子真正骑在 π 电子云上 → 真 π-cation；夹角大（>30°）表示环横在胺伞面一侧、
被三个取代基隔开 → PLIP 要排除的"穿过配体（through the ligand）"假阳性。

**作用范围**：仅当正电基团属于**配体**且 `fgroup == 'tertamine'` 时检查
（源码 `type(p).__name__ == 'lcharge'`，调用点 `pication(rings_bs, ligand.get_pos_charged(), False)`）；
**蛋白侧正电基团（含蛋白 N 末端 NH₃⁺）从不施加此检查**。

**决策（不实现）依据**：
1. **PLIP 该实现自身有两处缺陷**：① `amine_normal` 取 `OBAtomAtomIter` 前 3 个邻居而不检查
   元素——质子化叔胺 R₃NH⁺（3C+1H）会取到 H，法向量不再是取代基伞面且依赖迭代顺序；
   ② 特例分支后 `break` 跳出 `pos_charged` 循环——同一芳香环对其余配体正电基团（胍基/季铵等）
   的评估被跳过，多正电基团配体可能漏报。照抄即继承缺陷。
2. **蛋白侧本不受约束**：本项目核心场景（蛋白-蛋白界面等）在该判据下与 PLIP 无偏差。
3. **现有数据零触发**：三套真实体系（Amber/GROMOS/CHARMM）无真叔胺，判据在任何已验证
   结果上均不生效。
4. **属性定位**：这是 PLIP 针对特定几何假阳性的启发式过滤，非 π-cation 的定义性化学判据；
   不做不违背化学正确性，但**属已声明的已知偏差**（见 `docs/reference/limitations.md`）。

**遗留已知偏差表述**：当阳离子配体含叔胺、且芳香环位于其取代基侧（而非孤对侧）时，本项目
不施加额外角度过滤 → 相对 PLIP 可能**多报**该 π-cation（PLIP 眼中的"穿过配体"假阳性我们不剔除）。

**实现位置**（若未来要做）：识别器在 `Group.metadata` 存真叔胺的 3 个重原子邻居全局索引 →
pi_cation 检测器（3 策略）按索引取坐标算法向量并与环法向量比较；需先定义"配体/蛋白"判定规则
（PLIP 用 receptor/ligand 角色，本项目仅有 molecule 名 + `--select`）。

**优先级**：已关闭（2026-11-05 决策：不实现，设计如此）

### 8. 金属配位几何构型匹配 ❌

**推迟原因**：第一版只做距离判据，几何构型匹配复杂度高。

**PLIP 的操作**：按金属分组后，计算所有配位原子的角度，与已知几何构型（linear/trigonal/tetrahedral/octahedral 等）比较，选择 RMS 最小的构型，移除不符合该构型的多余配位原子。

**实现位置**：`interaction_detectors/metal_coordination_detector_*.py`

### 9. 金属配位排除纯水配位 ✅

**推迟原因**：当前 `metal_binding` 基团已排除水分子，不需要额外处理。但如果未来改为不排除水，则需要此逻辑。

**PLIP 的操作**：如果所有配位原子都来自水分子，不报告该金属配位。

**实现状态**（2026-08-28 更新）：`_find_metal_binding()` 中已通过 `if res.residue_name in WATER_RESIDUES: return` 排除水分子。

### 10. 水桥性能优化 ✅（部分）

**问题**：当前水桥检测候选三元组数量太大（24.9 万个），每个三元组需要逐帧遍历轨迹（~0.95 秒），总耗时约 65 小时，不可接受。

**根因**：`_process_tuple` 对每个三元组独立遍历轨迹，与其他检测器（候选对少）不同，水桥的三元组数量级太大。

**待优化**：减少候选三元组数量或优化轨迹加载方式。

**实现状态**（2026-08-28 更新）：
- PerFrame 版：用 KDTree 预筛选（水到供/受体 < 8.2Å），耗时从 65h 降至 ~5s
- TwoPass 版：Pass1 用 KDTree 逐帧发现，Pass2 只对 active 三元组补全
- PerTuple 版：仍为原始实现（参考用）

### 11. 水桥水分子去重 ✅（2026-11-06 改为显式-H 判据 + 三层配额）

**原记录（错误）**："PLIP 的操作：一个水分子最多参与 2 个氢键，保留 H-O-H 角度最接近 110° 的两个。"

**对抗性审查更正（2026-11-05）**：
1. PLIP 的 `refine_water_bridges()`（`preparation.py`）是**三步**——供体-氢键互斥 + (水,受体)去重 + 每水最多 2 桥，择优用 `|110° − omega|`（`omega` = 受体–水O–供体H，顶点**水氧**，不是水自身的 H-O-H 角）
2. PLIP 只用水的氧，是因为其输入为 **PDB**（结晶水通常无 H）；本项目 tpr 为**全原子拓扑**，水有显式 H（Amber 37021 个水全含 OW/HW1/HW2）
3. 现有实现在"每水桥数"上无任何收敛（实测单帧最大 **13** 桥，26.5% 的水实例 >2 桥）

**实现（2026-11-06）**：改用**显式水氢的氢键标准**（Baker & Hubbard 1984）+ 三层配额：
- 供体侧 `D–H···Ow`：H···Ow ≤ 2.5 Å 且 D-H···Ow ≥ 120°
- 受体侧 `Ow–Hw···A`：Hw···A ≤ 2.5 Å 且 Ow-Hw···A ≥ 120°
- 水在中间：`∠(供体H–Ow–受体A)` ∈ [71°, 140°]（保留 PLIP omega 语义）
- 配额：① 每水 H ≤1 受体；② 每水 O ≤2 供体；③ 每水 ≤2 桥（逐帧）
- 设计文档：`doc/water_bridge_design.md`

**实测效果（D927）**：
| 指标 | 旧（PLIP 宽松） | 新（显式-H） |
|---|---|---|
| 单帧桥数 | 249 | 53 |
| per_frame 候选对 | 1698 | 86（候选仅第一帧，设计原则） |
| two_pass 候选 / 有存在帧 | 14598 / 14562 | 3130 / 2913 |
| 每水最大桥数 | 13 | 2（配额违反 0） |
| 判据边界违反 | — | 0（max H···O=2.50、min 角=120.0°） |
| 帧 0 两策略一致 | — | 53 = 53（全轨迹对数因候选策略不同而不同，属设计） |

CHARMM 全原子体系（SMO-BST）：per_frame 1329→109、two_pass 1678→414（有存在帧 367）。

**遗留**：`per_tuple` 策略保留旧 PLIP 判据（逐候选组架构无法实现跨组配额），该策略按 AGENTS.md 已不被维护、测试不跑。

**优先级**：已关闭（2026-11-06）

### 12. 水桥单元测试结果断言 ✅（2026-11-03 补全 TwoPass 断言）

**原跟踪**：性能问题曾导致无法跑完全量检测、无法获取真实数据结果；性能解决后
（2026-08-28 更新）可补充断言。

**补齐（2026-11-03）**：
- **PerFrame**：已有数值断言（`n_pairs == 1698`、`occupancy` 范围、`top_pair` 组成与
  `occ > 0.9`）
- **TwoPass**（本次补，注意与 PerFrame 数值不同——策略差异）：
  - `test_n_pairs_candidates`：`n_pairs == 14598`（候选三元组总数，TwoPass 保留全部含
    从未存在者）
  - `test_active_pairs`：有存在帧的对数 `== 14562`（排除 36 个从未存在的候选）
  - `test_top_pair`：`occ > 0.9` 且第二元为 water
  - `test_top_pair_identity`：top pair 与 PerFrame 一致（LEU95 H_donor → SOL water →
    PHE104 H_acceptor，占有率 `0.970`）

**验证**：TwoPass 19 passed；全量回归 685 passed / 0 failed。

**优先级**：已关闭（2026-11-03）

### 13. PerFrame 检测器长轨迹内存优化 ❌

**问题**：PerFrame 检测器预分配 (n_pairs, n_frames) 的 distance/angle 矩阵。D927 体系氢键候选对 42,603 个，101 帧时矩阵约 71 MB，可接受。但外推到 1μs 轨迹（500,000 帧），矩阵将达到 342 GB，不可接受。

**根因**：当前实现对全部候选对存储每一帧的 distance/angle，即使大部分 pair 最终被阈值淘汰（42,603 候选 → 119 有效）。

**优化方案**：只存 existence 矩阵（bool，1 byte），不存 distance/angle。最终过滤后，对存活的 ~100 个 pair 重新遍迹计算 distance/angle。

| | 当前方案 | 优化方案 |
|:--|:---------|:---------|
| 氢键 101 帧 | 71 MB | 4.1 MB (existence only) |
| 氢键 500,000 帧 | 342 GB | 20.3 GB (existence only) |
| 最终结果 | 同 | 同（重算存活 pair） |

**进一步优化**：existence 矩阵也可改为稀疏存储或增量累积，将 20.3 GB 降到更低。

**实现位置**：`interaction_detectors/*_per_frame.py` 的 `detect()` 方法。

**优先级**：低（当前 101 帧测试场景无压力，1μs 轨迹时再做）。

### 14. Interaction 结果序列化（保存/加载） ✅（已被 HDF5 替代实现）

**需求**：将 Interaction 数据结构保存到文件，支持后续分析、可视化、跨工具共享。

> **实现状态（2026-10-09）**：本条的"JSON + .npz 双文件方案"未采用——已由
> `DuIvyInteractions/io/h5.py`（HDF5 单文件，format v1.0：existence/metrics/groups/
> metadata 无损存取 + gzip 压缩）**替代实现**。本条关闭。

**方案**：JSON + .npz 双文件格式。

| 文件 | 内容 | 格式 |
|:-----|:-----|:-----|
| `<name>.json` | 元数据（interaction_type, groups 信息：group_id, group_type, molecule, residue_name, residue_id 等） | JSON（人类可读） |
| `<name>.npz` | 数组数据（existence, distance, angle 等 metrics） | numpy compressed（高效压缩） |

**接口设计**：
- `save_interaction(interaction, path)` → 写出 json + npz
- `load_interaction(path)` → 读入 json + npz → 重建 Interaction 对象

**优势**：
- JSON 人类可读，方便调试和检查
- npz 压缩后通常 10:1（bool/float 数组压缩率高）
- numpy + json 无额外依赖

**替代方案**：
- pickle：最快但跨版本脆弱，不可读
- HDF5：单文件，支持切片读取，但需 h5py 依赖

**实现位置**：`utils/output.py`（待创建）

**优先级**：中（当前无紧迫需求，Pipeline 完成后实现）

### 15. 水桥 TwoPass 未添加 WATER_BRIDGE_MINDIST 下界 ✅（已修复 2026-10-30）

**原现状**：策略三（TwoPass）的水桥检测器在 `apply_threshold` 中未添加
`dist > WATER_BRIDGE_MINDIST`（2.5Å）的下界检查，仅检查 `dist < 4.1Å` 上界。
策略一（PerTuple）和策略二（PerFrame）均有此下界。

**原影响**：TwoPass 会保留 Ow-A 距离 < 2.5Å 的三元组（如 2.48Å），而策略一/二会排除。

**修复（2026-10-30）**：TwoPass `apply_threshold` 补上下界（`dist_dw`/`dist_wa` 均 > 2.5Å），
与策略一/二对齐；新增常量 `WATER_BRIDGE_MINDIST = 2.5`。

**验证**：D927 全轨迹（101 帧）修复前存在帧中 `dist_dw < 2.5` 有 2 个、`dist_wa < 2.5` 有
162 个对-帧组合，修复后**均为 0**；真实体系水桥数与策略二一致。新增回归断言
（`test_active_dist_dw_range`/`test_active_dist_wa_range` 补下界）。


### 16. 周期性问题 ✅（已实现并合并：v0.0.4）

**原问题**：相互作用（尤其水桥/氢键）可能受周期边界影响；原建议替换到 MDAnalysis 的周期接口。

**实现状态（2026-10-14，✅ 已实现于 feature/pbc-handling）**：按设计 `doc/pbc_design.md`
完成**全部 8 类相互作用 × two_pass/per_frame** 的最小镜像化——统一改用 MDAnalysis
`capped_distance`/`calc_bonds`/`calc_angles`/`minimize_vectors`（传完整 `ts.dimensions`，
正交/三斜自动处理）；聚合型基团（环/电荷中心）先以参考原子做 MIC 回绕再聚合；
未包裹坐标亦正确（实测 + 单测锁定，见 `Tests/unittests/test_pbc.py`）。对应 7 个
feat 提交 + 1 个 DRY 重构提交 + 29 个 PBC 单测。

**合并状态（2026-10-29，✅ 已合并）**：feature/pbc-handling → dev（`870872e`）→
master（`a9dca03`，tag `v0.0.4`）；dev 与 master 已同步至 `2a6c824`。合入时全量回归通过
（276 tests），三斜盒/未包裹坐标已有单测锁定。

**仍遗留**：
- `per_tuple` 策略保持非 PBC（将被舍弃，跨盒对有差异，limitations 注明）
- 水桥 two_pass D–A 距离下界（2.5 Å）：既有 TODO #15，另案处理

### 17. 正电残基的铵/胍氮被误判为 H 键受体 ✅（已修复 2026-09-30）

**现状**：受体判定逻辑（`_find_acceptors`）为"类型 ∈ 受体表 且 q<0"。各力场受体表均含氮类（如 OPLS `opls_287`/CHARMM `NH3`/GROMOS `NL` 等铵氮类型）。带负电荷的铵氮（如 OPLS LYSH NZ q=-0.30、CHARMM 质子化 His 咪唑氮 q=-0.51）因此被识别为 H_acceptor。

**化学问题**：铵氮（N⁺，无孤对电子）与胍氮不应作为 H 键受体；正电残基的氮被判为受体可能导致氢键检测假阳性。

**来源**：继承自 `amber_ff_identifier.py` 的 `ACCEPTOR_TYPES`（amber 同样含 N3 铵类），属跨力场既有行为，非 OPLS 独有。

**待处理**：
- 参考 PLIP/化学定义，明确"受体 = 有孤对电子的 N/O/S"判定
- 在 `_find_acceptors` 或受体表中排除正电残基的铵/胍氮（如 LYSH NZ、ARG NE/NH1/NH2、质子化 His 咪唑氮）
- 注意：需与 H_donor 判定（N-H 键）协调，避免剔除参与供体的氮

**实现位置**：`group_identifiers/*_ff_identifier.py` 的 `_find_acceptors`

> **修复状态（2026-09-30，✅ 已修复）**：四力场（amber/gromos/charmm/opls）已按
> "力场类型 + 化学事实 + 文献"锚定完成修复——剔除普通酰胺/铵/带 H 吡咯/胍基的受体资格，
> 保留 Pro N / His 无 H 吡啶 / 中性胺 / 核酸氨基（含三待查项查证：GAFF n2=亚胺受体、
> CGenFF NG* 剔 10 留 6、核酸糖苷 N* 剔）。
> 逐项判定依据见 `doc/acceptor_identification_evidence.md`；实证效果：
> GROMOS 真实体系 acceptor 336→206（=减 130 主链 N）、H 键 pair 108→86、水桥 2541→1698。
> 本条关闭。

**优先级**：已关闭（2026-09-30 修复；2026-10-29 状态更新——本次核查（`doc/group_identification_review.md`）实测
Amber/CHARMM/GROMOS 真实体系无带 H 的 N 被判受体，进一步佐证修复有效）

### 18. 识别器 H 选择非确定性（含等价氢残基）✅（2026-11-03 对抗性审查关闭）

**原现象（2026-10-14 记录）**：含等价氢的残基（如 ARG 胍基 NH2 的 HH21/HH22、LYS NH3 的多 H）
在基团鉴定时，供体基团 atoms 列表的 H 选择曾被认为受 Python set/dict 迭代顺序（哈希种子
随机）影响——同一体系两次运行可能选到不同等价 H（当时实测：ARG NH2 基团一次产出
`[NH2, HH22]`、另一次 `[NH2, HH21]`），导致水桥/氢键角度（theta/omega）跨运行不稳定。

**对抗性审查（2026-11-03）→ 关闭**：**当前代码下该现象不可复现**。证据：

1. **代码机制**：`_find_donors` 遍历 `res.bonds`（有序列表），**每个 D-H 键独立建组**
   （`atoms=[D, H]`）——ARG NH2 的两个等价 H（HH21/HH22）各自成为**独立供体组、同时存在**，
   不存在"从多个等价 H 中选择一个"的逻辑。`467d3c1`（2026-08-28，H_donor 原子结构优化：
   `metadata={"h_atom": h_atom}` → `atoms=[d_atom, h_atom]`）已结构性消除"选择一个 H"的路径。
2. **哈希依赖**：PYTHONHASHSEED 只影响 str/bytes 哈希，不影响 int；全部决定输出顺序的
   set/dict 遍历（去重 key、水桥 `set(dw_by_w)` 等）键均为 **int**（全局索引/group_id），
   迭代顺序稳定。
3. **实测（多种子 × 多策略 × 多运行）输出签名全部一致**：
   - 识别器：4 力场（Amber/GROMOS/CHARMM/OPLS）× PYTHONHASHSEED∈{0,42,123,random}，
     151985/1454/5677/113 基团的 SHA256 签名完全一致
   - 完整 pipeline（识别→检测→h5 含 exists 矩阵）：two_pass 6 种子 + per_frame 5 种子 +
     10 次随机 seed，water_bridge/hydrogen_bond 签名全部一致
   - per_frame 水桥 5 次独立进程：n_pairs=1698、存在帧签名 `d7973715...` 完全一致

**原"发现过程"说明**：PBC 改造对比时的 62 个"存在帧差异"三元组，经"固定 groups 后
old/new 判定差异 0/0.000%"反推为 groups 层差异；但该差异在当前代码下无法复现（见上），
应为当时中间版本（H 存于 metadata 时期）的产物，或 PBC 改造过程临时代码所致。

**遗留建议（非阻塞）**：可补一个"确定性守护测试"（固定 seed 双跑 assert 完全一致）防止
未来改动重新引入；已登记为可选加固项。

**实现位置**：`group_identifiers/*_ff_identifier.py` 的供体基团（H_donor）构建处

**优先级**：已关闭（2026-11-03，对抗性审查确认当前代码无此问题）

### 19. GROMOS 氢键 two_pass 测试基线过期（287→292）✅（2026-10-29 核对并更新）

**现象**：`Tests/unittests/test_gromos_interactions_two_pass.py::TestHydrogenBondTwoPass::test_n_pairs`
断言 `n_pairs == 287`，实际产生 **292**（+5），测试失败。

**发现过程**：核查 GROMOS 中性胺受体判据（3c，见 `acceptor_identification_evidence.md` §2/§3）
时跑全量回归发现。经 `git stash` 回退修改 + 切到 master（v0.0.4）验证：**master 上同样失败
（292 != 287）**——属**既有问题**，与 3c 修改无关（该体系 NT 类型仅 ASN 酰胺，判据不变）。

**根因推测**：基线 287 记录于 `f2b5796`（A1 受体修复），此后 v0.0.4 的 PBC 改造（MIC 距离
计算）改变了跨周期边界氢键对的判定，但该测试基线未同步更新。

**核对结论（2026-10-29）**：**292 为正确值**。逐帧对比 MIC 距离与欧氏距离：满足 MIC 判据
（D-A ≤ 4.1 Å）而欧氏距离不满足的「对-帧」组合共 6 个，涉及 5 个**跨周期边界真实氢键对**，
例如 SER35→5ZIN134（MIC 3.8 Å / 欧氏 78.6 Å）、ASN36→5ZIN134（MIC 3.2 Å / 欧氏 79.3 Å）。
即 PBC 修正正确检出了旧基线（287，非 PBC 时代）漏检的跨盒氢键——**287 已过期**。

**处理**：`test_n_pairs` 断言更新 287→292 并注释说明；文件头部注释氢键数同步 379→292。
其余力场（Amber/CHARMM）的氢键基线经全量回归验证**未过期**。

**优先级**：已关闭（2026-10-29）

---

*文档结束*
