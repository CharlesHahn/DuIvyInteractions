# TODO 清单

> 记录推迟到后续阶段实现的任务，避免遗忘。

---

## 相互作用检测阶段需要实现的任务

以下任务在基团识别阶段有意推迟，需要在实现 `interaction_detectors/` 时一并完成。

### 1. 芳香环平面性检验

**推迟原因**：初始帧构象可能未弛豫，平面性可能很差；热涨落导致瞬时平面性波动；力场原子类型已编码芳香性，不需要几何验证。

**实现位置**：`interaction_detectors/pi_stacking.py`

**判据**：
- 对环内每个原子，计算相邻两键的法向量
- 所有法向量两两夹角 ≤ 5.0°（PLIP `AROMATIC_PLANARITY`）
- 不满足 → 该帧此环不参与 π-π 堆积

**参考**：`doc/aromatic_ring_definition.md` §3.3

### 2. 疏水-芳香去重

**推迟原因**：π-π 堆积已包含疏水接触（芳香环之间的 van der Waals 力），如果同时报告 π-π 堆积和疏水相互作用，会重复计数同一物理接触。

**实现位置**：`interaction_detectors/` 的主流程编排

**规则**：
- 先检测 π-π 堆积
- 再检测疏水相互作用
- 移除与 π-π 堆积重叠的疏水接触（两个芳香环的原子之间的疏水接触）

**参考**：`doc/hydrophobic_definition.md` §4


### 3. 基团识别结果需要人为审查！

目前基于真实数据识别得到的基团的数据，还需要人去手工一个一个核对一下确认一下。单元测试的真实数据测试是基于现有代码输出结果做的，不一定对。以人工审查过的结果为准！

### 4. 可选的候选对预过滤优化

**检测算法**：遍历基团对，每对向量化计算全部帧的距离。预加载基团原子坐标到内存（F×G×3），每对直接得到一行 existence 和 distance 向量。

**问题**：候选对数量可能很大（如氢键 ~74000×38000），虽然每对计算很快，但总对数多时仍耗时。

**优化方案**：提供可选的预过滤参数（如 `--prefilter`），用户在确认体系构象变化不大的情况下，可在第一帧用粗略距离 cutoff 去除明显不可能的候选对。

**注意**：
- 此优化**默认关闭**，因为构象变化大的体系会遗漏真实相互作用（科学错误）
- 是否启用由用户通过命令行参数自行判定
- cutoff 应远大于相互作用距离阈值（如氢键 4.1 Å → 预过滤用 15 Å）

### 5. 基团组用户自定义过滤（tuple_filter）

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

### 6. 官能团分类粗糙问题

**问题**：PLIP 的 `is_functional_group` 对氮基团分类粗糙。"tertamine" 实际包含所有 sp3 氮（NR₃、NR₃H⁺、NH₃⁺ 等），不是化学意义上的"叔胺"。当前代码复现了 PLIP 的定义，导致蛋白 N 末端 NH₃⁺ 被误判为 tertamine。

**待改进**：需要更精细的氮基团分类，区分叔胺、伯胺、仲胺等不同化学基团。

**影响**：pi-cation 检测的叔胺角度检查会误触发。

### 7. PLIP 对 tertamine 的 pi-cation 特殊处理

**现象**：PLIP 在 `pication()` 中对 tertamine 做了额外角度检查，其他正电基团（quartamine、guanidine、sulfonium）没有此检查。

**PLIP 的操作**：计算胺的三个邻居平面的法向量（`amine_normal`），与环法向量（`ring.normal`）取夹角。若夹角 ≤ 30°，排除该 pi-cation 对。

**PLIP 注释**：`"Special case here if the ligand has a tertiary amine, check an additional angle. Otherwise, we might have a pi-cation interaction 'through' the ligand."`

**待理解**：为什么只有 tertamine 需要这个检查？这个角度条件（≤ 30°）的物理含义是什么？为什么它能防止"穿过配体"的假阳性？

### 8. 金属配位几何构型匹配

**推迟原因**：第一版只做距离判据，几何构型匹配复杂度高。

**PLIP 的操作**：按金属分组后，计算所有配位原子的角度，与已知几何构型（linear/trigonal/tetrahedral/octahedral 等）比较，选择 RMS 最小的构型，移除不符合该构型的多余配位原子。

**实现位置**：`interaction_detectors/metal_coordination_detector.py`

### 9. 金属配位排除纯水配位

**推迟原因**：当前 `metal_binding` 基团已排除水分子，不需要额外处理。但如果未来改为不排除水，则需要此逻辑。

**PLIP 的操作**：如果所有配位原子都来自水分子，不报告该金属配位。

### 10. 水桥性能优化

**问题**：当前水桥检测候选三元组数量太大（24.9 万个），每个三元组需要逐帧遍历轨迹（~0.95 秒），总耗时约 65 小时，不可接受。

**根因**：`_process_tuple` 对每个三元组独立遍历轨迹，与其他检测器（候选对少）不同，水桥的三元组数量级太大。

**待优化**：减少候选三元组数量或优化轨迹加载方式。

### 11. 水桥水分子去重

**推迟原因**：第一版暂不做。

**PLIP 的操作**：一个水分子最多参与 2 个氢键（两个 H 各做一个供体）。如果超过 2 个候选，保留水分子 H-O-H 角度最接近 110° 的两个。

### 12. 水桥单元测试结果断言

**推迟原因**：性能问题导致无法跑完全量检测，无法获取真实数据结果。

**待补充**：解决性能问题后，补充具体的结果数据断言（水桥数量、top pair 等）。

### 13. PerFrame 检测器长轨迹内存优化

**问题**：PerFrame 检测器预分配 (n_pairs, n_frames) 的 distance/angle 矩阵。D927 体系氢键候选对 42,603 个，101 帧时矩阵约 71 MB，可接受。但外推到 1μs 轨迹（500,000 帧），矩阵将达到 342 GB，不可接受。

**根因**：当前实现对全部候选对存储每一帧的 distance/angle，即使大部分 pair 最终被阈值淘汰（42,603 候选 → 119 有效）。

**优化方案**：只存 existence 矩阵（bool，1 byte），不存 distance/angle。最终过滤后，对存活的 ~100 个 pair 重新遍历轨迹计算 distance/angle。

| | 当前方案 | 优化方案 |
|:--|:---------|:---------|
| 氢键 101 帧 | 71 MB | 4.1 MB (existence only) |
| 氢键 500,000 帧 | 342 GB | 20.3 GB (existence only) |
| 最终结果 | 同 | 同（重算存活 pair） |

**进一步优化**：existence 矩阵也可改为稀疏存储或增量累积，将 20.3 GB 降到更低。

**实现位置**：`interaction_detectors/*_per_frame.py` 的 `detect()` 方法。

**优先级**：低（当前 101 帧测试场景无压力，1μs 轨迹时再做）。

### 14. Interaction 结果序列化（保存/加载）

**需求**：将 Interaction 数据结构保存到文件，支持后续分析、可视化、跨工具共享。

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

### 15. 水桥 TwoPass 未添加 WATER_BRIDGE_MINDIST 下界

**现状**：策略三（TwoPass）的水桥检测器在 `apply_threshold` 中未添加 `dist > WATER_BRIDGE_MINDIST`（2.5Å）的下界检查，仅检查 `dist < 4.1Å` 上界。策略一（PerTuple）和策略二（PerFrame）均有此下界。

**影响**：TwoPass 会保留 Ow-A 距离 < 2.5Å 的三元组（如 2.48Å），而策略一/二会排除这些。D927 体系单帧测试中，TwoPass 多出 3 对（381 vs 378）。

**待决定**：是否需要在 TwoPass 的 `apply_threshold` 中添加此下界。

---

*文档结束*
