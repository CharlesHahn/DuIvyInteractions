# 水桥判定设计（基于显式水氢的氢键标准）

> 状态：设计定稿（待实现）
> 分支：feature/water-bridge-hbond
> 日期：2026-11-05

## 1. 动机

### 1.1 现有实现的问题

现有水桥检测（`water_bridge_detector_*.py`）复现 PLIP 的判据：仅使用**水的氧（OW）**，
配合 PLIP 在**水氧处**的替代角（`omega`，受体–水O–供体H，71–140°）来近似水的取向，
再用 `refine_water_bridges` 式的事后裁剪补救"一个水被重复计入多个桥"的问题。

三个具体缺陷：

| # | 缺陷 | 后果 |
|---|---|---|
| 1 | 不使用水的两个氢（HW1/HW2），取向信息被丢弃 | 无法判断"水是否真的向受体供出氢键"，只能靠 `omega` 近似 |
| 2 | 重原子距离判据（D–A 2.5–4.1 Å）+ 宽松角度（≥100°） | 纳入 H 实际很远（如 H···O = 3.8 Å）、角度很差（如 45°）的"假桥" |
| 3 | 无配额收敛 | 同一水被 2 供体 × 5 受体组合成 8 个桥（实测单帧最大 13 桥），物理上不可能 |

### 1.2 我们的条件

本项目输入为 **GROMACS tpr（全原子拓扑）**，水的显式氢（HW1/HW2）**始终存在**：

| 力场 | 水原子 | 实测 |
|---|---|---|
| Amber（TIP3P） | `OW, HW1, HW2` | 37021 个水全部含 2 个 H，零异常 |
| GROMOS（SPC） | `OW, HW1, HW2` | 全原子水 |
| CHARMM（TIP3） | `OH2, H1, H2` | 全原子水 |

PLIP 只用 OW 的原因是其输入为 **PDB 文件**，结晶水（HOH 残基）通常只记录 O 原子。
这是输入格式约束，不是化学理想——本项目有条件采用更物理的判据。

## 2. 判定标准

### 2.1 水桥的物理定义

水桥 = 同一水分子同时作为**受体**（接受 D–H）和**供体**（向 A 供出 O–H），
桥接两个极性基团：

```
       D – H ··· O_w        （水作受体，② 侧）
                   \
                    H_w ··· A   （水作供体，① 侧）
```

### 2.2 两个氢键的几何判据

#### ① 供体侧：D–H ··· O_w

| 指标 | 定义 | 阈值 | 常量 |
|---|---|---|---|
| `dist_h_ow` | d(H, O_w) | ≤ 2.5 Å | `WATER_BRIDGE_HB_DIST_MAX` |
| `angle_dh_ow` | ∠(D–H···O_w)，顶点 **H** | ≥ 120° | `WATER_BRIDGE_HB_ANGLE_MIN` |
| 下界 | d(H, O_w) | > 0.5 Å | `WATER_BRIDGE_MIN_DIST` |
| `dist_d_ow` | d(D, O_w) | **仅诊断输出，不参与判定** | — |

#### ② 受体侧：O_w–H_w ··· A

| 指标 | 定义 | 阈值 | 常量 |
|---|---|---|---|
| `dist_hw_a` | d(H_w, A) | ≤ 2.5 Å | `WATER_BRIDGE_HB_DIST_MAX` |
| `angle_ow_hw_a` | ∠(O_w–H_w···A)，顶点 **H_w** | ≥ 120° | `WATER_BRIDGE_HB_ANGLE_MIN` |
| 下界 | d(H_w, A) | > 0.5 Å | `WATER_BRIDGE_MIN_DIST` |
| `dist_ow_a` | d(O_w, A) | **仅诊断输出，不参与判定** | — |

**为什么不用重原子 D···A ≤ 3.5 Å 判据**（对抗性审查修正）：

`D···A ≤ 3.5 Å` 看似与 `H···A ≤ 2.5 + D–H ≈ 1.0` 自洽，但**力场中 D–H 键长并非常数**：

- 实测 D–H 键长范围 **0.95–1.34 Å**（中位 1.01，含长键的 X–H）
- 当 D–H = 1.34 Å 且 H···A = 2.5 Å 时，D···A 可达 **3.84 Å > 3.5** → 合法桥被误杀
- 实测误杀率 **0.4%**（243 个 H···A 达标的供体侧对中 1 个 D···A = 3.54 Å 被判否）

因此重原子距离仅作**诊断指标**输出，判定只用 `H···A` 与角度——这也是 Baker-Hubbard 的核心判据形式。

**两个 H 的处理**：水的 HW1、HW2 分别判定；对同一 `(供体, 水, 受体)` 三元组，
取"满足受体侧条件且 `angle_ow_hw_a` 更优"的那个 H，并记录其选择（`hw_id`，见 §3.1 第 1 层配额）。

#### ③ "水在中间"约束（保留 PLIP 的 omega 语义）

| 指标 | 定义 | 阈值 |
|---|---|---|
| `angle_h_ow_a` | ∠(供体 H – O_w – 受体 A)，顶点 **O_w** | **71° ≤ angle ≤ 140°** |

**为什么必需**：仅两侧各有氢键不足以构成"桥"——若 D 与 A 位于水的**同一侧**，
水并未桥接两端。实测缺此约束时 **27.7%** 的桥 `∠D-Ow-A < 90°`（D/A 同侧）。
加入后保留 72.3%（D927：8987 → 6501 桥-帧）。

**与 PLIP 的对应**：PLIP 的 `omega` 即 `∠(受体 A – 水O – 供体H)`，范围 [71°, 140°]；
本设计的 `angle_h_ow_a` 定义与阈值完全相同——这是"水桥定义符合既有判定规则"的
几何保证（保留 PLIP 对水在中间的位置约束）。

## 3. 配额（去重）规则

### 3.1 三层配额（逐帧独立执行）

| 层 | 规则 | 依据 |
|---|---|---|
| **1. 每 H_w 唯一** | 同一水的每个 H_w 最多对应 1 个受体 A；择优取 `angle_ow_hw_a` 最接近 180° | 一个 O–H 的 σ* 反键轨道只能接受一个孤对 |
| **2. 每水 O 最多 2 个供体** | 同一水的 O_w 最多对应 2 个 D–H；择优取 `dist_h_ow` 最小 | 水有 2 对孤对，每对接受 1 个 D–H |
| **3. 每水最多 2 个桥** | 每水每帧最多输出 2 个三元组；按供体侧/受体侧各自排名配对（第 i 名配第 i 名） | 每桥消耗 1 对孤对 + 1 个 H；水有 2 对孤对和 2 个 H → 最多 2 桥 |

**配对说明**：`(D, W, A)` 三元组中 D 与 A 的物理关联较弱（二者仅共享同一水），
"第 i 名配第 i 名"是**输出表示法的约定**，不代表物理配对；配额本身保证每水 ≤ 2 桥。

### 3.2 配额在架构中的位置

配额是**跨对操作**（同一水的多个候选对之间择优），无法在逐对阈值判定中完成。
采用"先几何阈值、后配额收敛"两步：

```
每帧：compute_metrics → apply_threshold（几何）→ _apply_water_quota（配额）→ existence
```

`_apply_water_quota` 为水桥检测器专有方法，不修改基类模板。

## 4. 文献依据

| 判据 | 文献/来源 | 强度 |
|---|---|---|
| **H···A ≤ 2.5 Å** | Baker & Hubbard 1984, *Proteins* 4:1–13（"Hydrogen bonding in globular proteins"，被引 2555）；IUCr 权威转述："The distance **H···A is less than 2.5 Å**" | 经典标准 |
| **D-H···A ≥ 120°** | MDTraj 的 "Baker-Hubbard Hydrogen Bond Identification" 实现：`θ > 120 and r_H···Acceptor < 2.5 Å`；Schrödinger 判据；HBAT 经典氢键默认值（120°） | 多工具一致 |
| **D···A ≤ 3.5 Å** | 与 H···A 2.5 + D–H 键长 ≈1.0 Å 自洽 | 逻辑自洽 |
| **每水 ≤ 2 桥** | 水的 2 个 H + 2 对孤对，每桥消耗 1 H + 1 孤对（教科书水桥结构；冰/液态水 4 配位） | 化学约束 |
| **水的氢键几何** | 水氢键 O···O 平均 2.8 Å、偏离线性平均 14°（即 ∠≈166°）→ 达标氢键通常远优于 2.5 Å/120° 边界 | 实测统计 |

**未采用的判据（诚实披露）**：Schrödinger 完整判据含第三条件 `H···A-AA ≥ 90°`（AA = 与受体相连的原子），
用于约束氢键方向性。本设计**不采用**该条件：① 本项目受体 A 的孤对方向未建模；② 水作为供体的
氢键方向性弱于蛋白内部氢键；③ 加入会使判据过严、显著漏检。若后续需要更严格模式，可作为可选开关。

**与 PLIP 判据的关系**（说明差异来源，不作为对标目标）：

| 项 | PLIP | 本设计 |
|---|---|---|
| 供体侧距离 | 重原子 D–O 2.5–4.1 Å | H···O ≤ 2.5 Å（重原子 ≤ 3.5 Å 校验） |
| 供体侧角度 | `theta`（水O–供体H–供体D）≥ 100° | `angle_dh_ow` ≥ 120° |
| 受体侧 | `omega`（受体–水O–供体H）71–140° | `angle_ow_hw_a`（顶点 H_w）≥ 120° |
| 收敛 | `refine_water_bridges`（供体-氢键互斥 + (水,受体)去重 + 每水2桥，择优用 \|110°−omega\|） | 三层配额，择优用真实氢键角度/距离 |

> 本设计**不追求与 PLIP 数值对齐**：判据更严格，水桥数将显著减少（实测 D927 单帧 249 → 53）。

## 5. 实现架构

### 5.1 指标定义（metric_names）

| 指标 | 说明 | 兼容旧字段 |
|---|---|---|
| `dist_h_ow` | 供体 H 到水 O 距离 | 新 |
| `angle_dh_ow` | ∠D–H···O_w | 新（几何等价于旧 `theta`，阈值 100°→120°） |
| `dist_d_ow` | 供体 D 到水 O 距离 | 新（对应旧 `dist_dw`） |
| `dist_hw_a` | 水 H 到受体距离 | 新 |
| `angle_ow_hw_a` | ∠O_w–H_w···A | 新（取代旧 `omega`） |
| `dist_ow_a` | 水 O 到受体距离 | 新（对应旧 `dist_wa`） |
| `hw_id` | 受体侧选择的 H（0=HW1，1=HW2） | 新（配额第 1 层用，不导出） |

**导出兼容**：旧字段 `dist_dw`/`dist_wa`/`theta`/`omega` 由新字段替代。
具体映射如下（导出器 `DuIvyInteractions/io/water_bridge_exporter.py` 同步）：

| 旧字段 | 旧标签 | 新字段 | 新标签 |
|---|---|---|---|
| `dist_dw` | D-Ow Distance | `dist_d_ow` | D···Ow Distance (Å) |
| `dist_wa` | Ow-A Distance | `dist_ow_a` | Ow···A Distance (Å) |
| `theta` | O-H···Ow Angle | `angle_dh_ow` | D-H···Ow Angle (°) |
| `omega` | H-Ow···A Angle | `angle_ow_hw_a` | Ow-Hw···A Angle (°) |
| — | — | `dist_h_ow`（新） | H···Ow Distance (Å) |
| — | — | `dist_hw_a`（新） | Hw···A Distance (Å) |

`hw_id` 仅内部配额使用，不列入导出。

### 5.2 三个策略的改动

| 策略 | 改动 |
|---|---|
| **per_frame** | ① `metric_names` 更新；② **第一帧**候选预筛改为 `H(供体)↔O_w` 与 `H_w↔A`（候选仅基于第一帧，沿用既有设计原则）；③ `_compute_metrics` 计算新指标；④ `apply_threshold` 用新阈值；⑤ `detect` 逐帧插入 `_apply_water_quota` |
| **two_pass** | ① `metric_names` 更新；② `run_pass1` 预筛改用 H/H_w；③ `compute_pair_metrics` 计算新指标；④ `apply_threshold` 用新阈值；⑤ 基类 `run_pass2` 需在阈值后调用 `_apply_water_quota`（通过新增可选钩子，默认 no-op，保持其他检测器不受影响） |
| **per_tuple** | **不采用新判据**（见下） |

**per_tuple 策略的架构限制**：该策略为"逐候选组遍历全部帧"，无法实现跨候选组的配额
（配额需要同一水的多个候选一起比较）。按 `AGENTS.md` 约定，该策略处于"可能被舍弃"状态、
测试不跑。因此 per_tuple 保留 PLIP 旧判据并加文件头注释说明，不参与新判据实现——
三个策略中 per_frame 与 two_pass 采用新判据；二者帧级结果一致（实测帧 0 均为 53 桥），
但全轨迹对数不同（per_frame 候选仅第一帧=86，two_pass 累积候选=3130）——这是两策略的既有设计差异。

### 5.3 基类钩子（最小改动）

基类 `InteractionDetectorTwoPass.run_pass2` 已有 3 个钩子（`_build_indices_from_sparse`、
`_init_metrics`、`_get_box`）。本设计新增第 4 个钩子，插入位置在
`apply_threshold` 之后、写入 `existence` 之前：

```python
existence[:, f] = self._post_threshold(
    self.apply_threshold(frame_metrics), frame_metrics, group_tuples)
```

```python
def _post_threshold(self, ok, frame_metrics, group_tuples):
    """阈值判定后的可选后处理（默认原样返回）。"""
    return ok
```

水桥 `two_pass` 覆写该钩子执行配额；其他 7 类检测器使用默认实现，行为不变。

**配额所需的两个信息**及其通道：

| 信息 | 来源 | 是否入 `metric_names` |
|---|---|---|
| 每个对的水分子身份 | `group_tuples[i][1].group_id`（水基团） | 否（从元组直接取） |
| 受体侧选中的 H（0/1） | `frame_metrics["hw_id"]`（`compute_pair_metrics` 返回） | **否**（仅配额内部使用，不导出） |

`metric_names` 只列 6 个几何指标；`hw_id` 作为 `frame_metrics` 的额外键传递，
不进入稠密存储与导出。

**两策略共享配额实现**：`_apply_water_quota(ok, hw_id, water_uids)` 实现为本模块的公共方法，
`per_frame` 在逐帧循环中直接调用，`two_pass` 通过 `_post_threshold` 调用——同一份逻辑，避免重复。

### 5.4 水自身基团的排除（必须）

识别器把水的 O 识别为 `H_acceptor`（37021 个）、水的 H 识别为 `H_donor`（74042 个）。
水桥检测必须**排除属于水的供体/受体基团**，否则会把水-水氢键误算为水桥
（实测：不排除则单帧桥数从 53 爆炸到 59701）。现有 `per_frame.detect` 已用
`water_atoms` 集合实现该排除；`two_pass` 的基类 `_filter_groups` 也已按
`WATER_RESIDUES` 排除——实现时须确认两个策略行为一致。

### 5.5 PBC

全部距离/角度沿用现有 PBC 机制（`calc_bonds(box)`、`calc_angles(box)`、
`capped_distance(box)`）；预筛与配额不受 PBC 影响（配额是索引分组）。

## 6. 退化策略

**触发条件**：水基团不含 2 个 H 原子（如非全原子力场、PDB 输入）。

**处理**：
- 检测水的 H 数量；不足 2 个时，**跳过水桥检测并输出警告**
  （`[WARN] water_bridge: water has no explicit H, skipped`）
- 理由：该判据完全依赖显式 H；退化到 PLIP 方式会引入两套语义并增加复杂度，不符合 KISS

> 现有测试体系中 Amber 水有完整 H；GROMOS/CHARMM 测试体系为复合物（无水），
> 故退化分支在本仓库现有数据上不触发。

## 7. 影响面

### 7.1 代码改动清单

| 文件 | 改动 | 必需性 |
|---|---|---|
| `core/interfaces.py` | 新增 `_post_threshold` 钩子（默认 no-op）+ `run_pass2` 调用点 | 必需 |
| `water_bridge_detector_per_frame.py` | 指标、预筛、阈值、配额调用 | 必需 |
| `water_bridge_detector_two_pass.py` | 指标、预筛、阈值、`_post_threshold` 覆写 | 必需 |
| `water_bridge_detector_per_tuple.py` | 指标与判据同步（测试不跑，但保持代码一致） | 一致必 |
| `io/water_bridge_exporter.py` | 4 个字段名 + 标签更新，新增 2 个字段 | 必需 |
| `Tests/unittests/test_water_bridge_per_frame.py` | 指标名、阈值断言、数量基线重建 | 必需 |
| `Tests/unittests/test_water_bridge_two_pass.py` | 同上 | 必需 |
| `Tests/unittests/test_water_bridge.py` | （per_tuple，不跑）指标名同步 | 一致可选 |
| `docs/reference/criteria.md` | 水桥判据章节更新（用户可见） | 必需 |
| `docs/guide/result.md` | 水桥指标说明更新 | 必需 |
| `docs/reference/limitations.md`（中英） | 删除"水桥水分子去重未做"条目 | 必需 |

### 7.2 行为影响

| 项 | 影响 |
|---|---|
| **水桥数量** | 显著减少（D927 单帧实测 249 → 53，-79%） |
| **测试基线** | 水桥相关数量断言（`test_n_pairs == 1698`、`14598`/`14562`）需全部重建 |
| **导出字段** | 4 个改名 + 2 个新增 |
| **与其他检测器** | 无影响（基类钩子默认 no-op） |
| **与 PLIP 对标** | 不再可比（判据不同） |
| **per_tuple 策略** | 代码同步但测试不跑（按 AGENTS 约定） |

## 8. 验证计划

| 项 | 方法 |
|---|---|
| 判据正确性 | 构造合成水桥（已知几何）验证阈值边界 |
| 每个桥满足全部判据 | 逐桥重算 6 个指标，断言无违反（实测 53 桥/帧零违反） |
| 配额有效性 | 断言每水每帧桥数 ≤ 2；每 (水, H_w) ≤ 1 个 A；每水 O ≤ 2 个 D |
| 水自身排除 | 断言结果中不含"水-水-水"三元组 |
| PBC | 复用 `test_pbc.py` 的水桥用例，确认 MIC 下结果一致 |
| 真实数据 | D927 全轨迹跑通，新基线记录（数量 + top pair + occupancy） |
| 退化分支 | 构造无水 H 的 SystemData，断言警告且返回空 |

## 9. 关联文档

- `doc/plip_interaction_rules.md`（水桥判据参考）
- `doc/pbc_design.md` §6B（水桥 PBC 处理）
- `doc/TODO.md` #11（本设计的来源）
- `doc/group_atom_convention.md`（水 Group 原子约定）
