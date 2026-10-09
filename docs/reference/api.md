# Python API

本文档说明 DuIvyInteractions 的 Python 接口，供脚本化使用与二次开发。命令行用法见[命令参考](../guide/command.md)。所有接口签名以 `DuIvyInteractions/core/interfaces.py`、`core/datas.py`、`io/h5.py` 与 `pipeline.py` 为准。

## 包结构

```
DuIvyInteractions/
├── pipeline.py              # Pipeline 主流程（编排 Reader→Identifier→Detector→h5）
├── DII.py                   # 命令行入口（dii run / dii export）
├── system_readers/          # 拓扑读取器（GmxTprReader / GmxTprDumpReader）
├── group_identifiers/       # 基团识别器（amber/gromos/charmm/opls 4 力场，IDENTIFIER_CLASSES 注册表）
├── interaction_detectors/   # 相互作用检测器（8 类型 × 3 策略，共 24 个类）
├── io/                      # h5 序列化 + 8 个导出器
└── core/                    # 数据类（datas）、接口（interfaces）、常量（constants）
```

## 顶层流程：Pipeline

`Pipeline` 编排"读取 tpr → 基团识别 → 按类型检测 → 写 h5"全流程。构造参数 `ff` 支持 4 个力场：

```python
from DuIvyInteractions.pipeline import Pipeline

# ff: "amber" | "gromos" | "charmm" | "opls"（与 dii run --ff 一致）
pipeline = Pipeline(ff="amber", strategy="two_pass")

# 运行：tpr + xtc → 识别 → 检测 → 在 output/ 下写 <类型>.h5
pipeline.run("md.tpr", "md.xtc", "out/", interactions=None)
# interactions=None 表示检测全部 8 类；也可传子集列表
```

| 参数 | 类型 | 说明 |
|:-----|:-----|:-----|
| `ff` | str | 力场名：`"amber"` / `"gromos"` / `"charmm"` / `"opls"`（`IDENTIFIER_CLASSES` 的键；未知力场抛 `ValueError`） |
| `strategy` | str | `"two_pass"`（默认，性能最优）/ `"per_frame"` / `"per_tuple"` |
| `interactions` | list[str] \| None | 检测类型子集（`pipeline.ALL_INTERACTIONS` 中的名称），None=全部 8 类 |

`run()` 的行为细节：

- 基团识别只做一次、与帧无关；轨迹只加载一次（`mda.Universe(tpr, xtc)`）。
- 水残基排除使用识别器自带的 `WATER_RESIDUES`（见下文）：`_filter_groups` 对除水桥外的类型剔除残基名 ∈ `WATER_RESIDUES` 的基团，保证跨力场正确（如 CHARMM TIP3、OPLS HO4/HO5）。
- 逐类型检测并保存 `<output>/<类型>.h5`；单个类型失败只打印 `[WARN]`，不中断其余类型。
- 全部 8 类名称见 `pipeline.ALL_INTERACTIONS`；类型→检测器注册表为 `DETECTOR_CLASSES`（8 类型 × 3 策略，24 个类，键值为 `(TwoPass类, PerFrame类, PerTuple类)` 三元组），`STRATEGY_INDEX = {"two_pass": 0, "per_frame": 1, "per_tuple": 2}` 负责策略切换。

## 数据读取（system_readers）

```python
from DuIvyInteractions.system_readers import GmxTprReader, GmxTprDumpReader

sd = GmxTprReader().read("md.tpr")        # 读二进制 tpr（MDAnalysis）
# sd: SystemData（残基/原子/键/残基间键）
```

- `GmxTprReader`：经 MDAnalysis 读取二进制 tpr（主线）。
- `GmxTprDumpReader`：解析 `gmx dump` 文本输出，作为二进制读取的对照/兜底。
- 两者实现同一 `Reader` 接口（`core/interfaces.py`）：`name` 属性 + `read(source) -> SystemData`。

## 基团识别（group_identifiers）

4 个力场识别器类（均在 `DuIvyInteractions/group_identifiers/`，后三者子类化 `AmberFFGroupIdentifier` 复用力场无关逻辑）：

| 类 | 力场 | 水残基名（`WATER_RESIDUES` 类属性） |
|:---|:-----|:-----|
| `AmberFFGroupIdentifier` | Amber 家族（amber03/94/96/99/99SB/99SB-ildn/GS/14SB）+ GAFF/GAFF2 | `{"SOL", "HOH", "WAT"}` |
| `GromosFFGroupIdentifier` | GROMOS 53A6/54A7（联合原子，极性 H 显式） | `{"SOL"}` |
| `CharmmFFGroupIdentifier` | CHARMM36/C36m + CGenFF | `{"TIP3", "HOH", "SOL", "WAT"}` |
| `OplsFFGroupIdentifier` | OPLS-AA/L | `{"HOH", "HO4", "HO5", "SOL", "WAT"}` |

```python
from DuIvyInteractions.group_identifiers import AmberFFGroupIdentifier

groups = AmberFFGroupIdentifier().identify(sd)   # -> List[Group]
```

- 接口见 `core/interfaces.py` 的 `GroupIdentifier`：`name` 属性 + `identify(SystemData) -> List[Group]`。
- **`WATER_RESIDUES` 是类属性**：基类默认为空 `frozenset()`，子类必须覆盖为非空集合。pipeline 的 `_filter_groups` 与识别器内部的 `_find_water`/`_find_metal_binding` 统一经 `self.WATER_RESIDUES` 取用，是跨力场水排除的唯一入口（如 CHARMM TIP3、OPLS HO4/HO5）。
- `IDENTIFIER_CLASSES`（`group_identifiers/__init__.py`）为注册表：`{"amber", "gromos", "charmm", "opls"}`，是 `dii run --ff` 的 choices 与 `Pipeline(ff=...)` 的选择来源（`DII.py` 的 `--ff` 参数直接取 `IDENTIFIER_CLASSES` 的键）。

## 相互作用检测（interaction_detectors）

检测器按 类型 × 策略 命名，如 `SaltBridgeDetectorTwoPass`、`HydrogenBondDetectorPerFrame`。统一接口：

```python
from DuIvyInteractions.interaction_detectors import SaltBridgeDetectorTwoPass
import MDAnalysis as mda

det = SaltBridgeDetectorTwoPass()
# 按检测器所需基团类型过滤（与 Pipeline 内部一致）
filtered = [g for g in groups
            if g.group_type in det.required_group_types]

u = mda.Universe("md.tpr", "md.xtc")
results = det.detect(filtered, trajectory=u.trajectory)  # -> List[Interaction]
```

`detect()` 完整签名（三个策略基类一致，见 `core/interfaces.py`）：

```python
detect(groups, trajectory=None, n_workers=1,
       topology_path=None, trajectory_path=None,
       tuple_filter=None) -> List[Interaction]
```

| 参数 | 说明 |
|:-----|:-----|
| `groups` | 已按 `required_group_types` 过滤的基团列表 |
| `trajectory` | MDAnalysis 轨迹对象（串行必需；`PerFrame`/`TwoPass` 忽略 `n_workers`，必须传） |
| `n_workers` | 并行 worker 数（>1 时需同时提供 `topology_path` 与 `trajectory_path`；并行实现仅在 `PerTuple` 策略） |
| `tuple_filter` | 可选回调 `(Tuple[Group, ...]) -> bool`，对候选基团组做用户自定义过滤（如只保留不同分子间的 pair）；在候选生成后、几何计算前应用 |

- 三个策略基类：`InteractionDetectorPerTuple` / `InteractionDetectorPerFrame` / `InteractionDetectorTwoPass`（模板方法模式），`detect()` 接口一致、结果统一为 `List[Interaction]`。
- 子类须实现：`name`、`required_group_types`、`metric_names` 三个属性，以及策略对应的检测方法——PerTuple 为 `get_candidate_tuples`/`compute_metrics`/`apply_threshold`；PerFrame 为 `get_candidate_tuples`/`compute_metrics_for_frame`/`apply_threshold`；TwoPass 为 `initialize_candidates`/`compute_pair_metrics`/`apply_threshold`（后三者基类提供返回空值的默认实现，须覆写）。三者均可选覆写 `filter_candidate_tuples`（用首帧坐标预筛）与 `_post_process`（跨对后处理钩子）。
- 8 类型 × 3 策略共 24 个检测器全部注册于 `pipeline.DETECTOR_CLASSES`。
- 策略差异说明：`per_tuple` 策略（策略一）当前处于"可能被舍弃"状态，不作为现行结果（测试不跑，见[已知限制](limitations.md)）；现行策略为 `per_frame` 与 `two_pass`。

> 提示：日常使用推荐直接调 `Pipeline.run()`。直接使用检测器适用于需要精细控制候选生成/过滤（如 `tuple_filter`）的场景。

## 结果序列化（io.h5）

```python
from DuIvyInteractions.io import save_interactions, load_interactions

save_interactions(results, "out/salt_bridge.h5")   # List[Interaction] -> h5
its = load_interactions("out/salt_bridge.h5")      # h5 -> List[Interaction]
```

- `save_interactions(interactions, path, compress=True)`：`compress` 启用 gzip 压缩（默认开）。
- `path` 接受 `str` 或 `pathlib.Path`；空路径抛 `ValueError`，非 str/Path 类型抛 `TypeError`。
- 读取时校验顶层 `format_version`，与 `FORMAT_VERSION = "1.0"` 不匹配抛 `ValueError`。
- 无损往返（见[数据格式](data_format.md)），项目单测 `test_io_h5.py` 覆盖。

## 结果导出（io.exporters）

```python
from DuIvyInteractions.io import SaltBridgeExporter

it = load_interactions("out/salt_bridge.h5")[0]
exp = SaltBridgeExporter()
exp.save_xvg_count(it, "sb_count.xvg")             # 每帧活跃 pair 数
exp.save_xvg(it, "distance", "sb_distance.xvg")    # 数值指标时间序列
exp.save_xpm(it, "sb_existence.xpm")               # existence 热力图
exp.to_csv_summary(it, "sb_summary.csv")           # 每对汇总（occupancy + avg/std）
```

- 8 个导出器（`HydrogenBondExporter` / `PiStackingExporter` / `SaltBridgeExporter` / `HydrophobicExporter` / `HalogenBondExporter` / `MetalCoordinationExporter` / `WaterBridgeExporter` / `PiCationExporter`）继承自 `io/interaction_exporter.py` 的 `InteractionExporter`，均在 `io/__init__.py` 导出。
- 子类必须实现 `name`、`metric_labels` 两个抽象属性；`get_pair_label` 有默认实现（`"残基名残基号-残基名残基号"`，如 `ARG73-D927`），可覆写。
- `PiStackingExporter` 额外提供 `save_xpm_stacking_type`（π 堆积类型 P/T/N 热力图）。
- `to_csv_summary` 仅输出数值指标（字符串指标如 `pistacking_type` 跳过）；`dii export` 按 h5 内 `interaction_type` 查找对应导出器，未知类型报错退出。

## 数据类（core.datas）

| 类 | 说明 |
|:---|:-----|
| `SystemData` | 体系数据：`system_name`、`residues`（残基列表）、`inter_residue_bonds`（残基间键，如肽键/二硫键）；构造时校验残基/原子全局索引唯一性 |
| `ResidueData` | 残基：残基名/全局索引/分子内编号、所属分子名、原子列表、残基内键列表 |
| `AtomData` | 原子：全局索引、残基内索引、原子名、力场类型、元素、电荷、质量 |
| `BondData` | 残基内键：两端残基内索引 + 键类型（`BOND_TYPES`） |
| `InterResidueBond` | 残基间共价键（两端残基全局索引 + 残基内原子索引 + 键类型） |
| `Group` | 化学基团：`group_id`/`group_type`/`molecule`/`residue_name`/`residue_id`/`atoms`/`metadata`；属性 `num_atoms`/`atom_indices`/`net_charge`；`group_type` 必须是 `GROUP_TYPES` 的值，`atoms` 非空 |
| `Interaction` | 一种类型的全部检测结果：`interaction_type`/`groups`/`existence`/`metrics`/`times`；属性 `n_pairs`/`n_frames`、方法 `occupancy()` |
| `InteractionSparse` | TwoPass Pass1 的稀疏中间结果（键为 `(group_id, ...)` 元组） |

- `Group.metadata`：键值字典，键必须为字符串，值限 JSON 可序列化类型（str/int/float/bool/None/list/dict；numpy 标量/数组自动兜底）。
- `Interaction` 矩阵式存储：`existence[i][j]` 为第 i 个 pair 在第 j 帧是否存在，`metrics["distance"][i][j]` 为对应几何指标；`occupancy()` 返回每对存在比例。详见[数据格式](data_format.md)。
