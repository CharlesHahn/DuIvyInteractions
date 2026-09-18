# Python API

本文档说明 DuIvyInteractions 的 Python 接口，供脚本化使用与二次开发。命令行用法见[命令参考](../guide/command.md)。

## 包结构

```
DuIvyInteractions/
├── pipeline.py              # Pipeline 主流程
├── system_readers/          # 拓扑读取器
├── group_identifiers/       # 基团识别器
├── interaction_detectors/   # 相互作用检测器（8 类型 × 3 策略）
├── io/                      # 序列化 + 导出
└── core/                    # 数据类、接口、常量
```

## 顶层流程：Pipeline

```python
from DuIvyInteractions.pipeline import Pipeline

# 构造：力场 + 策略
pipeline = Pipeline(ff="amber", strategy="two_pass")

# 运行：tpr + xtc → 识别 → 检测 → 存 h5
pipeline.run("md.tpr", "md.xtc", "out/", interactions=None)
# interactions=None 表示检测全部 8 类；也可传子集列表
```

参数：

| 参数 | 类型 | 说明 |
|:-----|:-----|:-----|
| `ff` | str | 力场名，当前仅 `"amber"` |
| `strategy` | str | `"two_pass"`（默认）/ `"per_frame"` / `"per_tuple"` |
| `interactions` | list[str] \| None | 检测类型子集，None=全部 |

## 数据读取（system_readers）

```python
from DuIvyInteractions.system_readers import GmxTprReader, GmxTprDumpReader

sd = GmxTprReader().read("md.tpr")        # 读二进制 tpr（MDAnalysis）
# sd: SystemData（原子/残基/键）
```

`GmxTprDumpReader` 从 `gmx dump` 文本解析，作为二进制读取的对照/兜底。

## 基团识别（group_identifiers）

```python
from DuIvyInteractions.group_identifiers import AmberFFGroupIdentifier

groups = AmberFFGroupIdentifier().identify(sd)   # -> List[Group]
```

`IDENTIFIER_CLASSES` 为识别器注册表（`--ff` 参数的选择来源）。

## 相互作用检测（interaction_detectors）

检测器按类型 × 策略命名，如 `SaltBridgeDetectorTwoPass`、`HydrogenBondDetectorPerFrame`。统一接口：

```python
from DuIvyInteractions.interaction_detectors import SaltBridgeDetectorTwoPass
import MDAnalysis as mda

# 检测前需按检测器所需基团类型过滤（与 Pipeline 内部一致）
det = SaltBridgeDetectorTwoPass()
filtered = [g for g in groups
            if g.group_type in det.required_group_types]

u = mda.Universe("md.tpr", "md.xtc")
results = det.detect(filtered, trajectory=u.trajectory)  # -> List[Interaction]
```

> 提示：日常使用推荐直接调 `Pipeline.run()`，其内部已完成"识别 → 按 required_group_types 过滤 → 检测 → 存 h5"全流程。直接使用检测器适用于需要精细控制检测参数的场景。

三策略基类（`core/interfaces.py`）：`InteractionDetectorPerTuple` / `PerFrame` / `TwoPass`，`detect()` 接口一致、结果统一为 `List[Interaction]`。`pipeline.py` 的 `DETECTOR_CLASSES` 注册全部 24 个检测器（8 类型 × 3 策略）。

## 结果序列化（io.h5）

```python
from DuIvyInteractions.io import save_interactions, load_interactions

save_interactions(results, "out/salt_bridge.h5")   # List[Interaction] -> h5
its = load_interactions("out/salt_bridge.h5")      # h5 -> List[Interaction]
```

- `path` 接受 `str` 或 `pathlib.Path`
- `compress=True`（默认）启用 gzip 压缩
- 无损往返（见[数据格式](data_format.md)）

## 结果导出（io.exporters）

```python
from DuIvyInteractions.io import SaltBridgeExporter

it = load_interactions("out/salt_bridge.h5")[0]
exp = SaltBridgeExporter()
exp.save_xvg_count(it, "sb_count.xvg")             # 每帧活跃数
exp.save_xvg(it, "distance", "sb_distance.xvg")    # 指标时间序列
exp.save_xpm(it, "sb_existence.xpm")               # 存在性热力图
exp.to_csv_summary(it, "sb_summary.csv")           # 每对汇总
```

8 个导出器（`HydrogenBondExporter` / `PiStackingExporter` / `SaltBridgeExporter` / `HydrophobicExporter` / `HalogenBondExporter` / `MetalCoordinationExporter` / `WaterBridgeExporter` / `PiCationExporter`）继承自 `InteractionExporter` 基类。

## 数据类（core.datas）

| 类 | 说明 |
|:---|:-----|
| `SystemData` | 体系数据（残基列表、原子、分子间键） |
| `ResidueData` | 残基（残基名/号、分子名、原子列表） |
| `AtomData` | 原子（全局索引、名称、力场类型、元素、电荷、质量） |
| `BondData` | 键（两端原子索引、键级） |
| `Group` | 化学基团（group_id/group_type/molecule/residue/atoms/metadata） |
| `Interaction` | 一种类型的全部检测结果（groups/existence/metrics/times） |
| `InteractionSparse` | TwoPass Pass1 的稀疏中间结果 |

`Interaction` 矩阵式存储：`existence[i][j]` 为第 i 个 pair 第 j 帧是否存在，`metrics["distance"][i][j]` 为对应距离。详见[数据格式](data_format.md)。