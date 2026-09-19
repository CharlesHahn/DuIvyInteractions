---
name: duivyinteractions-skill
description: 操作 DuIvyInteractions（dii）命令行工具：对 GROMACS tpr+xtc 轨迹运行分子相互作用检测（氢键、π-π 堆积、盐桥、疏水、卤键、金属配位、水桥、π-阳离子共 8 类），并将 h5 结果导出为 xvg/xpm/csv。当用户要求运行 dii、做分子相互作用分析、处理 MD 轨迹的相互作用，或导出/查看 .h5 相互作用结果时使用本技能。即使用户没有明确说出 "dii"，只要任务涉及 MD 轨迹中蛋白-配体/分子间相互作用（氢键、堆积、盐桥等）的检测或结果导出，也应主动使用本技能。
compatibility: Requires Python 3.9+ and GROMACS gmx dump (for text-format tpr parsing); installs via pip install duivyinteractions
---

# DuIvyInteractions（dii）操作指南

DuIvyInteractions（命令 `dii`）是分子动力学（MD）相互作用分析工具：直接读 GROMACS tpr 中的力场原子类型确定性识别化学基团，逐帧检测相互作用，结果存 h5，可导出 xvg/xpm/csv。本技能说明如何正确使用 `dii run`（检测）与 `dii export`（导出），以及常见坑。具体参数细节以 `dii run --help` / `dii export --help` 为准，本技能只讲要点。

## 环境

- `dii` 是已安装的命令行工具：PyPI 安装 `pip install duivyinteractions`，或源码目录 `pip install -e .`
- 两种调用方式等价：
  - `dii ...`（CLI 入口）
  - `python -m DuIvyInteractions.DII ...`（入口不可用时兜底）
- 依赖：Python >= 3.9，numpy / MDAnalysis / h5py / scipy / DuIvyTools；文本 tpr 解析需要系统里有 GROMACS 的 `gmx dump`
- 动手前先确认可用：`dii --help`；需要完整参数时让 agent 自己跑 `dii run --help` / `dii export --help`

## 快速开始示例

**输入/输出对示例：**

> **输入**：用户提供 `md.tpr` 和 `md.xtc`，要求"分析这个体系里有哪些氢键和盐桥"
> **输出**：
> 1. `dii run -t md.tpr -f md.xtc -o out/ --ff amber --interactions hydrogen_bond,salt_bridge`
> 2. 确认 `out/hydrogen_bond.h5`、`out/salt_bridge.h5` 生成
> 3. `dii export -i out/salt_bridge.h5 -o out_export/`
> 4. 向用户汇报 `out_export/salt_bridge_summary.csv` 中 occupancy 最高的前几对

## dii run —— 检测相互作用

```
dii run -t <topology.tpr> -f <trajectory.xtc> -o <out_dir> --ff amber [选项]
```

- 必填：`-t/--tpr`、`-f/--xtc`、`-o/--output`、`--ff`（当前仅支持 `amber`；已验证 amber03/94/96/99/99sb/99sb-ildn/GS/14sb + GAFF）
- 可选：`--interactions`（逗号分隔，默认 all）；`--strategy`（`two_pass` 默认 / `per_frame` / `per_tuple`）
- 8 种类型：`hydrogen_bond, pi_stacking, salt_bridge, hydrophobic, halogen_bond, metal_coordination, water_bridge, pi_cation`

要点：
- 每个类型输出独立 `<type>.h5` 到输出目录
- 单个类型失败不中断其余，打印 `[WARN] <type> detection failed: <原因>`
- 策略：`two_pass` 默认且最快（KDTree 优化）；长轨迹（数万帧以上）优先 `two_pass`，避免 `per_frame` 预分配大矩阵；`per_tuple` 是参考实现、最慢
- 先跑少量类型验证命令正确，再全量

## dii export —— 导出结果

```
dii export -i <result.h5> -o <export_dir>
```

- 一个 h5 可含多个 Interaction；同类型重复自动加 `_2` 序号
- 每个 Interaction 的产物：
  - `<type>_count.xvg`：每帧活跃相互作用数
  - `<type>_<metric>.xvg`：数值指标时间序列（字符串指标如 pistacking_type 跳过）
  - `<type>_existence.xpm`：存在性热力图
  - `<type>_summary.csv`：每对汇总（occupancy、均值/标准差；occupancy=0 或全 NaN 的对跳过）
  - π-stacking 额外 `<type>_type.xpm`（0=无、1=T 型、2=P 型）
- 终端打印概览：`Type / Pairs / Frames / Time range / Top 5 occupancies`

## 能力边界（避免误用）

1. **CLI 不支持按分子/残基过滤**相互作用来源（如"只看分子 A 与 B 之间"）。需要时用 Python API 的 `tuple_filter` 参数（Pipeline 检测器的 `detect(..., tuple_filter=...)`，接受 `(Tuple[Group, ...]) -> bool` 函数）
2. `--ff` 仅支持 amber；**联合原子力场**（无显式 H，如 GROMOS）不支持——供体鉴定依赖显式 H，会失效

## 常见错误与处理

| 现象 | 原因 | 处理 |
|---|---|---|
| `cannot read h5 file: <path>` | h5 损坏或路径不存在 | 检查路径；重新 `dii run` |
| `no interaction results in h5: <path>` | h5 为空（0 个 Interaction） | 确认 run 是否成功产出 |
| `unknown interaction type: '<name>'` | 类型名拼错 | 用 8 类之一 |
| `unknown force field: '<name>'` | `--ff` 仅支持 amber | 改用 `--ff amber` |
| `'all' cannot be mixed with other types` | `all` 与具体类型混用 | 要么只用 `all`，要么只列具体类型 |
| `[skip] <type>: nothing to export` | 该类型 0 对或 0 帧 | 正常提示，跳过导出 |
| `[WARN] <type> detection failed: ...` | 该类型检测异常 | 看错误原因；其余类型不受影响 |

## 工作流建议

1. 先用少量类型/小轨迹验证可跑通，再全量
2. run 完成后用 `dii export` 出可视化产物
3. 分析 `_summary.csv` 按 occupancy 排序，找高占位相互作用对
4. 需要看图时用 DuIvyTools：`dit xvg_show -f <file>.xvg`、`dit xpm_show -f <file>.xpm`
