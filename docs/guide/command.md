# dii 命令参考

`dii` 是 DuIvyInteractions 的命令行入口（`pyproject.toml` 的 `[project.scripts]` 注册，指向 `DuIvyInteractions.DII:main`），提供两个子命令：`run`（检测）与 `export`（导出）。

## 全局用法

```bash
dii [-h] {run,export} ...
```

`dii --help` 输出即上示用法；`dii run --help` 与 `dii export --help` 可查看各子命令的完整参数。

## 子命令：run —— 运行相互作用检测

检测指定轨迹中的相互作用，并将结果保存为 HDF5 文件（格式版本 1.0）：

```bash
dii run -t TPR -f XTC -o OUTPUT --ff FF [--interactions LIST] [--strategy STRATEGY]
```

### 参数

| 参数 | 必填 | 默认 | 说明 |
|:-----|:----|:----|:-----|
| `-t, --tpr` | ✅ | — | GROMACS 拓扑文件（基团鉴定的输入） |
| `-f, --xtc` | ✅ | — | 轨迹文件（MDAnalysis 可读格式，如 xtc） |
| `-o, --output` | ✅ | — | 输出目录（自动创建） |
| `--ff` | ✅ | — | 力场类型，取值 `amber` / `gromos` / `charmm` / `opls`（argparse choices 由 `IDENTIFIER_CLASSES` 注册表生成） |
| `--interactions` | 否 | `all` | 要检测的相互作用类型，逗号分隔（大小写不敏感），如 `hydrogen_bond,pi_stacking` |
| `--strategy` | 否 | `two_pass` | 检测策略：`two_pass` / `per_frame` / `per_tuple` |

`--ff` 四个取值对应的识别器（`DuIvyInteractions/group_identifiers/__init__.py`）：

| 取值 | 识别器 | 覆盖力场 |
|:-----|:-----|:-----|
| `amber` | `AmberFFGroupIdentifier` | Amber 03/94/96/99/99sb/99sb-ildn/GS/14sb 蛋白 + GAFF/GAFF2 配体 |
| `gromos` | `GromosFFGroupIdentifier` | GROMOS 53A6 / 54A7（联合原子力场，配 SPC/SPC-E 水） |
| `charmm` | `CharmmFFGroupIdentifier` | CHARMM36 / C36m（含 CGenFF 配体） |
| `opls` | `OplsFFGroupIdentifier` | OPLS-AA / L |

传入未知 `--ff` 会报错并列出全部可用值；传入未知 `--interactions` 类型会报错退出并列出全部可用类型（`ALL_INTERACTIONS`）。

### 支持的 8 种相互作用类型

`hydrogen_bond`（氢键）、`pi_stacking`（π-π 堆积）、`salt_bridge`（盐桥）、`hydrophobic`（疏水）、`halogen_bond`（卤键）、`metal_coordination`（金属配位）、`water_bridge`（水桥）、`pi_cation`（π-阳离子）。

- `--interactions` 值先统一转小写、去除首尾空白，再按逗号切分（`DII.py` 的 `main`）。
- `all` 是默认值，等价于"全部 8 类"；**`all` 不得与其他类型混用**，混用时报错 `'all' cannot be mixed with other types`。
- 单类型检测失败（如体系不含金属、无卤素）只打印 `[WARN] <类型> detection failed: <原因>`，不中断其余类型。

### 内部流程

1. 用 `GmxTprReader` 读取 tpr，识别器做一次基团鉴定（与帧无关，只做一次）；
2. 用 `mda.Universe(tpr, xtc)` 加载轨迹；水残基名取识别器的 `WATER_RESIDUES`（跨力场正确排除水）；
3. 对每个类型：按检测器的 `required_group_types` 过滤基团（除水桥外排除水分子），逐帧检测，结果保存为 `<output>/<类型>.h5`。

### 输出

每个类型生成一个 h5 文件：`<output>/<interaction_type>.h5`。例如：

```bash
dii run -t md.tpr -f md.xtc -o out/ --ff amber
# 生成 out/hydrogen_bond.h5, out/pi_stacking.h5, out/salt_bridge.h5,
#      out/hydrophobic.h5, out/halogen_bond.h5, out/metal_coordination.h5,
#      out/water_bridge.h5, out/pi_cation.h5
```

h5 格式版本为 1.0（`io/h5.py` 的 `FORMAT_VERSION`），gzip 压缩，详见{doc}`/reference/data_format`。

### 策略说明

| 策略 | 机制 | 适用场景 |
|:-----|:-----|:---------|
| `two_pass`（默认） | Pass1 逐帧精确发现活跃基团组（KDTree 预筛 + 精确判据，稀疏存储）→ Pass2 对发现的基团组补全全部帧的指标 | 大体系、长轨迹，性能最优（水桥 KDTree 预筛后 65 h → ~5 s） |
| `per_frame` | 对全部候选基团组逐帧向量化计算，预分配稠密矩阵 | 候选基数极大的类型（如水桥）；长轨迹时预分配矩阵内存占用高 |
| `per_tuple` | 逐候选基团组遍历全部帧并向量化 | 对照实现，处于"可能被舍弃"状态（测试不跑），候选多时极慢 |

## 子命令：export —— 导出结果

读取 h5 文件，导出为 xvg/xpm/csv，并在终端打印概览：

```bash
dii export -i INPUT -o OUTPUT
```

### 参数

| 参数 | 必填 | 说明 |
|:-----|:----|:-----|
| `-i, --input` | ✅ | h5 文件路径 |
| `-o, --output` | ✅ | 输出目录（自动创建） |

### 行为

- 支持包含多个 Interaction 的 h5（`dii run` 生成的文件每类型一个，符合此场景；`io/h5.py` 的 `load_interactions` 返回 Interaction 列表）。
- 若 h5 中同类型含多个 Interaction，文件名追加序号避免覆盖（如 `salt_bridge_2_count.xvg` 等）。
- 空数据（0 对或 0 帧）打印 `[skip] <类型>: nothing to export (pairs=0, frames=0)` 并跳过该 Interaction。
- h5 缺失或损坏报错：`cannot read h5 file: <路径> (file may be corrupted or missing)`（`SystemExit`）。
- h5 中无任何结果（空列表）报错：`no interaction results in h5: <路径>`（`SystemExit`）。
- h5 含未知类型报错：`unknown interaction type: '<类型>'. Available: ...`（`SystemExit`）。
- 导出的 xvg 仅包含数值指标（字符串指标如 π-堆积的 `pistacking_type` 不导出为 xvg，仅在 type.xpm 中体现）。

### 概览输出

每个被导出的 Interaction 先打印概览（`DII.py` 的 `_print_overview`）：

```
===== Salt Bridge overview =====
Type:     salt_bridge
Pairs:    47
Frames:   101
Time range: 0.0 ~ 1000.0 ps

Top 5 occupancies:
  1. ARG210(3443-3451)···ASP211(3461-3463)  100.0%
  ...
```

`Pairs` 为该类型基团对数（`Interaction.n_pairs`），`Frames` 为帧数，`Time range` 为首尾帧时间（ps，保留 1 位小数），`Top 5` 按占位率降序（`OVERVIEW_TOP_N = 5`）。

### 输出文件

对每个 Interaction，在输出目录生成：

| 文件 | 内容 |
|:-----|:-----|
| `<类型>_count.xvg` | 每帧活跃相互作用数量（X 轴：时间 ps） |
| `<类型>_<指标>.xvg` | 每个数值指标的时间序列，每列一条曲线（对应一个基团对） |
| `<类型>_existence.xpm` | 存在性热力图（行=pair，列=帧；白=No，蓝=Yes） |
| `<类型>_type.xpm` | 仅 π-堆积：堆积类型图（0=无 / 1=T 型 / 2=平行） |
| `<类型>_summary.csv` | 每对的汇总统计（`pair_label, occupancy, avg_<指标>, std_<指标>`） |

各类型的数值指标（`<类型>_<指标>.xvg` 的文件名由此决定，数据来自各导出器的 `metric_labels`）：

| 类型 | 指标（对应 xvg 文件名） |
|:-----|:-----|
| `hydrogen_bond` | `distance`（D-A 距离）、`angle`（D-H···A 角） |
| `salt_bridge` | `distance`（电荷中心距离） |
| `pi_stacking` | `distance`（环心距离）、`angle`（法向量夹角）、`offset`（偏移） |
| `hydrophobic` | `distance` |
| `halogen_bond` | `distance`（X···A 距离）、`don_angle`（C-X···A 角）、`acc_angle`（X···A-R 角） |
| `metal_coordination` | `distance`（金属-配位原子距离） |
| `water_bridge` | `dist_dw`（D-Ow）、`dist_wa`（Ow-A）、`theta`（O-H···Ow 角）、`omega`（H-Ow···A 角） |
| `pi_cation` | `distance`（环-电荷距离）、`offset`（偏移） |

文件名不以 Interaction 在 h5 中的序号命名，而以类型名命名——`dii run` 输出的原始 h5 文件名（`<类型>.h5`）与导出文件名（`<类型>_xxx`）一一对应。

## 退出码

- `0`：成功
- 非零：参数错误、文件不存在/损坏、未知类型、空结果等，均通过 `SystemExit` 抛出（如 `unknown interaction type: 'xxx'. Available: ...`、`'all' cannot be mixed with other types`、`cannot read h5 file: ...`）。

## 下一步

- 每种输出文件的具体格式与解读，见[结果解读](result)。
- 真实数据上的完整运行示例，见[快速上手](quickstart)。