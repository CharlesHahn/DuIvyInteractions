# 结果解读

`dii export` 将 h5 结果导出为 xvg / xpm / csv 三类文件，并在终端打印概览。本文档说明每类文件的格式与含义，以及如何解读。h5 本身的数据结构见{doc}`/reference/data_format`。

## 概览（打印到终端）

运行 `dii export` 时，对每个被导出的 Interaction 先打印概览（`DII.py` 的 `_print_overview`）：

```
===== Salt Bridge overview =====
Type:     salt_bridge
Pairs:    47
Frames:   101
Time range: 0.0 ~ 1000.0 ps

Top 5 occupancies:
  1. ARG210(3443-3451)···ASP211(3461-3463)  100.0%
  2. LYS70(1157-1160)···ASP180(2983-2985)  100.0%
  3. ARG291(4737-4745)···ASP295(4795-4797)  100.0%
  4. ARG244(4009-4017)···ASP211(3461-3463)  100.0%
  ...
```

- **Pairs**：该类型下基团对的数量（`Interaction.n_pairs`）。检测器只保留至少在一帧被判定存在的基团对，从未出现的候选对不会进入结果与导出文件。
- **Frames / Time range**：轨迹总帧数与首尾帧时间（ps，保留 1 位小数）。
- **Top N 占位率**：占位率最高的 N 对（`OVERVIEW_TOP_N = 5`）。占位率（occupancy）= 该对存在的帧数 / 总帧数（`Interaction.occupancy()`），即该相互作用在整个模拟中出现的时间比例，取值 0–1。

## 输出文件总览

对每个 Interaction 生成 5 类文件（π-堆积额外生成 type.xpm）：

| 文件 | 内容 | 格式 |
|:-----|:-----|:-----|
| `<类型>_count.xvg` | 每帧活跃的相互作用数量 | GROMACS xvg 时间序列 |
| `<类型>_<指标>.xvg` | 每个数值指标的时间序列（每列一个基团对） | GROMACS xvg 时间序列 |
| `<类型>_existence.xpm` | 存在性热力图（行=pair，列=帧） | GROMACS xpm 矩阵 |
| `<类型>_type.xpm` | 仅 π-堆积：堆积类型图（无 / T 型 / 平行） | GROMACS xpm 矩阵 |
| `<类型>_summary.csv` | 每对的汇总统计 | CSV |

文件名以类型名而非 h5 内序号命名，与 `dii run` 的 h5 文件名（`<类型>.h5`）一一对应。

## XVG（时间序列）

xvg 是 GROMACS 标准时间序列格式，可用 DuIvyTools 或 Xmgrace 绘图（首列为时间，单位 ps）。

**count.xvg**：`<类型>_count.xvg` 的 Y 轴为每帧活跃的相互作用数量（`existence` 逐帧求和），用于观察整体活性随时间的变化（如盐桥总数是否在某段模拟中减少）。

**metric.xvg**：`<类型>_<指标>.xvg` 的第一列是时间，后续每列是一条曲线，对应一个基团对（图例为 pair 标签），值为该指标在**全部帧**上的取值。注意：

- 指标数值在全部帧上都会写出（`two_pass` 对 Pass1 发现的基团对补全全帧指标；`per_frame` 对至少一帧存在的候选对计算全帧指标）；个别几何量在特定帧可能因几何未定义而为 NaN（绘图中表现为断点）。
- **xvg 本身不标注该对在哪帧"存在"**——帧级存在性以 `existence.xpm` 与 `summary.csv` 的 occupancy 为准。判定"某对在何时成立"请结合 existence 热力图，或直接看 CSV 的汇总。

## CSV 汇总（每对统计）

`<类型>_summary.csv` 以行-列表格记录每对基团的汇总统计，**仅包含占位率 > 0 的基团对**（`to_csv_summary` 跳过 occupancy = 0 的行）：

| 列 | 含义 |
|:---|:-----|
| `pair_label` | 基团对标签（格式见下文） |
| `occupancy` | 占位率（0–1，保留 4 位小数） |
| `avg_<metric>` | 该指标在**活跃帧**（该对存在的帧）上的均值 |
| `std_<metric>` | 该指标在活跃帧上的标准差 |

统计只基于活跃帧（代码取 `metrics[i][existence[i]]` 并用 nanmean/nanstd）；某对活跃帧的指标全部为 NaN 时，对应单元格留空。字符串指标（如 π-堆积的 `pistacking_type`）不进入 CSV。

示例（数值为演示）：

```
pair_label,occupancy,avg_distance,std_distance
ARG210(3443-3451)···ASP211(3461-3463),1.0000,3.6200,0.0548
```

解读：`ARG210···ASP211` 盐桥占位率 100%，平均距离 3.62 Å，标准差 0.055 Å——该盐桥全程稳定存在且构象紧凑。反之若某对 occupancy 低而 std 大，说明该相互作用时断时续或几何摆动明显。

## XPM（热力图）

xpm 是 GROMACS 标准矩阵热力图格式，可直观查看"哪些基团对在哪些帧存在/处于何种构型"。

### existence 热力图

`<类型>_existence.xpm`：

- **行**：基团对（pair），**列**：帧（X 轴时间）
- **颜色**：白 = 该帧无此相互作用（No），蓝（#38A7D0）= 存在（Yes）
- **用途**：快速查看各基团对在时间上的存在模式——连续稳定（整行蓝）、间歇出现（蓝白交替）、只在局部时间出现。

### π-堆积类型热力图

`<类型>_type.xpm`（仅 π-堆积）：

- 三值颜色：白 = 无堆积（None）、粉（#F67088）= T 型堆积、蓝（#38A7D0）= 平行堆积
- **用途**：区分 π 堆积的几何构型（T 型/平行）随时间的演化；存在性以 existence.xpm 为准，type.xpm 只在存在帧上着色。

## 各类型的指标

`<类型>_<指标>.xvg` 与 `summary.csv` 的指标列由各类型定义（单位见各导出器 `metric_labels`）：

| 类型 | 指标（含义 / 单位） |
|:-----|:-----|
| `hydrogen_bond` | `distance`：D-A 距离（Å）；`angle`：D-H···A 角（°） |
| `salt_bridge` | `distance`：电荷中心距离（Å） |
| `pi_stacking` | `distance`：环心距离（Å）；`angle`：两环法向量夹角（°）；`offset`：环心间偏移（Å） |
| `hydrophobic` | `distance`：原子间距（Å） |
| `halogen_bond` | `distance`：X···A 距离（Å）；`don_angle`：C-X···A 角（°）；`acc_angle`：X···A-R 角（°） |
| `metal_coordination` | `distance`：金属-配位原子距离（Å） |
| `water_bridge` | `dist_dw`：D-Ow 距离（Å）；`dist_wa`：Ow-A 距离（Å）；`theta`：O-H···Ow 角（°）；`omega`：H-Ow···A 角（°） |
| `pi_cation` | `distance`：环心-正电中心距离（Å）；`offset`：偏移（Å） |

## 基团对标签格式

pair 标签用于在结构上定位相互作用的原子（各导出器的 `get_pair_label` 实现）。原子号是 tpr 中的**全局原子索引**：

| 类型 | 标签格式 | 示例 |
|:-----|:---------|:-----|
| 盐桥 / π-阳离子 / π-堆积 | `残基名残基号(起原子-止原子)···残基名残基号(起原子-止原子)`（括弧内为基团所含原子的最小-最大全局索引） | `ARG210(3443-3451)···ASP211(3461-3463)` |
| 氢键 | `残基名:供体原子(索引)-氢(索引)···残基名:受体原子(索引)` | `ARG210:NE(3443)-HE(3444)···ASP211:OD1(3461)` |
| 水桥 | `供体残基:供体原子(索引)-氢(索引)···水残基:OW(索引)···受体残基:原子(索引)` | `ARG210:NE(3443)-HE(3444)···SOL:OW(5000)···ASP211:OD1(3461)` |
| 卤键 | `供体残基:碳(索引)-卤素(索引)···受体残基:原子(索引)`（供体展示 C-X；受体侧 R 逐帧动态选择，不体现在标签中） | `D927:C(9)-CL(10)···ASP211:OD1(3461)` |
| 疏水 | `残基名:原子名(索引)···残基名:原子名(索引)`（取基团首原子） | `D927:C7(100)···PHE232:CG(3700)` |
| 金属配位 | `金属残基:金属原子(索引)···配位残基:配位原子(索引)` | `MG:MG(100)···ASP211:OD1(3461)` |

## 解读建议

- **找稳定相互作用**：按 `summary.csv` 的 occupancy 降序排列，高占位率（如 > 0.8）且 std 小的对是全程稳定的候选热点；低占位率对需结合 existence.xpm 判断是间歇出现还是模拟后期才形成。
- **看时间演化**：用 existence.xpm 观察对在时间轴上的连续/间歇模式；用 count.xvg 看整体活性变化；用 type.xpm（π-堆积）看构型切换。
- **核对几何合理性**：将 `<指标>.xvg` 的数值与判据阈值对照（如氢键 D-A ≤ 4.1 Å、D-H···A ≥ 100°——各类型判据见{doc}`/reference/criteria`），确认检出对确实落在判据内。

## 下一步

- 理解 8 类相互作用的检测判据与几何量定义，见{doc}`/reference/criteria`。
- 理解结果数据格式（h5），见{doc}`/reference/data_format`。