# 快速上手

本指南用仓库自带的真实 GROMACS 测试数据跑通 `dii run → dii export` 完整流程，并给出 4 个力场的可运行命令。所有示例命令均需在**项目仓库根目录**（含 `Tests/`、`pyproject.toml` 的目录）下执行；示例输出对齐当前代码（默认 `two_pass` 策略、HDF5 v1.0 存储）。

## 1. 准备数据

仓库 `Tests/` 下提供 3 套真实 GROMACS 测试数据（与 `Tests/unittests/` 集成测试共用）：

| 目录 | 力场 | 体系 | 文件 |
|:-----|:-----|:-----|:-----|
| `test_MD_case_amber/` | Amber 家族（amber14sb 蛋白 + GAFF 配体） | KRAS-RBD D927 复合物（含 GNP/Mg²⁺ 金属中心），116,383 原子 | `md.tpr` + `md1ns.xtc`（1 ns，101 帧，0–1000 ps） |
| `test_MD_case_gromos/` | GROMOS 53A6 | 130 残基蛋白 + 6 配体（ZIN1–6） | `gromos53a6_md.tpr` + `gromos53a6_md10ns.xtc` |
| `test_MD_case_charmm36/` | CHARMM36 | SMO-BST（Smoothened–β-谷甾醇）复合物 | `SMO-BST_md_complex.tpr` + `SMO-BST_md_complex.xtc` |

> **OPLS-AA**：基团识别已由 `Tests/unittests/test_opls_identifier.py` 验证（带电残基、芳香环等），但仓库尚未配套真实轨迹测试数据；`dii run --ff opls` 的用法与其余力场完全一致。
>
> CHARMM36 数据还提供全原子版 `SMO-BST_md_fullatom.tpr` + `SMO-BST_md_fullatom_100ps.xtc`（用于水桥验证，见 `Tests/unittests/test_charmm_water_bridge_real.py`）。

下面以 Amber 测试数据为例演示完整流程。

## 2. 运行相互作用检测

用 `dii run` 检测全部 8 类相互作用并保存为 h5 文件：

```bash
dii run -t Tests/test_MD_case_amber/md.tpr -f Tests/test_MD_case_amber/md1ns.xtc -o out_amber/ --ff amber
```

参数说明：

- `-t` / `--tpr`：GROMACS 拓扑文件（必填）
- `-f` / `--xtc`：轨迹文件（必填，MDAnalysis 可读格式）
- `-o` / `--output`：输出目录（必填，自动创建）
- `--ff`：力场类型，取值 `amber` / `gromos` / `charmm` / `opls`（必填，与 `DuIvyInteractions/group_identifiers/__init__.py` 的 `IDENTIFIER_CLASSES` 注册表键一一对应）
- 默认检测全部 8 类，每类存为一个 h5 文件：`<输出目录>/<类型>.h5`

检测成功后输出目录内容：

```
$ ls out_amber/
hydrogen_bond.h5          metal_coordination.h5     salt_bridge.h5
hydrophobic.h5            pi_cation.h5              water_bridge.h5
halogen_bond.h5           pi_stacking.h5
```

> `dii run` 成功时不打印中间过程；若某一类型检测失败，打印 `[WARN] <类型> detection failed: <原因>` 但不中断其余类型（`Pipeline.run` 逐类型 try/except，见 `DuIvyInteractions/DII.py`）。
>
> 在该测试集（101 帧）上，`two_pass` 策略检出盐桥 **47 对**（该数值即 `Tests/unittests/test_saltbridge_two_pass.py` 的断言基线，且全帧 distance 无缺失值）、π-堆积 **10 对**（38 个芳香环两两组合中实际成立的堆积对；已在当前代码上复核）。

### 其他力场

GROMOS 53A6：

```bash
dii run -t Tests/test_MD_case_gromos/gromos53a6_md.tpr \
        -f Tests/test_MD_case_gromos/gromos53a6_md10ns.xtc \
        -o out_gromos/ --ff gromos
```

CHARMM36：

```bash
dii run -t Tests/test_MD_case_charmm36/SMO-BST_md_complex.tpr \
        -f Tests/test_MD_case_charmm36/SMO-BST_md_complex.xtc \
        -o out_charmm/ --ff charmm
```

> 水的排除（水桥/金属配位检测不把水当作蛋白供体/配位点）由各识别器的 `WATER_RESIDUES` 类属性统一决定（`DuIvyInteractions/core/interfaces.py`）：Amber `{SOL, HOH, WAT}`、GROMOS `{SOL}`、CHARMM `{TIP3, HOH, SOL, WAT}`、OPLS `{HOH, HO4, HO5, SOL, WAT}`。用户无需手动指定，只需保证 `--ff` 与体系实际力场一致。

### 用 Python API 等价运行

```python
from DuIvyInteractions.pipeline import Pipeline

Pipeline(ff="amber", strategy="two_pass").run(
    "Tests/test_MD_case_amber/md.tpr",
    "Tests/test_MD_case_amber/md1ns.xtc",
    "out_amber/",
    interactions=None,          # None = 全部 8 类；也可传类型名列表
)
```

## 3. 导出结果

用 `dii export` 将 h5 结果导出为 xvg/xpm/csv，并在终端打印概览：

```bash
dii export -i out_amber/salt_bridge.h5 -o out_export/
```

终端输出概览（示例；具体残基与占位率取决于体系与版本）：

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

> 概览含义：`Pairs` 为至少在一帧被判定存在的基团对数（本示例 47，与单元测试基线一致）；`Frames` 为轨迹总帧数（101）；`Time range` 为首尾帧时间；`Top 5 occupancies` 按占位率降序排列，占位率 = 该对存在的帧数 / 总帧数。

输出目录中生成的文件：

```
$ ls out_export/
salt_bridge_count.xvg        # 每帧活跃盐桥数量
salt_bridge_distance.xvg     # 每对的距离时间序列（电荷中心距离，Å）
salt_bridge_existence.xpm    # 存在性热力图（行=pair，列=帧）
salt_bridge_summary.csv      # 每对的汇总统计（pair_label, occupancy, avg/std）
```

> 若导出 π-堆积（`pi_stacking.h5`），会额外生成 `<类型>_type.xpm`（堆积类型图：无 / T 型 / 平行）；若 h5 数据为空（0 对或 0 帧），打印 `[skip] <类型>: nothing to export` 并跳过导出。

## 4. 常用参数

只检测部分类型（逗号分隔；`all` 不得与其他类型混用，混用会报错退出）：

```bash
dii run -t md.tpr -f md.xtc -o out/ --ff amber \
    --interactions hydrogen_bond,pi_stacking
```

切换检测策略（默认 `two_pass`，大体系推荐；`per_frame` 为逐帧向量化实现；`per_tuple` 为逐候选组遍历的对照实现，仅供参考）：

```bash
dii run -t md.tpr -f md.xtc -o out/ --ff amber --strategy two_pass
```

## 5. 性能参考

以下为 Amber 测试体系（KRAS-RBD，116,383 原子，101 帧）上水桥检测的参考耗时（数字来自项目设计文档与 `Tests/unittests/`，环境为仓库开发环境；实际耗时随体系规模与机器而异）：

| 检测策略 | 水桥耗时 | 说明 |
|:---------|:---------|:-----|
| `two_pass`（默认） | ~5 s | Pass1 逐帧 KDTree（预筛半径 8.2 Å）发现活跃三元组 + Pass2 补全全帧指标；稀疏存储 |
| `per_frame` | ~5 s | 第一帧 KDTree 预筛候选三元组 + 逐帧向量化计算（约 56 ms/帧）；为全部候选对预分配矩阵 |
| `per_tuple` | 约 65 h | 逐候选三元组遍历全部帧；候选多时极慢，仅作对照参考 |

水桥候选三元组可达 24.9 万个，`per_tuple` 逐对遍历轨迹是其慢的根因。**大体系/长轨迹/候选组合多时请使用默认 `two_pass`**。

## 6. 已知限制与后续阅读

- 完整已知限制清单（PBC 未处理、水桥 TwoPass 缺 Ow-A 距离下界、疏水-芳香未去重、PerFrame 长轨迹内存占用等）见{doc}`/reference/limitations`。
- 8 类相互作用的检测判据与结果解读，见[结果解读](result)。
- `dii` 全部命令与参数，见[命令参考](command)。
- 为什么直接读 tpr 原子类型（核心原理），见{doc}`/reference/concepts`。
