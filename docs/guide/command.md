# dii 命令参考

`dii` 是 DuIvyInteractions 的命令行入口，提供两个子命令：`run`（检测）与 `export`（导出）。

## 全局用法

```bash
dii [-h] {run,export} ...
```

## 子命令：run —— 运行相互作用检测

检测指定轨迹中的相互作用，并将结果保存为 h5 文件。

```bash
dii run -t TPR -f XTC -o OUTPUT --ff FF [--interactions LIST] [--strategy STRATEGY]
```

### 参数

| 参数 | 必填 | 默认 | 说明 |
|:-----|:----|:----|:-----|
| `-t, --tpr` | ✅ | — | GROMACS 拓扑文件 |
| `-f, --xtc` | ✅ | — | 轨迹文件 |
| `-o, --output` | ✅ | — | 输出目录（自动创建） |
| `--ff` | ✅ | — | 力场类型，当前仅支持 `amber` |
| `--interactions` | 否 | `all` | 要检测的相互作用类型，逗号分隔（如 `hydrogen_bond,pi_stacking`） |
| `--strategy` | 否 | `two_pass` | 检测策略：`two_pass` / `per_frame` / `per_tuple` |

### 支持的 8 种相互作用类型

`hydrogen_bond`（氢键）、`pi_stacking`（π-π 堆积）、`salt_bridge`（盐桥）、`hydrophobic`（疏水）、`halogen_bond`（卤键）、`metal_coordination`（金属配位）、`water_bridge`（水桥）、`pi_cation`（π-阳离子）

### 输出

每个类型生成一个 h5 文件：`<output>/<interaction_type>.h5`。例如：

```bash
dii run -t md.tpr -f md.xtc -o out/ --ff amber
# 生成 out/hydrogen_bond.h5, out/pi_stacking.h5, out/salt_bridge.h5, ...
```

### 策略说明

| 策略 | 机制 | 适用场景 |
|:-----|:-----|:---------|
| `two_pass`（默认） | Pass1 逐帧发现活跃对（稀疏）→ Pass2 补全全帧 | 大体系，性能最优（水桥用 KDTree 预筛后 65h→~5s） |
| `per_frame` | 逐帧处理全部候选对，向量化 | 候选对数量极大的类型（如水桥） |
| `per_tuple` | 逐候选对加载全部帧 | 对照实现，参考用 |

## 子命令：export —— 导出结果

读取 h5 文件，导出为 xvg/xpm/csv，并打印概览。

```bash
dii export -i INPUT -o OUTPUT
```

### 参数

| 参数 | 必填 | 说明 |
|:-----|:----|:-----|
| `-i, --input` | ✅ | h5 文件路径 |
| `-o, --output` | ✅ | 输出目录（自动创建） |

### 行为

- 支持包含多个 Interaction 的 h5（`dii run` 生成的文件每种类型一个，符合此场景）
- 若 h5 含同类型多个 Interaction，文件名追加序号避免覆盖（如 `salt_bridge_2_*`）
- 空数据（0 对或 0 帧）自动跳过并提示
- 损坏或缺失的 h5 给出友好错误

### 输出文件

对每个 Interaction，在输出目录生成：

| 文件 | 内容 |
|:-----|:-----|
| `<类型>_count.xvg` | 每帧活跃相互作用数量 |
| `<类型>_<metric>.xvg` | 每个数值指标的时间序列（如 distance、angle） |
| `<类型>_existence.xpm` | 存在性热力图（行=pair，列=帧） |
| `<类型>_type.xpm` | π-堆积专属：堆积类型图（0=无/1=T 型/2=P 型） |
| `<类型>_summary.csv` | 每对的汇总统计 |

## 退出码

- `0`：成功
- 非零：参数错误、文件不存在/损坏、未知类型等（通过 `SystemExit` 抛出）