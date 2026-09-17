# 快速上手

本指南用仓库自带的测试数据（`Tests/test_MD_case/`）跑一遍完整流程，展示 `dii` 命令的真实用法与输出。数据为一个 KRAS-RBD 体系的 1 ns 轨迹（101 帧）。

## 1. 准备数据

测试数据在仓库 `Tests/test_MD_case/` 下，包含：

- `md.tpr` — GROMACS 拓扑文件（Amber 系力场：amber14sb 蛋白 + GAFF 配体）
- `md1ns.xtc` — 1 ns 轨迹（101 帧）

## 2. 运行相互作用检测

用 `dii run` 检测全部 8 类相互作用并保存为 h5 文件：

```bash
dii run -t Tests/test_MD_case/md.tpr -f Tests/test_MD_case/md1ns.xtc -o out/ --ff amber
```

- `-t` / `-f`：拓扑与轨迹文件
- `-o`：输出目录（自动创建）
- `--ff`：力场类型，当前仅支持 `amber`
- 默认检测全部 8 类，存为 `out/<类型>.h5`

**真实输出**（截取自测试数据，盐桥与 π-堆积）：

```
=== salt_bridge ===
基团对数: 47, 帧数: 101, 时间: 0~1000 ps
metrics: ['distance']
saved out/salt_bridge.h5

=== pi_stacking ===
基团对数: 10, 帧数: 101, 时间: 0~1000 ps
metrics: ['distance', 'angle', 'offset', 'pistacking_type']
saved out/pi_stacking.h5
```

## 3. 导出结果

用 `dii export` 将 h5 结果导出为 xvg/xpm/csv，并打印概览：

```bash
dii export -i out/salt_bridge.h5 -o out_export/
```

**真实输出**（概览 + 生成的文件）：

```
===== Salt Bridge 概览 =====
类型:     salt_bridge
基团对数: 47
帧数:     101
时间范围: 0.0 ~ 1000.0 ps

Top 5 占位率:
  1. ARG210(3443-3451)···ASP211(3461-3463)  100.0%
  2. LYS70(1157-1160)···ASP180(2983-2985)  100.0%
  3. ARG291(4737-4745)···ASP295(4795-4797)  100.0%
  4. ARG244(4009-4017)···ASP211(3461-3463)  100.0%
  5. ARG5(82-90)···PRO141(2352-2354)        100.0%

输出文件：
salt_bridge_count.xvg         # 每帧活跃相互作用数量
salt_bridge_distance.xvg      # 每对的距离时间序列
salt_bridge_existence.xpm     # 存在性热力图（行=pair，列=帧）
salt_bridge_summary.csv       # 每对的汇总统计
```

> 上表 Top 1 显示 `ARG210···ASP211` 占位率 100%，说明该盐桥在全部 101 帧都存在——即 RBD 与配体/受体间稳定形成的盐桥接触。

## 4. 更多参数

只检测部分类型：

```bash
dii run -t md.tpr -f md.xtc -o out/ --ff amber \
    --interactions hydrogen_bond,pi_stacking
```

切换检测策略（默认 `two_pass`，大体系推荐；`per_tuple`/`per_frame` 为对照实现）：

```bash
dii run -t md.tpr -f md.xtc -o out/ --ff amber --strategy two_pass
```

## 下一步

- 了解 8 类相互作用的检测判据与结果解读，见[结果解读](result.html)
- 查看 `dii` 全部命令与参数，见[命令参考](command.html)
- 理解核心原理（为什么直接读 tpr 类型），见[核心概念](concepts.html)（位于参考手册）