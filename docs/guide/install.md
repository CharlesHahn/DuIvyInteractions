# 安装与环境配置

DuIvyInteraction 是一个基于 MD 拓扑力场参数的分子间相互作用分析工具。本文档说明如何安装并验证环境。

## 系统要求

- **Python** >= 3.9（推荐 3.10+）
- **操作系统**：Linux / macOS / Windows（建议 Linux 或 macOS，完整支持 GROMACS tpr 解析）

## 依赖

安装 DuIvyInteractions 会自动安装以下 Python 依赖：

| 依赖 | 版本要求 | 用途 |
|:-----|:---------|:-----|
| numpy | >= 1.20 | 数值计算 |
| MDAnalysis | >= 2.0 | tpr 二进制读取、轨迹加载 |
| h5py | >= 3.0 | 结果 HDF5 序列化 |
| scipy | >= 1.7 | 空间索引（KDTree，TwoPass 检测器预筛） |
| DuIvyTools | >= 0.6.0 | xvg/xpm 结果文件解析与可视化 |

此外，如需用文本格式解析 tpr（`gmx_tpr_dump_reader`），需安装 **GROMACS** 并提供 `gmx dump` 命令。日常使用推荐直接用 `gmx_tpr_reader`（MDAnalysis 读二进制 tpr），无需 GROMACS。

## 安装

在项目根目录执行：

```bash
pip install -e .
```

`-e` 为可编辑安装，适合开发调试；发布环境可省略 `-e`。

## 验证安装

检查命令行工具是否可用：

```bash
dii --help
```

应输出 `usage: dii [-h] {run,export} ...`，表示 `dii` 命令已安装。

在 Python 中验证模块可导入：

```python
from DuIvyInteractions.pipeline import Pipeline
from DuIvyInteractions.io import load_interactions
print("DuIvyInteractions OK")
```

## 下一步

安装完成后，进入[快速上手](quickstart)用真实数据跑一遍完整流程。