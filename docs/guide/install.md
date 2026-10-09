# 安装与环境配置

本文档说明如何安装 DuIvyInteractions 并验证环境。包发布版本为 v0.0.1（`pyproject.toml`），结果文件采用 HDF5 格式版本 1.0。

## 系统要求

- **Python** >= 3.9（`pyproject.toml` 声明 `requires-python = ">=3.9"`，分类器覆盖 3.9–3.12，推荐 3.10+）
- **操作系统**：Linux / macOS / Windows（建议 Linux 或 macOS；tpr/xtc 解析为纯 Python + MDAnalysis 实现，与操作系统无关）
- **输入数据**：一个 GROMACS 拓扑文件（tpr）与一个轨迹文件（xtc，或其他 MDAnalysis 可读格式）。`dii run` 同时需要两者：tpr 提供力场原子类型与键合图（基团鉴定的唯一依据），轨迹提供逐帧坐标（几何判定）

> **无需安装 GROMACS**：默认读取路径 `GmxTprReader` 直接用 MDAnalysis 解析 tpr 二进制文件，不调用 `gmx`。仅在需要使用 `GmxTprDumpReader`（解析 `gmx dump` 文本的独立 Reader 实现）时才需要 GROMACS，且当前 `Pipeline` 固定使用二进制直读路径，不会自动回退。

## 依赖

安装 DuIvyInteractions 时会自动安装以下 Python 依赖（与 `pyproject.toml` 的 `dependencies` 完全一致）：

| 依赖 | 版本要求 | 用途 |
|:-----|:---------|:-----|
| numpy | >= 1.20 | 数值计算（矩阵式检测、统计） |
| MDAnalysis | >= 2.0 | tpr 二进制解析（`gmx_tpr_reader`）、轨迹加载（`Pipeline.run` 中的 `mda.Universe(tpr, xtc)`） |
| h5py | >= 3.0 | 结果 HDF5 序列化（格式版本 1.0，gzip 压缩） |
| scipy | >= 1.7 | 空间索引 KDTree（氢键/水桥 TwoPass 与 PerFrame 检测器的预筛） |
| DuIvyTools | >= 0.6.0 | xvg/xpm 结果文件的解析、构建与后续绘图 |

可选依赖：

- **开发/测试**：`pip install -e ".[dev]"` 额外安装 pytest（`[project.optional-dependencies].dev`）。

### MDAnalysis 版本与 tpr 兼容性

tpr 二进制解析依赖 MDAnalysis 对对应 GROMACS 版本（tpx 版本号）的支持。较新 GROMACS 生成的 tpr 需要较新版本的 MDAnalysis；若解析 tpr 时报 `Your tpx version is XXX, which this parser does not support`，请升级 MDAnalysis（本项目开发环境使用 2.7–2.10，实测 2.10 可解析 GROMACS 2018.1 的 tpr）。

## 安装

### 方式一：从 PyPI 安装（推荐）

```bash
pip install duivyinteractions
```

包名为 `duivyinteractions`（PyPI 已发布 v0.0.1）。安装时自动安装全部依赖，并在 PATH 注册 `dii` 命令（入口：`pyproject.toml` 的 `[project.scripts]`，指向 `DuIvyInteractions.DII:main`）。

### 方式二：从源码安装（开发）

在项目根目录（包含 `pyproject.toml` 的目录）执行：

```bash
pip install -e .
```

`-e` 为可编辑安装，源码改动即时生效，适合开发调试；仅需使用命令行工具时可用 `pip install .`。

## 验证安装

检查命令行工具是否可用：

```bash
dii --help
```

应输出（`run` 与 `export` 两个子命令已注册）：

```
usage: dii [-h] {run,export} ...
```

查看子命令参数（注意 `--ff` 的可选值与 `--strategy` 的可选值）：

```bash
dii run --help
dii export --help
```

在 Python 中验证模块可导入：

```python
from DuIvyInteractions.pipeline import Pipeline
from DuIvyInteractions.io import load_interactions

print("DuIvyInteractions OK")
```

对已安装的包查看版本与依赖：

```bash
pip show duivyinteractions
```

## 下一步

安装完成后，进入[快速上手](quickstart)，用仓库自带的真实测试数据（Amber / GROMOS / CHARMM36）跑一遍完整流程。
