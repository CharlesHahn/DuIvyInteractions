# DuIvyInteraction

基于 MD 拓扑力场参数的分子间相互作用判定工具。

## 解决什么问题

现有工具（PLIP / ProLIF）分析 MD 轨迹时，通过 OpenBabel / RDKit 从坐标重建化学信息（键序、芳香性、加氢），丢弃了 MD 力场拓扑中原有的化学语义。这导致：
- 对 trjconv 导出的 PDB（无 CONECT / 键序）芳香性判定失败
- 每帧重复推断，效率低
- 推断结果与力场参数不自洽

## 核心思路

**直接从 GROMACS tpr 拓扑中读取力场原子类型**，确定性识别化学基团，与模拟力场完全自洽。

力场原子类型（如 GAFF 的 `ca` = 芳香碳、`na` = 吡咯氮）是参数化时由 antechamber / sobtop 做出的化学判决的留存记录。直接读取 = 零损失、零歧义，不需要从坐标反推。

## 工具优势

- **确定性**：基团鉴定基于力场原子类型，不依赖几何推断
- **与力场自洽**：结果与模拟使用的力场参数同源
- **全原子显式 H**：H 键供体（D–H 键）、水桥（SOL 残基）、金属（元素+电荷）全部零推断
- **Amber 全家族兼容**：已验证 amber03/94/96/99/99sb/99sb-ildn/GS/14sb + GAFF，类型映射零冲突
- **支持 8 种相互作用类型**：氢键、π-π 堆积、盐桥、π-阳离子、卤键、疏水、金属配位、水桥

## 依赖

- Python >= 3.9
- NumPy >= 1.20
- SciPy >= 1.7
- MDAnalysis >= 2.0
- h5py >= 3.0
- DuIvyTools >= 0.6.0
- GROMACS（`gmx dump`，用于文本格式 tpr 解析）

## 安装

```bash
pip install -e .
```

## 项目状态

基团鉴定、相互作用检测、HDF5 结果存储、xvg/xpm/CSV 导出、Pipeline 编排和 DII 命令行工具均已完成。详见 `doc/` 目录下的设计文档。

## 使用

```bash
# 安装
pip install -e .

# 运行相互作用检测并保存 h5（--ff 必选，当前支持 amber）
dii run -t md.tpr -f md.xtc -o out/ --ff amber
# 可选参数：
#   --interactions hydrogen_bond,pi_stacking   只检测部分类型（默认 all=8类）
#   --strategy two_pass|per_frame|per_tuple    检测策略（默认 two_pass）

# 导出 h5 结果为 xvg/xpm/csv 并打印概览（支持多 Interaction h5）
dii export -i out/salt_bridge.h5 -o out_export/
```

支持的 8 种相互作用类型：氢键、π-π 堆积、盐桥、π-阳离子、卤键、疏水、金属配位、水桥。

在 Python 中调用：

```python
from DuIvyInteractions.pipeline import Pipeline

# 配置：力场 + 策略
pipeline = Pipeline(ff="amber", strategy="two_pass")
pipeline.run("md.tpr", "md.xtc", "out/", interactions=None)  # None=全部8类
```
