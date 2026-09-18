# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## 编程哲学（项目规范）

1. **KISS**：函数 ≤ 20 行，类 ≤ 100 行；不做"可能有用"的功能；简单方案优先。
2. **显式优于隐式**：数据流清晰可见；无魔法数字（全部常量）；类型注解必加。
3. **可读性即文档**：docstring 只写"是什么"，不写"为什么"（"为什么"属 doc/ 设计文档）。

## 项目：基于 MD 拓扑力场参数的相互作用判定工具

**输入** GROMACS tpr + xtc；**输出** PLIP 级精细分类的相互作用（π-π 堆积、H 键、盐桥、疏水、π-阳离子、卤键、水桥、金属配位）。

**核心创新**：直接读 tpr 力场原子类型做确定性基团鉴定，不用 OpenBabel（PLIP）/RDKit（ProLIF）重建化学。

### 两段式架构

```
① 基团鉴定（难点）：tpr 原子类型 + 键合图 + 显式 H + 电荷 → 确定性识别基团，只做一次、与帧无关
② 几何判定：逐帧 PLIP 式距离/角度/平面临近判据 → 每帧相互作用列表 → 时间统计
```

### 动机（详见 doc/project_background.md）

PLIP/ProLIF 通过 OpenBabel/RDKit 重建化学（键序/芳香/加氢），丢弃 MD 拓扑化学信息，对 trjconv 导出 PDB 判芳香失败。本项目差异化：**直接读 tpr 类型语义**——这是现有工具无一做到的空白。

## 类型映射关键结论

类型→化学特征的映射来自两条独立证据链（GAFF 命名规则 + rtp 残基反推），已验证 Amber 全家族零冲突。详细论证见 `doc/project_background.md`。

⚠️ **类型名坑**（务必记住）：
- `C` = 主链羰基碳 AND Tyr CZ 芳环碳（双身份，需"环内≥4强芳香邻居"升级）
- `N3`（蛋白正电氨基）vs `n3`（GAFF 中性氨基）——大小写不同义
- `CA` 跨力场不同义（amber=芳香碳，GROMOS=α碳）——必须特征映射表
- N-H 键在 `Constraint:` 段，不在 `Bond:` 段！供体识别必须合并两者

## 基团鉴定可行性矩阵（检测器开发核心参考）

| 相互作用 | 所需基团特征 | tpr 中证据 | 确定性 |
|:----|:----|:----|:----|
| π-π / π-阳 / 卤-π | 芳香环 | 类型 `ca/na/nb/cp/cg` → 图论环检测（SSSR）→ 环内全芳香；可选键序 func=4/5 与平面性交叉验证 | 高 |
| H 键供体 | D–H（D=N/O/S/F） | D–H 键条目显式存在（全原子）＋ q(H)>0 | 100%（零推断） |
| H 键受体 | 孤对可用 | 类型规则表（`o/oh/os/nb/n/f`…）＋ H 计数（`nb` 无 H 必受体）＋ 电荷负验证 | 高 |
| 盐桥 | 形式电荷对 | 蛋白残基名字典（LYS/ARG/ASP/GLU/HIP）；配体靠类型＋H 计数 | 蛋白 100%，配体高 |
| 疏水 | 非极性 C/S/X | 类型集合（`c3/c2`…）＋无极性取代 | 高 |
| 卤键 | σ-hole 卤素 | 卤素类型（`f/cl/br/i`）＋邻接碳区分芳香/烷基卤 | 高 |
| 金属配位 | 金属中心 | 元素/类型＋电荷（本项目 Mg²⁺） | 100% |
| 水桥 | 水分子 | 残基名 SOL/HOH ＋ OW/HW 原子名 | 100% |

## 两个工程难点与对策

1. **力场类型语义映射**：类型名跨力场不同义（GAFF `ca` vs CHARMM `CG2R61` vs OPLS `CA`）。对策：建**特征空间映射**（类型→{杂化, 芳香性, 极性, 带H, 孤对}），基团由特征组合；新力场入库 = 填特征表。首版支持 Amber 系。**边界**：全原子力场（GROMOS 联合原子无显式 H → 供体鉴定失效，声明不支持）。
2. **tpr 二进制读取**：MDAnalysis 垫底 → 自研解析 → `gmx dump` 文本兜底（验收基准）。

## 目录结构

```
DuIvyInteraction/
├── DuIvyInteractions/            # 主包
│   ├── core/                     # datas.py(数据类), interfaces.py(ABC), constants.py
│   ├── system_readers/           # gmx_tpr_reader.py, gmx_tpr_dump_reader.py
│   ├── group_identifiers/        # amber_ff_identifier.py
│   ├── interaction_detectors/    # 8 类型 × 3 策略（*_detector_{per_tuple,per_frame,two_pass}.py）
│   ├── io/                       # h5.py(序列化), interaction_exporter.py + 8 导出子类
│   ├── pipeline.py               # 编排：Reader→Identifier→Detector→h5
│   ├── DII.py                    # 命令行入口（dii run / dii export）
│   ├── utils/                    # （空，待实现）
│   └── visualizers/              # （空，待实现）
├── Tests/                        # 单元测试
│   ├── unittests/                # 单元测试用例
│   ├── test_MD_case/             # 测试数据（KRAS-RBD 体系，已 gitignore）
│   └── original_draft/           # 早期验证脚本与调研（历史参考，不随包发布）
├── doc/                          # 中文设计文档 + project_background.md(论证归档)
├── docs/  docs_en/               # 文档站点（中/英，ReadTheDocs）
└── dist/                         # PyPI 构建产物（发布用）
```

**架构原则**：单一职责；依赖高层→低层不反向（detectors/identifiers → core）；策略模式可插拔。

**检测器接口**：3 个基类（PerTuple/PerFrame/TwoPass，模板方法模式），`detect()` 接口一致、结果统一为矩阵式 `List[Interaction]`，由 `pipeline.py` 的 `STRATEGY_INDEX` 切换。详见 `core/interfaces.py`。

## 使用

```bash
dii run -t md.tpr -f md.xtc -o out/ --ff amber     # 检测→h5（--ff 必选，当前仅 amber）
#   --interactions hydrogen_bond,pi_stacking      # 只检部分类型（默认 all=8）
#   --strategy two_pass|per_frame|per_tuple       # 策略（默认 two_pass）
dii export -i out/salt_bridge.h5 -o out_export/   # 导出 xvg/xpm/csv + 概览（支持多 Interaction h5）
```

```python
from DuIvyInteractions.pipeline import Pipeline
Pipeline(ff="amber", strategy="two_pass").run("md.tpr", "md.xtc", "out/", interactions=None)
```

## 当前代码状态（2026-09-16）

- ✅ 阶段一：基团鉴定（D927 验证完成，已迁移进新架构；`Tests/original_draft/` 保留历史脚本与调研，不随包发布）
- ✅ 阶段二：相互作用检测（8 类型 × 3 策略，TwoPass 水桥 KDTree 优化 65h→~5s）
- ✅ 阶段三：结果存储与导出（HDF5 序列化 + xvg/xpm/CSV 导出器 + XPM 手动构建修复）
- ✅ 阶段四：命令行（dii run / dii export，多 Interaction 遍历 + 空数据/损坏 h5/0 帧防护 + 索引校验）

**未实现**（详见 `doc/TODO.md`）：可视化、`utils/output.py`、基团识别结果人工审查、长轨迹 PerFrame 内存优化、PBC 处理等。

## 常用命令

```bash
# tpr dump（⚠️ 必须用 DIP 环境的 gmx，MDsoftware 的 gmx 跑不起来）
/home/hanyl/.micromamba/envs/DIP/bin/gmx dump -s <md.tpr> > dump.txt 2> dump.log
```