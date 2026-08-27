# DuIvyInteraction

基于 MD 拓扑力场参数的分子间相互作用判定工具。

## 核心思路

现有工具（PLIP/ProLIF）在分析 MD 轨迹时，通过 OpenBabel/RDKit 重建化学信息（键序、芳香性、加氢），丢弃了 MD 力场拓扑中原有的化学语义。本项目直接从 **GROMACS tpr 拓扑中读取力场原子类型**，确定性识别官能团，与模拟力场完全自洽。

## 两段式架构

1. **基团鉴定**：tpr 力场参数（原子类型 + 键合图 + 显式 H + 电荷）→ 确定性识别化学基团（芳香环、H 键供/受体、带电基团…）。只做一次，与帧无关。
2. **几何判定**：逐帧用 PLIP 式距离/角度/平面临近判据 → 相互作用列表。

## 目录结构

```
DuIvyInteraction/
│
├── DuIvyInteractions/                  # 主包
│   ├── core/                           # 核心基础：数据类 + 接口 + 常量
│   │   ├── datas.py                    # Group, AtomData, SystemData, Interaction, ...
│   │   ├── interfaces.py               # Reader, GroupIdentifier, InteractionDetector ABC
│   │   └── constants.py                # GROUP_TYPES, BOND_TYPES, 元素周期表, ...
│   │
│   ├── input_readers/                  # 文件 → SystemData
│   │   ├── gmx_tpr_dump_reader.py      # gmx dump 文本解析
│   │   └── gmx_tpr_reader.py           # tpr 二进制解析（MDAnalysis）
│   │
│   ├── group_identifiers/              # SystemData → List[Group]
│   │   └── amber_ff_identifier.py      # Amber 力场基团识别
│   │
│   ├── interaction_detectors/          # Group[] + 坐标 → Interaction[]（待实现）
│   ├── visualizers/                    # 结果可视化（待实现）
│   ├── utils/                          # 无状态工具函数
│   └── pipeline.py                     # 主流程编排（待实现）
│
├── docs/                               # 设计文档
├── pyproject.toml                      # 包管理配置
├── CLAUDE.md                           # 项目规范
└── README.md
```

## 依赖

- Python >= 3.9
- NumPy >= 1.20
- MDAnalysis >= 2.0（用于 tpr 二进制读取）
- GROMACS（`gmx dump` 用于文本格式 tpr 解析）

## 安装

```bash
pip install -e .
```

## 力场支持

已验证 Amber 家族全部力场（amber03/94/96/99/99sb/99sb-ildn/GS/14sb + GAFF 配体），类型名→化学特征映射**零冲突**。