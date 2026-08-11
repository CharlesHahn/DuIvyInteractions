# DuIvyInteraction

基于 MD 拓扑力场参数的分子间相互作用判定工具。

## 核心思路

现有工具（PLIP/ProLIF）在分析 MD 轨迹时，通过 OpenBabel/RDKit 重建化学信息（键序、芳香性、加氢），丢弃了 MD 力场拓扑中原有的化学语义。本项目直接从 **GROMACS tpr 拓扑中读取力场原子类型**，确定性识别官能团，与模拟力场完全自洽。

## 两段式架构

1. **基团鉴定**：tpr 力场参数（原子类型 + 键合图 + 显式 H + 电荷）→ 确定性识别化学基团（芳香环、H 键供/受体、带电基团…）。只做一次，与帧无关。
2. **几何判定**：逐帧用 PLIP 式距离/角度/平面临近判据 → 相互作用列表。

## 当前状态

已完成第一阶段（基团鉴定）的可行性验证（D927 体系）：

- ✅ 解析器：从 `gmx dump` 提取原子/类型/键/约束/残基
- ✅ 特征映射表：GAFF + Amber 系列力场（amber03~14sb）的类型→化学特征覆盖
- ✅ 环检测：SSSR 最小环 + 稠合冗余去除，D927 3 芳香环鉴定正确
- ✅ H 键供体/受体、卤素、带电基团鉴定
- ✅ 映射表自动验证器（`verify_type_mapping.py`）

## 目录结构

```
├── CLAUDE.md                    # 项目文档
├── tpr_to_functional_groups.md  # 验证记录
├── plip_md_research_survey.md   # 竞品调研
├── parse_tpr_dump.py            # tpr dump 解析器
├── functional_groups.py         # 特征映射+环检测+官能团鉴定
├── verify_type_mapping.py       # 映射表自动验证器
├── fg_report.txt                # 官能团报告
├── dump_md_D927.tpr.txt/.log    # 原始 tpr dump 数据
└── README.md
```

## 依赖

- GROMACS（`gmx dump` 用于提取 tpr 拓扑）
- Python 3.8+（标准库，无外部依赖）

## 力场支持

已验证 Amber 家族全部力场（amber03/94/96/99/99sb/99sb-ildn/GS/14sb + GAFF 配体），类型名→化学特征映射**零冲突**。详见 `verify_type_mapping.py`。