# DuIvyInteractions 文档

基于 MD 拓扑力场参数的分子间相互作用判定工具。

**核心思路**：直接从 GROMACS tpr 拓扑读取力场原子类型，确定性识别化学基团，与模拟力场完全自洽——不用 OpenBabel/RDKit 从坐标重建化学。

```{toctree}
:maxdepth: 1
:caption: 使用指南

guide/index
```

```{toctree}
:maxdepth: 1
:caption: 参考手册

reference/index
```

```{toctree}
:maxdepth: 1
:caption: 其他

changelog
```