# 参考手册

本部分面向**想理解原理与自定义的科研用户**：软件核心概念、力场类型映射、基团识别规则、相互作用判据、数据格式、Python API 与扩展方法等技术细节。

## 当前状态概览

DuIvyInteractions v0.0.1 的参考实现状态如下（与代码——`DuIvyInteractions/`——逐项对应）：

- **两段式架构**：① 基团鉴定——只做一次、与帧无关，从 tpr 原子类型 + 键合图 + 显式 H + 电荷确定性地识别基团；② 几何判定——逐帧按距离/角度/平面临近判据判定 8 类相互作用。详见[核心概念](concepts)、[基团识别规则](group_rules)与[相互作用判据](criteria)。
- **4 力场识别器**：Amber 家族（蛋白 + GAFF/GAFF2）、GROMOS 53A6/54A7、CHARMM36/C36m（含 CGenFF）、OPLS-AA/L，经 `IDENTIFIER_CLASSES` 注册表统一管理；水残基排除用各识别器自带的 `WATER_RESIDUES` 类属性（跨力场正确排除 CHARMM TIP3、OPLS HO4/HO5 等）。详见[力场类型映射](force_field)与[Python API](api)。
- **8 类相互作用 × 3 策略**：`two_pass`（默认，性能最优）/ `per_frame`（现行）/ `per_tuple`（策略一，处于"可能被舍弃"状态）；H 键受体判定已完成 A1 修复（剔除铵/胍基/带 H 吡咯等非受体 N，证据见 `doc/acceptor_identification_evidence.md`）。
- **存储与导出**：HDF5（格式版本 1.0）无损序列化 + 8 个导出器输出 xvg/xpm/csv；命令行 `dii run` / `dii export`。详见[数据格式](data_format)与[命令参考](../guide/command.md)。
- **真实测试数据**：Amber KRAS-RBD D927、GROMOS 53A6（蛋白 + 6 配体）、CHARMM36 SMO-BST 三套真实 MD 数据（`Tests/test_MD_case_*`），8 类相互作用全覆盖。
- **范围与限制**：PBC、金属配位几何构型、水桥去重、疏水-芳香去重、可视化等均未实现，详见[已知限制](limitations)。

```{toctree}
:maxdepth: 1
:caption: 核心概念

concepts
glossary
```

```{toctree}
:maxdepth: 1
:caption: 科学细节

force_field
group_rules
criteria
```

```{toctree}
:maxdepth: 1
:caption: 数据与扩展

data_format
api
extension
```

```{toctree}
:maxdepth: 1
:caption: 范围与限制

limitations
```