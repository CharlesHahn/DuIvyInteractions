# 未来待优化事项清单

> 创建日期：2026-08-25

---

## 1. bond_type 精度丢失

GmxTprReader 无法区分 Bond/Constraint/Settle 来源，全部标记为 `"bond"`。GmxTprDumpReader 无此问题。

## 2. 键级缺失

tpr 文件不存储键级（单键/双键/芳香键）。两个 Reader 的 bond_type 只能标记为 `"bond"`（未知键级）。

## 3. segid 命名不一致

GmxTprReader 输出 `seg_0_RBD_pro`，GmxTprDumpReader 输出 `RBD_pro`。MDAnalysis 自动加前缀。

## 4. resid 编号不一致

GmxTprReader 输出 MDA 连续编号（1~142），GmxTprDumpReader 输出 PDB 原始编号（157~298）。MDAnalysis 默认 `tpr_resid_from_one=True`。

## 5. SOL 键数不一致

GmxTprDumpReader 从 SETTLE 生成 3 个键（含 H1-H2），GmxTprReader 只有 2 个（O-H1, O-H2）。MDAnalysis 不把 H1-H2 算作键。

---

*文档结束*
