# 疏水基团识别定义

> 状态：定义定稿

---

## 1. 化学定义

疏水性 = 不能与水形成氢键 = 周围没有电负性原子。

疏水相互作用源于熵变（疏水效应），不是吸引力。观察到的疏水接触随距离指数衰减。

---

## 2. 疏水原子判定规则

一个原子是疏水原子，当且仅当同时满足以下两个条件：

**条件 1**：原子是碳（element == 'C'）。

**条件 2**：所有邻居原子都是碳或氢（neighbor element ∈ {'C', 'H'}）。

```
对每个原子 a：
  if a.element != 'C' → 非疏水
  if 任何邻居的 element ∉ {'C', 'H'} → 非疏水
  else → 疏水原子
```

### 2.1 判定示例

| 原子 | 邻居 | 疏水? | 原因 |
|------|------|-------|------|
| 甲基 C | 3H + 1C | ✅ | 邻居全是 C/H |
| 亚甲基 C | 2H + 2C | ✅ | 邻居全是 C/H |
| 苯环 C | 2C + 1H | ✅ | 邻居全是 C/H |
| 羰基 C | 1O + ... | ❌ | 有 O 邻居 |
| 醇 C | 1O + ... | ❌ | 有 O 邻居 |
| 酰胺 C | 1N + 1O | ❌ | 有 N/O 邻居 |
| 腈 C | 1N + ... | ❌ | 有 N 邻居 |
| 氟代 C | 3F + 1C | ❌ | 有 F 邻居 |

### 2.2 与芳香环的关系

芳香环中的碳原子（如苯环 C）也满足疏水条件（C + 邻居全是 C/H）。

**在基团识别阶段不排除芳香环原子**。疏水基团和芳香环基团是独立识别的，允许原子重叠。

**去重在相互作用检测阶段处理**：如果两个芳香环已经形成 π-π 堆积，则移除它们之间的疏水相互作用（避免重复计数同一物理接触）。

---

## 3. 基团数据结构

每个疏水原子生成一个 Group：

| 字段 | 说明 |
|------|------|
| group_type | `"hydrophobic"` |
| atoms | 单个碳原子的 AtomData |
| metadata | 空（无额外属性） |

---

## 4. 相互作用检测阶段的重叠标记

疏水基团识别不排除芳香环原子。若两个芳香环已经形成 π-π 堆积，同一物理接触可能同时被
报告为 π-π 与疏水；默认行为保持完整（不删除），由 `dii run --mark-redundant` 在检测后
以原子级、逐帧方式标记：

**规则**：如果两个芳香环已经形成 π-π 堆积，则把这两个环的原子之间的疏水对-帧标记为
`redundant_pi_stacking`（`Interaction.flags`），原始疏水记录完整保留。

**原因**：π-π 堆积已经包含了芳香环之间的疏水接触；但直接删除会丢失环-脂肪碳接触、
未形成 π-π 的环-环接触以及原始记录的可追溯性。

**实现**：`DuIvyInteractions/redundancy.py` 的 `pi_overlap_flags()`；编排见
`pipeline.py`；设计见 `doc/redundancy_marking_design.md`。

---

## 5. 参考文献

1. PLIP 源码: https://github.com/pharmai/plip — `structure/preparation.py` (`hydrophobic_atoms`)
2. PLIP 文档: "An atom is classified as hydrophobic if it is a carbon and has only carbon or hydrogen atoms as neighbors."

---

*文档结束*
