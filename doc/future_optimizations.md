# 未来待优化事项清单

> 创建日期：2026-08-25
> 状态：待处理

---

## 1. bond_type 精度问题

### 1.1 问题描述

当前两个 Reader 的 bond_type 精度不同：

| Reader | bond | constrained | settle | 信息来源 |
|--------|------|-------------|--------|---------|
| GmxTprDumpReader | ✅ | ✅ | ✅ | dump 文本的 Bond/Constraint/Settle 段 |
| GmxTprReader | ✅ 全部标记为 bond | ❌ | ❌ | MDAnalysis 扁平化输出，丢失来源 |

### 1.2 影响

- **基团鉴定**：H 键供体识别需要区分 D-H 键是真正的共价键还是约束（constrained）。当前 MDA Reader 无法精确区分。
- **水分子分析**：无法从 MDA Reader 识别哪些键来自 SETTLE。

### 1.3 优化方案

**方案 A：启发式规则推断**

```python
# 伪代码
if resname == 'SOL':
    bond_type = 'settle'
elif is_NH_bond(atom1, atom2) and is_protein(resname):
    bond_type = 'constrained'  # 可能误判
else:
    bond_type = 'bond'
```

风险：N-H 键在非蛋白分子中可能是真正的共价键，不能简单标记为 constrained。

**方案 B：MDAnalysis 底层扩展**

修改 MDAnalysis 的 TPRParser，保留 ilist 来源信息。需要：
- Fork MDAnalysis 或提交 PR
- 在 Bond 对象中添加 `source` 属性
- 工作量大，维护成本高

**方案 C：混合方案**

优先使用 GmxTprDumpReader（精度最高），MDA Reader 作为降级方案。在文档中明确说明精度差异。

### 1.4 优先级

**P2（低优先级）**：当前两个 Reader 的差异不影响核心功能。基团鉴定主要依赖原子类型和键连接图，不依赖 bond_type。

---

## 2. 键级（bond order）缺失

### 2.1 问题描述

tpr 文件的 Bond 段全部是 func=1（ harmonic potential），不直接存储键级（单键/双键/芳香键）。当前两个 Reader 的 bond_type 只能标记为 `"bond"`（未知键级）。

### 2.2 影响

- **芳香环检测**：不能依赖键级判断芳香性。需要从原子类型推断（`ca`=芳香碳、`nb`=吡啶氮等）。
- **键类型分析**：无法区分单键、双键、芳香键。

### 2.3 优化方案

**方案 A：从原子类型推断（推荐）**

```python
# 伪代码
AROMATIC_TYPES = {'ca', 'nb', 'na', 'cp', 'cg'}
if atom1.type in AROMATIC_TYPES and atom2.type in AROMATIC_TYPES:
    bond_type = 'aromatic'
```

风险：某些原子类型可能有歧义（如 `CA` 在 amber 中是芳香碳，在 GROMOS 中是 α 碳）。

**方案 B：从键长推断**

```python
# 伪代码
b0 = bond_length(atom1, atom2)
if 1.34 < b0 < 1.41:
    bond_type = 'aromatic'
elif b0 < 1.34:
    bond_type = 'double'
elif b0 > 1.50:
    bond_type = 'single'
```

风险：b0 是平衡键长，不是实际键长。用户明确禁止此方案（"一定不允许 b0 来做键级推断！"）。

**方案 C：从 tpr 的 func 参数推断**

tpr 文件中不同 func 类型对应不同键：
- func=1: harmonic (可能是单键/双键/芳香键)
- func=2: constraint (LINCS 约束)
- func=4/5: 离散/芳香二面角（可间接推断芳香性）

需要解析 tpr 的力场参数段，获取 func 类型。

### 2.4 优先级

**P1（中优先级）**：基团鉴定需要芳香性信息，但可以从原子类型推断，不依赖键级。

---

## 3. segid 命名差异

### 3.1 问题描述

两个 Reader 的分子名（molecule_name）不同：

| dump Reader | MDA Reader |
|-------------|------------|
| `RBD_pro` | `seg_0_RBD_pro` |
| `D927` | `seg_1_D927` |
| `SOL` | `seg_5_SOL` |

### 3.2 影响

- **下游代码**：需要处理两种命名格式。
- **用户交互**：输出结果中的分子名不一致。

### 3.3 优化方案

**方案 A：GmxTprReader 剥离前缀**

```python
# 在 _build_residues 中
molecule_name = seg.segid
if molecule_name.startswith('seg_') and '_' in molecule_name[4:]:
    # seg_0_RBD_pro → RBD_pro
    molecule_name = '_'.join(molecule_name.split('_')[2:])
```

风险：某些 segid 本身以 `seg_` 开头，剥离会误伤。

**方案 B：统一格式**

两个 Reader 都输出带前缀的格式，或都不带前缀。需要设计文档明确规定。

**方案 C：保持现状**

在文档中说明差异，下游代码自行处理。

### 3.4 优先级

**P2（低优先级）**：不影响核心功能，下游代码可以适配。

---

## 4. resid 编号差异

### 4.1 问题描述

两个 Reader 的 `residue_idx_in_molecule` 含义不同：

| Reader | 含义 | RBD_pro 第一个残基 |
|--------|------|-------------------|
| GmxTprDumpReader | PDB 原始编号 | 157 |
| GmxTprReader | MDA 连续编号 | 1 |

### 4.2 影响

- **残基引用**：无法从 MDA Reader 获取 PDB 原始编号。
- **输出可读性**：MDA Reader 的编号对用户不直观。

### 4.3 优化方案

**方案 A：使用 MDA 的 resnum 属性**

MDAnalysis 有 `resnum` 属性，可能保留原始编号。需要验证。

**方案 B：从 tpr 文件直接读取**

tpr 文件的 moltype 段包含 `nr` 字段（PDB 编号）。需要解析 tpr 的残基信息段。

**方案 C：保持现状**

在文档中说明差异，用户按需选择 Reader。

### 4.4 优先级

**P2（低优先级）**：不影响核心功能。如果需要 PDB 编号，使用 GmxTprDumpReader。

---

## 5. SOL 键数差异

### 5.1 问题描述

SETTLE 约束了 3 个距离（O-H1, O-H2, H1-H2），但：

| Reader | SOL 键数 | 说明 |
|--------|---------|------|
| GmxTprDumpReader | 3 | 包含 H1-H2 |
| GmxTprReader | 2 | 只有 O-H1, O-H2 |

### 5.2 影响

- **水分子分析**：H1-H2 键对水桥分析无影响。
- **环检测**：水分子不参与环检测，无影响。

### 5.3 优化方案

**保持现状**：2 个键已足够。H1-H2 不是真正的共价键，不需要作为"键"存储。

### 5.4 优先级

**P3（无需处理）**：当前行为正确。

---

## 6. 总结

| 事项 | 优先级 | 建议方案 |
|------|--------|---------|
| bond_type 精度 | P2 | 启发式规则推断或保持现状 |
| 键级缺失 | P1 | 从原子类型推断芳香性 |
| segid 命名差异 | P2 | 剥离前缀或保持现状 |
| resid 编号差异 | P2 | 使用 resnum 或保持现状 |
| SOL 键数差异 | P3 | 无需处理 |

---

## 7. 决策记录

| 日期 | 决策 | 理由 |
|------|------|------|
| 2026-08-25 | bond_type 全部标记为 bond（MDA Reader） | MDA 不区分来源，无法精确还原 |
| 2026-08-25 | SOL 保留 3 个键（dump Reader） | 尊重原始 dump 数据，SETTLE 约束 3 个距离 |
| 2026-08-25 | 不从 b0 推断键级 | 用户明确禁止 |
| 2026-08-25 | 两个 Reader 并存 | 精度 vs 便利性，用户按需选择 |

---

*文档结束*
