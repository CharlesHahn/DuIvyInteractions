# Group 基团数据结构设计文档

> 创建日期：2026-08-24
> 状态：已定稿
> 作者：dsh-tui+hanyl

---

## 1. 设计背景

### 1.1 项目目标

构建通用的分子相互作用分析工具，核心流程：

```
拓扑(tpr) → 基团识别 → 基团列表
轨迹(xtc) → 逐帧坐标 → 基团坐标 → 几何判定 → 相互作用矩阵
```

### 1.2 Group 的角色

Group 是整个系统的**核心数据单元**，是"识别器"和"判定器"之间的桥梁：

```
identifiers/amber.py  →  产出 Group 对象
                              ↓
detectors/hbond.py    ←  消费 Group 对象
```

所有模块必须使用统一的 Group 数据结构，否则无法对接。

---

## 2. 第一性原理分析

### 2.1 Group 是什么？

**Group = 一组能参与相互作用的原子 + 它们的化学语义**

### 2.2 Group 需要回答的问题

| 问题 | 用途 | 必要性 |
|:----|:----|:----|
| **我是谁？** | 唯一标识，去重 | 必须 |
| **我是什么类型？** | 决定能参与哪种相互作用 | 必须 |
| **我包含哪些原子？** | 查坐标、算距离 | 必须 |
| **我属于哪个分子/残基？** | 排除同分子内相互作用 | 必须 |
| **每个原子的属性？** | 化学判定（电荷、元素、类型） | 必须 |
| **我的几何中心在哪？** | 算距离 | 可动态计算 |
| **我的法向量？** | 算 π-π 角度 | 可动态计算 |
| **我的键连接？** | 调试、验证、扩展 | 可选（存 metadata） |
| **谁识别了我？** | 调试、可复现 | 可选（存 metadata） |

### 2.3 关键决策：几何信息不存 Group

**问题**：几何信息（中心、法向量）要存在 Group 里吗？

**答案**：**不存，动态计算**

理由：
1. 中心/法向量依赖坐标，每帧都变
2. Group 是拓扑信息，与帧无关
3. 让 Detector 自己算，保持 Group 简单

```
Group（静态）      坐标（动态）      Detector（计算）
atom_indices  +   coordinates   →   center, normal, distance, angle
```

### 2.4 关键决策：bond 信息不存 Group

**问题**：bond 信息需要存在 Group 里吗？

**逐类型分析**：

| 基团类型 | 识别时需要 bonds？ | 判定时需要 bonds？ | 结论 |
|:----|:----|:----|:----|
| **H_donor** | ✓ 找 D-H 对 | ✗ 已知 D 和 H 原子 | 不需要 |
| **aromatic_ring** | ✓ 图论找环 | ✗ 用环原子算中心/法向量 | 不需要 |
| **H_acceptor** | ✗ 看原子类型 | ✗ 只用受体坐标 | 不需要 |
| **charged_*** | ✗ 看电荷 | ✗ 只用原子坐标 | 不需要 |
| **halogen** | ✗ 看元素 | ✗ 只用原子坐标 | 不需要 |
| **metal** | ✗ 看元素 | ✗ 只用原子坐标 | 不需要 |
| **water** | ✗ 看残基名 | ✗ 用 O 坐标 | 不需要 |
| **hydrophobic** | ✗ 看原子类型 | ✗ 只用原子坐标 | 不需要 |

**结论**：bond 信息只在 Identifier 内部使用，不传递到 Group。

**数据流**：
```
拓扑（含 bonds）
    ↓
Identifier（内部使用 bonds 找 Group）
    ↓
Group（不含 bonds，只有 atom_indices）
    ↓
Detector（用 atom_indices + 坐标算距离/角度）
```

**备选方案**：如需保留键信息，可存入 `metadata["bonds"]`（见第 5 节）。

---

## 3. 最终设计

### 3.1 数据结构定义

```python
from dataclasses import dataclass, field
from typing import List, Optional, Dict

@dataclass
class Group:
    """一个可参与相互作用的基团"""
    
    # === 身份信息 ===
    group_id: int                           # 唯一标识
    group_type: str                         # 基团类型（见常量）
    
    # === 位置信息 ===
    molecule: str                           # 所属分子名（如 "D927"）
    residue_name: str                       # 残基名（如 "TYR"）
    residue_id: int                         # 全局残基号
    atom_indices: List[int]                 # 全局原子索引列表
    
    # === 原子属性（与 atom_indices 等长）===
    atom_types: List[str]                   # 力场原子类型（如 ["ca", "ca", "ca"]）
    elements: List[str]                     # 元素符号（如 ["C", "C", "C"]）
    charges: List[float]                    # 原子电荷（如 [-0.1, -0.1, -0.1]）
    
    # === 附加信息 ===
    metadata: Dict = field(default_factory=dict)  # 识别器、键信息等
```

### 3.2 基团类型常量

```python
GROUP_TYPES = {
    # H键相关
    "H_donor",              # H键供体（D-H 对）
    "H_acceptor",           # H键受体
    
    # π 相关
    "aromatic_ring",        # 芳香环（π-π, π-阳离子, 卤-π）
    
    # 电荷相关
    "charged_positive",     # 正电基团（盐桥）
    "charged_negative",     # 负电基团（盐桥）
    
    # 其他
    "halogen",              # 卤素（卤键）
    "metal",                # 金属中心（金属配位）
    "water",                # 水分子（水桥）
    "hydrophobic",          # 疏水基团
}
```

### 3.3 设计要点

| 设计决策 | 理由 |
|:----|:----|
| `charges` 用列表 | 每个原子电荷不同，供体需要 D 和 H 的电荷 |
| `elements` 用列表 | 与 `atom_indices` 对齐，方便遍历 |
| `atom_types` 必填 | 保留力场类型信息，方便调试和扩展 |
| `metadata` 替代 `identifier_name` | 更灵活，可存任意附加信息 |
| `residue_id` 全局 | 跨分子唯一标识 |
| `H_donor` / `H_acceptor` | 更明确，避免歧义 |
| 不存几何信息 | 依赖坐标，动态计算 |
| 不存 bond 信息 | 只在识别阶段使用，判定阶段不需要 |

---

## 4. 使用示例

### 4.1 芳香环

```python
Group(
    group_id=1,
    group_type="aromatic_ring",
    molecule="D927",
    residue_name="D927",
    residue_id=200,
    atom_indices=[2359, 2360, 2361, 2362, 2363, 2364],
    atom_types=["ca", "ca", "ca", "ca", "ca", "ca"],
    elements=["C", "C", "C", "C", "C", "C"],
    charges=[0.34, -0.36, -0.40, -0.33, -0.33, -0.36],
    metadata={"identifier": "amber", "ring_size": 6}
)
```

### 4.2 H键供体

```python
Group(
    group_id=2,
    group_type="H_donor",
    molecule="RBD_pro",
    residue_name="TYR",
    residue_id=10,
    atom_indices=[133, 134],  # O-H 对
    atom_types=["OH", "HO"],
    elements=["O", "H"],
    charges=[-0.558, 0.399],
    metadata={"identifier": "amber", "donor_element": "O"}
)
```

### 4.3 金属中心

```python
Group(
    group_id=3,
    group_type="metal",
    molecule="Mg",
    residue_name="Mg",
    residue_id=500,
    atom_indices=[5056],
    atom_types=["MG"],
    elements=["Mg"],
    charges=[2.0],
    metadata={"identifier": "amber"}
)
```

### 4.4 带键信息的芳香环（可选）

```python
Group(
    group_id=1,
    group_type="aromatic_ring",
    molecule="D927",
    residue_name="D927",
    residue_id=200,
    atom_indices=[2359, 2360, 2361, 2362, 2363, 2364],
    atom_types=["ca", "ca", "ca", "ca", "ca", "ca"],
    elements=["C", "C", "C", "C", "C", "C"],
    charges=[0.34, -0.36, -0.40, -0.33, -0.33, -0.36],
    metadata={
        "identifier": "amber",
        "ring_size": 6,
        "bonds": {
            "2359_2360": "aromatic",
            "2360_2361": "aromatic",
            "2361_2362": "aromatic",
            "2362_2363": "aromatic",
            "2363_2364": "aromatic",
            "2359_2364": "aromatic",
        },
        "external_bonds": {
            "2355_2359": "single",  # 连接到侧链
        },
    }
)
```

---

## 5. Bond 信息存储方案（可选）

### 5.1 设计原则

1. **可选**：metadata 是 dict，bond 信息按需添加
2. **高性能**：用字典，O(1) 查找
3. **规范化**：键格式保证唯一性

### 5.2 存储格式

```python
metadata["bonds"] = {
    "{较小索引}_{较大索引}": "键级",
    ...
}

metadata["external_bonds"] = {
    "{较小索引}_{较大索引}": "键级",
    ...
}
```

### 5.3 键生成函数

```python
def make_bond_key(idx1: int, idx2: int) -> str:
    """生成规范化的键字符串（较小索引在前）"""
    return f"{min(idx1, idx2)}_{max(idx1, idx2)}"
```

### 5.4 键级常量

```python
BOND_ORDERS = {
    "single",      # 单键
    "double",      # 双键
    "triple",      # 三键
    "aromatic",    # 芳香键
    "constrained", # LINCS 约束
    "virtual",     # 虚拟位点
}
```

### 5.5 访问示例

```python
# 检查两个原子之间是否有键
bond_key = make_bond_key(2359, 2360)
if bond_key in group.metadata.get("bonds", {}):
    bond_order = group.metadata["bonds"][bond_key]
    print(f"键级: {bond_order}")

# 获取某个原子的所有连接
def get_bonds_for_atom(group, atom_idx):
    """获取指定原子的所有键"""
    bonds = group.metadata.get("bonds", {})
    return {
        k: v for k, v in bonds.items()
        if str(atom_idx) in k.split("_")
    }
```

### 5.6 性能对比

| 方式 | 查找复杂度 | 示例 |
|:----|:----|:----|
| **tuple 列表** | O(n) | `for (a, b, order) in bonds: if a == 2359 ...` |
| **字典** | O(1) | `bonds["2359_2360"]` |

---

## 6. 设计原则总结

| 原则 | 体现 |
|:----|:----|
| **最小化** | 只存必要信息，几何信息动态算 |
| **与帧无关** | Group 是拓扑信息，不依赖轨迹 |
| **自描述** | 包含分子、残基、来源等上下文 |
| **可扩展** | `metadata` 字典可存任意附加信息 |
| **高性能** | bond 信息用字典存储，O(1) 查找 |
| **类型安全** | 列表等长约束（atom_indices, atom_types, elements, charges） |

---

## 7. 依赖关系

```
core/data.py          ← 定义 Group 数据类
      ↑
identifiers/amber.py  ← 产出 Group 对象
      ↓
detectors/*.py        ← 消费 Group 对象
```

core/ 是最底层，谁都不依赖，所有人都依赖它。

---

## 8. 后续扩展

| 扩展场景 | 实现方式 |
|:----|:----|
| 新基团类型 | 在 `GROUP_TYPES` 常量中添加 |
| 新原子属性 | 在 Group 中添加新列表字段 |
| 新附加信息 | 在 `metadata` 字典中添加 |
| 键信息 | 存入 `metadata["bonds"]`（见第 5 节） |

---

*文档结束*
