# SystemData 体系分子数据层设计文档

> 创建日期：2026-08-24
> 状态：已定稿
> 作者：dsh-tui+hanyl

---

## 1. 设计背景

### 1.1 项目目标

构建通用的分子相互作用分析工具，核心流程：

```
输入文件（tpr/prmtop/pdb） → Reader → SystemData → Identifier → List[Group]
```

### 1.2 SystemData 的角色

SystemData 是 Reader 和 Identifier 之间的**统一数据层**：

```
Reader 层                    数据层                Identifier 层
(可插拔)                     (统一)                (可插拔)

TprReader    ──┐                                ┌── AmberIdentifier
AmberReader  ──┼──→   SystemData   ──→  ──┼── RdkitIdentifier
PDBReader    ──┘                                └── SmartsIdentifier
```

**核心价值**：将不同引擎的文件格式抽象为统一的分子数据，解耦输入和识别。

---

## 2. 第一性原理分析

### 2.1 体系是什么？

**体系 = 一组残基 + 残基间的连接关系**

关键特征：
- **残基是基本单位**：相互作用发生在残基级别（如 TYR165 ↔ D927）
- **原子属于残基**：每个原子都有明确的残基归属
- **键有两种**：残基内键 + 残基间键（如肽键、二硫键）

### 2.2 为什么以残基为中心？

**问题**：蛋白是一个"分子"，但有 142 个残基。相互作用发生在残基级别。

**分析**：

| 粒度 | 优点 | 缺点 |
|:----|:----|:----|
| 以分子为中心 | 简单 | 粒度太粗，不利于交互分析 |
| 以原子为中心 | 粒度最细 | 粒度太细，缺乏语义 |
| **以残基为中心** | 符合化学直觉，粒度适中 | - |

**结论**：残基是相互作用分析的自然单位。

### 2.3 需要什么信息？

**Identifier 需要的信息**：

| 信息 | 用途 | 必要性 |
|:----|:----|:----|
| 原子类型名 | 特征映射（`ca` → 芳香碳） | 必须 |
| 键连接 | 环检测、供体识别 | 必须 |
| 残基信息 | 输出定位 | 必须 |
| 电荷 | 供受体/带电基团判定 | 必须 |
| 元素符号 | 元素识别 | 必须 |
| 原子质量 | 质心计算、质量加权分析 | 建议有 |
| 残基间键 | 骨架连接（未来扩展） | 建议有 |

### 2.4 索引设计

**核心问题**：索引是全局的还是局部的？

**分析**：
- 原子需要两个索引：全局唯一（用于距离计算）+ 残基内（用于键连接）
- 残基需要两个索引：全局唯一（用于标识）+ 分子内（PDB 编号）

**原则**：索引命名必须清晰体现是全局还是局部。

---

## 3. 最终设计

### 3.1 AtomData

```python
@dataclass
class AtomData:
    """单个原子的信息。"""
    atom_global_idx: int        # 全局原子索引（整个体系唯一）
    atom_idx_in_residue: int    # 残基内索引（用于残基内键连接）
    atom_name: str              # 原子名（如 "CG"）
    atom_type: str              # 力场类型（如 "ca"）
    atom_element: str           # 元素符号（如 "C"）
    atom_charge: float          # 电荷
    atom_mass: float            # 原子质量（原子质量单位，-1.0 表示未设置）
```

**设计决策**：

| 决策 | 选择 | 理由 |
|:----|:----|:----|
| 两个索引 | `atom_global_idx` + `atom_idx_in_residue` | 全局用于距离计算，局部用于键连接 |
| 无 atomic_number | 删除 | 与 element 冗余，可推导 |
| 字段前缀 `atom_` | 统一 | 与残基字段区分，避免混淆 |
| mass 默认值 | `-1.0` | 区分"未设置"与虚拟粒子（mass=0） |

### 3.2 BondData

```python
@dataclass
class BondData:
    """残基内的一个键。"""
    atom1_idx_in_residue: int           # 残基内索引
    atom2_idx_in_residue: int           # 残基内索引
    bond_type: str                      # 键类型，见 BOND_TYPES
```

**设计决策**：

| 决策 | 选择 | 理由 |
|:----|:----|:----|
| 合并 constraints | 是 | 约束也是一种"键"，统一存储 |
| 使用局部索引 | `atom_idx_in_residue` | 残基内键，局部索引足够 |
| 键类型字段 | `bond_type` | 区分来源：常规键、约束、SETTLE |

**BOND_TYPES 取值**：

| bond_type | 说明 | 来源 |
|:----|:----|:----|
| `"bond"` | 常规键（键级未知） | dump Bond 段 |
| `"constrained"` | LINCS 约束 | dump Constraint 段（如 N-H 键） |
| `"settle"` | SETTLE 水约束 | dump Settle 段（如 O-H、H-H 距离） |
| `"single"` | 单键 | tpr func=1, b0≈1.5Å |
| `"double"` | 双键 | tpr func=1, b0≈1.3Å |
| `"triple"` | 三键 | tpr func=1, b0≈1.2Å |
| `"aromatic"` | 芳香键 | tpr func=1, b0≈1.4Å + 芳香类型 |
| `"virtual"` | 虚拟位点 | tpr func=3+ |

### 3.3 ResidueData

```python
@dataclass
class ResidueData:
    """一个残基的数据。"""
    residue_name: str                   # 残基名（如 "TYR"）
    residue_global_idx: int             # 全局残基索引（整个体系唯一）
    residue_idx_in_molecule: int                # 分子内残基编号（PDB 编号）
    molecule_name: str                  # 所属分子名（如 "RBD_pro"）
    atoms: List[AtomData]               # 残基内的原子
    bonds: List[BondData]               # 残基内的键
```

**设计决策**：

| 决策 | 选择 | 理由 |
|:----|:----|:----|
| 两个残基索引 | `residue_global_idx` + `residue_idx_in_molecule` | 全局唯一标识，分子内 PDB 编号 |

### 3.4 InterResidueBond

```python
@dataclass
class InterResidueBond:
    """残基间的共价键（如肽键、二硫键）。"""
    residue1_global_idx: int            # 残基 1 的全局索引
    atom_idx_in_residue1: int                # 原子在残基 1 内的索引
    residue2_global_idx: int            # 残基 2 的全局索引
    atom_idx_in_residue2: int                # 原子在残基 2 内的索引
    bond_type: str                      # 键类型
```

**设计决策**：

| 决策 | 选择 | 理由 |
|:----|:----|:----|
| 命名方式 | `residue1/2` + `atom_idx_in_residue1/2` | 清晰标识两个端点 |
| 索引类型 | 残基用全局，原子用局部 | 残基跨残基唯一，原子在残基内 |

### 3.5 SystemData

```python
@dataclass
class SystemData:
    """体系数据，连接 Reader 和 Identifier。"""
    system_name: str
    residues: List[ResidueData]
    inter_residue_bonds: List[InterResidueBond]
```

**设计决策**：

| 决策 | 选择 | 理由 |
|:----|:----|:----|
| 无分子列表 | 删除 | 冗余，可从 `residue.molecule_name` 推导 |
| 有体系名称 | `system_name` | 用于输出/日志 |
| 有残基间键 | `inter_residue_bonds` | 骨架连接，未来扩展 |

---

## 4. 索引规则

### 4.1 索引清单

| 字段 | 所属数据类 | 范围 | 用途 |
|:----|:----|:----|:----|
| `atom_global_idx` | AtomData | 整个体系 | 创建 Group、计算距离 |
| `atom_idx_in_residue` | AtomData | 残基内 | 残基内键连接 |
| `atom1_idx_in_residue` | BondData | 残基内 | 键的原子 1 |
| `atom2_idx_in_residue` | BondData | 残基内 | 键的原子 2 |
| `residue_global_idx` | ResidueData | 整个体系 | 标识残基、跨残基键 |
| `residue_idx_in_molecule` | ResidueData | 分子内 | PDB 残基编号 |
| `atom_idx_in_residue1` | InterResidueBond | 残基内 | 跨残基键的原子 1 |
| `atom_idx_in_residue2` | InterResidueBond | 残基内 | 跨残基键的原子 2 |

### 4.2 命名原则

- **全局索引**：带 `global_idx` 后缀，明确标识范围
- **局部索引**：带 `atom_idx_in_residue` 或 `residue_idx_in_molecule`，明确标识相对位置
- **字段前缀**：`atom_` 或 `residue_`，避免混淆

---

## 5. 使用示例

### 5.1 完整示例

```python
SystemData(
    system_name="RBD_D927_KRAS",
    residues=[
        # TYR165（蛋白残基）
        ResidueData(
            residue_name="TYR",
            residue_global_idx=0,
            residue_idx_in_molecule=165,
            molecule_name="RBD_pro",
            atoms=[
                AtomData(atom_global_idx=0, atom_idx_in_residue=0,
                         atom_name="N", atom_type="N",
                         atom_element="N", atom_charge=-0.416,
                         atom_mass=14.007),
                AtomData(atom_global_idx=1, atom_idx_in_residue=1,
                         atom_name="CA", atom_type="CX",
                         atom_element="C", atom_charge=0.023,
                         atom_mass=12.011),
                # ...
            ],
            bonds=[
                BondData(atom1_idx_in_residue=0, atom2_idx_in_residue=1, bond_type="single"),
                BondData(atom1_idx_in_residue=0, atom2_idx_in_residue=20, bond_type="constrained"),
                # ...
            ],
        ),
        # D927（配体残基）
        ResidueData(
            residue_name="D927",
            residue_global_idx=1,
            residue_idx_in_molecule=901,
            molecule_name="D927",
            atoms=[
                AtomData(atom_global_idx=2355, atom_idx_in_residue=0,
                         atom_name="C3", atom_type="c3",
                         atom_element="C", atom_charge=-0.46,
                         atom_mass=12.011),
                # ...
            ],
            bonds=[
                BondData(atom1_idx_in_residue=0, atom2_idx_in_residue=1, bond_type="single"),
                # ...
            ],
        ),
        # ...
    ],
    inter_residue_bonds=[
        # 肽键：TYR165 的 C 与 ASP166 的 N 相连
        InterResidueBond(
            residue1_global_idx=0,
            atom_idx_in_residue1=5,
            residue2_global_idx=1,
            atom_idx_in_residue2=0,
            bond_type="single"
        ),
        # 二硫键：CYS42 的 SG 与 CYS58 的 SG 相连
        InterResidueBond(
            residue1_global_idx=10,
            atom_idx_in_residue1=6,
            residue2_global_idx=25,
            atom_idx_in_residue2=6,
            bond_type="single"
        ),
    ],
)
```

---

## 6. 数据层次图

```
SystemData
├── system_name: "RBD_D927_KRAS"
├── residues: List[ResidueData]
│     ├── ResidueData (TYR165)
│     │     ├── residue_name: "TYR"
│     │     ├── residue_global_idx: 0
│     │     ├── residue_idx_in_molecule: 165
│     │     ├── molecule_name: "RBD_pro"
│     │     ├── atoms: List[AtomData]
│     │     │     ├── AtomData (N)
│     │     │     │     ├── atom_global_idx: 0
│     │     │     │     ├── atom_idx_in_residue: 0
│     │     │     │     ├── atom_name: "N"
│     │     │     │     ├── atom_type: "N"
│     │     │     │     ├── atom_element: "N"
│     │     │     │     ├── atom_charge: -0.416
│     │     │     │     └── atom_mass: 14.007
│     │     │     ├── AtomData (CA)
│     │     │     └── ...
│     │     └── bonds: List[BondData]
│     │           ├── BondData (atom1_idx_in_residue=0, atom2_idx_in_residue=1, "single")
│     │           ├── BondData (atom1_idx_in_residue=0, atom2_idx_in_residue=20, "constrained")
│     │           └── ...
│     ├── ResidueData (D927)
│     └── ...
└── inter_residue_bonds: List[InterResidueBond]
      ├── InterResidueBond (TYR165.C → ASP166.N)
      └── ...
```

---

## 7. 设计原则总结

| 原则 | 体现 |
|:----|:----|
| **以残基为中心** | 残基是相互作用分析的自然单位 |
| **索引清晰** | 全局/局部索引命名明确，避免混淆 |
| **自包含** | 残基内信息完整，不依赖外部数据 |
| **可扩展** | 支持残基间键，未来可扩展分子列表 |
| **KISS** | 只存必要信息，不冗余 |

---

## 8. 依赖关系

```
core/data.py          ← 定义 AtomData, BondData, ResidueData, InterResidueBond, SystemData
      ↑
identifiers/amber.py  ← 实现 TprReader，产出 SystemData
      ↓
identifiers/amber.py  ← 实现 AmberGroupIdentifier，消费 SystemData
```

---

## 9. 扩展性设计

| 扩展场景 | 实现方式 |
|:----|:----|
| **新输入格式** | 实现新的 Reader，产出 SystemData |
| **新力场** | 实现新的 Identifier，消费 SystemData |
| **新键类型** | 在 `bond_type` 中添加新值 |
| **分子级别操作** | 从 `residue.molecule_name` 推导分子列表 |

---

*文档结束*
