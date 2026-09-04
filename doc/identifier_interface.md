# GroupIdentifier 接口设计文档

> 创建日期：2026-08-24
> 状态：已定稿
> 作者：dsh-tui+hanyl

---

## 1. 设计背景

### 1.1 项目目标

构建通用的分子相互作用分析工具，核心流程：

```
拓扑(tpr) → 基团识别 → 基团列表
轨迹(xtc) → 逐帧坐标 → 几何判定 → 相互作用列表
```

### 1.2 GroupIdentifier 的角色

GroupIdentifier 是"识别器"的接口，定义了如何从拓扑中识别基团：

```
拓扑数据（tpr/MDAnalysis Universe）
       ↓
GroupIdentifier.identify()  →  List[Group]
       ↓
detectors/*.py 消费 Group 对象
```

---

## 2. 第一性原理分析

### 2.1 识别器是什么？

**识别器 = 从拓扑数据中提取基团的算法**

关键特征：
- **输入**：拓扑数据（tpr 文件、MDAnalysis Universe、SMILES 等）
- **输出**：基团列表（List[Group]）
- **可插拔**：不同力场/分子需要不同的识别器

### 2.2 需要回答的问题

| 问题 | 用途 | 必要性 |
|:----|:----|:----|
| **识别器叫什么？** | 日志、调试 | 必须 |
| **识别器能识别哪些基团？** | 前置检查 | 必须 |
| **如何识别？** | 核心算法 | 必须 |

### 2.3 输入数据类型

| 识别器 | 输入数据 | 说明 |
|:----|:----|:----|
| AmberGroupIdentifier | tpr 文件路径 | 需要解析 tpr |
| RDKitGroupIdentifier | SMILES/Mol 对象 | 需要 RDKit |
| SMARTSGroupIdentifier | SMILES + SMARTS | 模式匹配 |

**结论**：输入类型不统一，需要泛化。

---

## 3. 最终设计

### 3.1 接口定义（初始版本）

```python
from abc import ABC, abstractmethod
from typing import List

class GroupIdentifier(ABC):
    """基团识别器接口。"""

    @property
    @abstractmethod
    def name(self) -> str:
        """识别器名称。"""
        ...

    @abstractmethod
    def identify(self, source: str) -> List[Group]:
        """从数据源识别基团。

        Args:
            source: 数据源（文件路径、SMILES 等）

        Returns:
            识别到的基团列表
        """
        ...
```

### 3.1.1 接口定义（当前版本，2026-08-28 更新）

> **演进说明**：初始版本中 `identify()` 直接接收文件路径字符串。后因引入 `Reader` 层
> （见 `system_data_design.md`），数据流变为 `文件 → Reader → SystemData → Identifier → List[Group]`，
> `identify()` 的输入从 `str` 改为 `SystemData`，实现读取与识别的解耦。
> 以下为 `core/interfaces.py` 中的实际定义。

```python
from abc import ABC, abstractmethod
from typing import List

class GroupIdentifier(ABC):
    """基团识别器接口。"""

    @property
    @abstractmethod
    def name(self) -> str:
        """识别器名称。"""
        ...

    @abstractmethod
    def identify(self, system_data: SystemData) -> List[Group]:
        """从 SystemData 识别基团。

        Args:
            system_data: 体系数据（由 Reader 产出）

        Returns:
            识别到的基团列表
        """
        ...
```

### 3.2 设计原则

| 原则 | 体现 |
|:----|:----|
| **KISS** | 只定义必要方法 |
| **显式** | 方法名清晰表达含义 |
| **灵活** | `source` 参数支持多种类型 |

### 3.3 设计决策

| 决策 | 选择 | 理由 |
|:----|:----|:----|
| **输入类型** | `str` | 统一为文件路径/SMILES |
| **返回类型** | `List[Group]` | 与 Group 数据结构对接 |
| **方法数量** | 2 个 | `name` + `identify` |

---

## 4. 子类设计

### 4.1 子类示例

```python
class AmberGroupIdentifier(GroupIdentifier):
    """Amber 力场基团识别器。"""

    @property
    def name(self) -> str:
        return "amber"

    def identify(self, source: str) -> List[Group]:
        """从 tpr 文件识别基团。

        Args:
            source: tpr 文件路径

        Returns:
            识别到的基团列表
        """
        # 1. 解析 tpr
        # 2. 特征映射
        # 3. 环检测
        # 4. 供体/受体鉴定
        # 5. 返回 Group 列表
        ...
```

### 4.2 子类定义位置

子类定义在 `identifiers/amber.py` 中，**不违反依赖方向**。

---

## 5. 使用示例

### 5.1 基本使用（初始版本，直接传文件路径）

```python
identifier = AmberGroupIdentifier()
groups = identifier.identify("md.tpr")
print(f"识别到 {len(groups)} 个基团")
```

### 5.1.1 基本使用（当前版本，通过 Reader 读取）

```python
from DuIvyInteractions.input_readers import GmxTprReader
from DuIvyInteractions.group_identifiers import AmberFFGroupIdentifier

reader = GmxTprReader()
system_data = reader.read("md.tpr")

identifier = AmberFFGroupIdentifier()
groups = identifier.identify(system_data)
print(f"识别到 {len(groups)} 个基团")
```

### 5.2 Pipeline 中使用

```python
def run_pipeline(identifier: GroupIdentifier, topology_path: str):
    groups = identifier.identify(topology_path)
    # 后续处理
    ...
```

---

## 6. 依赖关系

### 初始版本
```
core/interfaces.py        ← 定义 GroupIdentifier 接口
      ↑
identifiers/amber.py      ← 实现 AmberGroupIdentifier
      ↓
pipeline.py               ← 消费 GroupIdentifier
```

### 当前版本（2026-08-28 更新）
```
core/interfaces.py                  ← 定义 GroupIdentifier 接口
      ↑
group_identifiers/amber_ff_identifier.py  ← 实现 AmberFFGroupIdentifier
      ↓
pipeline.py（待实现）               ← 消费 GroupIdentifier
```

---

## 7. 扩展性设计

| 扩展场景 | 实现方式 |
|:----|:----|
| **新力场识别器** | 继承 GroupIdentifier，实现 `identify()` |
| **新输入类型** | 在 `identify()` 中处理 |
| **新输出格式** | 返回 List[Group]，格式统一 |

---

*文档结束*
