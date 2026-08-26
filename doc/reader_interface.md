# Reader 接口设计文档

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

### 1.2 Reader 的角色

Reader 是"数据解析器"的接口，定义了如何从文件中读取数据：

```
输入文件（tpr/prmtop/pdb）
       ↓
Reader.read()  →  SystemData
       ↓
GroupIdentifier 消费 SystemData
```

---

## 2. 第一性原理分析

### 2.1 Reader 是什么？

**Reader = 从文件中读取数据，转换为 SystemData**

关键特征：
- **输入**：文件路径（字符串）
- **输出**：SystemData（体系分子数据）
- **可插拔**：不同文件格式需要不同的 Reader

### 2.2 需要回答的问题

| 问题 | 用途 | 必要性 |
|:----|:----|:----|
| **Reader 叫什么？** | 日志、调试 | 必须 |
| **能读什么格式？** | 前置检查 | 必须 |
| **如何读取？** | 核心算法 | 必须 |

### 2.3 当前需求

| Reader | 输入格式 | 说明 | 优先级 |
|:----|:----|:----|:----|
| TprDumpReader | gmx dump 文本 | 从文本解析 | **P0** |
| TprReader | .tpr 二进制 | 直接解析 tpr | P1 |

---

## 3. 最终设计

### 3.1 接口定义

```python
from abc import ABC, abstractmethod

class Reader(ABC):
    """数据读取器接口。"""

    @property
    @abstractmethod
    def name(self) -> str:
        """读取器名称。"""
        ...

    @abstractmethod
    def read(self, source: str) -> SystemData:
        """从文件读取数据。

        Args:
            source: 文件路径

        Returns:
            SystemData 实例
        """
        ...
```

### 3.2 设计原则

| 原则 | 体现 |
|:----|:----|
| **KISS** | 只定义必要方法 |
| **显式** | 方法名清晰表达含义 |
| **灵活** | 支持多种文件格式 |

### 3.3 设计决策

| 决策 | 选择 | 理由 |
|:----|:----|:----|
| **输入类型** | `str` | 文件路径 |
| **返回类型** | `SystemData` | 统一数据层 |
| **方法数量** | 2 个 | `name` + `read` |

---

## 4. 子类设计

### 4.1 子类示例

```python
class TprDumpReader(Reader):
    """从 gmx dump 文本读取数据。"""

    @property
    def name(self) -> str:
        return "tpr_dump"

    def read(self, source: str) -> SystemData:
        """从 gmx dump 文本读取数据。

        Args:
            source: gmx dump 输出的文本文件路径

        Returns:
            SystemData 实例
        """
        # 1. 解析文本
        # 2. 构建 SystemData
        ...
```

### 4.2 子类定义位置

子类定义在 `identifiers/tpr_dump_reader.py` 中。

---

## 5. 使用示例

### 5.1 基本使用

```python
reader = TprDumpReader()
system_data = reader.read("dump_md_D927.tpr.txt")
print(f"读取到 {system_data.n_residues} 个残基")
```

### 5.2 Pipeline 中使用

```python
def run_pipeline(reader: Reader, identifier: GroupIdentifier,
                 file_path: str):
    # 1. 读取数据
    system_data = reader.read(file_path)
    
    # 2. 识别基团
    groups = identifier.identify(system_data)
    
    # 3. 后续处理
    ...
```

---

## 6. 依赖关系

```
core/interfaces.py          ← 定义 Reader 接口
      ↑
identifiers/tpr_dump_reader.py  ← 实现 TprDumpReader
      ↓
pipeline.py                 ← 消费 Reader
```

---

## 7. 扩展性设计

| 扩展场景 | 实现方式 |
|:----|:----|
| **新文件格式** | 继承 Reader，实现 `read()` |
| **新解析逻辑** | 在子类中实现 |
| **新输出格式** | 返回 SystemData，格式统一 |

---

*文档结束*
