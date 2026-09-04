# InteractionDetector 接口设计文档

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

### 1.2 InteractionDetector 的角色

InteractionDetector 是"检测器"的接口，定义了如何检测相互作用：

```
Group 对象（静态）
       ↓
InteractionDetector.detect_frame()  →  List[Interaction]
       ↓
统计/可视化/输出
```

---

## 2. 第一性原理分析

### 2.1 检测器是什么？

**检测器 = 从基团和坐标中检测相互作用的算法**

关键特征：
- **输入**：基团列表 + 坐标 + 帧信息
- **输出**：相互作用列表（List[Interaction]）
- **逐帧**：每帧独立检测
- **可插拔**：不同相互作用类型需要不同的检测器

### 2.2 需要回答的问题

| 问题 | 用途 | 必要性 |
|:----|:----|:----|
| **检测器叫什么？** | 日志、调试 | 必须 |
| **需要哪些基团类型？** | 前置过滤 | 必须 |
| **如何检测？** | 核心算法 | 必须 |

### 2.3 输入输出关系

```
输入：
  - groups: List[Group]           # 所有基团
  - coordinates: np.ndarray       # 坐标数组
  - frame: int                    # 帧号
  - time_ps: float                # 时间（皮秒）

输出：
  - List[Interaction]             # 检测到的相互作用
```

---

## 3. 最终设计

### 3.1 接口定义（初始版本）

```python
from abc import ABC, abstractmethod
from typing import List

class InteractionDetector(ABC):
    """相互作用检测器接口。"""

    @property
    @abstractmethod
    def name(self) -> str:
        """检测器名称。"""
        ...

    @property
    @abstractmethod
    def required_group_types(self) -> List[str]:
        """需要的基团类型列表。"""
        ...

    @abstractmethod
    def detect_frame(self, groups: List[Group], coordinates: np.ndarray,
                     frame: int, time_ps: float) -> List[Interaction]:
        """检测单帧的相互作用。

        Args:
            groups: 所有基团列表
            coordinates: 坐标数组
            frame: 帧号
            time_ps: 时间（皮秒）

        Returns:
            检测到的相互作用列表
        """
        ...
```

### 3.1.1 接口定义（当前版本，2026-08-28 更新）

> **演进说明**：初始版本采用逐帧调用模式（`detect_frame`），由 Pipeline 外部控制帧遍历。
> 后因性能需求（候选对数量大时逐帧调用无法做跨帧向量化、候选预筛选），改为三种策略 ABC，
> 检测器内部控制帧遍历，一次调用处理全部帧。详见 `perframe_detection_strategies.md`。
>
> 以下为 `core/interfaces.py` 中的实际定义。三种 ABC 共享相同的公共接口 `detect()`，
> 便于 Pipeline 统一调用，子类只需实现不同的抽象方法。

**公共接口（三种 ABC 共享）**：

```python
def detect(self, groups: List[Group], trajectory=None,
           n_workers: int = 1,
           topology_path: str = None,
           trajectory_path: str = None,
           tuple_filter=None) -> List[Interaction]:
    """检测全部帧的相互作用。
    Args:
        groups: 基团列表（已按 required_group_types 过滤）
        trajectory: MDAnalysis 轨迹对象
        n_workers: 并行 worker 数（1=串行）
        topology_path: 拓扑文件路径（并行时必填）
        trajectory_path: 轨迹文件路径（并行时必填）
        tuple_filter: 可选的基团组过滤函数 (Tuple[Group,...]) -> bool
    """
```

**策略一：InteractionDetectorPerTuple**

```python
class InteractionDetectorPerTuple(ABC):
    """逐组遍历轨迹，向量化跨帧。适用于候选 tuple 数量少的场景。"""

    @property
    @abstractmethod
    def name(self) -> str: ...

    @property
    @abstractmethod
    def required_group_types(self) -> List[str]: ...

    @property
    @abstractmethod
    def metric_names(self) -> List[str]: ...

    @abstractmethod
    def get_candidate_tuples(self, groups, coordinates=None) -> List[Tuple[Group, ...]]: ...

    @abstractmethod
    def compute_metrics(self, group_tuple, coords) -> Dict[str, np.ndarray]: ...

    @abstractmethod
    def apply_threshold(self, metrics) -> np.ndarray: ...

    def filter_candidate_tuples(self, tuples, coordinates) -> List[Tuple[Group, ...]]:
        """可选覆盖：用一帧坐标预过滤候选组。"""
        return tuples
```

**策略二：InteractionDetectorPerFrame**

```python
class InteractionDetectorPerFrame(ABC):
    """逐帧向量化计算全部候选组。适用于候选组数量极大的场景（如水桥）。"""

    @property
    @abstractmethod
    def name(self) -> str: ...

    @property
    @abstractmethod
    def required_group_types(self) -> List[str]: ...

    @property
    @abstractmethod
    def metric_names(self) -> List[str]: ...

    @abstractmethod
    def get_candidate_tuples(self, groups, coordinates=None) -> List[Tuple[Group, ...]]: ...

    @abstractmethod
    def compute_metrics_for_frame(self, tuples, all_positions, frame) -> Dict[str, np.ndarray]: ...

    @abstractmethod
    def apply_threshold(self, metrics) -> np.ndarray: ...
```

**策略三：InteractionDetectorTwoPass**

```python
class InteractionDetectorTwoPass(ABC):
    """两轮遍历 + 稀疏存储。Pass1 逐帧发现 active pairs，Pass2 补全全帧数据。"""

    @property
    @abstractmethod
    def name(self) -> str: ...

    @property
    @abstractmethod
    def required_group_types(self) -> List[str]: ...

    @property
    @abstractmethod
    def metric_names(self) -> List[str]: ...

    def initialize_candidates(self, groups, trajectory, tuple_filter=None) -> list: ...

    def compute_pair_metrics(self, group_tuples, all_positions) -> Dict[str, np.ndarray]: ...

    def apply_threshold(self, metrics) -> np.ndarray: ...

    def run_pass1(self, groups, trajectory, tuple_filter=None) -> InteractionSparse: ...

    def run_pass2(self, sparse, trajectory) -> List[Interaction]: ...
```

### 3.2 设计原则

| 原则 | 体现 |
|:----|:----|
| **KISS** | 只定义必要方法 |
| **显式** | 方法名清晰表达含义 |
| **灵活** | 支持任意类型的相互作用 |

### 3.3 设计决策

| 决策 | 选择 | 理由 |
|:----|:----|:----|
| **输入类型** | List[Group] + np.ndarray | 统一接口 |
| **返回类型** | List[Interaction] | 与 Interaction 数据结构对接 |
| **方法数量** | 3 个 | `name` + `required_group_types` + `detect_frame` |

---

## 4. 子类设计

### 4.1 子类示例

```python
class HydrogenBondDetector(InteractionDetector):
    """氢键检测器。"""

    def __init__(self, distance_cutoff: float = 0.35, angle_cutoff: float = 150.0):
        self.distance_cutoff = distance_cutoff
        self.angle_cutoff = angle_cutoff

    @property
    def name(self) -> str:
        return "hydrogen_bond"

    @property
    def required_group_types(self) -> List[str]:
        return ["H_donor", "H_acceptor"]

    def detect_frame(self, groups: List[Group], coordinates: np.ndarray,
                     frame: int, time_ps: float) -> List[Interaction]:
        """检测单帧的氢键。"""
        # 1. 过滤出 donor 和 acceptor
        # 2. 计算距离和角度
        # 3. 判断是否满足条件
        # 4. 返回 Interaction 列表
        ...
```

### 4.2 子类定义位置

子类定义在 `detectors/hydrogen_bond.py` 中，**不违反依赖方向**。

---

## 5. 使用示例

### 5.1 基本使用（初始版本，逐帧调用）

```python
detector = HydrogenBondDetector(distance_cutoff=0.35, angle_cutoff=150.0)
interactions = detector.detect_frame(groups, coordinates, frame=100, time_ps=2000.0)
print(f"检测到 {len(interactions)} 个氢键")
```

### 5.1.1 基本使用（当前版本，一次调用处理全部帧）

```python
from DuIvyInteractions.interaction_detectors import HydrogenBondDetectorTwoPass

detector = HydrogenBondDetectorTwoPass()
interactions = detector.detect(groups, trajectory)
# interactions 是 List[Interaction]，每个 Interaction 包含该类型下所有基团对的全帧数据
```

### 5.2 Pipeline 中使用（初始版本）

```python
def run_pipeline(detectors: List[InteractionDetector], groups: List[Group],
                 trajectory):
    all_interactions = []
    for ts in trajectory:
        for detector in detectors:
            interactions = detector.detect_frame(groups, ts.positions, ts.frame, ts.time)
            all_interactions.extend(interactions)
    return all_interactions
```

### 5.2.1 Pipeline 中使用（当前版本）

```python
def run_pipeline(detectors, groups, trajectory):
    all_interactions = []
    for detector in detectors:
        interactions = detector.detect(groups, trajectory)
        all_interactions.extend(interactions)
    return all_interactions
```

---

## 6. 依赖关系

### 初始版本
```
core/interfaces.py        ← 定义 InteractionDetector 接口
      ↑
detectors/hydrogen_bond.py ← 实现 HydrogenBondDetector
      ↓
pipeline.py               ← 消费 InteractionDetector
```

### 当前版本（2026-08-28 更新）
```
core/interfaces.py                            ← 定义三种 ABC
      ↑
interaction_detectors/hydrogen_bond_detector_*.py  ← 实现（每种策略一个文件）
      ↓
pipeline.py（待实现）                         ← 消费检测器
```

---

## 7. 扩展性设计

| 扩展场景 | 实现方式 |
|:----|:----|
| **新相互作用类型** | 继承三种 ABC 之一，实现对应抽象方法 |
| **新检测参数** | 在子类的 `__init__` 中添加 |
| **新输出格式** | 返回 List[Interaction]，格式统一 |
| **同一类型不同策略** | 实现多个子类（如 HydrogenBondDetectorPerTuple / PerFrame / TwoPass） |

---

*文档结束*
