# Interaction 相互作用数据结构设计文档

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

### 1.2 Interaction 的角色

Interaction 是"判定器"的输出，是整个系统的最终数据产品：

```
Group 对象（静态）
       ↓
detectors/hydrogen_bond.py  →  产出 Interaction 对象
       ↓
统计/可视化/输出
```

---

## 2. 第一性原理分析

### 2.1 相互作用是什么？

**相互作用 = 一组基团在 MD 轨迹中满足某种条件**

关键特征：
- **基团组**：不一定是 2 个（水桥涉及 3 个）
- **时间序列**：每帧的判定结果
- **几何指标**：距离、角度、偏移量等

### 2.2 数据规模问题

**问题**：百万帧 × 千个相互作用 = 数据爆炸

**错误方案**：每帧一个对象
```
1000 个相互作用 × 1000000 帧 = 10 亿个对象  ❌ 内存爆炸
```

**正确方案**：一个相互作用 = 一个对象 + 时间序列数组
```
1000 个相互作用 × 1 个对象 = 1000 个对象  ✅
每个对象包含 1000000 帧的 numpy 数组
```

### 2.3 灵活性需求

不同相互作用类型需要不同的几何指标：

| 类型 | 基团数量 | 几何指标 |
|:----|:----|:----|
| H-bond | 2 | distance, angle |
| π-π 堆积 | 2 | distance, angle, offset |
| 盐桥 | 2 | distance |
| 疏水 | 2 | distance |
| 水桥 | 3 | distance_donor_water, distance_water_acceptor, angle |

**结论**：基类需要灵活支持任意数量的基团和指标。

---

## 3. 最终设计

### 3.1 基类定义（初始版本）

```python
from dataclasses import dataclass, field
from typing import Tuple, Dict
import numpy as np

@dataclass
class Interaction:
    """一组基团之间的相互作用记录（全帧）。

    Attributes:
        interaction_type: 相互作用类型
        groups: 参与的基团元组（2个或更多）
        existence: bool 数组，每帧是否存在
        metrics: 几何指标字典，键为指标名，值为 numpy 数组
    """

    # === 身份信息 ===
    interaction_type: str
    groups: Tuple[Group, ...]

    # === 时间序列数据 ===
    existence: np.ndarray
    metrics: Dict[str, np.ndarray]

    def __post_init__(self):
        """验证数据完整性。"""
        n_frames = len(self.existence)
        for name, arr in self.metrics.items():
            if len(arr) != n_frames:
                raise ValueError(
                    f"{name} length {len(arr)} != existence length {n_frames}"
                )

    @property
    def n_frames(self) -> int:
        """帧数。"""
        return len(self.existence)

    @property
    def occupancy(self) -> float:
        """存在比例。"""
        return np.sum(self.existence) / self.n_frames

    def __repr__(self) -> str:
        groups_str = ", ".join(g.group_type for g in self.groups)
        return (f"Interaction(type='{self.interaction_type}', "
                f"groups=({groups_str}), frames={self.n_frames})")
```

### 3.1.1 基类定义（当前版本，2026-08-28 更新）

> **演进说明**：初始版本中一个 Interaction 对象 = **一对**基团的全帧数据（1D 数组）。
> 后改为矩阵存储：一个 Interaction 对象 = 该类型下**所有**基团对的全帧数据（2D 数组）。
>
> **改动原因**：
> 1. numpy 批量操作——`np.sum(existence, axis=1)` 一次算出所有对的 occupancy
> 2. 内存连续——一个大矩阵比 N 个小数组布局更紧凑
> 3. 构建简单——检测器内部用 numpy 预分配矩阵直接填入
>
> **附带影响**：初始版本设想的 `HydrogenBond`、`PiStacking`、`WaterBridge` 子类
> （有 `donor`、`acceptor` 等便捷属性）在矩阵存储方案下不再适用——一个 Interaction
> 包含多对，无法说 `.donor` 是谁。这些子类未实现，当前也不计划实现。

```python
from dataclasses import dataclass
from typing import List, Tuple, Dict
import numpy as np

@dataclass
class Interaction:
    """一种相互作用类型的全部检测结果。

    按类型组织：一个 Interaction 对象包含该类型下所有基团对在全部帧上的结果。
    groups[i] 对应 existence[i] 和 metrics 中各数组的第 i 行。

    Attributes:
        interaction_type: 相互作用类型（如 "salt_bridge", "hydrogen_bond"）
        groups: 基团对列表，每对是一个 tuple
        existence: (n_pairs, n_frames) bool 数组
        metrics: 几何指标字典，值为 (n_pairs, n_frames) 数组
    """

    interaction_type: str
    groups: List[Tuple[Group, ...]]
    existence: np.ndarray
    metrics: Dict[str, np.ndarray]

    @property
    def n_pairs(self) -> int:
        """基团对数量。"""
        return len(self.groups)

    @property
    def n_frames(self) -> int:
        """帧数。"""
        return self.existence.shape[1]

    def occupancy(self) -> np.ndarray:
        """每对基团的存在比例，shape=(n_pairs,)。"""
        return np.sum(self.existence, axis=1) / self.n_frames
```

**两版本对比**：

| 维度 | 初始版本 | 当前版本 |
|:---|:---|:---|
| groups 类型 | `Tuple[Group, ...]` | `List[Tuple[Group, ...]]` |
| existence shape | `(n_frames,)` | `(n_pairs, n_frames)` |
| metrics shape | `(n_frames,)` | `(n_pairs, n_frames)` |
| occupancy | `float` 属性 | `np.ndarray` 方法 |
| 一个对象包含 | 1 对基团 | 同类型所有基团对 |
| 子类（HydrogenBond 等） | 有意义 | 不适用（未实现） |

### 3.2 设计原则

| 原则 | 体现 |
|:----|:----|
| **KISS** | 只存必要信息，统计量计算得到 |
| **显式** | 指标名清晰表达含义 |
| **灵活** | 支持任意数量的基团和指标 |

### 3.3 设计决策

| 决策 | 选择 | 理由 |
|:----|:----|:----|
| **groups 用 Tuple** | 不可变 | 交互定义后不应改变 |
| **metrics 用 Dict** | 灵活 | 不同类型需要不同指标 |
| **预存 occupancy** | 计算属性 | 高频使用，避免重复计算 |
| **existence 为 bool** | 节省内存 | numpy bool 数组紧凑 |

---

## 4. 子类设计

### 4.1 子类的作用（初始版本）

子类为特定相互作用类型提供：
1. **便捷构造**：明确的参数名（`donor`, `acceptor`）
2. **便捷访问**：`hbond.donor` 比 `hbond.groups[0]` 更清晰

### 4.2 子类定义位置

子类定义在 `core/data.py` 中，与基类一起，**不违反依赖方向**。

### 4.3 子类示例（初始版本，未实现）

> **2026-08-28 说明**：以下子类在矩阵存储方案下不再适用——一个 Interaction 包含
> 多对基团，无法定义 `.donor`、`.acceptor` 等单对属性。这些子类未实现，当前也不
> 计划实现。保留此节作为设计演进的历史记录。

```python
@dataclass
class HydrogenBond(Interaction):
    """氢键相互作用。"""

    @classmethod
    def create(cls, donor: Group, acceptor: Group,
               existence: np.ndarray, distance: np.ndarray,
               angle: np.ndarray) -> "HydrogenBond":
        """创建氢键实例。"""
        return cls(
            interaction_type="hydrogen_bond",
            groups=(donor, acceptor),
            existence=existence,
            metrics={"distance": distance, "angle": angle}
        )

    @property
    def donor(self) -> Group:
        return self.groups[0]

    @property
    def acceptor(self) -> Group:
        return self.groups[1]


@dataclass
class PiStacking(Interaction):
    """π-π 堆积相互作用。"""

    @classmethod
    def create(cls, ring1: Group, ring2: Group,
               existence: np.ndarray, distance: np.ndarray,
               angle: np.ndarray, offset: np.ndarray) -> "PiStacking":
        """创建 π-π 堆积实例。"""
        return cls(
            interaction_type="pi_stacking",
            groups=(ring1, ring2),
            existence=existence,
            metrics={"distance": distance, "angle": angle, "offset": offset}
        )


@dataclass
class WaterBridge(Interaction):
    """水桥相互作用。"""

    @classmethod
    def create(cls, water: Group, donor: Group, acceptor: Group,
               existence: np.ndarray, distance_donor_water: np.ndarray,
               distance_water_acceptor: np.ndarray,
               angle: np.ndarray) -> "WaterBridge":
        """创建水桥实例。"""
        return cls(
            interaction_type="water_bridge",
            groups=(water, donor, acceptor),
            existence=existence,
            metrics={
                "distance_donor_water": distance_donor_water,
                "distance_water_acceptor": distance_water_acceptor,
                "angle": angle
            }
        )
```

---

## 5. 使用示例

### 5.1 氢键（初始版本，单对子类）

```python
hbond = HydrogenBond.create(
    donor=donor_group,
    acceptor=acceptor_group,
    existence=np.array([True, False, True, ...]),
    distance=np.array([0.28, 0.35, 0.27, ...]),
    angle=np.array([165.0, 0.0, 168.0, ...]),
)

# 便捷访问
print(hbond.donor)
print(hbond.acceptor)
print(hbond.occupancy)  # 预计算属性

# 计算其他统计量
dist_mean = np.mean(hbond.metrics["distance"][hbond.existence])
```

### 5.1.1 氢键（当前版本，矩阵存储）

```python
# 由检测器直接产出
from DuIvyInteractions.interaction_detectors import HydrogenBondDetectorTwoPass

detector = HydrogenBondDetectorTwoPass()
results = detector.detect(groups, trajectory)
hbond = results[0]  # 一个 Interaction 对象包含所有氢键对

# 访问
print(hbond.interaction_type)  # "hydrogen_bond"
print(hbond.n_pairs)           # 基团对数量
print(hbond.n_frames)          # 帧数

# 第 i 个基团对
pair = hbond.groups[i]         # (donor_group, acceptor_group)
occ = hbond.occupancy()        # (n_pairs,) 每对的存在比例

# 第 i 个基团对的全帧距离
dist_i = hbond.metrics["distance"][i]  # (n_frames,)

# 第 i 个基团对在第 f 帧的距离
dist_if = hbond.metrics["distance"][i, f]
```

### 5.2 π-π 堆积（初始版本）

```python
pi_stack = PiStacking.create(
    ring1=ring1,
    ring2=ring2,
    existence=np.array([True, True, False, ...]),
    distance=np.array([0.38, 0.39, 0.45, ...]),
    angle=np.array([15.0, 12.0, 25.0, ...]),
    offset=np.array([0.12, 0.10, 0.18, ...]),
)
```

### 5.3 水桥（初始版本）

```python
water_bridge = WaterBridge.create(
    water=water_group,
    donor=donor_group,
    acceptor=acceptor_group,
    existence=np.array([True, False, True, ...]),
    distance_donor_water=np.array([0.28, 0.35, 0.27, ...]),
    distance_water_acceptor=np.array([0.29, 0.36, 0.28, ...]),
    angle=np.array([160.0, 0.0, 165.0, ...]),
)
```

---

## 6. 内存估算

假设 1000 个相互作用，100 万帧（1μs，2fs 步长）：

| 数据 | 类型 | 每个相互作用 | 1000 个总计 |
|:----|:----|:----|:----|
| existence | bool | 1MB | 1GB |
| distance | float64 | 8MB | 8GB |
| angle | float64 | 8MB | 8GB |

**总计约 17GB**（2 个指标）

可通过以下方式优化：
1. **稀疏存储**：只存储 existence=True 的帧
2. **降精度**：float32 替代 float64
3. **分批处理**：不一次性加载全部帧

---

## 7. 依赖关系

### 初始版本
```
core/data.py          ← 定义 Interaction 基类 + 子类
      ↑
detectors/*.py        ← 产出 Interaction 对象
      ↓
utils/output.py       ← 输出 Interaction 对象
visualize/plotter.py  ← 可视化 Interaction 对象
```

### 当前版本（2026-08-28 更新）
```
core/datas.py                              ← 定义 Interaction + InteractionSparse
      ↑
interaction_detectors/*_per_tuple.py       ← 产出 Interaction（策略一）
interaction_detectors/*_per_frame.py       ← 产出 Interaction（策略二）
interaction_detectors/*_two_pass.py        ← 产出 InteractionSparse → Interaction（策略三）
      ↓
utils/output.py（待实现）                  ← 输出 Interaction 对象
visualizers/plotter.py（待实现）           ← 可视化 Interaction 对象
```

> **InteractionSparse**：TwoPass 策略的中间产物，Pass1 输出、Pass2 消费。
> 以 `(group_id1, group_id2, ...)` 为键，只存储 existence=True 的帧，
> 长轨迹场景下比稠密存储节省 99%+ 内存。定义在 `core/datas.py` 中。

---

## 8. 扩展性设计

| 扩展场景 | 实现方式 |
|:----|:----|
| **新相互作用类型** | 继承 Interaction，定义子类 |
| **新几何指标** | 在子类的 `create` 方法中添加参数 |
| **新基团组合** | 在子类的 `create` 方法中定义参数 |
| **新输出格式** | 在 utils/output.py 添加新函数 |

---

## 9. 文件组织

```
DuIvyInteractions/
├── core/
│   ├── data.py              # Group + Interaction 基类 + 子类
│   ├── constants.py         # GROUP_TYPES 等常量
│   └── interfaces.py        # ABC 接口
├── detectors/
│   ├── hydrogen_bond.py     # HydrogenBondDetector
│   ├── pi_stacking.py       # PiStackingDetector
│   └── ...
└── ...
```

---

*文档结束*
