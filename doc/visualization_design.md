# 可视化与数据分析设计文档

> 创建日期：2026-09-05
> 更新日期：2026-09-05
> 状态：**初步设计（待讨论）**
> 作者：dsh-tui+hanyl

---

## 1. 设计背景

### 1.1 项目现状

已完成的核心功能：
- ✅ 基团识别（AmberFFGroupIdentifier）
- ✅ 相互作用检测（8 种类型 × 3 种策略）
- ✅ 结果存储（HDF5 格式）

### 1.2 设计目标

基于已检测的相互作用数据，提供：
1. 数据查看与筛选
2. 基础统计分析
3. 数据可视化
4. 结构可视化
5. 数据导出

---

## 2. 数据结构回顾

### 2.1 核心数据结构

```python
@dataclass
class Interaction:
    interaction_type: str          # 相互作用类型
    groups: List[Tuple[Group, ...]]  # 基团对列表
    existence: np.ndarray          # (n_pairs, n_frames) bool
    metrics: Dict[str, np.ndarray] # {name: (n_pairs, n_frames)}
```

### 2.2 相互作用类型

| 类型 | metrics | 说明 |
|:-----|:--------|:-----|
| hydrogen_bond | distance, angle | 氢键 |
| pi_stacking | distance, angle, offset, pistacking_type | π-π 堆积 |
| salt_bridge | distance | 盐桥 |
| hydrophobic | distance | 疏水 |
| halogen_bond | distance, don_angle, acc_angle | 卤键 |
| metal_coordination | distance | 金属配位 |
| water_bridge | dist_dw, dist_wa, theta, omega | 水桥 |
| pi_cation | distance, offset | π-阳离子 |

### 2.3 Group 结构

```python
@dataclass
class Group:
    group_id: int
    group_type: str        # "aromatic_ring", "H_donor", ...
    molecule: str          # 所属分子名
    residue_name: str      # 残基名
    residue_id: int        # 残基号
    atoms: List[AtomData]  # 原子列表
    metadata: Dict         # 附加信息
```

---

## 3. 功能设计

### 3.1 数据查看与筛选

#### 筛选维度

| 维度 | 筛选项 | 说明 |
|:-----|:-------|:-----|
| **相互作用类型** | hydrogen_bond, pi_stacking, ... | 按类型筛选 |
| **Group 属性** | group_id, group_type | 按基团筛选 |
| **分子** | molecule | 按分子筛选（如 D927、RBD） |
| **残基** | residue_name, residue_id | 按残基筛选（如 ARG73） |
| **占位率** | min, max | 按占位率范围筛选 |
| **距离** | min, max | 按距离范围筛选（均值/最小值/最大值） |
| **角度** | min, max | 按角度范围筛选 |
| **帧范围** | start, end | 按时间范围筛选 |
| **存在帧数** | min, max | 按存在帧数筛选 |

#### 排序选项

| 排序项 | 方向 | 说明 |
|:-------|:-----|:-----|
| 占位率 | 升序/降序 | 按占位率排序 |
| 平均距离 | 升序/降序 | 按平均距离排序 |
| 存在帧数 | 升序/降序 | 按存在帧数排序 |
| 残基 ID | 升序/降序 | 按残基排序 |
| 分子名 | 升序/降序 | 按分子排序 |

#### 分组选项

| 分组项 | 说明 |
|:-------|:-----|
| 按相互作用类型 | 同类型的放在一起 |
| 按分子 | 同分子的放在一起 |
| 按残基 | 同残基的放在一起 |

---

### 3.2 基础统计

#### 单对统计

| 统计项 | 说明 |
|:-------|:-----|
| 占位率 | existence 中 True 的比例 |
| 存在帧数 | existence 中 True 的数量 |
| 总帧数 | existence 的总帧数 |
| 距离均值 | 活跃帧的距离均值 |
| 距离标准差 | 活跃帧的距离标准差 |
| 距离最小值 | 活跃帧的距离最小值 |
| 距离最大值 | 活跃帧的距离最大值 |
| 距离中位数 | 活跃帧的距离中位数 |
| 角度均值 | 活跃帧的角度均值（如有） |
| 角度标准差 | 活跃帧的角度标准差（如有） |
| 首次出现帧 | 首次存在的帧号 |
| 最后出现帧 | 最后存在的帧号 |
| 连续存在最长 | 最长连续存在帧数 |

#### 按类型统计

| 统计项 | 说明 |
|:-------|:-----|
| 相互作用对数 | 该类型有多少对 |
| 平均占位率 | 所有对的平均占位率 |
| 占位率分布 | 占位率的直方图 |
| 高占位率对数 | 占位率 > 阈值的对数 |
| 平均距离 | 所有对的平均距离 |

#### 按残基统计

| 统计项 | 说明 |
|:-------|:-----|
| 参与相互作用数 | 该残基参与的相互作用对数 |
| 参与类型数 | 参与几种类型 |
| 平均占位率 | 该残基所有相互作用的平均占位率 |
| 作为 donor/acceptor | 该残基在氢键中的角色 |

#### 按分子统计

| 统计项 | 说明 |
|:-------|:-----|
| 参与相互作用数 | 该分子参与的相互作用对数 |
| 参与残基数 | 涉及几个残基 |
| 平均占位率 | 该分子所有相互作用的平均占位率 |

#### 界面分析

| 统计项 | 说明 |
|:-------|:-----|
| 界面残基对 | 两个分子之间接触的残基对 |
| 界面面积 | 接触面积（近似） |
| 界面相互作用数 | 界面上的相互作用数量 |

---

### 3.3 数据可视化

**方案**：使用 [DuIvyTools](https://github.com/CharlesHahn/DuIvyTools) 进行可视化，本项目只需实现 **Interaction → xvg/xpm 格式转换**。

#### 设计思路

```
Interaction 数据 → 格式转换器 → xvg/xpm 文件 → DuIvyTools 可视化
```

**优点**：
1. 复用成熟工具，不需要重新实现可视化代码
2. xvg/xpm 是 GROMACS 标准格式，用户熟悉
3. DuIvyTools 有 30+ 命令，覆盖大部分需求

#### 需要实现的转换器

| 转换器 | 输入 | 输出 | 说明 |
|:-------|:-----|:-----|:-----|
| `to_xvg_timeline` | Interaction | .xvg | 时间线图数据（existence 随时间变化） |
| `to_xvg_distance` | Interaction | .xvg | 距离数据（metrics["distance"]） |
| `to_xvg_angle` | Interaction | .xvg | 角度数据（metrics["angle"]） |
| `to_xvg_occupancy` | Interaction | .xvg | 占位率数据 |
| `to_xpm_existence` | Interaction | .xpm | 存在矩阵（existence → 热力图） |
| `to_xpm_occupancy_matrix` | Interaction | .xpm | 占位率矩阵（残基×残基） |

#### DuIvyTools 可用命令

转换后可使用 DuIvyTools 的命令进行可视化：

```bash
# 时间线图
dit xvg_show -f timeline.xvg

# 距离分布
dit xvg_show_distribution -f distance.xvg

# 散点图
dit xvg_show_scatter -f distance_angle.xvg

# 热力图
dit xpm_show -f existence.xpm

# 比较图
dit xvg_compare -f file1.xvg file2.xvg

# 柱状图
dit xvg_ave_bar -f occupancy.xvg
```

#### 转换规则设计（待实现）

**xvg 格式**：
- 第一行：标题（@ title）
- 第二行：X 轴标签（@ xaxis label）
- 第三行：Y 轴标签（@ yaxis label）
- 后续行：数据（x y1 y2 ...）

**xpm 格式**：
- 头部：矩阵维度、标签等
- 数据：字符矩阵（每个字符代表一个值）

---

### 3.4 结构可视化

**方案**：实现 Interaction → PyMOL 脚本转换。

#### 需要实现的转换器

| 转换器 | 输入 | 输出 | 说明 |
|:-------|:-----|:-----|:-----|
| `to_pymol_interactions` | Interaction | .pml | 显示相互作用（虚线） |
| `to_pymol_distances` | Interaction | .pml | 显示距离标签 |
| `to_pymol_highlight` | Interaction | .pml | 高亮残基 |

#### PyMOL 脚本示例

```pymol
# 显示氢键
distance hbond_1, /RBD_pro/ARG/NE, /D927/O1, 4.1
color green, hbond_1

# 高亮残基
select interface, /RBD_pro/ARG+LYS+ASP
color yellow, interface
```

---

### 3.5 数据导出

#### 表格格式

| 格式 | 说明 |
|:-----|:-----|
| CSV | 逐对导出（interaction_type, group1, group2, occupancy, avg_distance, ...） |
| CSV | 逐帧导出（interaction_type, group1, group2, frame, existence, distance, ...） |

#### 图片格式

| 格式 | 说明 |
|:-----|:-----|
| PNG | 栅格图 |
| SVG | 矢量图 |
| PDF | 矢量图 |

#### 报告格式

| 格式 | 说明 |
|:-----|:-----|
| HTML | 交互式报告 |
| Markdown | 文本报告 |

---

## 4. 功能优先级

| 优先级 | 功能 | 理由 |
|:-------|:-----|:-----|
| **P0** | 占位率统计 | 最基本的分析 |
| **P0** | xvg 转换器（时间线、距离） | DuIvyTools 可视化 |
| **P0** | xpm 转换器（existence 矩阵） | DuIvyTools 热力图 |
| **P0** | CSV 导出 | 数据交换 |
| **P1** | PyMOL 脚本转换 | 结构展示 |
| **P1** | 残基接触频率统计 | 热点分析 |
| **P2** | 其他转换器 | 按需添加 |

---

## 5. 待讨论事项

### 5.1 转换器设计

- 转换器放在哪个目录？（建议：`DuIvyInteractions/converters/`）
- 转换器的接口设计（函数式 vs 面向对象）
- xvg/xpm 格式的详细规范

### 5.2 与 DuIvyTools 的集成

- 是否需要在 DuIvyTools 中添加新的命令？
- 如何处理 DuIvyTools 不支持的可视化类型？

### 5.3 数据处理

- 数据筛选在转换之前还是之后？
- 如何处理大数据量？

---

## 6. 后续步骤

1. 讨论并确认转换器设计
2. 实现 xvg 转换器
3. 实现 xpm 转换器
4. 实现 PyMOL 脚本转换器
5. 测试验证

---

*文档结束（初步设计，待讨论）*
