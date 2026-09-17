# 扩展指南

本文档说明如何为 DuIvyInteractions 添加新力场、新相互作用类型或新判据。面向集成者/开发者。

## 架构总览

工具采用**策略模式**：Reader → GroupIdentifier → InteractionDetector → io（序列化/导出）四层，各层可插拔。

```
输入(tpr+xtc) → Reader → SystemData → GroupIdentifier → List[Group]
             → Detector（per_tuple/per_frame/two_pass）→ List[Interaction]
             → io/h5 序列化 → io/exporter 导出 xvg/xpm/csv
```

三策略检测器接口一致（`detect()` 返回统一矩阵式 `List[Interaction]`），由 `Pipeline` 的 `STRATEGY_INDEX` 切换，见 `core/interfaces.py`。

## 添加新力场

**原理**：基团判定基于特征空间映射（类型→{杂化, 芳香性, 极性, 带H, 孤对}），而非类型名硬编码。新力场只需填特征表，无需重写逻辑。

### 步骤

1. 在 `group_identifiers/amber_ff_identifier.py` 顶部的特征集合中补充该力场的类型：
   - `STRONG_AROMATIC`：芳香类型
   - `COMPATIBLE_TYPES`：兼容类型（n-1 个芳香原子强制下参与共轭）
   - `ACCEPTOR_TYPES`：受体类型（有孤对）
   - `METAL_IONS`：金属元素（若新增金属）
   - `WATER_RESIDUES`：水分子残基名（若不同）
2. 在 `group_identifiers/__init__.py` 的 `IDENTIFIER_CLASSES` 注册新识别器类
3. `dii run` 的 `--ff` 参数即可选择新力场

**已验证兼容**：amber03/94/96/99/99sb/99sb-ildn/GS/14sb + GAFF/GAFF2。类型名跨版本差异（如 CX vs CT）已处理。

## 添加新相互作用类型

### 步骤

1. **实现检测器**：继承 `InteractionDetectorTwoPass`（推荐，性能最优）或 `PerFrame`/`PerTuple`，实现必需抽象方法：
   - `name`：类型名（如 `"pi_cation"`）
   - `required_group_types`：需要的基团类型（Pipeline 据此过滤）
   - `metric_names`：指标名列表
   - `initialize_candidates` / `compute_pair_metrics` / `apply_threshold`（TwoPass）
2. **注册到 Pipeline**：在 `pipeline.py` 的 `DETECTOR_CLASSES` 添加 `(TwoPass, PerFrame, PerTuple)` 三元组
3. **实现导出器**：继承 `io/interaction_exporter.py` 的 `InteractionExporter`，实现 `name`/`metric_labels`/`get_pair_label`，在 `io/__init__.py` 注册
4. **注册到 DII**：在 `DII.py` 的 `ALL_INTERACTIONS` 与 `exporter_classes` 添加
5. **补充测试**：新建 `Tests/unittests/test_<type>_*.py` 单测

## 添加新判据（同一类型多判据）

同一类型可有多个 Detector（如 `HBondStrict`, `HBondLoose`），只需继承基类并覆盖 `apply_threshold`。判据阈值是模块级常量，调参只需改常量。

## 数据格式扩展

### 新增指标（metric）

- 检测器 `metrics` 字典加入新键
- 导出器 `metric_labels` 加对应标签（数值型自动导出 xvg；字符串型如 `pistacking_type` 自动跳过 xvg，但可扩展专属 xpm）

### 基团 metadata

`Group.metadata` 为键值字典，键须为字符串，值限 JSON 可序列化类型（numpy 标量/数组自动兜底）。识别器可用它携带额外化学信息。

## 测试与验证

- 运行单测：`pytest Tests/unittests/`
- 真实数据验证：`Tests/test_MD_case/`（KRAS-RBD 体系 1 ns 轨迹）跑 `dii run` + `dii export` 端到端
- 结果校验：h5 往返无损（`test_io_h5.py`）、XPM 值语义一致（`test_exporters.py`）