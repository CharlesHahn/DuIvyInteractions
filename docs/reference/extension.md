# 扩展指南

本文档说明如何为 DuIvyInteractions 添加新力场、新相互作用类型或新判据。面向集成者/开发者。所有接口以 `DuIvyInteractions/core/interfaces.py` 与 `pipeline.py` 的注册表为准。

## 架构总览

工具采用**策略模式**：Reader → GroupIdentifier → InteractionDetector → io（序列化/导出）四层，各层可插拔：

```
输入(tpr+xtc) → Reader → SystemData → GroupIdentifier → List[Group]
             → Detector（per_tuple/per_frame/two_pass）→ List[Interaction]
             → io/h5 序列化 → io/exporter 导出 xvg/xpm/csv
```

三策略检测器接口一致（`detect()` 返回统一矩阵式 `List[Interaction]`），由 `Pipeline` 的 `STRATEGY_INDEX = {"two_pass": 0, "per_frame": 1, "per_tuple": 2}` 切换，检测器注册表为 `pipeline.DETECTOR_CLASSES`。

## 添加新力场

**原理**：基团判定基于特征空间映射（类型→{杂化, 芳香性, 极性, 带H, 孤对}）+ 结构判据，而非类型名硬编码。现有 `GromosFFGroupIdentifier` / `CharmmFFGroupIdentifier` / `OplsFFGroupIdentifier` 均**子类化 `AmberFFGroupIdentifier`** 复用力场无关逻辑（供体、卤键、金属、水、金属配位、带电官能团），只替换力场相关部分。

### 步骤

1. **子类化基类**：新识别器继承 `GroupIdentifier`（全量重写）或 `AmberFFGroupIdentifier`（推荐，复用 A 类逻辑），实现 `name` 属性与（如继承时）必要的覆写。
2. **在新识别器模块顶部定义类型特征表**（模块级 `frozenset`/`dict` 常量，**不要**改动 `amber_ff_identifier.py` 内的表），参考现有实现：
   - `XXX_ACCEPTOR_TYPES`：H 键受体类型（有孤对；各力场已按化学事实逐类型判定并剔除铵/酰胺/带 H 吡咯/胍基 N）
   - `XXX_STRONG_AROMATIC` / `XXX_COMPATIBLE_TYPES`：芳香类型 / 在 n-1 个芳香原子强制下参与共轭的兼容类型
   - 芳香环强信号/残基白名单（如 GROMOS 的 `GROMOS_RING_TYPES`/`GROMOS_AROMATIC_RESIDUES` 分级判定）
   - 疏水类型（如 GROMOS 的 `GROMOS_HYDROPHOBIC_TYPES` 反列举 + `GROMOS_HYDROPHOBIC_EXCLUDED` 极性排除）
   - `METAL_IONS`（若新增金属）、蛋白正/负电残基字典（残基名 → 原子名列表）
3. **声明 `WATER_RESIDUES` 类属性**：必须覆盖为非空集合（基类默认空 `frozenset()`），pipeline 与识别器内部 `_find_water`/`_find_metal_binding` 统一经 `self.WATER_RESIDUES` 排除水。参照现有实现：amber=`{SOL,HOH,WAT}`、gromos=`{SOL}`、charmm=`{TIP3,HOH,SOL,WAT}`、opls=`{HOH,HO4,HO5,SOL,WAT}`。
4. **类型粒度粗时覆写结构判据**：若类型名无法区分化学环境（如 GROMOS 的 `N` 通吃主链酰胺/Pro/芳香 N），需覆写 `_find_acceptors` 增加结构判据——GROMOS 判据为：键数 ≥ 4（铵）或带 H 邻居（普通酰胺/带 H 吡咯/胍基）→ 非受体；无 H（Pro N / His 吡啶型 N）→ 受体（设计依据见 `doc/force_field_compatibility_survey.md` 与 `doc/gromos_identifier_design.md`）。CHARMM/OPLS 受体表已按类型逐一判定（无二义 N 类型），仅换表即可，无需结构判据。供体（`_classify_dh_pair`：D=N/O/S/F 且 q(H)>0）与疏水判定依赖显式 H：**联合原子力场**（无脂肪族显式 H）须按类型白名单/反列举处理（参照 GROMOS）。
5. **注册到注册表**：在 `group_identifiers/__init__.py` 的 `IDENTIFIER_CLASSES` 添加 `"力场名": 识别器类`。
6. **CLI 自动生效**：`dii run --ff` 的 choices 直接来自 `IDENTIFIER_CLASSES`（`DII.py`），无需改动 CLI。

### 已验证力场

- **Amber 家族**：amber03/94/96/99/99SB/99SB-ildn/GS/14SB + GAFF/GAFF2（类型名跨版本差异已处理，映射验证零冲突）
- **GROMOS**：53A6/54A7（联合原子，极性 H 显式；芳香分级判定 + 疏水反列举，见 `doc/gromos_identifier_design.md`）。**边界**：GROMOS 配体不支持（需 ATB 自动拓扑参数化，已声明）；SPC/SPC-E 水模型（残基名 `SOL`）
- **CHARMM**：CHARMM36/C36m + CGenFF 配体（与 CHARMM-GUI 官方 rtf 逐字一致，见 `doc/charmm_identifier_design.md`）；TIP3/HOH 水模型
- **OPLS**：OPLS-AA/L（2001，见 `doc/opls_identifier_design.md`）；HOH/SPC、HO4/TIP4P、HO5/TIP5P 水模型

新力场入库前建议先跑一遍对应力场的真实数据端到端验证（见文末"测试与验证"一节）。

## 添加新相互作用类型

### 步骤

1. **实现检测器**：继承 `InteractionDetectorTwoPass`（推荐，性能最优）或 `PerFrame`/`PerTuple`，实现必需方法：
   - `name`：类型名（如 `"pi_cation"`）
   - `required_group_types`：需要的基团类型（Pipeline 据此过滤）
   - `metric_names`：指标名列表
   - TwoPass：`initialize_candidates` / `compute_pair_metrics` / `apply_threshold`；PerTuple/PerFrame 对应 `get_candidate_tuples` / `compute_metrics`（PerFrame 为 `compute_metrics_for_frame`）/ `apply_threshold`
2. **注册到 Pipeline**：在 `pipeline.py` 的 `DETECTOR_CLASSES` 添加 `"类型名": (TwoPass类, PerFrame类, PerTuple类)` 三元组，并把类型名加入 `ALL_INTERACTIONS` 元组（`dii run --interactions` 的合法值来源）
3. **实现导出器**：继承 `io/interaction_exporter.py` 的 `InteractionExporter`，实现抽象属性 `name`/`metric_labels`（`get_pair_label` 已有默认实现，可覆写），并在 `io/__init__.py` 导出
4. **注册到 DII**：在 `DII.py` 的 `_run_export` 中 `exporter_classes` 字典添加 `"类型名": 导出器类`
5. **补充测试**：新建 `Tests/unittests/test_<type>_per_frame.py` 与 `test_<type>_two_pass.py`（按项目测试约定，策略一 per_tuple 的 `test_<type>.py` 处于"可能被舍弃"状态，不跑，见下）

#### 检测器代码骨架（TwoPass）

```python
from ..core.interfaces import InteractionDetectorTwoPass

# 判据阈值：模块级常量（与 PLIP 一致）
MY_DIST_MAX = 4.0  # Å


class MyInteractionDetectorTwoPass(InteractionDetectorTwoPass):
    """新相互作用检测器（TwoPass 策略）。"""

    @property
    def name(self) -> str:
        return "my_interaction"

    @property
    def required_group_types(self) -> list[str]:
        return ["H_donor", "H_acceptor"]

    @property
    def metric_names(self) -> list[str]:
        return ["distance"]

    def initialize_candidates(self, groups, trajectory, tuple_filter=None):
        # Pass1 前：生成候选对（可在此预筛或直接覆写 run_pass1 实现 KDTree 等）
        return super().initialize_candidates(groups, trajectory, tuple_filter)

    def compute_pair_metrics(self, group_tuples, all_positions):
        # 对候选对计算当前帧指标 -> {"distance": (n_groups,)}
        ...

    def apply_threshold(self, metrics):
        # 按判据判定本帧哪些 pair 存在 -> (n_groups,) bool
        return metrics["distance"] <= MY_DIST_MAX
```

#### 导出器代码骨架

```python
from typing import Dict
from .interaction_exporter import InteractionExporter


class MyInteractionExporter(InteractionExporter):
    """新相互作用导出器。"""

    @property
    def name(self) -> str:
        return "My Interaction"

    @property
    def metric_labels(self) -> Dict[str, str]:
        return {"distance": "Distance (Å)"}

    def get_pair_label(self, interaction, pair_idx: int) -> str:
        # 生成基团对标签（可选覆写；基类默认 "残基名残基号-残基名残基号"）
        g1, g2 = interaction.groups[pair_idx]
        return f"{g1.residue_name}{g1.residue_id}···{g2.residue_name}{g2.residue_id}"
```

## 添加新判据（同一类型多判据）

同一类型可有多个 Detector（如 `HBondStrict`, `HBondLoose`），只需继承基类并覆盖 `apply_threshold`。判据阈值是模块级常量，调参只需改常量；需重跑受影响类型的单测。

## 数据格式扩展

### 新增指标（metric）

- 检测器 `metrics` 字典（形状 `(n_pairs, n_frames)`）加入新键，并在 `metric_names` 列出
- 导出器 `metric_labels` 加对应标签。数值指标自动导出 xvg；字符串指标（dtype kind 为 `U`/`S`/`O`，如 `pistacking_type`）自动跳过 xvg，h5 中以"展平 1D 字符串 + attrs['shape']"存储（见[数据格式](data_format.md)）
- 如需专属可视化，可在导出器内新增方法（参照 `PiStackingExporter.save_xpm_stacking_type`），并在 `DII.py` 的导出流程中按类型调用

### 基团 metadata

`Group.metadata` 为键值字典，键须为字符串，值限 JSON 可序列化类型（numpy 标量/数组自动兜底）。识别器可用它携带额外化学信息（如 `_find_metal_binding` 的 `{"source": "element"}`）。

## 测试与验证

- 运行单测：按项目测试约定**显式文件列表**跑 `Tests/unittests/test_<type>_per_frame.py` 与 `test_<type>_two_pass.py`（策略一 per_tuple 的 `test_<type>.py` 处于"可能被舍弃"状态，不跑；不要 `pytest Tests/unittests/` 无差别跑）
- 真实数据验证（3 力场端到端）：
  - `Tests/test_MD_case_amber/`（Amber KRAS-RBD D927）
  - `Tests/test_MD_case_gromos/`（GROMOS 53A6 蛋白 + 6 配体）
  - `Tests/test_MD_case_charmm36/`（CHARMM36 SMO-BST）
- 结果校验：h5 往返无损（`test_io_h5.py`）、XPM 值语义一致（`test_exporters.py`）