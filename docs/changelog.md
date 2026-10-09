# 变更日志

## v0.1.0（开发中）

### 2026-10-09

- **CHARMM36 真实测试案例（SMO-BST）**：新增 CHARMM36 真实测试数据（Smoothened–β-谷甾醇复合物，CHARMM36 + CGenFF + SPC/E 水，GROMACS 2018.1，来源 Mendeley v94vzbwzf3，Omar et al. 2020, *Data in Brief*）与 64 项真实数据测试（识别器 28 项 + 策略二/三各 16 项 + 水桥 4 项）；含两套 tpr（删水复合物 / 全原子 63055 个 SOL 水）；真实数据测试覆盖全部 8 类相互作用（SMO-BST 体系无卤素/金属中心，卤键与金属配位以零结果断言覆盖）
- **测试输出 gitignore**：`Tests/interaction_h5data/` 整目录加入 `.gitignore`（测试运行生成的 h5/csv/xpm 等产物，不进入版本库）

### 2026-10-08

- **跨力场水残基处理**：`GroupIdentifier.WATER_RESIDUES` 类属性成为统一入口（amber={SOL,HOH,WAT}、gromos={SOL}、charmm 含 TIP3、opls 含 HO4/HO5）；`_find_water` / `_find_metal_binding` 改用 `self.WATER_RESIDUES`（修复金属配位水漏排）；water_bridge 检测器按 water 组原子集合数据驱动排除水供体/受体（不依赖残基名，跨力场通用）
- **测试约定（策略一不跑）**：per_tuple 策略（无策略后缀的 `test_*.py`）标注"可能被舍弃"；全量测试改用显式文件列表，只跑策略二（`*_per_frame.py`）与策略三（`*_two_pass.py`）

### 2026-10-06 ~ 10-07

- **GROMOS 53A6 真实数据集成测试**：基于 `gromos53a6_md.tpr`（130 残基蛋白 + 6 个配体 ZIN1-6）的识别器 + 策略二/三集成测试；断言基线标注"尚未经人工核验"（除芳香环/带电残基与蛋白序列交叉对照外）
- **GROMOS 真实测试数据**：新增 `Tests/test_MD_case_gromos/`；amber 真实测试目录重命名 `test_MD_case` → `test_MD_case_amber`

### 2026-10-04

- **4 力场识别器可用**：GROMOS 53A6/54A7、CHARMM36/C36m、OPLS-AA/L 基团识别器合入 dev（各含合成 SystemData 单元测试）；`IDENTIFIER_CLASSES` 注册表 amber/gromos/charmm/opls 齐备，`dii run --ff` 可选 4 力场

### 2026-09-30

- **H 键受体误判修复（A1，4 力场）**：剔除带 H 的非受体 N 类型（普通酰胺/铵/带 H 吡咯/胍基），保留 Pro N、His 无 H 吡啶、中性胺、核酸氨基；按力场分别以类型 + H 邻居计数（Amber）、rtf 类型（CHARMM）、类型表（OPLS）、结构判据（GROMOS）实现；逐项判定依据见 `doc/acceptor_identification_evidence.md`
- **三待查项查证落地**：逐一核实 GAFF `n2`（亚胺型 sp2 N，按官方类型表判为受体）、CGenFF `NG*`（按官方 MASS 注释逐项判定留/剔）、核酸糖苷 `N*`（嘌呤 N9 / 嘧啶 N1，连糖 3 键吡咯型 → 剔除非受体）；对抗性审查修正（剔核酸糖苷 `N*`、锁定 CGenFF 配体测试），判定证据见 `doc/acceptor_identification_evidence.md`
- **调研与设计文档**：多力场兼容性调研（GROMOS/CHARMM/OPLS 对应评估）、GROMOS 识别器设计方案、竞品调研系列 + 论文框架 v1

### 2026-09-28

- **英文文档站点**：新增 `docs_en/`（独立 ReadTheDocs 项目）并由 CLI 同步英文化

### 2026-09-18 ~ 09-19

- **国际化与发布**：CLI 消息英文化 + 双语 README；补齐 PyPI 元数据（许可证 GPL-3.0）并发布 v0.0.1 初版包
- **测试补强**：GmxTprDumpReader / Pipeline / DII CLI 单元测试、h5 加载错误路径、读取器一致性测试
- **中文文档专业化修订**：清除营销话术、补全细节、修正技术准确性问题

### 2026-09-17

- **文档对抗性修订**：新建中文用户文档站点（guide + reference）；优化主 index 目录结构；修复 MyST 内部交叉引用；修正结果解读中疏水/卤键标签格式、h5 字符串 metric 存储机制描述

### 2026-09-16

- **XPM bug 修复**：绕过 DuIvyTools `refresh_by_value_matrix` 重映射，新增 `_build_discrete_xpm` 手动构建 Discrete XPM（value_matrix 即颜色索引，colors/notes 按索引对齐），修复值集合不完整时热力图颜色错位
- **dii export 增强**：遍历全部 Interaction（同类型重复加序号），空数据/0 帧跳过导出，损坏 h5 友好报错，pair_indices 拒绝 bool + 保序去重
- **h5 序列化加固**：metadata 支持 numpy 类型兜底序列化；path 参数类型校验（拒绝 BytesIO/None 垃圾文件）
- **依赖补全**：pyproject 新增 `h5py` / `scipy` / `DuIvyTools`
- **对抗性测试**：新增 XPM 手动构建（9 例）、DII export（4 例）、path 校验（4 例）测试

### 2026-09-14

- **命令行工具**：新增 `dii run`（相互作用检测 + h5 保存）与 `dii export`（导出 xvg/xpm/csv + 概览）
- **Pipeline 编排**：串联 Reader → Identifier → Detector → h5 保存

### 2026-09-05 ~ 09-13

- **结果序列化**：HDF5 无损存储/加载（Interaction/Group/Atom 全量 roundtrip）
- **导出器**：InteractionExporter 基类 + 8 个子类（xvg/xpm/CSV）；CSV 汇总、π堆积类型三值 XPM
- **结果时间序列**：Interaction 增加 times 字段

### 2026-08-12 ~ 09-03

- **相互作用检测**：8 种类型 × 3 策略（PerTuple / PerFrame / TwoPass）全部实现
- **TwoPass 性能优化**：水桥用 KDTree 预筛后耗时从 65h 降至 ~5s
- **目录重构**：`input_readers` → `system_readers`，新增 `io/` 目录

### 2026-08-11（基团识别完成）

- 基团识别模块完成（Amber 力场）
- 从 tpr 到官能团全链路打通（D927 体系验证）
- 类型映射表对 Amber 全家族（amber03/94/96/99/99sb/99sb-ildn/GS/14sb + GAFF）零冲突