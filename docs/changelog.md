# 变更日志

## v0.1.0（开发中）

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