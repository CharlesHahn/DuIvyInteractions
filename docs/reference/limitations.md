# 已知限制与版本信息

本文档说明工具的版本信息、支持范围与已知限制，帮助用户判断适用性。

## 版本信息

| 项 | 值 |
|:---|:---|
| 软件版本 | v0.0.1 |
| HDF5 结果格式版本 | 1.0 |
| 支持 Python | >= 3.9 |
| 包名（PyPI） | `duivyinteractions` |
| 命令行入口 | `dii` |

## 支持范围

- **力场（4 家族）**：
  - **Amber 家族**：amber03/94/96/99/99SB/99SB-ildn/GS/14SB 蛋白 + GAFF/GAFF2 配体（全原子显式 H）
  - **GROMOS 53A6/54A7**：联合原子力场，极性 H 显式；配 SPC/SPC-E 水模型（`SOL`）；类型粒度粗，由结构判据补偿
  - **CHARMM36/C36m**：蛋白 + CGenFF 配体；默认 TIP3/HOH 水模型
  - **OPLS-AA/L**：含 `opls_XXX` 编号映射；HOH/SPC、HO4/TIP4P、HO5/TIP5P 水模型
- **拓扑**：GROMACS tpr（经 MDAnalysis 读取）
- **轨迹**：MDAnalysis 可读取的格式（xtc 等）
- **结构**：依赖显式 H 做供体判定（D–H 键且 q(H)>0）；联合原子力场（GROMOS）无脂肪族显式 H，供体/疏水邻接判定相应受限（疏水已用类型白名单 + 极性邻居排除补偿）

## 边界

- **GROMOS 配体不支持**：GROMOS 配体需 ATB（Automatic Topology Builder）自动拓扑参数化，本项目不覆盖此流程；识别器验证范围为 GROMOS 蛋白残基。
- **per_tuple 策略（策略一）处于"可能被舍弃"状态**：该策略为逐候选组遍历的对照实现，测试不跑，不作为现行结果；现行策略为 `per_frame` 与 `two_pass`。

## 已知限制

### 检测判据相关

- **PBC 未处理**：周期边界条件下的相互作用（尤其水桥、氢键）未做周期性处理，长轨迹或跨边界相互作用可能不准确
- **金属配位几何构型未匹配**：仅用距离判据，未做 linear/trigonal/tetrahedral/octahedral 等配位构型匹配；纯水配位已排除（`metal_binding` 基团排除 `WATER_RESIDUES`）
- **水桥水分子去重未做**：一个水分子参与多个氢键时，未按 H-O-H 角度择优保留两个
- **水桥 TwoPass 策略缺距离下界**：TwoPass 策略未加 Ow-A 距离下界（2.5 Å），可能与 PerTuple/PerFrame 策略结果有少量差异
- **疏水-芳香未去重**：π-堆积与疏水可能对同一物理接触重复计数
- **跨策略差异**：疏水相互作用的"同残基去重"（同一 `(group_id, residue_id)` 对仅保留平均距离最近的 pair）仅在 per_frame 策略实现，two_pass 策略缺失（检测器代码内留有 TODO）；水桥 TwoPass 缺距离下界亦属此类
- **氮基团分类粗糙**：正电氮基团统一归为"叔胺"类（复现 PLIP 定义，实际涵盖所有 sp3 N），未细分伯/仲/叔胺；π-阳离子检测的叔胺角度检查可能误触发
- **基团识别结果未经人工全面审查**：单元测试基于代码输出，尚需人工核对"标准答案"（`doc/TODO.md` 待办）

### 性能相关

- **长轨迹内存占用**：PerFrame 策略为全部候选对预分配 `(n_pairs, n_frames)` 矩阵，1μs 级长轨迹（数十万帧）内存占用高（如 42,603 候选 × 500,000 帧的 distance 矩阵可达数百 GB 量级）；TwoPass 策略无此问题（稀疏存储 + 仅对活跃对补全）

### 功能范围

- **可视化**：结果导出为 xvg/xpm（可用 DuIvyTools/Xmgrace 绘图），工具内置绘图功能未实现
- **H 键受体判定（A1 修复，2026-09-30）**：已修复"铵/胍基/带 H 吡咯等非受体 N 被误判为受体"的跨力场问题——四力场按"力场类型 + 化学事实 + 文献"逐类型剔除普通酰胺/铵/带 H 吡咯/胍基 N，保留 Pro N / His 无 H 吡啶 / 中性胺 / 核酸氨基；实证效果（GROMOS 体系）：acceptor 336→206、H 键 pair 108→86、水桥 2541→1698。逐项判定依据见 `doc/acceptor_identification_evidence.md`

## 注意事项

- 阈值与 PLIP 一致，便于结果对照，但具体体系可能需按需调整（常量在检测器文件顶部，修改后需重跑受影响类型单测）
- 基团鉴定只做一次、与帧无关；几何判定逐帧进行，计算量随帧数线性增长
- 水残基排除按力场区分（`WATER_RESIDUES` 类属性：Amber SOL/HOH/WAT、GROMOS SOL、CHARMM TIP3/HOH/SOL/WAT、OPLS HOH/HO4/HO5/SOL/WAT），使用 `--ff` 时务必选择与体系一致的力场