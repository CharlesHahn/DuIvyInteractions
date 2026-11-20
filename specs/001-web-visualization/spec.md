# Feature Specification: Web 可视化（本地 Web 分析界面）

**Feature Branch**: `feature/visualization`

**Created**: 2026-11-20

**Status**: Draft

**Input**: User description: "构建一个基于 web 的本地工具：读入 h5，展示各种结果；甚至可以直接读入 dii 的输入文件（tpr/xtc），通过 web 界面用户可以做计算和可视化。"

> 本 spec 采用 GitHub Spec Kit 的模板骨架（轻量用法，未 bootstrap `.specify/`）；
> 简要方向见 `doc/visualization_web_design.md`（W1/W2/W3 三阶段 + 已锁定的架构决定）。

## User Scenarios & Testing *(mandatory)*

### User Story 1 - 浏览与下钻 h5 检测结果 (Priority: P1)

用户在本机运行 `dii web`（仅端口/host 参数），浏览器自动打开本地页面，在界面内输入数据源路径（如 `dii run` 的输出目录 `out/`，含 8 类 h5）：先看到**总览**（8 个类型卡片：对数、帧数、occupancy 概况）；
点进某类型看到**全对表格**（pair_label、occupancy、avg/std 指标、flags 列，可排序/筛选/搜索/分页）；
点进某一对看到**下钻页**（existence 时间线、各指标逐帧曲线、基团元信息、原始数值表），
并可下载当前视图的 CSV。

**Why this priority**: 这是全部价值的地基——h5 里已有的一切（8 类检测结果、flags 标记、
元数据）都通过它变成"可看、可查、可导出"；不依赖 tpr/xtc 与任何计算，单独交付即完整可用。

**Independent Test**: 用仓库自带的真实数据（`dii run` 生成 D927 两类 h5）即可完整测试：
`dii web` 启动 → 页面内输入路径 → 总览渲染 → 类型页排序筛选 → 下钻页曲线 → CSV 下载，全程无需轨迹文件。

**Acceptance Scenarios**:

1. **Given** 输出目录含有效的 `hydrogen_bond.h5`（2366 对 × 101 帧），**When** 打开总览页，**Then** hydrogen_bond 卡片显示对数 2366、帧数 101，页面加载 < 3 秒
2. **Given** 类型页表格，**When** 按 occupancy 降序排序并搜索 "ARG"，**Then** 表格实时更新，仅显示匹配行
3. **Given** 某对的下钻页，**When** 查看，**Then** existence 时间线与全部数值指标曲线正确渲染（活跃帧与 `dii export` 的 CSV 数值一致）
4. **Given** h5 中存在 `flags/redundant_pi_stacking`，**When** 打开 hydrophobic 类型页，**Then** 表格含该 flags 列（冗余比例），无 flags 的 h5 则**不出现**该列
5. **Given** 下钻页，**When** 点击"下载 CSV"，**Then** 下载的文件可在表格软件中打开且数值与页面一致

---

### User Story 2 - 结构可视化：3D 分子上叠加相互作用 (Priority: P2)

用户在下钻页点击"3D 视图"，提供 tpr/xtc 路径后，浏览器内渲染任意一帧的 3D 分子结构，
并把该对（或该类型全部对）的相互作用画成**原子间的虚线/高亮**；
拖动帧滑块可切换帧，观察相互作用随时间的出现/消失。

**Why this priority**: 这是相对 PLIP 静态报告的核心差异化（全轨迹、8 类、力场驱动），
但依赖 3D 前端库选型与坐标导出管线，放在 W1 验证架构之后。

**Independent Test**: 用 D927 tpr/xtc + 任一 h5 测试：渲染 1 帧结构 → 叠加该对的原子对
虚线 → 滑块切帧 → 虚线随 existence 变化，不依赖 Story 1 之外的任何功能。

**Acceptance Scenarios**:

1. **Given** 有效的 tpr/xtc 与 h5，**When** 在下钻页打开 3D 视图，**Then** 结构渲染成功且该对的原子间虚线可见
2. **Given** 帧滑块，**When** 切换到 existence=False 的帧，**Then** 该对的虚线消失
3. **Given** 1001 帧轨迹，**When** 悔块切换，**Then** 无可感知卡顿（坐标按需加载，不整包下发）

---

### User Story 3 - Web 界面发起计算 (Priority: P3)

用户在 Web 表单里填入 tpr/xtc 路径、`--ff`、`--strategy`、`--interactions` 等
参数（等价 `dii run`），提交后任务后台执行并显示进度（复用 CLI 日志行），
完成后结果自动出现在总览页——形成"计算 → 浏览"闭环。

**Why this priority**: 编排性质（subprocess 复用 `dii run` 全部防护），价值在闭环但
不引入新分析能力；且全 8 类长轨迹计算是分钟~小时级，任务管理的复杂度必须最后处理。

**Independent Test**: 用短轨迹（CHARMM 100ps，11 帧）测试：表单提交 → 进度出现 →
完成 → 总览页出现新结果，计算结果与命令行 `dii run` 逐位一致。

**Acceptance Scenarios**:

1. **Given** 表单参数与命令行等价，**When** 提交，**Then** 后台任务启动且页面不阻塞
2. **Given** 计算进行中，**When** 刷新任务页，**Then** 可见最新日志行（[INFO]/[WARN]）
3. **Given** 任务完成，**When** 返回总览页，**Then** 新生成的 h5 出现在卡片中
4. **Given** 参数非法（路径不存在、--ff 拼错），**When** 提交，**Then** 友好报错（与 CLI 的 `error: ...` 一致），不产生半成品输出

---

### Edge Cases

- - 数据源选择：Web 内输入路径（单个 h5 文件或输出目录均可——目录扫描其中全部有效 h5，跳过损坏文件并告警，与 `dii export` 防护一致）；不做文件系统树浏览器（范围控制），仅路径输入 + 最近使用记忆
- **空数据 h5**（0 对或 0 帧）：类型卡片显示"无数据"，不报错（与 `dii export` 的 `[skip]` 语义一致）
- **损坏 h5**：加载失败只在总览页告警（`load_interactions` 抛 OSError/ValueError），不拖垮整个页面
- **旧版 h5 无 `flags/` 组**：正常显示，不出现 flags 列（向后兼容已验证）
- **大矩阵**：单类型 > 10⁴ 对或轨迹 > 10³ 帧时，服务端分页 + 降采样（existence 热力图按行聚合），绝不整包下发
- **端口被占用**：默认 8888，被占用自动 +1 顺延（最多尝试 10 个），页面提示实际地址
- **浏览器兼容**：仅保证最新两版 Chrome/Firefox/Edge（本地工具，不做兼容回退）

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: 系统 MUST 提供 CLI 入口 `dii web`，仅绑定 `127.0.0.1`（`--host` 可覆盖），参数仅 `--port`（默认 8888，被占用自动顺延并在页面提示实际地址），启动后自动打开默认浏览器；数据源（h5 文件/输出目录）不作为命令行参数，MUST 在 Web 界面内通过路径输入选择，并记住最近使用的数据源
- **FR-002**: 系统 MUST 通过 `load_interactions()` 读取 h5（existence、times、metrics、groups、可选 flags），不重新实现任何解析
- **FR-003**: 总览页 MUST 展示每个类型的对数、帧数、occupancy 分布概览
- **FR-004**: 类型页 MUST 提供表格视图，支持按 pair_label 搜索、按 occupancy/avg 指标排序、分页
- **FR-005**: 下钻页 MUST 展示单对的 existence 时间线与各指标逐帧曲线，数值与 h5 原始数据一致
- **FR-006**: 系统 MUST 在 h5 含 `flags` 时展示对应列（如 `avg_redundant_pi_stacking`），不含时不展示
- **FR-007**: 系统 MUST 提供 CSV 下载，字段与 `dii export` 的 summary 约定一致
- **FR-008**: 服务端 MUST 对大数据做聚合/降采样/分页，任何页面响应不整体下发 (n_pairs × n_frames) 全矩阵
- **FR-009**: 所有筛选/统计/聚合/绘图数据组装逻辑 MUST 位于 `visualizers/` 纯函数层并可单元测试；路由与模板仅做胶水
- **FR-010**: Web 相关依赖 MUST 全部声明在可选 extra `[viz]` 中，核心包安装不引入任何 web/绘图依赖
- **FR-011**: 绘图 MUST 采用"服务端组装图表 JSON + 前端库纯渲染"的方式；前端图表库定为 **Plotly.js**（用户拍板 2026-11-20；JS 内嵌进包保证离线可用，不引入 plotly.py）
- **FR-012**: 3D 视图 MUST 支持逐帧切换并把 h5 的 `atom_global_idx` 映射到坐标虚线；3D 前端库定为 **Mol\***（选型依据：我们的全原子体系达 10⁵ 原子级，Mol* 的 GPU 渲染与空间分区原生支撑该规模；原生支持 MD 轨迹格式与逐帧切换；PDB 官方 viewer，NAR 2021，长期维护有保障。3Dmol.js 轻量、API 简单，但在 10⁵+ 原子体系有内存吃紧的社区报告，不采用）
- **FR-013**: 计算发起 MUST 通过 subprocess 调用 `dii run`（复用 CLI 全部参数校验与日志语义），且 MUST 支持任务进行中查看日志、完成或失败后有明确终态
- **FR-014**: 任何页面对损坏/空数据的 h5 MUST 只告警不崩溃，行为与 `dii export` 的防护语义一致

### Key Entities *(include if feature involves data)*

- **Interaction**（既有，只读）：interaction_type、groups（含 Group 的 atom_global_idx/residue 元数据）、existence (n_pairs, n_frames)、metrics、times、flags
- **h5 输出目录**（既有）：`out/` 下每类型一个 `<type>.h5`，格式版本 1.0（可选 flags/ 组）
- **PairView**（新，纯函数层）：单对的下钻视图数据——元信息 + 逐帧 existence/metrics 序列
- **TypeTable**（新，纯函数层）：类型页的表格数据——行（对）× 列（occupancy、avg/std 指标、flags 派生列），支持排序/筛选/分页参数
- **VizJob**（新，W3）：一次计算任务的参数集、状态（pending/running/done/failed）、日志缓冲、产出目录

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: D927 真实数据（4196 对 × 101 帧的 two_pass 水桥）总览与类型页首次加载 < 3 秒（本机 127.0.0.1）
- **SC-002**: 单对下钻页 1001 帧曲线渲染 < 1 秒
- **SC-003**: 零数据孤岛——`dii run` 产出的任何 h5 字段（含 flags）都能在 UI 中找到对应展示，覆盖率 100%
- **SC-004**: 核心包零污染——`pip install duivyinteractions` 后不出现任何 web/绘图依赖；`pip install duivyinteractions[viz]` 后 W1 全功能可用
- **SC-005**: W3 发起的计算与等价命令行 `dii run` 的 h5 输出逐位一致
- **SC-006**: 纯函数层单测覆盖筛选/聚合/图表 JSON 组装的核心路径（≥ 20 个用例），传输层薄胶水不承担逻辑

## Assumptions

- 本机单用户工具：无鉴权、无并发压力设计，仅绑定 127.0.0.1；Web 内路径输入只做本机文件读取，不提供文件系统树浏览（与本地定位一致，无安全增量）
- h5 全部由 `dii run` 生成（格式版本 1.0），不处理手工编辑的 h5
- 主力环境 Linux/MacOS；Windows 尽力支持（不阻塞交付）
- tpr/xtc 与 h5 位于本机文件系统，**不做文件上传**（W3 用路径输入）
- 现有 `dii export` 保持不变，与 viz 并存不互替
- Python >= 3.9 的支持边界与 [viz] extra 的依赖版本兼容性在实现期核实（如冲突则在 spec 修订中明确放宽到 >= 3.10）
