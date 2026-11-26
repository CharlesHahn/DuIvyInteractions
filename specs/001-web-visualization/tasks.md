# Tasks: Web 可视化 W1（h5 查看器）

**Input**: Design documents from `/specs/001-web-visualization/`

**Prerequisites**: spec.md, plan.md

**Organization**: 按 User Story 分组；W2/W3（US2/US3）不在本清单，US1 完成并合并后另开。

## Format: `[ID] [P?] [Story] Description`

- **[P]**: 可并行（不同文件，无依赖）
- **[Story]**: 所属用户故事（本清单全部为 US1）

## Phase 1: Setup

- [x] **T001** [US1] `specs/001-web-visualization/plan.md`：技术上下文与结构（本分支已完成）
- [x] **T002** [US1] `pyproject.toml` 新增 `[viz]` extra：fastapi / uvicorn / jinja2 / httpx（TestClient 用）

## Phase 2: Pure Function Layer（`DuIvyInteractions/visualizers/`，零 web 依赖）

- [x] **T003** [P] [US1] `visualizers/data.py::parse_datasource(path)`：路径 → `[(type_name, h5_path)]`；目录扫全部 h5、单文件取一；损坏/空数据跳过并告警（spec Edge Case）
- [x] **T004** [P] [US1] `visualizers/data.py::overview(cards)`：每类型卡片数据（n_pairs、n_frames、flags 名单、occupancy 数组）
- [x] **T005** [P] [US1] `visualizers/data.py::type_table(interaction, ...)`：行 = 全部 pair（pair_label、occupancy、avg/std × 每指标、flags 派生列），支持 search / sort_by / sort_dir / page（常量 `PAGE_SIZE=50`）
- [x] **T006** [P] [US1] `visualizers/data.py::pair_drilldown(interaction, idx)`：元信息 + times + existence(0/1) + 指标序列（> `MAX_SERIES_POINTS=2000` 时降采样，保留首末）
- [x] **T007** [P] [US1] `visualizers/data.py::recent_datasource`：读/写 `~/.cache/duivyinteractions/web_recent.json`（最多 10 条，损坏文件静默重建）
- [x] **T008** [P] [US1] `visualizers/figures.py`：`fig_existence_timeline` / `fig_metric_timeseries` / `fig_occupancy_hist`——返回纯 dict `{"data": [...], "layout": {...}}`（常量：标题/轴标签，plotly.js-basic 的 step/scatter/histogram trace）
- [x] **T009** [US1] `DII.py`：type→Exporter 映射提为模块级 `_EXPORTER_CLASSES`（`_run_export` 复用，`visualizers/` 经 io 层导入调用 `get_pair_label`）

## Phase 3: Transport Layer（`DuIvyInteractions/web/`，薄胶水）

- [x] **T010** [P] [US1] 内嵌 `web/static/plotly-basic.min.js`（plotly.js-basic-dist-min，MIT，~1MB）+ `web/static/app.css`（极简）
- [x] **T011** [P] [US1] `web/templates/`：base / index（路径表单 + 最近数据源 + 卡片）/ type（表格 + 搜索排序分页表单）/ pair（元信息 + 图容器 + CSV 链接）
- [x] **T012** [US1] `web/app.py::create_app()`：`GET /`（数据源表单 + 卡片）；`POST /select`（校验路径 → 服务端状态 → redirect `/`）；`GET /type/{name}`（表格）；`GET /pair/{name}/{idx}`（下钻）；`GET /api/figures/{name}/{idx}.json`（figures JSON）；`GET /download/csv/{name}`（CSV 响应）；损坏 h5 只告警（FR-014）
- [x] **T013** [US1] `DII.py` + `web` 子命令：仅 `--port`/`--host`；socket 试连顺延端口（D7）；惰性导入 `web` 包；启动后线程延迟打开浏览器

## Phase 4: Verification

- [x] **T014** [US1] `Tests/unittests/test_web_viz.py`：纯函数层 ≥ 20 用例（SC-006）——数据源解析（文件/目录/损坏/空）、overview、type_table（分页/搜索/排序 asc+desc/flags 列）、drilldown（长度/值域/元信息）、downsample（上限/首末保留/小输入不变）、figures（结构/长度/trace 数）、CSV 字段、recent roundtrip + 路由冒烟（TestClient：/ 200、/select redirect、/type 200、figures JSON 结构、CSV 200）
- [x] **T015** [US1] 仓库测试约定回归：显式文件列表全量跑（含 test_web_viz）
- [x] **T016** [US1] 冒烟：`dii web` 启动 + 真实 D927 数据页面可达 + SC-001/SC-002 计时

## Notes

- 全部任务已实现并通过 T015 全量回归（728+ 基线 + test_web_viz）。
- W2（Mol\* 结构视图）与 W3（Web 发起计算）待 US1 合并后另立 spec/tasks。
