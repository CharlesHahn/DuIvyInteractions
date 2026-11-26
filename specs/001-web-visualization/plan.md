# Implementation Plan: Web 可视化（本地 Web 分析界面）

**Branch**: `feature/visualization` | **Date**: 2026-11-20 | **Spec**: [spec.md](spec.md)

**Input**: Feature specification from `/specs/001-web-visualization/spec.md`

## Summary

实现 W1（h5 查看器，P1）：`dii web` 启动 FastAPI 本地服务（仅 127.0.0.1），
Web 内输入路径选择数据源 → 总览卡片 → 类型表格（搜索/排序/分页）→ 单对下钻
（existence 时间线 + 指标曲线）+ CSV 下载。逻辑全部下沉 `visualizers/` 纯函数层
（可单测），FastAPI 路由与 Jinja 模板仅做胶水。W2（Mol\* 3D）与 W3（Web 发起计算）
按 spec 的 P2/P3 顺延，本 plan 只规划 W1 的落地结构。

## Technical Context

**Language/Version**: Python >= 3.9（[viz] extra 依赖版本兼容性在实现期核实）

**Primary Dependencies**: FastAPI（路由）+ Uvicorn（服务）+ Jinja2（模板）；
图表为**内嵌 Plotly.js**（plotly.js-basic-dist-min，MIT，打包进静态资源，离线可用）；
**不引入 plotly.py**——服务端直接组装图表 JSON（dict）。

**Storage**: 只读 h5（`io.h5.load_interactions`）；"最近使用数据源"存
`~/.cache/duivyinteractions/web_recent.json`（无数据库）。

**Testing**: pytest（仓库约定：显式文件列表）；纯函数层用 tmp 目录构造 h5 fixture，
不依赖 MDAnalysis 轨迹；路由用 FastAPI TestClient（httpx）。

**Target Platform**: 本机 127.0.0.1（--host 可覆盖）；Linux/MacOS 优先，Windows 尽力。

**Project Type**: 既有 CLI 库的子命令（`dii web`）+ 本地 web-service（同包内 `DuIvyInteractions/web/`）。

**Performance Goals**: SC-001 总览/类型页首次加载 < 3s（D927 4196 对 × 101 帧）；
SC-002 单对下钻 1001 帧曲线 < 1s。

**Constraints**: 只绑 127.0.0.1；核心包零 web 依赖（全部在 `[viz]` extra）；
任何响应不整体下发 (n_pairs × n_frames) 全矩阵（分页 + 按需下钻）；
不替代 `dii export`（并存）。

**Scale/Scope**: 单用户；8 类 h5 × 每类可达 10⁴ 对 / 10³ 帧。

## Constitution Check

*GATE: 通过（对照 AGENTS.md 编程哲学）*

| 门 | 结论 |
|---|---|
| KISS（函数 ≤ 20 行、类 ≤ 100 行） | ✅ 纯函数层按"一个视图数据一个函数"切分；模板/路由薄胶水 |
| 显式优于隐式（无魔法数字、类型注解） | ✅ 分页大小、降采样上限、端口顺延次数全部常量；全函数带注解 |
| 可读性即文档（docstring 写"是什么"） | ✅ "为什么"留在 specs/ 与 doc/ 设计文档 |
| 单一职责、依赖高层→低层 | ✅ `visualizers/`（逻辑）与 `web/`（传输）分层，均只依赖 core/io；core 零改动 |
| 不做"可能有用"的功能 | ✅ 明确非目标：树浏览器、鉴权、上传、SPA、浏览器重计算 |

## Project Structure

### Documentation (this feature)

```text
specs/001-web-visualization/
├── spec.md              # 功能规格（已完成，决策全部定案）
├── plan.md              # 本文件
└── tasks.md             # 任务拆解（W1，按 US1 分组）
```

### Source Code (repository root)

```text
DuIvyInteractions/
├── visualizers/                 # 纯函数层（可单测，零 web 依赖）
│   ├── __init__.py
│   ├── data.py                  # 数据源解析/总览/类型表格/下钻/降采样/最近数据源
│   └── figures.py               # Plotly JSON 组装（dict，不依赖 plotly.py）
├── web/                         # 传输层（薄胶水）
│   ├── __init__.py
│   ├── app.py                   # create_app()：路由 + Jinja + 静态挂载
│   ├── static/
│   │   ├── plotly-basic.min.js  # 内嵌 plotly.js-basic-dist-min（MIT）
│   │   └── app.css              # 极简样式
│   └── templates/
│       ├── base.html
│       ├── index.html           # 数据源选择 + 总览卡片
│       ├── type.html            # 类型表格（搜索/排序/分页）
│       └── pair.html            # 下钻页（时间线 + 曲线 + 元信息）
├── DII.py                       # + `web` 子命令（惰性导入 web 包，端口顺延，自动开浏览器）
└── pyproject.toml               # + [viz] extra：fastapi / uvicorn / jinja2 / httpx(TestClient)
Tests/unittests/
└── test_web_viz.py              # 纯函数层 ≥ 20 用例 + 路由冒烟（SC-006）
```

## 关键技术决策

| # | 决策 | 依据 |
|---|---|---|
| D1 | **图表 JSON 服务端组装**：`visualizers/figures.py` 返回 `{"data": [...], "layout": {...}}` 纯 dict | 不装 plotly.py（拖 pandas/kaleido）；组装逻辑可单测（spec FR-011） |
| D2 | **Plotly.js 内嵌打包资源**：`web/static/plotly-basic.min.js`（basic bundle ~1MB，覆盖 line/scatter/bar/heatmap） | 本地工具必须离线可用；basic bundle 已含 W1 全部图型 |
| D3 | **表格用服务端渲染 + GET 参数**（search/sort_by/sort_dir/page），不引入 htmx/HTMX/Alpine | W1 的表格交互 GET 参数即可表达；少一个前端依赖（KISS） |
| D4 | **数据源状态存服务端内存**（单用户本地） + 最近数据源持久化 json | spec FR-001/Edge Case；`/select` POST 后 redirect，无会话体系 |
| D5 | **pair 标签复用现有导出器**：把 `DII.py` 的 type→Exporter 映射提为模块级常量 `_EXPORTER_CLASSES`，`visualizers/` 直接实例化调 `get_pair_label` | 单一来源，避免复制 8 个标签格式（water_bridge 三元组等） |
| D6 | **降采样常量**：`MAX_SERIES_POINTS = 2000`（下钻曲线超过则 index 抽稀，保留首末）；表格分页 `PAGE_SIZE = 50` | spec Edge Case"大矩阵"；SC-002 |
| D7 | **端口顺延**：`--port` 起 socket 试连，占用则 +1 最多 10 次，实际地址在启动日志与页面显示 | spec FR-001/Edge Case |
| D8 | **`dii web` 惰性导入**：web/fastapi/uvicorn 全部在子命令内 import | 核心包不装 [viz] extra 时 `dii run/export` 完全不受影响（spec FR-010/SC-004） |

## Evaluation (SC 对照)

| SC | 验证方式 |
|---|---|
| SC-001（< 3s） | 本机 D927 two_pass 水桥数据冒烟计时 |
| SC-002（< 1s） | 下钻页冒烟计时（GROMOS 1001 帧数据） |
| SC-003（覆盖率 100%） | 测试断言：existence/metrics/flags 在 UI 数据组装中全部出现 |
| SC-004（核心零污染） | `pip install duivyinteractions`（无 extra）时 import 核心包不触碰 fastapi |
| SC-005 | W3 范围，本 plan 不涉及 |
| SC-006（纯函数 ≥ 20 用例） | `Tests/unittests/test_web_viz.py` |
