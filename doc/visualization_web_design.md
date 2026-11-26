# Web 可视化设计（简要版）

> 创建日期：2026-11-17
> 状态：**简要设计（方向已定，细节后续迭代）**
> 分支：feature/visualization
> 前置文档：`doc/visualization_design.md`（2026-09-05 初版，其中 P0 导出项已由 `dii export` 实现，本文档聚焦 Web 形态）

## 1. 目的

给 DII 一个**本地 Web 界面**，把检测结果变成可交互探索的东西：

1. 读入 `dii run` 产出的 h5 文件，在浏览器里浏览、筛选、下钻、画图——替代翻 xvg/csv
2. 读入 tpr/xtc，直接在界面上发起计算（跑 `dii run`），完成后自动展示结果——形成"计算 → 分析"闭环
3. （远期）3D 结构上叠加相互作用标注

## 2. 简要目标（三阶段，每阶段独立可用）

| 阶段 | 目标 | 说明 |
|---|---|---|
| **W1：h5 查看器** | `dii viz -i out/` 起本地 Web | 总览页（类型卡片）→ 类型页（表格：排序/筛选/搜索，含 flags 列）→ 下钻页（单对 existence 时间线 + 指标曲线 + 元信息）+ CSV 下载 |
| **W2：结构可视化** | 3D 分子 + 相互作用叠加 | 从 tpr/xtc 抽帧导出坐标，浏览器 3D 渲染（Mol\*/3Dmol.js），h5 的 `atom_global_idx` 生成虚线/高亮；帧滑块看逐帧状态 |
| **W3：Web 发起计算** | 表单参数 → 跑 `dii run` → 展示 | 本地路径输入（不搞大文件上传）；subprocess 复用 CLI 全部防护；异步任务 + 进度；完成后跳转 W1 |

## 3. 关键设计决定（本次锁定）

| # | 决定 | 理由 |
|---|---|---|
| 1 | **本地 Web，只绑 127.0.0.1** | 定位是本机分析工具，非多人服务；顺带规避安全面 |
| 2 | **输入以 h5 为主**（`load_interactions()`），W3 另接受 tpr/xtc 路径 | h5 是唯一持久化载体，flags/metrics/groups 全在里面 |
| 3 | **技术栈：FastAPI + Jinja + 轻量前端（htmx/Alpine）**，3D 用浏览器库（CDN） | 纯 Python、无 node 构建链、KISS；3D 局部引入，不搞 SPA |
| 4 | **与 `dii export` 并存不替代** | export = 静态文件交付（论文/归档）；viz = 交互探索 |
| 5 | **核心包不背 web 依赖**：`pip install duivyinteractions[viz]` 可选 extra | 检测核心保持轻量 |
| 6 | **大数据不整包下发**：服务端聚合 + 按需下钻 | μs 级轨迹的 (n_pairs × n_frames) 矩阵不能塞浏览器 |

## 4. 架构草图

```
浏览器（127.0.0.1:PORT）
   │  HTTP
   ▼
FastAPI（薄传输层：路由 + 模板）
   │
   ▼
visualizers/（纯函数层：筛选/统计/聚合/画图数据组装——可单测）
   │
   ├── load_interactions()      # io 层现成
   ├── Pipeline / subprocess    # W3：跑 dii run，复用 CLI 防护与日志
   └── MDAnalysis 抽帧          # W2：坐标导出
```

原则：**逻辑全部下沉 `visualizers/` 纯函数，传输层薄到只剩胶水**（Web UI 难测，纯函数层必须可单测）。

## 5. W1 页面清单（第一个可交付）

| 页面 | 内容 |
|---|---|
| 总览 `/` | 每类型卡片：对数、帧数、occupancy 分布、flags 概况 |
| 类型页 `/type/<name>` | 全对表格：pair_label、occupancy、avg/std 指标、flags；排序/筛选/搜索/分页 |
| 下钻页 `/pair/<type>/<idx>` | existence 时间线、各指标曲线、基团元信息、原始数值表 |
| 下载 | 每页 CSV（复用现有导出逻辑） |

## 6. 非目标（明确不做）

- 不做多用户/鉴权/部署到公网
- 不在浏览器里做重计算（聚合全在服务端）
- 不改检测核心与 h5 格式
- 不替代 `dii export`

## 7. 迭代方式

每阶段完成即合并；细节（页面布局、筛选维度、3D 库选型）在实现中逐步优化，不预先过度设计。
