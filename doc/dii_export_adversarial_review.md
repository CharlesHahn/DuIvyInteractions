# dii export 对抗性审查报告

> 创建日期：2026-09-16
> 更新日期：2026-09-16
> 状态：**审查完成；2.1/2.2 已修复（2026-09-16）**
> 作者：hanyl

---

## 1. 审查背景

### 1.1 审查对象

对 `dii export` 命令相关代码做对抗性审查（adversarial review），以"找茬"心态检查正确性、边界情况、错误处理与性能问题。

### 1.2 审查范围

| 文件 | 职责 |
|:-----|:-----|
| `DuIvyInteractions/DII.py` | CLI 入口（`dii run` / `dii export`） |
| `DuIvyInteractions/io/interaction_exporter.py` | 导出基类（xvg/xpm/csv） |
| `DuIvyInteractions/io/pi_stacking_exporter.py` | π-堆积导出器（type.xpm） |
| `DuIvyInteractions/io/h5.py` | HDF5 序列化/反序列化 |
| `DuIvyInteractions/core/datas.py` | Interaction/Group 数据类 |
| `DuIvyInteractions/io/*_exporter.py`（8 子类） | 各类型 pair 标签 |
| `DuIvyInteractions/pipeline.py` | 检测流水线 |
| `Tests/unittests/test_exporters.py` | 现有测试 |

### 1.3 审查方法

所有疑点均**实测验证**：
- 构造对抗场景（全活跃矩阵、空帧、n_pairs=0、多 Interaction h5、metadata 含 numpy 类型等）
- 真实数据端到端：`Tests/test_MD_case/md.tpr + md1ns.xtc` → 盐桥检测（47 对 × 101 帧）→ save h5 → `dii export` 全流程

---

## 2. 🚨 严重 Bug（静默错误，输出与事实相反/数据丢失）

### 2.1 XPM 值被 `refresh_by_value_matrix` 静默重映射 → 热力图颜色错位

**位置**：`interaction_exporter.py:211`（`to_xpm_existence`）、`pi_stacking_exporter.py:98`（`to_xpm_stacking_type`）

**根因**：DuIvyTools 的 `refresh_by_value_matrix(is_Continuous=False)` 会把 value_matrix 从"原始值"改写为"去重排序后集合中的索引"：

```python
out_value_list = sorted(list(set(chain(*self.value_matrix))))
value_to_index = {value: idx for idx, value in enumerate(out_value_list)}
...
if not is_Continuous:
    self.value_matrix = indices.tolist()   # 值被重映射！
```

而 exporter 是在 refresh **之后**才设置 `colors`/`notes`（固定按 0/1/2 语义写死），两处语义错位。**当值集合不完整时**（矩阵缺某个值类别），重映射不再是恒等映射，数据被静默改写。

**实测复现**：

| 场景 | 输入 | 期望 value_matrix | 实际 value_matrix |
|:-----|:-----|:-----|:-----|
| existence 全活跃 | `[[True,True,True]]` | `[[1,1,1]]` | `[[0,0,0]]`（全白"No"！） |
| stacking 型全活跃 | `[["T","P"]]` | `[[1,2]]` | `[[0,1]]`（P 型显示成 T 型） |

- existence 全活跃：值集合 `{1}` → 重映射 `{1:0}` → 全部变 0
- 全活跃 stacking：值集合 `{1,2}` → 重映射 `{1:0, 2:1}` → 颜色与 notes 错位

**触发条件**：只要某个值类别在整个矩阵中缺位（如全部帧都活跃、或数据中没有某种堆积类型）。真实数据下完全可能发生。

**偶然性说明**：仓库已提交的真实 `Tests/interaction_h5data/pi_stacking_type.xpm` 恰好值集合 {0,1,2} 完整、映射恒等，**只是运气好没暴露**。这是一个数据依赖的隐蔽 bug。

**✅ 已修复（2026-09-16）**：绕过 `refresh_by_value_matrix`，新增 `_build_discrete_xpm` 手动构建 Discrete XPM（value_matrix 即颜色索引，colors/notes 按索引对齐），并新增 `_validate_xpm_indices` 校验（拒绝越界/负/float 索引）。修复原理与代码详见 `doc/dii_xpm_fix.md`。

### 2.2 `dii export` 只导出第一个 Interaction，其余静默丢弃

**位置**：`DII.py:63-66`

```python
interactions = load_interactions(args.input)
...
it = interactions[0]    # ← 只取第一个！
```

`save_interactions` 明确支持 Interaction **列表**，但 export 只处理 `interactions[0]`，其余无任何警告地丢弃。

**实测复现**：构造含 salt_bridge + pi_stacking 两种类型的 h5，`dii export` 只导出 salt_bridge 的 xvg/xpm/csv，pi_stacking 数据全部丢失。

---

## 3. ⚠️ 中等 Bug（崩溃/异常）

### 3.1 n_pairs=0 的 h5 → 导出中途崩溃，留下部分输出

- `to_xvg_count` 先成功写文件 → `to_xpm_existence` 抛 `ValueError: No pairs found`
- `dii export` 无 try/except 兜底 → traceback + 不完整输出目录

### 3.2 空帧（0 帧）h5 → overview IndexError

- `_print_overview` 中 `it.times[0]` / `it.times[-1]` 对 0 帧数组直接 `IndexError`（实测确认）

### 3.3 metadata 含 numpy 数组 → save_interactions 崩溃

- `h5.py:194`：`json.dumps(g.metadata)` 对 ndarray 抛 `TypeError: Object of type ndarray is not JSON serializable`（实测确认）
- 当前真实数据 metadata 只有 str/None 未触发，但 API 脆弱——检测器/未来代码一旦写入 numpy 值即崩溃
- Group.metadata 的 docstring 声称"值只能是 JSON 支持的类型"，但**无运行时校验**

### 3.4 `save_interactions` 对 file-like 输入产生垃圾文件

- `h5.py:28`：`path = str(path)` 把 `io.BytesIO` 转成字符串当文件名，直接在工作目录创建 `<_io.BytesIO object at ...>` 垃圾文件（测试中实踩）
- 应校验 path 类型（str/Path）或支持 file object

---

## 4. 🔸 低风险/健壮性问题

| # | 问题 | 位置 | 说明 |
|:--|:-----|:-----|:-----|
| 1 | `bool` 是 `int` 子类，`[True]` 静默当索引 1 通过 `_validate_pair_indices` | `interaction_exporter.py:366` | `isinstance(True, int)` 为 True |
| 2 | 重复 pair_indices `[0,0]` 未去重 | `interaction_exporter.py:371` | xvg 产生重复列/legend |
| 3 | `to_csv_summary` 对全 NaN 活跃帧抛 `RuntimeWarning: Mean of empty slice` | `interaction_exporter.py:298-299` | `np.nanmean`/`np.nanstd` 对空 slice 警告（实测 2 条） |
| 4 | CSV 不过滤 occupancy=0 的 pair | `interaction_exporter.py:292` | 氢键 4 万对时输出 4 万行几乎全空行，文件膨胀 |
| 5 | **pyproject.toml 缺 h5py 依赖** | `pyproject.toml:12-15` | `dependencies` 只有 numpy+MDAnalysis；README 声称 SciPy>=1.7 也未列入。干净环境 `dii` 直接跑不了 |
| 6 | 损坏/缺失 h5 报原始 `OSError: file signature not found` | `h5.py:52` | 无友好 CLI 错误 |
| 7 | 子类 `get_pair_label` 强假设 groups 元组结构与原子数 | 各 exporter 子类 | 如氢键必须 donor 有 2 原子，不符则 IndexError |
| 8 | 测试盲区 | `test_exporters.py` | `TestRealData` 因无 h5 全部 skip；无多-interaction/空帧/n_pairs=0/全活跃矩阵用例——**正是 2.1/2.2/3.1/3.2 全漏过的原因** |

---

## 5. ✅ 通过项（排除误报）

以下疑点经实测**不是 bug**，记录以免后续重复排查：

| 疑点 | 实测结论 |
|:-----|:---------|
| 字符串 metric（pistacking_type）h5 roundtrip | dtype=`<U1` 无损，shape 正确 |
| 检测器非活跃帧字符串填充 | 检测器用 `'N'` 填充，`to_xpm_stacking_type` 的 T/P mask 逻辑本身正确（2.1 纯粹是 refresh 重映射造成） |
| 负索引 `[-1]` 校验 | 正确拒绝 |
| `pair_indices=[]` 空列表 | 合法，只产生时间列 |
| 空 interaction 列表 roundtrip | 正常；`dii export` 正确报"h5 中无相互作用结果" |
| PRO 出现在盐桥结果 | C 端 COO⁻，非 bug（用户确认） |
| 真实数据端到端 | 101 帧、47 对盐桥：detect → save h5 → `dii export` 全流程跑通，xvg/xpm/csv 全部生成 |

---

## 6. 🎯 修复优先级建议

| 优先级 | 问题 | 建议方案 |
|:-------|:-----|:---------|
| **P0** | 2.1 XPM 值重映射 | ① refresh 后重建 value_matrix；② colors/notes 依据实际值集合动态生成（不写死 0/1/2 语义）；③ 绕开 DuIvyTools `is_Continuous=False` 的重映射（如直接构造 dot_matrix，不依赖其改写） |
| **P0** | 2.2 多 Interaction 丢弃 | `_run_export` 遍历全部 interactions 逐个导出；或检测到多个时警告/要求用户指定 type |
| **P1** | 3.1/3.2 export 前校验 | export 前校验 `n_pairs>0`、`n_frames>0`，给出友好错误 |
| **P1** | 低风险 #5 依赖 | pyproject 补 `h5py`（README 的 SciPy 也补进 dependencies） |
| **P2** | 3.3/3.4 输入防御 | metadata 序列化前校验/转换 numpy 类型；path 参数类型校验 |
| **P2** | 低风险 #8 测试盲区 | 补对抗性测试：全活跃矩阵、空帧、n_pairs=0、多 Interaction、metadata 含 numpy |

---

## 7. 复现脚本要点

```python
# 2.1 复现：全活跃 existence → XPM 全白
it = Interaction("salt_bridge", groups,
                 existence=np.array([[True, True, True]]),
                 metrics={"distance": np.array([[3.0, 3.1, 3.2]])},
                 times=np.array([0.0, 10.0, 20.0]))
xpm = SaltBridgeExporter().to_xpm_existence(it)
assert np.array(xpm.value_matrix).tolist() == [[1, 1, 1]]  # 实际得到 [[0, 0, 0]]

# 2.2 复现：多 Interaction h5 只导出第一个
save_interactions([salt_bridge_it, pi_stacking_it], "multi.h5")
# dii export -i multi.h5 -o out/  → 只生成 salt_bridge_*.xvg/xpm/csv
```

---

*文档结束（审查完成，待修复）*