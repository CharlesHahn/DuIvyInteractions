# XPM 值重映射 Bug 修复方案

> 创建日期：2026-09-16
> 更新日期：2026-09-16
> 状态：**已实施并验证（2026-09-16）**
> 作者：hanyl

---

## 1. Bug 根因

### 1.1 症状

`to_xpm_existence` / `to_xpm_stacking_type` 导出的 XPM 在值集合不完整时颜色与数据错位：
- existence 全活跃 → XPM 全显示白色"No"（活跃帧被画成非活跃）
- stacking 全活跃 `[T,P]` → P 型显示成 T 型，蓝色丢失

### 1.2 根因链

1. **DuIvyTools Discrete 契约**：`parse_xpm` 明示 "for Discrete, value store the index of chars|notes|colors"——Discrete 型 XPM 的 `value_matrix` 存的是**索引**，不是任意语义值。
2. **`refresh_by_value_matrix(is_Continuous=False)`** 会把 value_matrix 从"任意值"重映射为"去重排序后的索引"（`value_to_index`），且**强制覆盖** notes/colors（源码 TODO 自证"不支持预定义 notes/chars/colors"）。
3. **exporter 误用**：在 refresh 之后又按原始语义覆盖 colors/notes，与 refresh 已定死的 `color_num`/`chars` 数量冲突 → `save()` 的 `zip(chars, colors, notes)` 截断 → 错位。
4. **数据依赖**：值集合完整 `{0,1,2}` 时重映射为恒等映射，碰巧正确；缺任一类别即错。

### 1.3 责任

- DuIvyTools `refresh_by_value_matrix`：行为自洽（规整器），无责
- matplotlib：忠实按索引取色（`ListedColormap`），无责
- **根因：本项目 exporter 对 DuIvyTools API 的误用**

---

## 2. 修复思路（第一性原理）

- 我们的 `vm`（0=无/1=T/2=P；或 0/1）**本身就是 0..k-1 连续索引**，且与 colors/notes 顺序天然对齐——正好是 DuIvyTools Discrete 语义，**不需要 refresh 的"任意值→索引"规整**。
- 因此正确做法：**绕过 refresh**，直接构建 XPM 的公开字段（官方构建模式：`XPM("", is_file=False, new_file=True)` + 设置字段 + `save()`，见 DuIvyTools `otherCommands.py` DCCM 示例）。
- 对已是指数的数据调用 refresh 属于误用，值集合不完整时破坏语义。

---

## 3. 修改内容

### 3.1 `DuIvyInteractions/io/interaction_exporter.py`

| 改动 | 说明 |
|:--|:--|
| 加 `import string` | XPM 字符表需要 |
| 加常量 `XPM_LETTERS` | 与 DuIvyTools refresh 的 letters 方案一致（82 字符） |
| 加模块级 `_validate_xpm_indices(value_matrix, n_colors)` | 校验整数类型 + 索引范围 [0, n_colors) |
| 加基类方法 `_build_discrete_xpm(...)` | 手动构建 Discrete XPM：value_matrix 即索引，colors/notes 按索引对齐 |
| 重写 `to_xpm_existence` | 改用 `_build_discrete_xpm`，删除 refresh 调用与事后覆盖 |

### 3.2 `DuIvyInteractions/io/pi_stacking_exporter.py`

| 改动 | 说明 |
|:--|:--|
| 重写 `to_xpm_stacking_type` | vm 构造不变（0/1/2），改用 `_build_discrete_xpm` |

### 3.3 `Tests/unittests/test_exporters.py`

新增 `TestDiscreteXpmBuild` 对抗性测试类（9 个用例）：
- 全活跃 existence 不被重映射（原 bug 回归）
- stacking `[T,P]` 全活跃不错位（原 bug 回归）
- 部分活跃输出与旧版一致（回归）
- 全活跃 roundtrip 读回一致
- 负索引/float dtype/越界正索引/>82 色/colors-notes 长度不一致 全部拒绝

---

## 4. 验证结果（2026-09-16 实测）

### 4.1 单元测试

`pytest Tests/unittests/test_exporters.py` → **57 passed, 6 skipped**（6 skip 为缺真实 h5 的 TestRealData，与本次改动无关）。

### 4.2 真实数据端到端（Tests/test_MD_case/md.tpr + md1ns.xtc）

| 验证项 | 结果 |
|:--|:--|
| 盐桥 47 对×101 帧 `dii export` | ✅ existence.xpm 读回 value_matrix 与源 existence 完全一致 |
| π堆积 10 对×101 帧 `dii export` | ✅ type.xpm 读回 T=1/P=2/无=0 三值与源 metrics 完全一致（T:213, P:67, 无:730） |
| 概览/Top 占位率 | ✅ 正常打印 |

### 4.3 边界防御（预验证脚本实测）

| 输入 | 行为 |
|:--|:--|
| `value_matrix` 越界正索引 | `ValueError: index out of range [0, 2): min=0, max=3` |
| 负索引 | `ValueError: index out of range [0, 2): min=-1, max=1`（原 refresh 静默错乱） |
| float dtype | `TypeError: must be integer, got float64` |
| colors > 82 | `ValueError: too many colors (90), max 80`（显式拒绝，不实现 char_per_pixel=2） |
| colors/notes 长度不一致 | `ValueError: colors(2) != notes(1)` |

---

## 5. 影响面

- 仅改 exporter 的 XPM 构建；h5 序列化、其他 exporter、pipeline、DII.py 均未触碰
- 值集合完整的数据输出字节级不变（恒等映射）；不完整的数据从错变对
- 下游 DIT `xpm_show`/`xpm2csv`/`xpm2dat` 无接口变化

---

*文档结束*