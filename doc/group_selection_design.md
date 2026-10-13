# CLI 用户基团选择与过滤设计方案（group_selection）

> 创建日期：2026-10-09
> 分支：feature/user-specified-groups
> 状态：方案定稿（待实现）
> 参考：MDAnalysis 官方 Atom selection language（userguide.mdanalysis.org/selections.html，本项目已依赖 MDAnalysis）
> 现有能力：`detector.detect(..., tuple_filter=)` 已在全部检测器实现（`core/interfaces.py` L141-172 等），
> 本项目标即把这一 Python API 能力以 CLI 参数形式暴露给用户，检测器内部几何/阈值**零改动**。

---

## 1. 背景与需求

用户需要限定相互作用检测的作用范围，典型场景：

- 蛋白-蛋白界面分析（如只检测 RBD 与 KRAS 之间的氢键/盐桥/π 堆积）；
- 蛋白-配体相互作用分析；
- 分子内分析（某链/某分子内部的全部相互作用）；
- 三聚体/多链复合物的两两界面。

现状：`tuple_filter` 仅 Python API 可用；`Pipeline.run()` 与 CLI（`dii run`）未暴露，
`skills/duivyinteractions-skill/SKILL.md` 明确记录"CLI 不支持按分子/残基过滤"。
本设计把该能力打通到 CLI，遵循项目哲学：**一切对用户显式；任何歧义、空、不清楚的状态一律报错，把问题丢给用户；不做隐式猜测、自动转换、静默忽略。**

## 2. 设计哲学

1. **显式优于隐式**：软件只做用户明确表达的事；歧义状态绝不静默降级。
2. **复用成熟语法**：原子选择采用 MDAnalysis 选择语言（我们已依赖 MDAnalysis，零新语法、零新依赖）。
3. **KISS**：只新增 2 个可选参数；配对规则二值枚举，不做通用 DSL。
4. **帧无关架构**：基团鉴定只做一次、与帧无关；select 求值用第一帧坐标一次完成。
5. **不侵入检测器**：复用既有 `tuple_filter` 钩子，检测器几何/阈值零改动。

## 3. 参数模型

```
dii run -t md.tpr -f md.xtc -o out/ --ff amber \
       [--select "SELEXP" ["SELEXP" ...]] [--pairing any|between]
```

| 参数 | 形式 | 默认 | 语义 |
|---|---|---|---|
| `--select` | `append + nargs='+'`，可重复出现或单次多值，共 0..N 个选择串 | 无 | 每个串是 MDAnalysis 选择语法，定义一个**独立**的原子集合 |
| `--pairing` | 二值枚举 | `any` | `any` = 池内全部候选不设配对约束；`between` = 仅保留**跨 select** 候选对 |

### 3.1 MDAnalysis 选择语法速查（用户可用能力）

- **属性选择**：`segid NAME`、`resname NAME`、`resid 1:30`、`resnum 1 to 5`、`name CA`、`type X`、`moltype M`、`index 0:100`、`bynum 1:100`
- **标准集**：`protein`、`backbone`、`nucleic`、`nucleicbase`、`nucleicsugar`
- **布尔组合**：`and` / `or` / `not`、括号分组（如 `protein and not (resname ALA or resname LYS)`）
- **通配符**：`GL*`、`GL[YN]`、`?`
- **几何关键字**（第一帧静态求值，见 §4.4）：`around 6 (resname LIG)`、`byres ...`、`sphzone 6 (protein)`、`point X Y Z D`、`prop z > 5`

## 4. 核心语义（全部显式规则）

### 4.1 参与池与基团保留

- 各 `--select` 串在 Universe 上求值为**独立的原子集合** S₁..Sₙ（不合并）；
- **基团保留 ⇔ ∃ i：基团全部原子 ⊆ Sᵢ**（基团必须完整属于某一个 select）；
- 横跨多个 select 的基团（原子散落于不同 select）→ **不保留**（规则确定，非歧义）；
- 仅属并集但不完整属于任何单一 select → **不保留**。

### 4.2 水（隐式恒定"水侧"）

水是**隐式 select（水侧）**，按 P0-1 拍板：
- 用户无需选择、不可关闭；water 基团天然属于水侧，始终可用；
- **between 配对判定跳过水侧**（水桥判 (供体, 受体) 两端，水不参与）；
- 水的角色由力场 `WATER_RESIDUES` 系统定义（Amber `SOL/HOH/WAT`、CHARMM `TIP3`…），
  不属于用户自定义范畴，故隐式不违背显式哲学。

**归属优先级（唯一规则）**：**group_type == "water" 的基团**恒属水侧、跳过 select
归属判定。水残基（SOL 等）上识别出的 H_donor/H_acceptor 等基团**不属于水侧**——
它们须满足 select 覆盖才保留（防止水残基供体污染参与池）。
即使用户显式 `--select "resname SOL"` 选中水原子，也不改变 water 基团的归属
（仍在水侧、不参与配对判定）——保证显式选水不会破坏水桥的配对语义。

**显式选水的参与（方案 A，2026-10-09 补）**：`_filter_groups` 对非水桥检测器
默认排除 `WATER_RESIDUES` 残基的基团；但当用户**显式** `--select` 选中水残基时，
被归属到**具体 select（非水侧）**的水残基 H_donor/H_acceptor 基团**豁免该排除**，
可参与氢键等检测（如"蛋白 × 水之间的氢键"）。
- 无 select（`assignment=None`）时行为与既有完全一致（水仍被排除）；
- 归属为水侧（`WATER_SIDE`）的水基团仍不参与非水桥检测。

### 4.3 配对规则（tuple 层）

| pairing | 条件 | 行为 |
|---|---|---|
| `any` | select 数 ≥ 0 | 池内全部候选对（含同 select 内部对 + 跨 select 对） |
| `between` | **select 数 ≥ 2**（否则报错） | 仅保留两端基团**分属不同 select** 的对 |

**约束位置表**（配对判定作用于 tuple 哪些成员，固定不开放）：

| 相互作用类型 | tuple | 参与配对判定的位置 |
|---|---|---|
| hydrogen_bond / salt_bridge / hydrophobic / pi_stacking / pi_cation / halogen_bond / metal_coordination | (A, B) | (0, 1) |
| **water_bridge** | (donor, **water**, acceptor) | **(0, 2)** —— 水分子中立 |

### 4.4 帧无关与几何选择

- 基团鉴定与帧无关（架构铁律）；select 求值用**第一帧坐标一次完成**；
- 几何关键字（`around`/`sphzone`/`point`/`prop`）为**静态近似**——第一帧构象决定选中集合；
- 运行日志醒目打印："几何选择按第一帧坐标静态求值，构象大幅变化时请用静态属性选择或 Python API"。

## 5. 校验与报错（一切歧义/空 → 报错）

| # | 情况 | 处理 |
|---|---|---|
| 1 | select 串语法错误 | 报错（含串原文 + MDAnalysis 原因） |
| 2 | 某 select 匹配 0 原子 | 报错（"select 未选中任何原子"） |
| 3 | select 间原子重叠 | 报错（列出重叠原子：resid/残基名），用户自行修正（如 `and not`） |
| 4 | `between` 且 select 数 < 2 | 报错（"between 需要至少 2 个 select"） |
| 5 | 基团池中**非水基团为 0**（用户选择未命中任何真实化学基团） | 报错（提示：若使用 `around`/`sphzone` 请考虑 `byres`，基团要求全部原子被选中） |
| 6 | 结果因过滤为 0 对 | 正常（有基团无对是合法结果），统计为空 |

> 第 5 条说明：**group_type=="water" 的桥接水基团恒真保留**，故空池检查以
> "非水基团为 0" 为准——避免用户选错（如 `name CA` 只选到裸原子）时被水基团
> 撑满池而静默产出空结果。
> 常见触发：`--select "around 6 (resname LIG)"` 未加 `byres` ——
> 接触残基部分原子在圈外 → 基团全被排除 → 报错，用户自己加 `byres`。

## 6. 执行流程（帧无关）

```
1. mda.Universe(tpr, xtc) 已存在（pipeline 现有）
2. 显式 u.trajectory[0] 加载第一帧坐标（不依赖 MDAnalysis 默认帧状态）
3. 每个 --select 串 → u.select_atoms()（第一帧坐标求值）→ frozenset Sᵢ
4. 校验：语法 / 空 / 两两重叠（§5）
5. 基团归属：非水基团 → 找 i 使 atoms ⊆ Sᵢ（select 两两不相交 → 归属唯一）
              水基团（group_type == "water"）→ 恒属水侧（跳过 select 判定，见 §4.2）
6. 基团池 = 归属成功（含水侧）的全部基团；非水基团为 0 → 报错（见 §5 第 5 条）
7. 组装 tuple_filter 谓词 → 传给 detector.detect(..., tuple_filter=)
   - any：不传 tuple_filter（None，等价恒真、零开销）
   - between：两端基团 select 归属不同（水侧不参与）
     - pair 型检测器的 water 基团已被 _filter_groups 排除，谓词无需特判水
     - water_bridge 3 元组由约束位置表 (0,2) 限定参与判定的成员
8. 检测器内部几何/阈值零改动
```

## 7. 使用示例（验收基线）

| 意图 | 命令 |
|---|---|
| 全量（现状） | `dii run -t md.tpr -f md.xtc -o out/ --ff amber` |
| 分子内分析 | `dii run ... --select "segid A"` |
| 蛋白-配体界面 | `dii run ... --select "protein" --select "resname LIG" --pairing between` |
| 蛋白-蛋白界面 | `dii run ... --select "segid RBD_pro" --select "segid KRAS_pro" --pairing between` |
| 三聚体两两界面 | `dii run ... --select "segid A" --select "segid B" --select "segid C" --pairing between` |
| 配体 6Å 口袋残基 | `dii run ... --select "byres around 6 (resname LIG)"` |
| 界面水桥（水自动桥接） | 同蛋白-蛋白界面命令，water_bridge 自动判两端、水中立 |

## 8. 实现落点（不碰检测器几何逻辑）

```
新增 DuIvyInteractions/selection.py（单文件；遵循 KISS：函数 ≤20 行、类 ≤100 行）：
  - resolve_selection(universe, sel_str) -> frozenset[int]        # MDAnalysis 求值
  - validate_selections(sets) -> None                             # 空/重叠校验（报错）
  - assign_groups(groups, sets) -> Dict[int, int|"water"]         # 基团→select 归属
  - filter_groups(groups, assignment) -> List[Group]              # 保留完整属于某 select 的基团
  - build_tuple_filter(pairing, assignment, interaction_type)     # tuple 层谓词（查 CONSTRAINT_POSITIONS）
改动 pipeline.py：run() 增 select/pairing 参数；池过滤 + tuple_filter 透传
改动 DII.py：CLI 参数 + 全部报错路径
新增 Tests/unittests/test_selection_cli.py：合成 SystemData + 三力场真实各 1 例
文档：docs/guide/command.md、docs/reference/api.md（语法 + 语义 + 错误表 + 示例）
```

## 9. 已锁定决策清单（经拍板/对抗性审查确认）

1. ✅ 选择语法 = MDAnalysis（复用依赖，不发明 DSL）
2. ✅ `--select` 可多给，各 select 独立成集合（不合并）
3. ✅ **基团保留 ⇔ 基团全部原子 ⊆ 某一个 select**
4. ✅ `--pairing` 二值 any/between，无 within（分子内 = 单 select 天然表达）
5. ✅ select 间原子重叠 → 报错，用户手写修正，绝不自动转换
6. ✅ 空/歧义/不清楚 → 一律报错丢给用户
7. ✅ between 需 ≥2 select，否则报错
8. ✅ 水 = 隐式恒定水侧，默认在池、不参与配对（water_bridge 判两端）
9. ✅ 帧无关：几何选择第一帧静态求值 + 日志警示
10. ✅ 检测器内部零改动，复用既有 `tuple_filter`

## 10. 边界与已知限制

- **几何选择是静态近似**（第一帧构象），非逐帧动态；配体漂移后选中集合不更新（见 §4.4）；
- **原子级 select**（如 `name N` + `name H`）会使基团横跨 → 基团不保留（规则确定，非报错）；
  残基/分子级选择（`segid`/`resname`/`byres`）永不触发该分支；
- **重叠必须由用户手写修正**（如 `protein and not (byres around 6 ...)`），不做自动去重/归属；
- **显式选水 + water_bridge**：water_bridge 的水是桥接体非端（检测器内部排除水的
  供/受体），用户显式 `--select "resname SOL"` 选水时水不做端 → 供/受两端若同侧则
  无跨侧对 → 空结果（设计语义，非 bug）；如需水桥分析用默认隐式水侧（不显式选水）；
- 复杂自定义过滤（如"只保留不同分子间的 pair"之外的特殊规则）仍回落 Python API `tuple_filter`；
- 空池报错基于**非水基团**：用户只选水（如 `--select "resname SOL"` + 非水桥类型）
  时，水的供/受体基团归属具体 select 而豁免——此时"非水基团"含水的显式基团，
  空池检查按归属判定（非 WATER_SIDE 即计入），不误报。若用户希望强制排除水，
  用 `--select "resname SOL and not name H*"` 等显式写法或回落 Python API；
- 与 `WATER_RESIDUES`/`needs_water` 既有机制正交：`_filter_groups` 的 WATER_RESIDUES
  排除在无 select（assignment=None）时保持既有行为完全不变；仅显式归属到具体
  select 索引（int）的水残基基团豁免。

---

*文档结束*
