# dii 命令参考

`dii` 是 DuIvyInteractions 的命令行入口（`pyproject.toml` 的 `[project.scripts]` 注册，指向 `DuIvyInteractions.DII:main`），提供两个子命令：`run`（检测）与 `export`（导出）。

## 全局用法

```bash
dii [-h] {run,export} ...
```

`dii --help` 输出即上示用法；`dii run --help` 与 `dii export --help` 可查看各子命令的完整参数。

## 子命令：run —— 运行相互作用检测

检测指定轨迹中的相互作用，并将结果保存为 HDF5 文件（格式版本 1.0）：

```bash
dii run -t TPR -f XTC -o OUTPUT --ff FF [--interactions LIST] [--strategy STRATEGY]
           [--select SELEXP ...] [--pairing any|between]
```

### 参数

| 参数 | 必填 | 默认 | 说明 |
|:-----|:----|:----|:-----|
| `-t, --tpr` | ✅ | — | GROMACS 拓扑文件（基团鉴定的输入） |
| `-f, --xtc` | ✅ | — | 轨迹文件（MDAnalysis 可读格式，如 xtc） |
| `-o, --output` | ✅ | — | 输出目录（自动创建） |
| `--ff` | ✅ | — | 力场类型，取值 `amber` / `gromos` / `charmm` / `opls`（argparse choices 由 `IDENTIFIER_CLASSES` 注册表生成） |
| `--interactions` | 否 | `all` | 要检测的相互作用类型，逗号分隔（大小写不敏感），如 `hydrogen_bond,pi_stacking` |
| `--strategy` | 否 | `two_pass` | 检测策略：`two_pass` / `per_frame` / `per_tuple` |
| `--select` | 否 | 无 | MDAnalysis 选择串，`action="append", nargs="+"`：可**重复出现**或**单次多值**（如 `--select "segid A" "segid B"`）；每个值定义一个**独立**的原子集合（不合并）。基团须**完整属于**某一个 select 才进入参与池（见下文"用户基团选择与过滤"） |
| `--pairing` | 否 | `any` | 配对规则：`any` = 参与池内全部候选不设配对约束（默认）；`between` = 仅保留两端基团分属不同 select 的候选对（需 ≥2 个 select，否则报错） |

`--ff` 四个取值对应的识别器（`DuIvyInteractions/group_identifiers/__init__.py`）：

| 取值 | 识别器 | 覆盖力场 |
|:-----|:-----|:-----|
| `amber` | `AmberFFGroupIdentifier` | Amber 03/94/96/99/99sb/99sb-ildn/GS/14sb 蛋白 + GAFF/GAFF2 配体 |
| `gromos` | `GromosFFGroupIdentifier` | GROMOS 53A6 / 54A7（联合原子力场，配 SPC/SPC-E 水） |
| `charmm` | `CharmmFFGroupIdentifier` | CHARMM36 / C36m（含 CGenFF 配体） |
| `opls` | `OplsFFGroupIdentifier` | OPLS-AA / L |

传入未知 `--ff` 会报错并列出全部可用值；传入未知 `--interactions` 类型会报错退出并列出全部可用类型（`ALL_INTERACTIONS`）；`--select` / `--pairing` 配置错误（语法错误、空 select、select 间重叠、`between` 数量不足、空池）一律以 `error: ...` 形式友好报错退出（见下文"用户基团选择与过滤"与"退出码"）。

### 支持的 8 种相互作用类型

`hydrogen_bond`（氢键）、`pi_stacking`（π-π 堆积）、`salt_bridge`（盐桥）、`hydrophobic`（疏水）、`halogen_bond`（卤键）、`metal_coordination`（金属配位）、`water_bridge`（水桥）、`pi_cation`（π-阳离子）。

- `--interactions` 值先统一转小写、去除首尾空白，再按逗号切分（`DII.py` 的 `main`）。
- `all` 是默认值，等价于"全部 8 类"；**`all` 不得与其他类型混用**，混用时报错 `'all' cannot be mixed with other types`。
- 单类型检测失败（如体系不含金属、无卤素）只打印 `[WARN] <类型> detection failed: <原因>`，不中断其余类型。

### 用户基团选择与过滤（--select / --pairing）

`--select` 与 `--pairing` 用于限定检测的作用范围：`--select` 定义参与检测的基团**参与池**，`--pairing` 定义候选对的配对规则。选择语法复用 **MDAnalysis Atom selection language**（本项目已依赖 MDAnalysis，不引入新语法）；选择与配对仅作用于基团池与候选对层面，**检测器内部的几何判据与阈值零改动**。核心用法：

```bash
dii run -t md.tpr -f md.xtc -o out/ --ff amber \
       [--select "SELEXP" ["SELEXP" ...]] [--pairing any|between]
```

- `--select` 为 `action="append", nargs="+"`：可**重复出现**，也可**单次多值**，最终拍平为一列选择串。每个值是一条 MDAnalysis 选择串，定义一个**独立**的原子集合，多个 select **不合并**。以下两种写法等价：

```bash
dii run -t md.tpr -f md.xtc -o out/ --ff amber --select "segid A" --select "segid B"
dii run -t md.tpr -f md.xtc -o out/ --ff amber --select "segid A" "segid B"
```

- `--pairing` 取 `any`（默认）或 `between`：`any` = 参与池内全部候选不设配对约束；`between` = 仅保留**跨 select** 的候选对（需 ≥2 个 select，否则报错）。

#### 语法速查

`--select` 的值是一条 MDAnalysis 选择语法字符串，常用写法：

| 类别 | 语法 | 示例 |
|:-----|:-----|:-----|
| 属性选择 | `segid` / `resname` / `resid` / `resnum` / `name` / `type` / `moltype` / `index` / `bynum` | `segid seg_0_RBD_pro`、`resname LIG`、`resid 1:30`、`resnum 1 to 5`、`name CA`、`type X`、`moltype M`、`index 0:100`、`bynum 1:100` |
| 标准集 | `protein` / `backbone` / `nucleic` / `nucleicbase` / `nucleicsugar` | `--select "protein"` |
| 布尔组合 | `and` / `or` / `not`、括号分组 | `protein and not (resname ALA or resname LYS)` |
| 通配符 | `*` / `[ ]` / `?` | `GL*`、`GL[YN]`、`?` |
| 几何关键字（第一帧静态求值） | `around` / `byres` / `sphzone` / `point` / `prop` | `around 6 (resname LIG)`、`byres around 6 (resname LIG)`、`sphzone 6 (protein)`、`point X Y Z D`、`prop z > 5` |

> 完整语法与各关键字定义以 MDAnalysis 官方 Atom selection language 文档为准。

#### 语义

**基团保留规则**（参与池的构成）：

- 基团保留 ⇔ 该基团**全部原子**都被**某一个** select 集合完整覆盖（全部原子 ∈ 某个 select）；
- 不接受"仅属多个 select 的并集"、不接受横跨（原子散落于不同 select）的基团；
- 因此**原子级** select（如 `name N` + `name H`）会使基团横跨 → 不保留，属确定规则而非报错；**残基/分子级**选择（`segid` / `resname` / `byres`）不会触发该分支。

**水侧（隐式恒定）**：

- `group_type == "water"` 的桥接水基团**恒属"水侧"**——隐式 select，用户不可关闭、无需选择；水基团始终保留，跳过 select 归属判定；
- 水残基（SOL 等）上识别出的 H_donor/H_acceptor **不属于水侧**，须满足 select 覆盖才保留（防止水残基供体污染参与池）；
- between 配对判定**跳过水侧**。

**配对规则（between）**：

- 仅保留两端基团分属**不同 select** 的候选对；任一端归属为水侧（或未归属）的对被排除；
- 约束位置表固定、不开放给用户：氢键 / π-堆积 / 盐桥 / 疏水 / 卤键 / 金属配位 / π-阳离子 → `(0, 1)`（两端）；水桥 `(donor, water, acceptor)` → `(0, 2)`，水分子中立（不参与配对判定）；
- 分子内分析 = 单个 `--select` + 默认 `any`（不设配对约束），无需额外参数。

**显式选水（方案 A）**：

- 未指定 `--select`（`assignment=None`）时，水残基基团按 `WATER_RESIDUES` 排除（与既有行为完全一致）；
- 显式选中水残基（如 `--select "resname SOL"`）时，被归属到**具体 select（非水侧）**的水残基 H_donor/H_acceptor 基团**豁免** `WATER_RESIDUES` 排除，可参与氢键等检测（如"蛋白 × 水之间的氢键"）。

#### 示例

```bash
# 蛋白-蛋白界面（Amber 真实数据：RBD × KRAS，between）
dii run -t Tests/test_MD_case_amber/md.tpr -f Tests/test_MD_case_amber/md1ns.xtc \
        -o out_rbd/ --ff amber \
        --select "segid seg_0_RBD_pro" --select "segid seg_2_KRAS_pro" --pairing between

# 蛋白-配体界面（CHARMM36 真实数据：SMO 蛋白 × BST 配体，between）
dii run -t Tests/test_MD_case_charmm36/SMO-BST_md_complex.tpr \
        -f Tests/test_MD_case_charmm36/SMO-BST_md_complex.xtc \
        -o out_smo/ --ff charmm \
        --select "protein" --select "resname BST" --pairing between

# 分子内分析（单 select，pairing 默认 any）
dii run -t Tests/test_MD_case_amber/md.tpr -f Tests/test_MD_case_amber/md1ns.xtc \
        -o out_rbd_only/ --ff amber --select "segid seg_0_RBD_pro"
```

通用写法（占位符替换）：蛋白-配体界面 `--select "protein" --select "resname LIG" --pairing between`；三聚体两两界面 `--select "segid A" --select "segid B" --select "segid C" --pairing between`（N 个 select 的两两之间均保留）；配体口袋残基 `--select "byres around 6 (resname LIG)"`；显式选水（蛋白 × 水氢键）`--select "segid seg_0_RBD_pro" --select "resname SOL" --pairing between`。

> **几何选择注意**：`around 6 (resname LIG)` 未加 `byres` 时，接触残基的部分原子可能在圈外，基团因"未完整覆盖"被排除，可能触发空池报错（见下文错误表）；加 `byres` 可保证整残基进入选中集合。

#### 错误信息

选择配置错误一律转换为**友好报错**（`SystemExit`，输出 `error: ...`，不暴露 traceback）：

| 情况 | 输出示例 |
|:-----|:-----|
| select 串语法错误 | `error: select 语法错误: '<串>' → <MDAnalysis原因>` |
| 某 select 匹配 0 原子 | `error: select #N 未选中任何原子` |
| select 间原子重叠 | `error: select #N 与 select #M 原子重叠 (K 个原子)，请用 and not 显式剔除后重试` |
| `between` 且 select 数 < 2 | `error: --pairing between 需要至少 2 个 --select` |
| 参与池内非水基团为 0（选择未命中任何真实化学基团） | `error: 参与池内无可用基团（非水基团 0 个）。若使用了 around/sphzone 几何选择，基团可能因部分原子在圈外被排除，请改用 byres around …` |

空池检查以**非水基团为 0** 为准：水侧桥接水基团恒真保留，因此用户选错（如 `--select "name CA"` 只选到裸原子）时不会因水基团撑满参与池而静默产出空结果。select 间原子重叠**必须由用户手写修正**（如 `protein and not ...`），工具不做自动去重/归属。

#### 限制

- **几何选择是静态近似**：`around` / `sphzone` / `point` / `prop` 按**第一帧坐标**静态求值，非逐帧动态；配体大幅漂移后选中集合不更新。select 求值是唯一的坐标依赖点，第一帧完成后即固定，与"基团鉴定一次、与帧无关"的架构不冲突；
- **显式选水 + 水桥**：water_bridge 中水是桥接体（非端），显式选水时水不做端；若供体/受体两端同属一个 select，则无跨侧候选对 → 空结果（设计语义，非缺陷）；如需水桥分析，使用默认隐式水侧（不显式选水）；
- 复杂自定义配对规则（超出"跨 select"语义）请回落 Python API 的 `tuple_filter`（见{doc}`/reference/api`）。

### 内部流程

1. 用 `GmxTprReader` 读取 tpr，识别器做一次基团鉴定（与帧无关，只做一次）；
2. 用 `mda.Universe(tpr, xtc)` 加载轨迹；水残基名取识别器的 `WATER_RESIDUES`（跨力场正确排除水）；
3. （若指定了 `--select`）先 `u.trajectory[0]` 定位第一帧，把每条选择串求值为独立原子集合（几何选择按第一帧坐标静态求值）；校验语法 / 空 / 两两重叠后，将基团归属到对应的 select（`group_type == "water"` 的桥接水基团恒属水侧），过滤出参与池（见上文"用户基团选择与过滤"）；
4. 对每个类型：按检测器的 `required_group_types` 过滤基团（除水桥外排除水分子；显式归属到具体 select 的水残基基团豁免），`between` 时经 tuple 层谓词只保留跨 select 候选对，逐帧检测，结果保存为 `<output>/<类型>.h5`。

### 输出

每个类型生成一个 h5 文件：`<output>/<interaction_type>.h5`。例如：

```bash
dii run -t md.tpr -f md.xtc -o out/ --ff amber
# 生成 out/hydrogen_bond.h5, out/pi_stacking.h5, out/salt_bridge.h5,
#      out/hydrophobic.h5, out/halogen_bond.h5, out/metal_coordination.h5,
#      out/water_bridge.h5, out/pi_cation.h5
```

h5 格式版本为 1.0（`io/h5.py` 的 `FORMAT_VERSION`），gzip 压缩，详见{doc}`/reference/data_format`。

### 策略说明

| 策略 | 机制 | 适用场景 |
|:-----|:-----|:---------|
| `two_pass`（默认） | Pass1 逐帧精确发现活跃基团组（KDTree 预筛 + 精确判据，稀疏存储）→ Pass2 对发现的基团组补全全部帧的指标 | 大体系、长轨迹，性能最优（水桥 KDTree 预筛后 65 h → ~5 s） |
| `per_frame` | 对全部候选基团组逐帧向量化计算，预分配稠密矩阵 | 候选基数极大的类型（如水桥）；长轨迹时预分配矩阵内存占用高 |
| `per_tuple` | 逐候选基团组遍历全部帧并向量化 | 对照实现，处于"可能被舍弃"状态（测试不跑），候选多时极慢 |

## 子命令：export —— 导出结果

读取 h5 文件，导出为 xvg/xpm/csv，并在终端打印概览：

```bash
dii export -i INPUT -o OUTPUT
```

### 参数

| 参数 | 必填 | 说明 |
|:-----|:----|:-----|
| `-i, --input` | ✅ | h5 文件路径 |
| `-o, --output` | ✅ | 输出目录（自动创建） |

### 行为

- 支持包含多个 Interaction 的 h5（`dii run` 生成的文件每类型一个，符合此场景；`io/h5.py` 的 `load_interactions` 返回 Interaction 列表）。
- 若 h5 中同类型含多个 Interaction，文件名追加序号避免覆盖（如 `salt_bridge_2_count.xvg` 等）。
- 空数据（0 对或 0 帧）打印 `[skip] <类型>: nothing to export (pairs=0, frames=0)` 并跳过该 Interaction。
- h5 缺失或损坏报错：`cannot read h5 file: <路径> (file may be corrupted or missing)`（`SystemExit`）。
- h5 中无任何结果（空列表）报错：`no interaction results in h5: <路径>`（`SystemExit`）。
- h5 含未知类型报错：`unknown interaction type: '<类型>'. Available: ...`（`SystemExit`）。
- 导出的 xvg 仅包含数值指标（字符串指标如 π-堆积的 `pistacking_type` 不导出为 xvg，仅在 type.xpm 中体现）。

### 概览输出

每个被导出的 Interaction 先打印概览（`DII.py` 的 `_print_overview`）：

```
===== Salt Bridge overview =====
Type:     salt_bridge
Pairs:    47
Frames:   101
Time range: 0.0 ~ 1000.0 ps

Top 5 occupancies:
  1. ARG210(3443-3451)···ASP211(3461-3463)  100.0%
  ...
```

`Pairs` 为该类型基团对数（`Interaction.n_pairs`），`Frames` 为帧数，`Time range` 为首尾帧时间（ps，保留 1 位小数），`Top 5` 按占位率降序（`OVERVIEW_TOP_N = 5`）。

### 输出文件

对每个 Interaction，在输出目录生成：

| 文件 | 内容 |
|:-----|:-----|
| `<类型>_count.xvg` | 每帧活跃相互作用数量（X 轴：时间 ps） |
| `<类型>_<指标>.xvg` | 每个数值指标的时间序列，每列一条曲线（对应一个基团对） |
| `<类型>_existence.xpm` | 存在性热力图（行=pair，列=帧；白=No，蓝=Yes） |
| `<类型>_type.xpm` | 仅 π-堆积：堆积类型图（0=无 / 1=T 型 / 2=平行） |
| `<类型>_summary.csv` | 每对的汇总统计（`pair_label, occupancy, avg_<指标>, std_<指标>`） |

各类型的数值指标（`<类型>_<指标>.xvg` 的文件名由此决定，数据来自各导出器的 `metric_labels`）：

| 类型 | 指标（对应 xvg 文件名） |
|:-----|:-----|
| `hydrogen_bond` | `distance`（D-A 距离）、`angle`（D-H···A 角） |
| `salt_bridge` | `distance`（电荷中心距离） |
| `pi_stacking` | `distance`（环心距离）、`angle`（法向量夹角）、`offset`（偏移） |
| `hydrophobic` | `distance` |
| `halogen_bond` | `distance`（X···A 距离）、`don_angle`（C-X···A 角）、`acc_angle`（X···A-R 角） |
| `metal_coordination` | `distance`（金属-配位原子距离） |
| `water_bridge` | `dist_dw`（D-Ow）、`dist_wa`（Ow-A）、`theta`（O-H···Ow 角）、`omega`（H-Ow···A 角） |
| `pi_cation` | `distance`（环-电荷距离）、`offset`（偏移） |

文件名不以 Interaction 在 h5 中的序号命名，而以类型名命名——`dii run` 输出的原始 h5 文件名（`<类型>.h5`）与导出文件名（`<类型>_xxx`）一一对应。

## 退出码

- `0`：成功
- 非零：参数错误、文件不存在/损坏、未知类型、空结果等，均通过 `SystemExit` 抛出（如 `unknown interaction type: 'xxx'. Available: ...`、`'all' cannot be mixed with other types`、`cannot read h5 file: ...`）。
- `--select` / `--pairing` 相关错误同样非零退出：输出带 `error:` 前缀的友好提示（如 `error: select 语法错误: ...`、`error: select #1 未选中任何原子`、`error: --pairing between 需要至少 2 个 --select`、`error: 参与池内无可用基团（非水基团 0 个）...`），不暴露 traceback。

## 下一步

- 每种输出文件的具体格式与解读，见[结果解读](result)。
- 真实数据上的完整运行示例，见[快速上手](quickstart)。