# OPLS 力场基团识别器设计方案（opls_ff_identifier）

> 创建日期：2026-10-08
> 分支：feature/opls-identifier
> 状态：方案定稿（待实现）
> 参考实现：`amber_ff_identifier.py`（OPLS 与 Amber 同为"类型名编码化学"的全原子力场，走"填特征表"路线）+ 复用 CHARMM 的 `_is_terminal_ammonium` 结构判据
> 实证依据：本地 `oplsaa.ff`（**OPLS-AA/L v2001**，forcefield.doc 明确）+ [Towhee OPLS-aa 类型语义表](https://towhee.sourceforge.net/forcefields/oplsaa.html)（Jorgensen 组参数权威实现）+ [Robertson 2015（OPLS-AA/L 综述）](https://pmc.ncbi.nlm.nih.gov/articles/PMC4504185/) + [GROMACS 手册（OPLS-AA/L 定义）](https://manual.gromacs.org/documentation/2022/reference-manual/topologies/force-field-organization.html) 交叉验证

---

## 1. 背景与核心差异

DII 的基团鉴定依赖四类证据源：类型名语义、键合图、元素+电荷、残基/原子名。

OPLS 是**全原子力场**（GROMACS 内置为 **OPLS-AA/L**，即 2001 年氨基酸二面角更新版，见 [GROMACS 手册](https://manual.gromacs.org/documentation/2022/reference-manual/topologies/force-field-organization.html)、[Robertson 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4504185/)），类型语义存在但与 Amber/CHARMM 有两个关键差异：

| 维度 | Amber/CHARMM | **OPLS** |
|---|---|---|
| 类型名 | 符号（CA/CT/N3…）| **数字编号 `opls_XXX`**（805 类型）|
| 类型→符号映射 | 类型名即符号 | **`ffnonbonded.itp` 第二列含原始符号**（CA/CT/NA/NB…）|
| 版本 | — | **OPLS-AA/L（2001）** + GROMACS 移植扩展（如 CX=HIP 专用）|
| 带电 LYS | LYS=NH3+（+1）| **LYS=中性(0.0)、LYSH=质子化(+1)**（同 GROMOS）|
| His 命名 | HSD/HSE/HSP（CHARMM）| **HISD/HISE（中性）+ HISH（=HISP/HIS+，质子化）** |

**关键结论**：
1. OPLS 是"填表"路线（与 CHARMM 同构），但**类型表 = opls 数字编号**，且可用符号映射驱动构建
2. **LYS 中性问题与 GROMOS 相同** → 带电字典用 **ARG/LYSH**（正电）；`_is_terminal_ammonium` 结构判据**直接复用**（LYSH NZ=opls_287 -0.30 结构同 CHARMM NH3）
3. 类型→化学语义以 **符号** 为准（Towhee/Robertson 权威），opls 编号是 GROMACS 外壳
4. ⚠️ **GROMACS oplsaa.ff 含 Towhee-OPLS-aa(1996) 没有的扩展类型**（如 `CX`=HIP 咪唑碳，ffbonded 注释 "jtr: HIP"）——**权威核对必须本地文件 + 论文双源，不能只依赖 Towhee**

---

## 2. 关键事实（本地 oplsaa.ff 实测 + Towhee 权威验证）

### 2.1 类型映射（opls_XXX ↔ 符号 ↔ 语义，双重来源）

| 符号 | opls_XXX | 语义（Towhee）| 本地证据 |
|---|---|---|---|
| `CA` | opls_145 | 中性芳香碳 | ffnonbonded 第2列 + rtp PHE 环碳全 opls_145 |
| `CT` | opls_135/136/137 | sp3 脂肪碳 | ffnonbonded + rtp 骨架 |
| `HC/HA` | opls_140/146 | 脂肪/芳香氢 | ffnonbonded |
| `C` | opls_235 | 羰基碳 | ffnonbonded + rtp 主链 C |
| `N` | opls_238 | 酰胺氮 | ffnonbonded + rtp 主链 N |
| `C_3/O2` | opls_271/272 | 羧酸碳/氧 | ffnonbonded + rtp ASP/GLU |
| `N3` | opls_287 | sp3 铵氮 | ffnonbonded + rtp LYSH NZ / N 端 n.tdb |
| `N2` | opls_300/303 | sp2 胍/芳香胺氮 | ffnonbonded + rtp ARG |
| `CA` | opls_302 | 胍碳 | ffnonbonded + rtp ARG CZ |
| `C*` | opls_500 | 5 元环芳香碳 | rtp TRP CG |
| `CB/CN/CW/CV/CR` | opls_501/502/506/507/508/514 | TRP/His 环碳 | rtp TRP/HIS |
| **`CX`** | **opls_510** | **HIP（质子化 His）咪唑碳 CG/CD2**（GROMACS 移植扩展，ffbonded 注释 "jtr: HIP"，仿 CV/CW）| ffbonded.itp + rtp HISH CG/CD2 |
| `CR` | opls_509 | HIP 咪唑碳 CE1（仿 C* 语义）| ffbonded + rtp HISH CE1 |
| `NA` | opls_503/512 | sp2 芳香 N 带 H | rtp HISD ND1 / TRP NE1 / HISH ND1+NE2 |
| `NB` | opls_511 | 5 元环脱质子 N | rtp HISD NE2 |
| `H` | opls_504/513 | N 上 H | rtp |
| `NT` | opls_900 | sp3 胺氮（中性）| rtp LYS NZ（中性）|

### 2.2 蛋白残基实测（rtp 净电荷核算）

| 残基 | 净电荷 | 关键原子 | 带电结论 |
|---|---|---|---|
| PHE | 0.000 | 环碳 6×opls_145 | 芳香（CA 语义）|
| TRP | 0.000 | CG=opls_500、CD1=opls_514、NE1=opls_503 | 双环（5+6）|
| **HISD** | 0.000 | ND1=opls_503(H)、NE2=opls_511 | 中性 His（Nd 质子化）|
| **HISE** | 0.000 | ND1=opls_511、NE2=opls_503(H) | 中性 His（Ne 质子化）|
| **HISH** | **+1.000** | ND1/NE2=opls_512(NA,-0.54)+H(opls_513)、咪唑碳 CG/CD2=opls_510(CX,+0.215)、CE1=opls_509(CR,+0.385) | **质子化 His（正电）** |
| **LYS** | **0.000** | NZ=opls_900(-0.90)+2H | **中性 ε-胺** |
| **LYSH** | **+1.000** | NZ=opls_287(-0.30)+3H(opls_290) | **质子化 ε-铵** |
| **ARG** | **+1.000** | CZ=opls_302、NH1/NH2=opls_300 | 胍基阳离子 |
| **ASP/GLU** | **-1.000** | CG/CD=opls_271、OD/OE=opls_272 | 羧酸根 |
| PRO | 0.000 | CA/CB/CG/CD | 五元环（非芳香）|

### 2.3 末端 patch（n.tdb/c.tdb 实测）

- N 端 `GLY-NH3+`：N=opls_287(-0.30) + 3H(opls_290 +0.33) → **N+3H=+0.69**（结构同 CHARMM/GROMOS ）
- N 端 `PRO-NH2+`：N=opls_309(-0.20) + 2H(opls_310 +0.31) → **N+2H=+0.42**
- C 端 `GLY-COO-`：C=opls_271(+0.70) + 2×opls_272(-0.80) → **净 -0.90**
- 封端：ACE（乙酰化 N 端）

### 2.4 水模型（rtp + ion）

- `HOH`（opls_116/117，SPC，-0.82/+0.41）
- `HO4`（TIP4P）、`HO5`（TIP5P）
- 离子：`F/CL/BR/LI/NA/K/RB/CS/MG/CA/SR/BA`（rtp 尾）

---

## 3. 逐项识别方案（对照 amber，改动分类）

### A. 零改动复用（力场无关）
供体（`_classify_dh_pair` 元素判断）、卤键、金属、金属配位、配体带电（官能团层）、疏水（元素+邻居 C/H）

### B. 仅换类型表（填表，用 opls 编号）
1. **`OPLS_ACCEPTOR_TYPES`**（新建）
   ⚠️ **符号来源必须是 GROMACS 实际符号（ffnonbonded 第 2 列 + rtp 实测），不能照抄 Towhee**（Towhee 的 OHa/OHm/OHp 在 GROMACS 中不存在，实际为 `OH`；S/SH 是元素非符号）。
   GROMACS OPLS 受体符号：氧 `O/O2/OH/OS/OW/ON/OU/OY/OL/O_2/O_3` + 氮 `N/N2/N3/NA/NB/NC/NO/NT/NY/NZ` + 硫（元素）+ 卤素（元素 F/Cl/Br/I）→ **映射为对应 opls 编号全集**（实现时从 ffnonbonded 按符号生成，勿手写）
2. **`OPLS_STRONG_AROMATIC`**（新建）
   符号：`CA`（opls_145）+ `C*/CB/CN/CW/CV/CR`（opls_500-514 系）+ `NA/NB`（opls_503/511）→ **opls 编号全集**（实测：opls_145/500/501/502/503/506/507/508/509/510/511/512/514）
3. **`WATER_RESIDUES`** → 加 `HOH`/`HO4`/`HO5`（+ 现有 SOL/WAT 兼容）

### C. 改字典（数据）
- **正电**：`ARG`（opls_302/300，胍基）、`LYSH`（opls_287 + 3H）、`HISH`（opls_512 + H + 咪唑碳）
- **负电**：`ASP`（opls_271/272）、`GLU`（opls_271/272）
- **中性不进字典**：LYS（0.0！）、ARGN、ASPH、GLUH、HISD、HISE
- ⚠️ **与 Amber 差异**：正电用 **LYSH**（非 LYS）、HISH（非 HIP）

### D. 结构判据（复用 CHARMM 的 `_is_terminal_ammonium`）
N 端正电/侧链铵（LYSH）识别：**N 邻居∈{C,H} + 重原子邻居不连 ≥2N** + N+H 净电荷验证——**跨力场直接复用**（已验证 CHARMM 场景）

---

## 4. 实现结构

```
DuIvyInteractions/group_identifiers/opls_ff_identifier.py   # 新文件
DuIvyInteractions/group_identifiers/__init__.py              # 注册表加 "opls"
```

```python
IDENTIFIER_CLASSES = {
    "amber": ...,  "gromos": ...,  "charmm": ...,
    "opls": OplsFFGroupIdentifier,   # 新增
}
```

类 `OplsFFGroupIdentifier(AmberFFGroupIdentifier)`，覆写：`_find_acceptors`（换 opls 表）、`_identify_protein_charged`（换字典 ARG/LYSH/HISH + ASP/GLU）、`_filter_aromatic_rings`（换 opls 表）、`_identify_functional_group_charged`（N 端结构判据，同 CHARMM）。

> **`_is_terminal_ammonium` 复用方式**：因它是 `@staticmethod`，OPLS 子类可直接调用（`CharmmFFGroupIdentifier._is_terminal_ammonium(...)` 或复制）。建议**提升到公共模块**（三力场共用），待实现时定夺。

---

## 5. 范围与里程碑

**第一版范围**：蛋白（标准残基 + His 三态 + 带电 + 末端 patch）；**不含**核酸/配体（OPLS 配体参数化靠 LigParGen，属后续）。

**里程碑**（合并 dev 判据）：
1. `opls_ff_identifier.py` + 注册表 → `dii run --ff opls` 可跑通
2. 非芳香类真实数据跑通（含 LYSH 带电、N 端/C 端）
3. 芳香类（π-π/π-阳离子）通过 opls CA/C*/NA/NB 表识别
4. 单元测试绿（`test_opls_identifier.py`：合成 SystemData，N 端/LYSH/HISH/ARG/ASP/GLU + PHE/TRP/HIS 环）

---

## 6. 风险与边界（对抗性预审）

| 项 | 说明 |
|---|---|
| **opls 编号 ≠ 语义** | 805 个类型，但语义由符号决定（ffnonbonded 第2列 + Towhee）。**必须用符号驱动建表**，不能手工猜编号 |
| **LYS 中性坑** | 若误用 LYS（0.0）进正电字典 → 盐桥漏判/误判。**必须用 LYSH**（同 GROMOS 教训）|
| **HISH 原子集** | HISH 咪唑 N（opls_512 -0.54）负电荷 → 带电基团须含咪唑碳（opls_510/509）才显正电（同 CHARMM HSP 教训）|
| **受体卤素** | OPLS 卤素是元素 `F/CL/BR/I`（非有符号类型）→ 受体表按元素加 |
| **N 端结构判据** | 复用 `_is_terminal_ammonium`；已验证 LYSH(287)/ARG(300) 场景区分 ✓ |
| **5 元环芳香** | TRP/HIS 5 元环类型（**C*/CB/CN/CW/CV/CR + CX**（HIP 专用）/NA/NB）必须全在 STRONG_AROMATIC；CX 仅 HISH 咪唑 CG/CD2 用，CE1 用 CR |

---

## 7. 证据来源（双重验证）

| 论断 | 本地证据 | 权威来源 |
|---|---|---|
| 版本 = OPLS-AA/L(2001) | forcefield.doc "OPLS-AA/L (2001 aminoacid dihedrals)" | [GROMACS 手册](https://manual.gromacs.org/documentation/2022/reference-manual/topologies/force-field-organization.html)、[Robertson 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4504185/) |
| opls 类型=数字 | oplsaa.ff 805 类型 atp | [Towhee OPLS-aa](https://towhee.sourceforge.net/forcefields/oplsaa.html) |
| 符号映射 | ffnonbonded.itp 第2列 | Towhee（OPLS-aa 用 Amber 命名约定）|
| 符号→语义 | —— | [Towhee 类型表](https://towhee.sourceforge.net/forcefields/oplsaa.html)（CA/C*/NA/NB/N2/N3/O2 逐个语义）|
| **CX = HIP 咪唑碳** | ffbonded.itp "jtr: HIP CB-CG" + rtp HISH CG/CD2=opls_510 | GROMACS 移植扩展（Towhee 无 HIP 蛋白实现故无 CX）|
| LYS 中性/LYSH +1 | rtp 净电荷核算 | Towhee `k+`=protonated lysine |
| HISH=HISP | rtp `[ HISH ] ; also known as HISP` | Towhee `h+`=both N protonated |
| N 端 NH3+ 电荷 | n.tdb（287/-0.30+3×0.33）| Towhee N3+/NT+1 ammonium |
| 配体参数化 | —— | [LigParGen（Jorgensen 组）](https://jorgensenresearch.com/ligpargen) |

本地力场路径：`D:\CharlieAPP\gmx2019_05_GPU\share\gromacs\top\oplsaa.ff\`

---

*文档结束*