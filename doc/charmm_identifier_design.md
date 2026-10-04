# CHARMM 力场基团识别器设计方案（charmm_ff_identifier）

> 创建日期：2026-10-07
> 分支：feature/charmm-identifier
> 状态：方案定稿（待实现）
> 参考实现：`amber_ff_identifier.py`（CHARMM 与 Amber 同为"类型名强编码化学"的全原子力场，走"填特征表"路线，无需 GROMOS 的逻辑重写）
> 实证依据：本地 `charmm36-feb2026_cgenff-5.0.ff`（MacKerell 官方 2026-02 版）+ CHARMM-GUI 官方 `top_all36_prot.rtf` 交叉验证一致

---

## 1. 背景

DII 的基团鉴定依赖四类证据源（见 `amber_ff_identifier.py`）：类型名语义、键合图、元素+电荷、残基/原子名。

CHARMM36/C36m 是**全原子力场**，原子类型**强编码化学环境**（`CA`=芳香碳、`CT1/2/3`=sp3 脂肪碳、`NR1/2/3`=咪唑氮、`NC2`=胍氮…），与 Amber 同构 → **"新力场入库 = 填特征表"成立**，识别器参考 amber 实现。

**与 GROMOS 的本质区别**（决定实现复杂度）：
| 维度 | CHARMM | GROMOS |
|---|---|---|
| 显式 H | ✅ 全原子 | ⚠️ 仅极性/芳香（UA）|
| 芳香碳类型 | ✅ CA（强编码）| ❌ C（通名）|
| 疏水判据 | ✅ 正列举（邻居∈{C,H}）| ⚠️ 需反列举 |
| 识别器逻辑 | **仅换表 + 改字典** | 需重写芳香/疏水 |

---

## 2. 关键事实（本地 charmm36 5.0 实测 + 官方 rtf 验证）

### 2.1 原子类型表（atomtypes.atp，蛋白 C36m 部分）

| 类型 | 语义 | 证据 |
|---|---|---|
| `CA` | **芳香碳** | atp + rtf `MASS -1 CA ! aromatic C`；PHE 环碳全 CA |
| `CAI` | TRP 桥头芳香碳 | rtf `aromatic C next to CPT in trp` |
| `CPH1/CPH2` | HIS 咪唑碳（CG/CD2 / CE1）| rtf `his CG and CD2 / CE1 carbon` |
| `CPT/CY` | TRP 环碳（桥 / 吡咯）| rtf `trp C between rings / TRP C in pyrrole ring` |
| `CT1/CT2/CT3/CT` | sp3 脂肪碳（CH/CH2/CH3/无H）| rtf `aliphatic sp3 C`；MATCH 论文 |
| `CP1/CP2/CP3` | Pro 四面体碳 | rtf |
| `CC/CD` | 羧基/酰胺碳（Asp/Glu 侧链）| rtf `carbonyl C, asn,asp,gln,glu` |
| `C` | 肽主链羰基碳 | rtf `carbonyl C, peptide backbone` |
| `NR1/NR2/NR3` | 中性 His 质子化/未质子化/正电咪唑氮 | rtf `neutral his protonated/... / charged his ring nitrogen` |
| `NY` | TRP 吡咯氮 | rtf `TRP N in pyrrole ring` |
| `NH1/NH2/NH3` | 肽胺/酰胺/铵氮 | rtf |
| `NC2` | 胍氮（Arg）| rtf `guanidinium nitrogen` |
| `N` | Proline N | rtf |
| `HP` | 芳香 H | rtf `aromatic H` |
| `H/HA/HB1/HA1-3/HR1-3/HC/HS` | 极性/烷基/主链/His 特化 H | rtf |
| `O/OB/OC/OH1/OS` | 羰基/醋酸/羧酸根/羟基/酯氧 | rtf |
| `S/SM/SS` | 硫/二硫/硫醇根 | rtf |

> 以上与 CHARMM-GUI `top_all36_prot.rtf` 的 MASS 表**逐字一致**（本地 GROMACS 移植版 = charmm2gmx 官方转换，来源可信）。

### 2.2 蛋白残基实测（aminoacids.rtp，逐原子）

| 残基 | 净电荷 | 关键原子 | 说明 |
|---|---|---|---|
| PHE | 0.000 | 环碳 6×`CA`、环氢 `HP` | 全 C 环 + CA 类型 |
| TRP | 0.000 | CG=`CY`、NE1=`NY`、CE2/CD2=`CPT`、苯环=`CA/CAI` | 双环（5+6）|
| HSD | 0.000 | ND1=`NR1`(H)、NE2=`NR2` | 中性 His（ND1 质子化）|
| HSE | 0.000 | ND1=`NR2`、NE2=`NR1`(H) | 中性 His（NE2 质子化）|
| **HSP** | **+1.000** | ND1/NE2 均=`NR3`(-0.51)、咪唑碳=`CPH1/CPH2`(CG/CE1/CD2 正) | **质子化 His（正电）**；⚠️ 咪唑 N 负电荷，带电基团须含咪唑碳（CG/CE1/CD2）方为正电荷中心 |
| **LYS** | **+1.000** | NZ=`NH3`(HZ1/HZ2/HZ3) | **ε-铵 NH3+（正电）** |
| LSN | 0.000 | NZ=`NH2` | 中性 Lys（变体）|
| **ARG** | **+1.000** | CZ=`C`、NE/NH1/NH2=`NC2` | 胍基阳离子 |
| ARGN | 0.000 | — | 中性 Arg |
| **ASP** | **-1.000** | CG=`CC`、OD1/OD2=`OC` | 羧酸根 |
| ASPP | 0.000 | — | 质子化 Asp |
| **GLU** | **-1.000** | CD=`CC/CD`、OE1/OE2=`OC` | 羧酸根 |
| GLUP | 0.000 | — | 质子化 Glu |
| **CYM** | **-1.000** | SG=`SS` | 硫醇根 Cys⁻ |
| PRO | 0.000 | CA/CB/CG/CD=`CP1/CP2/CP3` | 五元环（非芳香）|

### 2.3 末端/封端（patches，n.tdb/c.tdb 实测）

- N 端 `NH3+`（标准）：N=`NH3`（-0.30）+ 3×H(`HC` +0.33) → **N+3H 净 +0.69**（正电）
- Pro N 端 `NH2+`：N=`NP`（-0.07）+ 2×H → 正电
- 中性 N 端 `NH2`：N=`NH2`（-0.96），中性
- C 端 `COO-`（标准）：C=`CC`（+0.34）+ 2×OT(`OC` -0.67) → **净 -1.0**（负电）；`COOH` 质子化变体中性
- 乙酰化 N 端：`ACE` 残基（CH3-CO-，中性）

### 2.4 水与离子（solvent.rtp + watermodels）

- 水：`TIP3`（OH2=`OT`, H=`HT`，-0.834/+0.417）、`HOH`（同 TIP3 原子集）；tip3p.itp 残基名 `SOL`；另有 TIP4P/5P/SPC/SPCE
- 离子：Metals/solvent.rtp 大量金属（`AG1P...YB3P`、`CAL/MG/NA/POT/CLA` 等）

### 2.5 CGenFF 配体类型（cgenff.rtp 9.3k 行，285+ 类型）

类型名**强编码**：`CG2R61`（环6元sp2芳香碳）、`CG2R51/52/53`（5元环）、`CG2DC1/2`（共轭烯）、`CG2O1-7`（羰基）、`NG2R*`（环氮）、`OG3*`（氧）等。
> 配体识别可扩展性强，但**第一版不做配体**（同 GROMOS 边界，蛋白优先）。

---

## 3. 逐项识别方案（对照 amber，改动分类）

### A. 零改动复用（力场无关）
- H 键供体（`_classify_dh_pair`：D∈{N,O,S,F} + q(H)>0）— 全原子显式 H ✓
- 卤键供/受体、金属、金属配位、配体带电（官能团层）
- 水判断逻辑（`WATER_RESIDUES` 需加 `TIP3`/`HOH`，见 B）

### B. 仅换类型表（填表）
1. **`CHARMM_ACCEPTOR_TYPES`**（新建）
   `O/OB/OC/OH1/OS/OW/OT` + `N/NH1/NH2/NH3/NC2/NR1/NR2/NR3/NY` + `S/SM/SS` + **CGenFF 卤素**（`FGA1-3/FGR1/CLGA1/CLGA3/CLGR1/BRGA1-3/BRGR1/IGR1`，卤素可作 H 键受体，Lin 2017）+ CGenFF 氧/氮（如做配体）
   ⚠️ **不可用 GAFF 裸名 `F/CL/BR/I`**：CHARMM 中 `CL` 是羰基碳、裸卤素不存在
2. **`CHARMM_STRONG_AROMATIC`**（新建，替代 amber ACCEPTOR 芳香表）
   `CA/CAI/CPH1/CPH2/CPT/CY`（蛋白）+ `NR1/NR2/NR3/NY`（芳香氮）+ CGenFF `CG2R51-67/NG2R50-67`（配体，扩展用）
3. **`WATER_RESIDUES`** → 加 `TIP3`、`HOH`
4. **疏水**：全原子显式 H → **复用 amber 正列举**（C + 邻居∈{C,H}），类型集合换 CHARMM 碳类（CA 需排除或保留？见 §6 边界）

### C. 改字典（数据）
- **正电字典**：
  ```python
  "LYS": ["NZ","HZ1","HZ2","HZ3"],   # ε-铵 NH3+（同 Amber）
  "ARG": ["CZ","NE","NH1","NH2","HE","HH11","HH12","HH21","HH22"],
  "HSP": ["ND1","NE2","HD1","HE2","CG","CE1","CD2"],  # 质子化 His；须含咪唑碳（N 负电荷，4原子集净 -0.14）
  ```
- **负电字典**：`"ASP": ["CG","OD1","OD2"]`、`"GLU": ["CD","OE1","OE2"]`、`"CYM": ["SG"]`
- **中性变体不进字典**：LSN/ARGN/ASPP/GLUP/HSD/HSE
- ⚠️ **与 Amber 的命名差异**：HIS 用 HSD/HSE（中性）+ HSP（正电）；Amber 用 HIP（正电）/HID/HIE。字典键按 CHARMM 名。
- 🔴 **N 端正电需覆写官能团层**：因 CHARMM N 端 N 带负电荷（NH3 -0.30），继承的 tertamine 单原子验证失效；实现时覆写为 **N + 键连 H 净电荷** 验证（见 §6）。

### D. 重写逻辑
**无需重写**：
- 芳香环：类型表直接可用（CA/CPH*/CPT/CY/NR* 强编码）→ 复用 amber `_filter_aromatic_rings` 逻辑，仅换表
- 疏水：全原子显式 H → 复用 amber 正列举
- 带电去重：继承父类即可（CHARMM 咪唑/铵 N 全为负电荷 ND1/NE2=-0.5、NZ(NH3)=-0.3 → tertamine 单原子净电荷验证不通过 → **无双识别问题**，无需 GROMOS 式子集去重覆写）

---

## 4. 实现结构

```
DuIvyInteractions/group_identifiers/charmm_ff_identifier.py   # 新文件
DuIvyInteractions/group_identifiers/__init__.py               # 注册表加 "charmm"
```

```python
IDENTIFIER_CLASSES = {
    "amber": AmberFFGroupIdentifier,
    "gromos": GromosFFGroupIdentifier,
    "charmm": CharmmFFGroupIdentifier,   # 新增
}
```

类 `CharmmFFGroupIdentifier(AmberFFGroupIdentifier)`，覆写：`_find_acceptors`（换表）、`_identify_protein_charged`（换字典）、`_filter_aromatic_rings`（换表，逻辑可复用 amber 式三条件）。

---

## 5. 范围与里程碑

**第一版范围**：蛋白（标准残基 + His 三态 + 带电 + 末端 patch 经官能团层）；**不含**核酸/脂质/糖/配体（力场文件已含，属后续扩展）。

**里程碑**（合并 dev 判据）：
1. `charmm_ff_identifier.py` + 注册表 → `dii run --ff charmm` 可跑通
2. 非芳香类（供体/受体/盐桥/卤键/金属/水桥/金属配位）真实数据跑通
3. 芳香类（π-π/π-阳离子）通过 CA/CPH/CPT/CY/NR 类型表识别（PHE/TYR/TRP/HSP 环）
4. 单元测试绿（合成 SystemData，复用 GROMOS 测试模式：`test_charmm_identifier.py`；**含 N 端 NH3+ 正电识别断言、C 端 COO- 负电断言**）

---

## 6. 风险与边界（对抗性预审）

| 项 | 说明 |
|---|---|
| **CA 双义坑** | 原子名 `CA`（α碳）≠ 类型名 `CA`（芳香碳）——DII 读 `atom_type` 字段，天然免疫（同 Amber `C` 双身份）|
| **疏水含 CA 环碳** | CHARMM 芳香环碳（CA）邻居为 C/H → 会命中疏水判据（与 Amber 行为一致）；TODO #2 去重范畴，测试锁定现状 |
| **🔴 N 端正电识别（P1 缺陷）** | CHARMM N 端 NH3+/NH2+ 的 **N 原子带负电荷**（`NH3` -0.30、`NP` -0.07，n.tdb 实测）→ 继承的 tertamine 单原子净电荷验证（>0.1）**失败 → N 端正电漏识别**。**Amber 无此问题**（N 端 N3 带正电荷 0.03~0.29）。**修正方案**：覆写官能团层验证逻辑为 **N + 键连 H 的净电荷**（NH3+: -0.30 + 3×0.33 = **+0.69** > 0.1 ✓，兼容 Amber）|
| **C 端识别正常** | COO-（C=+0.34 + 2×OT -0.67）→ carboxylate 验证三原子净电荷 **-1.0** < -0.1 ✓（多原子验证天然通过，不受影响）|
| **双识别** | VERIFIED 无风险：CHARMM 正电残基 N 均为负电荷（NH3 -0.3、NC2 -0.7、NR3 -0.51）→ tertamine 单原子验证必不通过 → 字典层唯一 |
| **C36 vs C36m** | 类型一致（par_all36m_prot 仅改参数不改类型），识别器兼容两者 |
| **水模型** | TIP3（OT/HT）为 CHARMM 默认；SOL/HOH/TIP3 残基名均需入库 |

## 7. 原子类型交叉验证（本地 + 官方）

设计依赖的全部原子类型均经 **本地 `charmm36-feb2026_cgenff-5.0.ff` 与 CHARMM-GUI 官方 `top_all36_prot.rtf` 逐条比对一致**：

| 类型 | 本地 atp 注释 | 官方 rtf MASS 注释 | 结论 |
|---|---|---|---|
| `CA` | aromatic C | aromatic C | ✅ |
| `CAI` | aromatic C next to CPT in trp | 同 | ✅ |
| `CPH1/CPH2` | his CG/CD2 / CE1 | 同 | ✅ |
| `CPT/CY` | trp C between rings / TRP C in pyrrole | 同 | ✅ |
| `CT1/CT2/CT3/CT` | aliphatic sp3 C for CH/CH2/CH3/no H | 同 | ✅ |
| `CC/CD` | carbonyl C, asn/asp/gln/glu | 同 | ✅ |
| `NR1/NR2/NR3` | neutral his prot/unprot / charged his ring N | 同 | ✅ |
| `NY` | TRP N in pyrrole ring | 同 | ✅ |
| `NH1/NH2/NH3` | peptide/amide/ammonium N | 同 | ✅ |
| `NC2` | guanidinium nitrogen | 同 | ✅ |
| `N` | proline N | 同 | ✅ |
| `HP` | aromatic H | 同 | ✅ |
| `O/OB/OC/OH1/OS` | carbonyl/acetic/carboxylate/hydroxyl/ester O | 同 | ✅ |
| `S/SM/SS` | sulphur/disulfide/thiolate S | 同 | ✅ |
| 末端 `NH3`(N端)/`NP`(Pro N端)/`CC`+`OT`(C端) | n.tdb/c.tdb patch 实测 | rtf PRES 定义 | ✅ |

---

## 7. 证据来源

| 论断 | 本地证据 | 权威来源 |
|---|---|---|
| CA=芳香 C、CT1-3=脂肪 C | `charmm36.../atomtypes.atp` | CHARMM-GUI `top_all36_prot.rtf` MASS 表；[MATCH 论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC3228871/) |
| NR1/NR2/NR3 咪唑氮 | rtp HSD/HSE/HSP | 同上 rtf |
| HSP=+1、LYS=+1、LSN=0 | rtp 净电荷核算 | 同上 rtf（RESI HSD/HSE/HSP/LYS/LSN）|
| ASP/GLU=-1、CYM=-1 | rtp 净电荷 | 同上 rtf（RESI ASP/GLU/CYM）|
| CGenFF 命名规则 | cgenff.rtp + 参数文件 | [Vanommeslaeghe 2012, JCIM](https://doi.org/10.1021/ci300363c)；[CGenFF 官网说明](https://docs.silcsbio.com/2022.2/cgenff/cgenff.html) |
| 下载来源 | 本地 4 包 | [MacKerell Lab 官网](http://mackerell.umaryland.edu/charmm_ff.shtml)（2026-02 版）|

本地力场路径：`D:\CharlieAPP\gmx2019_05_GPU\share\gromacs\top\charmm36-feb2026_cgenff-5.0.ff\`

---

*文档结束*