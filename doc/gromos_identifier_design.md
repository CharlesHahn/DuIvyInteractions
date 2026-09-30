# GROMOS 力场基团识别器设计方案

> 创建日期：2026-10-07
> 分支：feature/gromos-identifier
> 状态：方案定稿（待实现）
> 用途：新力场入库——在现有 "新力场入库 = 填特征表" 架构承诺下，评估并设计 GROMOS 54A7 的基团识别器

---

## 1. 背景

DII 的基团鉴定依赖四类证据源（见 `amber_ff_identifier.py`）：

1. 力场类型名语义（STRONG_AROMATIC / ACCEPTOR_TYPES）
2. 键合图（环检测、官能团邻居、卤键连接、疏水邻居、D–H 对）
3. 元素 + 电荷（供体 D 元素、q(H)>0、q<0、净电荷、金属元素）
4. 残基/原子名（带电残基字典、水残基名）

当前映射表基于 Amber 家族（GAFF + ff14SB 等）。本设计评估 GROMOS 54A7 能否按"填特征表"接入，并给出逐项识别方案。

**前置调研**：`doc/force_field_compatibility_survey.md`（三力场兼容性调研）结论——GROMOS 是唯一"类型名不编码化学环境"的力场，对应需降级证据源；本设计为该结论的实现落地。

**实证依据**：本地 GROMACS 2019（`D:\CharlieAPP\gmx2019_05_GPU\share\gromacs\top\gromos54a7.ff\`）的 `atomtypes.atp`、`aminoacids.rtp`、`spc.itp` 实际文件内容，与 GitHub GROMACS master 逐字一致。

---

## 2. GROMOS 54A7 关键事实（本地 rtp/atp 实测）

### 2.1 原子类型表（atomtypes.atp）

| 类型 | 语义 | 备注 |
|---|---|---|
| O / OM / OA / OE / OW | 羰基 / 羧基 / 羟基·糖·酯 / 醚·酯 / 水 氧 | |
| N / NT / NL / NR / NZ / NE | 肽胺 / 终端 NH2 / 终端 NH3 / **芳香氮** / Arg NH2 / Arg NE | NR 是芳香环氮 |
| C | **裸碳（bare carbon）** | 羰基碳、羧基碳、**芳香环碳全部同名** |
| CH0 / CH1 / CH2 / CH3 / CH4 / CH2r | sp3 脂肪碳（UA） | CH2r = 环内 CH2 |
| CR1 | 芳香 CH 组 | 类型表存在，但蛋白 rtp 实际未用 |
| HC / H | 碳上氢 / 非碳氢（极性 H） | 极性/芳香 H 显式 |
| S / P / F / CL / BR | 硫 / 磷 / 卤素 | |
| CU1+ / CU2+ / FE / ZN2+ / MG2+ / CA2+ / NA+ / CL- | 金属/离子 | 元素可解析 |

### 2.2 蛋白芳香残基实测（aminoacids.rtp）

| 残基 | 环碳类型 | 环氮类型 | 环氢类型 |
|---|---|---|---|
| PHE | **C ×6** | — | HC |
| TYR | **C ×6** | — | HC |
| TRP | **C ×8**（含五元环）| **NR ×1**（NE1）| HC + H(HE1 极性) |
| HISA/HISB/HISH/HIS1/HIS2 | **C ×3**（咪唑）| **NR ×2**（ND1/NE2）| HC + H(HD1 极性) |

**结论**：GROMOS 蛋白芳香碳类型全为 `C`，与羰基/羧基碳同名；芳香性**不在类型中编码**，靠 improper 二面角几何维持。类型语义证据源（Amber 的核心）在 GROMOS 下失效。

### 2.3 其他命名

- 带电残基：**ARG / LYS / ASP / GLU**（标准名）+ **HISH**（质子化 His）+ 中性变体（ARGN/LYSH/ASPH/GLUH）
- 水：**SOL**（spc.itp 实测，OW/HW1/HW2；另有 spce/tip3p/tip4p itp）
- 离子残基：`NA+ / CL- / MG2+ / ZN2+ / CA2+` 等
- 显式 H：极性 H（N-H、O-H）与芳香 C-H 全部显式；仅脂肪族非极性 H 并入 CH1/CH2/CH3（UA）

---

## 3. 逐项基团识别方案

| # | 基团 | Amber 判据 | GROMOS 判据 | 做法 |
|---|---|---|---|---|
| 1 | H_donor | D–H 键 + q(H)>0，D∈{N,O,S,F} | **同 amber**（极性 H 显式） | 复用 `_classify_dh_pair`，零改动 |
| 2 | H_acceptor | 类型 ∈ ACCEPTOR_TYPES + q<0 | 类型 ∈ GROMOS 受体表 + q<0 | 新建 `GROMOS_ACCEPTOR_TYPES`（见 §4），逻辑同 amber |
| 3a | 蛋白带电 | 残基名字典 + 原子名 | 字典键改 GROMOS 名（ARG/LYS/ASP/GLU/HISH）| 需扩展字典 |
| 3b | 配体带电 | 官能团模式（元素+邻居） | **同 amber**（力场无关） | 复用，零改动 |
| 4 | 卤键供/受体 | 卤素+邻居 | **同 amber**（力场无关） | 复用，零改动 |
| 5 | 金属 | 元素 ∈ METAL_IONS | **同 amber**（GNOMOS 离子类型元素可解析） | 复用，零改动 |
| 6 | 水 | 残基名 ∈ WATER_RESIDUES | **SOL**（已在集合内） | 零改动 |
| 7 | 芳香环 | 类型 ∈ STRONG_AROMATIC + 环检测 | **降级：残基名 ∈ 芳香四类 + 环检测 + 环原子全 C/NR**（见 §5） | 需新实现 |
| 8 | 疏水 | C + 邻居 ∈ {C,H} | **重写：C/CH* + 邻居 ∉ {O,N,S}**（反列举） | 需新实现 |
| 9 | 金属配位 | 元素 ∈ {O,N,S} 非水 | **同 amber**（力场无关） | 复用，零改动 |

---

## 4. GROMOS 受体类型表（新建）

```python
GROMOS_ACCEPTOR_TYPES = frozenset({
    "O", "OM", "OA", "OE", "OW",       # 氧：羰基/羧基/羟基/醚/水
    "N", "NT", "NL", "NR", "NZ", "NE", # 氮：肽胺/终端/芳香/胍基
    "S",                                # 硫
    "F", "CL", "BR",                   # 卤素
})
```
判定逻辑与 amber `_find_acceptors` 相同（类型 ∈ 表 且 q < 0）。实测电荷支持：主链 O(-0.38)、ASP OD(-0.635)、ARG NH1(-0.26)、TRP NE1(-0.10) 均负。

---

## 5. 芳香环降级方案（核心难点）

**问题**：GROMOS 芳香碳类型 = `C`（与羰基/羧基碳同名），STRONG_AROMATIC 类型表失效。

**方案（第一版）**：
1. **环检测**：复用 amber `_detect_rings`（残基内 BFS）。
2. **环成分过滤**：环内原子类型全 ∈ {`C`, `NR`}（GROMOS 环碳/环氮）。
3. **残基名前置过滤**：仅当残基 ∈ {PHE, TYR, TRP, HISA, HISB, HISH, HIS1, HIS2} 时，该残基内的 `C`/`NR` 环判定为芳香环。

> GROMOS 蛋白只有这 4 类残基含芳香环，残基名是强信号；第一版不做非蛋白芳香环（GROMOS 配体参数化靠 ATB，第一版不支持配体）。

**不采用**（第二版候选）：几何平面性兜底（PLIP 风格，环原子到拟合平面最大偏差 ≤ 阈值）——用于非蛋白芳香环，如核酸碱基。

**环序约定**：`Group.atoms` 仍按 BFS 路径序输出，检测器（PiStacking/π-cation 的 ring normal/center）与力场无关，零改动。

---

## 6. 疏水重写方案

- amber 判据："C 且所有邻居 ∈ {C, H}"——GROMOS UA 无碳上 H 邻居，失效。
- **GROMOS 判据**："C（含 CH0/CH1/CH2/CH3/CH4/CH2r）且所有邻居 ∉ {O, N, S}"（反列举排除极性取代）。
- 原子粒度按 UA 碳原子计，检测器距离判据不变。

---

## 7. 架构与注册

```
DuIvyInteractions/group_identifiers/gromos_identifier.py   # 新文件
DuIvyInteractions/group_identifiers/__init__.py             # 注册表加一项
```

```python
IDENTIFIER_CLASSES = {
    "amber": AmberFFGroupIdentifier,
    "gromos": GromosFFGroupIdentifier,   # 新增
}
```

`dii run --ff gromos` 即可用；pipeline.py / DII.py / 全部检测器零改动。

**功能边界（明确声明）**：
- ✅ 可做：供体/受体/带电（蛋白）/卤键/金属/水/金属配位 + 降级芳香环 + 重写疏水
- ❌ 第一版不支持：GROMOS 配体（无 GAFF 等价物，需 ATB 工具链，属未来工作）

---

## 8. 可靠性与风险

| 项 | 置信度 | 说明 |
|---|---|---|
| 供体/受体/金属/水 | 高 | 显式 H + 电荷，与 amber 同证据强度 |
| 带电 | 高 | GROMOS 残基命名已实测（ARG/LYS/ASP/GLU/HISH）|
| 芳香环 | 中 | PHE/TYR 仅靠残基名+环检测，弱于 amber 类型验证；非天然环不可识别 |
| 疏水 | 中 | 反列举语义，边缘情况（如 C 连卤素）需测试 |
| 配体 | — | 第一版不支持 |

---

## 9. 里程碑（合并 dev 判据）

1. `gromos_identifier.py` + 注册表 → `dii run --ff gromos` 可跑通不报错
2. 非芳香类（供体/受体/盐桥/卤键/金属/水桥/金属配位）真实数据跑通
3. 芳香类（π-π / π-阳离子）降级方案跑通（PHE/TYR/TRP/HIS 环检出）
4. 单元测试绿（`Tests/unittests/test_gromos_identifier.py`：GROMOS 54A7 体系基团计数断言 + 与 Amber 结果交叉验证）→ 合并 dev

---

## 10. 关联文档

- `doc/force_field_compatibility_survey.md` —— 三力场调研（§5 结论为本设计依据）
- `doc/paper_framework_v1.md` / `doc/why_force_field_group_identification.md` —— 论文框架与差异化论证
- `CLAUDE.md` —— "新力场入库 = 填特征表"架构承诺

---

*文档结束*