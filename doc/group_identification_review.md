# 基团识别结果自动核验（TODO #3）

> 状态：核查完成（自动化 diff + 8 个子 agent 文献核验 + 对抗性审查）
> 日期：2026-10-14
> 分支：feature/group-identification-review

## 1. 目标

TODO #3：基团识别结果需要人工核验——单元测试基于代码输出（自举），需独立的
化学事实基准验证识别器输出是否与真实体系化学一致。

**核查目标**（拆分为三个子目标）：
1. **规则符合化学**：识别器的类型→基团映射规则是否化学正确
2. **实现忠实规则**：代码是否忠实执行规则（无漏检/误报）
3. **标准答案锚定**：核验通过的清单固化为 ground truth，回归测试断言独立于代码

## 2. 方法：独立化学基准 + 自动 diff

### 2.1 独立基准

以"残基名 + 原子名 + 键合图"推导期望基团（**不读识别器任何代码**），依据标准
氨基酸化学事实：

| 基团 | 期望规则（化学事实） |
|:-----|:-----|
| `hydrophobic` | C 且所有邻居 ∈ {C, H} |
| `H_donor` | N/O/S 与 ≥1 个 H 键连 |
| `H_acceptor` | O 一律（含羟基 OH，有孤对电子）；N 无 H 邻居且为吡啶/亚胺型（`nb/NB/n/n2/na/NA`）——含 Pro 亚胺 N、His 无 H 吡啶 N |
| `charged_positive` | LYS NZ、ARG NE（基团中心原子）+ 蛋白 N 末端主链 N（≥2 H 邻居） |
| `charged_negative` | ASP CG、GLU CD（中心原子）+ 蛋白 C 末端羧基（C 连 ≥2 O） |

- **代表原子约定**：与 `group_atom_convention.md` 对齐——H_donor/H_acceptor
  代表原子为 atoms[0]（D/A）；charged 为基团中心原子（如 LYS=NZ、ASP=CG）
- **范围**：20 标准氨基酸（含 HIS 质子化态 HID/HIE/HIP）；水/离子/配体跳过

### 2.2 自动 diff

```
识别器输出 G_identified vs 化学基准 G_expected
  缺失 missing = G_expected - G_identified   （识别器漏检 → 错误）
  多余 extra  = G_identified - G_expected    （识别器多检 → 需甄别）
```

- **缺失 = 0** → 识别器无漏检（核心正确性证据）
- **多余 ≠ 0** → 逐类甄别：基准偏窄（识别器正确）或识别器误报（代码 bug）

### 2.3 补充核查（diff 之外）

- **直接核查**：water 组成、aromatic_ring 环大小、metal_binding 候选池元素、
  halogen_acceptor PLIP 定义符合性
- **文献核验**（8 个子 agent 并行）：每个"多余"类别用 anysearch 文献 + 本地
  力场文件裁决"识别器正确 or 误报"
- **对抗性审查**：对每个"疑似缺陷"做反向验证（是否文献支持、是否忠实类型语义）

## 3. 工具

`Tests/tools/verify_identification.py`：

```bash
python Tests/tools/verify_identification.py \
    --tpr <md.tpr> --ff <amber|gromos|charmm|opls> \
    [--report <path>]
```

- 输出 markdown 报告：核验残基数、缺失明细、多余明细（按原子类型分布 + 全量）
- 退出码：`0`（无缺失）/ `1`（有缺失——识别器漏检）

## 4. 结果（三力场真实体系 + OPLS 合成体系）

### 4.1 自动 diff

| 力场 | 体系 | 残基数 | 缺失 | 多余 | 多余构成 |
|:-----|:-----|:---:|:---:|:---:|:-----|
| Amber | KRAS-RBD D927 | 309 | 0 | 32 | PRO N 14 + MET SD 10 + CYS SG 8 |
| GROMOS | 53A6 蛋白+6配体 | 130 | 0 | 5 | MET SD 5 |
| CHARMM | SMO-BST | 495 | 0 | 65 | PRO N 27 + His ND1 7 + CYS/MET S 22 + 其他 9 |
| OPLS | 20 AA 合成（rtp） | 20 | 0 | 5 | PRO N / MET SD / CYS SG / His ND1 |

**缺失 = 0（4 力场）**：识别器无漏检——所有化学上应有的基团（疏水 C、供体
N/O/S、受体、带电）均被识别。

### 4.2 "多余"项文献核验结论（全部为识别器正确）

| 类别 | 文献证据 | 结论 |
|:-----|:-----|:-----|
| PRO 亚胺 N 受体 | Deepak 2016 Biophys J（>15% Pro 参与 N-H···N）；Adhikary 2014 JACS | **正确**（弱受体；GROMOS rtp 实测 N 无 H） |
| MET SD 硫醚受体 | Zhou 2009 Proteins；van Bergen 2016 Sci Rep；Allen 1997 Acta Cryst B | **正确**（弱受体；三力场 SD 均 q<0） |
| CYS SG 硫醇受体 | Mazmanian 2016 JPCB；Gregoret 1991 Proteins；Adhav 2023 QRB Discov | **正确**（-SH 弱受体、Cys⁻ 强受体） |
| HIS ND1 咪唑受体 | Li 2011 JACS（无 H 咪唑 N 可作受体）；Deepak 2016 Biochemistry | **正确**（HIE 型；4 力场已按质子化态区分） |

### 4.3 直接核查通过项

| 项 | 核查方法 | 结果 |
|:-----|:-----|:-----|
| water | 37021 个水基团组成 | 全部 3 原子（O+2H）✓ |
| aromatic_ring | 38/10/86 环大小 | 全部 5/6 元环 ✓ |
| metal_binding | 候选池元素 + 水排除 | 元素 ⊆ {O,N,S}、零水污染（3 力场）✓ |
| halogen_acceptor | PLIP 定义符合性 | 881/1332/337 全部 C/P/S 连 O/N/PS ✓ |

## 5. 对抗性审查：两个"疑似缺陷"被推翻

### 5.1 F1 卤素作 H 键受体 —— 不是缺陷

- **原判断**："有机卤素不作 H 键受体"（引 Dunitz & Taylor 1997）
- **审查后**：该文献标题实为《Organic Fluorine **Hardly Ever** Accepts Hydrogen
  Bonds》（几乎从不，非绝不）；且 30 余年来是争议话题，近年多篇文献明确报告
  F 作受体（Scerba 2011 JOC "NH⁺–F Hydrogen Bonding"；Wang 2024 JOC
  "CF₂–F···H–N⁺ Says Yes"；Champagne 2015 Synthesis 综述）
- **结论**：工具按类型 + 电荷纳入卤素，**忠实于类型语义**；与 PLIP 的差异
  （PLIP 显式剔除）是设计选择差异，**非缺陷**。文献支持尤其是强供体
  （N⁺–H）与 F 的弱氢键。

### 5.2 F2 GNP 嘌呤环漏检 —— 不是缺陷

- **原判断**："GNP 嘌呤环漏检是识别器缺陷，应补 cc/cd/nc 类型"
- **审查后**：GAFF 官方定义明确区分——`ca` = "Sp2 C in **pure aromatic**
  systems"；`cc/cd` = "Sp2 C in **non-pure aromatic** systems"（原始论文：
  "inner sp2 carbon in conjugated ring systems"）；`nc/nd` = "inner sp2 nitrogen
  in conjugated ring systems"（GAFF 论文 Wang 2004 + AMBER 官方档案）
- **工具现状**：cc/cd 在 COMPATIBLE（非强芳香）、nc 不收录——**忠实 GAFF 语义**
- **结论**：按"尊重力场类型编码 + 尊重 Sobtop 输出"原则，工具行为正确；
  改 STRONG 表反而违背类型语义。**非缺陷**（设计边界：`aromatic_ring`
  仅覆盖纯芳香环，不含非纯芳香共轭环——建议文档声明）

## 6. 确认的修复：GROMOS 中性胺受体判据（3c）

**问题**：GROMOS `NT` 类型同时覆盖中性胺（LYS NZ，无羰基 C 邻居，胺 N 孤对
自由）与侧链酰胺（ASN ND2 / GLN NE2，有羰基 C 邻居，孤对离域至羰基）。
原"带 H 一律排除"会连带排除中性胺受体。

**化学依据**：酰胺的化学定义 = N 与羰基 C 相连（孤对离域→非受体）；胺 N
孤对自由→可作受体（教科书 + Eildal 2013 JACS）。

**修复**：对 `NT` 且邻居无羰基 C（C 邻接 O）者保留为受体；其余 N 判据不变。

**验证**：GROMOS 全 20 残基逐一测试——仅 LYS NZ 变为受体，ASN/GLN 酰胺、
ARG 胍基、HIS 芳香 N、LYSH 质子化胺全部不变；真实数据零影响（体系 NT 类型
仅 ASN 酰胺）。新增 2 个单测（LYS NZ 受体 / ASN ND2 非受体）。

## 7. 核查中发现的既有问题

**GROMOS 氢键 two_pass 测试基线过期（287→292）**：`test_n_pairs` 断言 287，
实际 292（+5）。经 `git stash` 回退 + 切 master（v0.0.4）验证——**master 上同样
失败**，属 PBC 改造后未同步更新基线。已登记 `doc/TODO.md` #19。

## 8. 结论

1. **4 力场基团识别忠实于力场类型语义**：无漏检（缺失=0）；"多余"项经文献
   核验全部为识别器正确（基准偏窄）
2. **无确认的识别器缺陷**（F1/F2 经对抗性审查推翻；同残基相互作用按设计保留）
3. **唯一实际修复**：GROMOS 中性胺受体判据（3c，化学定义驱动的结构补偿）
4. **跨力场一致性**：4 力场受体判定方向一致（铵/酰胺/吡咯排除，Pro/吡啶型 N
   保留）；GROMOS 中性胺修复后与其余 3 力场一致
5. **核查工具可复用**：`Tests/tools/verify_identification.py` 支持 4 力场自动 diff

## 9. 关联文档

- `group_atom_convention.md`（代表原子约定）
- `charged_group_design.md`（带电基团设计）
- `acceptor_identification_evidence.md`（受体判定证据，含 GROMOS 中性胺登记）
- `aromatic_ring_definition.md`（芳香环定义；纯芳香 vs 非纯芳香边界）
- `doc/TODO.md`（#3 关联、#19 既有基线过期）
