# PLIP 相互作用鉴定规则

> 来源：PLIP v3.0.1 config.py + DOCUMENTATION.md
> 用途：作为 DuIvyInteraction InteractionDetector 的判据参考

---

## 1. 疏水相互作用（Hydrophobic）

| 参数 | 值 | 说明 |
|:-----|:---|:-----|
| `HYDROPH_DIST_MAX` | 4.0 Å | 疏水原子间最大距离 |

**规则**：两个疏水原子（C + 邻居 ∈ {C,H}）之间距离 ≤ 4.0 Å。

**后处理**：
1. 移除已形成 π-π 堆积的环之间的疏水接触
2. 一个配体原子与同残基多个蛋白原子的接触，只保留最近的
3. 一个蛋白原子与多个配体原子的接触，只保留最近的

**参考文献**：无特定文献，经验值

---

## 2. 氢键（Hydrogen Bond）

| 参数 | 值 | 说明 |
|:-----|:---|:-----|
| `HBOND_DIST_MAX` | 4.1 Å | D-A 最大距离 |
| `HBOND_DON_ANGLE_MIN` | 100° | D-H···A 最小角度 |

**规则**：
- D-A 距离 ≤ 4.1 Å
- D-H···A 角度 ≥ 100°

**约束**：
- 一个供体 D 只能参与一个氢键（多个候选时保留角度最接近 180° 的）
- 受体 A 可参与多个氢键（如分叉氢键）
- 已形成盐桥的原子对不再报告氢键

**参考文献**：Hubbard & Haider, 2001（距离 +0.6 Å，角度 +10°）

---

## 3. π-π 堆积（Pi Stacking）

| 参数 | 值 | 说明 |
|:-----|:---|:-----|
| `PISTACK_DIST_MAX` | 5.5 Å | 环心最大距离 |
| `PISTACK_ANG_DEV` | 30° | 角度偏差上限 |
| `PISTACK_OFFSET_MAX` | 2.0 Å | 最大偏移量 |

**规则**：
- 环心距离 ≤ 5.5 Å
- 检测两环法向量夹角 θ：
  - **平行堆积（P-stacking）**：θ 接近 0° 或 180°，偏差 ≤ 30°
  - **T 型堆积（T-stacking）**：θ 接近 90°，偏差 ≤ 30°
- 环心投影偏移 ≤ 2.0 Å（一环环心投影到另一环平面，投影点到环心的距离）

**参考文献**：McGaughey, 1998

---

## 4. π-阳离子 / 阳离子-π（Pi-Cation / Cation-Pi）

| 参数 | 值 | 说明 |
|:-----|:---|:-----|
| `PICATION_DIST_MAX` | 6.0 Å | 电荷中心到环心最大距离 |
| `PISTACK_OFFSET_MAX` | 2.0 Å | 最大偏移量 |

**规则**：
- 正电荷中心到芳香环环心距离 ≤ 6.0 Å
- 偏移量 ≤ 2.0 Å
- 若正电基团是配体的叔胺，附加角度判据

**参考文献**：Gallivan and Dougherty, 1999

---

## 5. 盐桥（Salt Bridge）

| 参数 | 值 | 说明 |
|:-----|:---|:-----|
| `SALTBRIDGE_DIST_MAX` | 5.5 Å | 相反电荷中心最大距离 |

**规则**：两个相反电荷中心距离 ≤ 5.5 Å。无额外角度限制。

**参考文献**：Barlow and Thornton, 1983（+1.5 Å 扩展）

---

## 6. 水桥（Water Bridge）

| 参数 | 值 | 说明 |
|:-----|:---|:-----|
| `WATER_BRIDGE_MINDIST` | 2.5 Å | 水 O 到极性原子最小距离 |
| `WATER_BRIDGE_MAXDIST` | 4.1 Å | 水 O 到极性原子最大距离 |
| `WATER_BRIDGE_OMEGA_MIN` | 71° | 受体-水O-供体H 最小角度 |
| `WATER_BRIDGE_OMEGA_MAX` | 140° | 受体-水O-供体H 最大角度 |
| `WATER_BRIDGE_THETA_MIN` | 100° | 水O-供体H-供体D 最小角度 |

**规则**：
- 水分子位于供体/受体之间
- 水 O 到供体/受体距离在 2.5~4.1 Å
- ω 角（受体-水O-供体H）在 71°~140°
- θ 角（水O-供体H-供体D）≥ 100°
- 一个水分子最多参与 2 个氢键（两个 H 各做一个供体）。若超过 2 个候选，保留水分子 H-O-H 角度最接近 110° 的两个

**参考文献**：Jiang et al., 2005

---

## 7. 卤键（Halogen Bond）

| 参数 | 值 | 说明 |
|:-----|:---|:-----|
| `HALOGEN_DIST_MAX` | 4.0 Å | 卤素到受体最大距离 |
| `HALOGEN_DON_ANGLE` | 165° | C-X···A 最优角度 |
| `HALOGEN_ACC_ANGLE` | 120° | X···A-R 最优角度 |
| `HALOGEN_ANGLE_DEV` | 30° | 角度偏差上限 |

**规则**：
- X-A 距离 ≤ 4.0 Å
- C-X···A 角度在 165° ± 30°（即 135°~195°）
- X···A-R 角度在 120° ± 30°（即 90°~150°）

**参考文献**：Auffinger et al.

---

## 8. 金属配位（Metal Complex）

| 参数 | 值 | 说明 |
|:-----|:---|:-----|
| `METAL_DIST_MAX` | 3.0 Å | 金属到配位原子最大距离 |

**规则**：
- 金属离子到蛋白侧链原子距离 ≤ 3.0 Å
- 蛋白侧链配位原子：Cys(S), His(N), Asn/Glu/Ser/Thr/Tyr(O), 主链 O
- 配体配位原子：醇/酚盐/羧酸盐/磷酸盐/硫醇盐/咪唑/吡咯/铁硫簇（特殊构型）
- 自动匹配配位几何构型（线性/三角/四面体/八面体等）
- 移除与最优几何构型不匹配的多余配位原子

**参考文献**：Harding, 2001

---

## 9. 全局参数

| 参数 | 值 | 说明 |
|:-----|:---|:-----|
| `BS_DIST` | 7.5 Å | 结合位点距离 cutoff |
| `MIN_DIST` | 0.5 Å | 所有距离的最小值 |
| `AROMATIC_PLANARITY` | 5.0° | 芳香环平面性偏差上限 |

---

## 10. 检测顺序（有依赖关系）

```
1. 盐桥          （charged groups → salt bridge）
2. 氢键          （移除已形成盐桥的）
3. π-π 堆积      （aromatic rings → pi stacking）
4. π-阳离子      （aromatic ring + charged_positive）
5. 卤键          （halogen_donor + halogen_acceptor）
6. 水桥          （water + H_donor/H_acceptor）
7. 金属配位      （metal + 配位原子）
8. 疏水          （移除已形成 π-π 堆积的，再聚类去重）
```

**顺序重要性**：后面的相互作用需要移除与前面重叠的部分（盐桥→氢键去重，π-π→疏水去重）。

---

## 11. PLIP 与我们的对应关系

| PLIP 相互作用 | 我们的基团 | 判据要点 |
|:-------------|:----------|:---------|
| Hydrophobic | hydrophobic | 距离 ≤ 4.0 Å |
| Hydrogen Bond | H_donor + H_acceptor | D-A ≤ 4.1 Å, DHA ≥ 100° |
| Pi Stacking | aromatic_ring | 环心 ≤ 5.5 Å, 角度偏差 ≤ 30°, 偏移 ≤ 2.0 Å |
| Pi-Cation | aromatic_ring + charged_positive | 距离 ≤ 6.0 Å, 偏移 ≤ 2.0 Å |
| Salt Bridge | charged_positive + charged_negative | 距离 ≤ 5.5 Å |
| Water Bridge | water + H_donor/H_acceptor | 距离 2.5~4.1 Å, 角度约束 |
| Halogen Bond | halogen_donor + halogen_acceptor | 距离 ≤ 4.0 Å, 双角度约束 |
| Metal Complex | metal | 距离 ≤ 3.0 Å, 配位几何匹配 |

---

*文档结束*
