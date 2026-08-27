# TODO 清单

> 记录推迟到后续阶段实现的任务，避免遗忘。

---

## 相互作用检测阶段需要实现的任务

以下任务在基团识别阶段有意推迟，需要在实现 `interaction_detectors/` 时一并完成。

### 1. 芳香环平面性检验

**推迟原因**：初始帧构象可能未弛豫，平面性可能很差；热涨落导致瞬时平面性波动；力场原子类型已编码芳香性，不需要几何验证。

**实现位置**：`interaction_detectors/pi_stacking.py`

**判据**：
- 对环内每个原子，计算相邻两键的法向量
- 所有法向量两两夹角 ≤ 5.0°（PLIP `AROMATIC_PLANARITY`）
- 不满足 → 该帧此环不参与 π-π 堆积

**参考**：`docs/aromatic_ring_definition.md` §3.3

### 2. 疏水-芳香去重

**推迟原因**：π-π 堆积已包含疏水接触（芳香环之间的 van der Waals 力），如果同时报告 π-π 堆积和疏水相互作用，会重复计数同一物理接触。

**实现位置**：`interaction_detectors/` 的主流程编排

**规则**：
- 先检测 π-π 堆积
- 再检测疏水相互作用
- 移除与 π-π 堆积重叠的疏水接触（两个芳香环的原子之间的疏水接触）

**参考**：`docs/hydrophobic_definition.md` §4

---

## 待实现的模块

### 3. interaction_detectors（相互作用检测器）

待实现的检测器：

| 检测器 | 文件 | 依赖基团类型 |
|:-------|:-----|:------------|
| 氢键 | `hydrogen_bond.py` | H_donor, H_acceptor |
| π-π 堆积 | `pi_stacking.py` | aromatic_ring |
| 盐桥 | `salt_bridge.py` | charged_positive, charged_negative |
| 疏水 | `hydrophobic.py` | hydrophobic |
| 卤键 | `halogen_bond.py` | halogen_donor, halogen_acceptor |
| 金属配位 | `metal_coordination.py` | metal |
| 水桥 | `water_bridge.py` | water, H_donor, H_acceptor |
| π-阳离子 | `pi_cation.py` | aromatic_ring, charged_positive |

### 4. utils（工具函数）

| 文件 | 职责 |
|:-----|:-----|
| `geometry.py` | 距离、角度、平面法向量计算 |
| `trajectory.py` | MDAnalysis 轨迹读取封装 |

### 5. pipeline.py（主流程编排）

编排 Reader → GroupIdentifier → InteractionDetector → Output 的完整流程。

### 6. visualizers（可视化）

结果可视化：时间线图、热图等。

---

## 项目收尾任务

### 7. original_draft 清理

项目完成后，处理 `DuIvyInteractions/original_draft/` 目录（移至 `archive/` 或删除）。

---

*文档结束*
