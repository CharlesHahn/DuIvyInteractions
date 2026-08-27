# 芳香环识别定义

> 状态：定义定稿

---

## 1. 化学定义

芳香性 = 环状 + 平面 + 共轭 + Hückel 规则（4n+2 π 电子）。

Hückel 规则不限制环大小。5 元、6 元、7 元环均可芳香。

---

## 2. 力场原子类型分类

### 2.1 芳香类型（STRONG_AROMATIC）

明确标记为芳香的原子类型。

| 来源 | 类型 | 说明 |
|------|------|------|
| GAFF | ca, cg, ch, cm, cn, cp, cq, c1 | 芳香碳 |
| GAFF | na, nb, nh, ni, nj, n1, n2 | 芳香氮 |
| GAFF | pb | 芳香磷 |
| Amber 蛋白 | CA, CB, CC, CK, CM, C5, C6, C7, C*, CW, CR, CN, CV, CQ | 芳香碳 |
| Amber 蛋白 | NA, NB, NC, N* | 芳香氮 |

### 2.2 兼容原子类型（COMPATIBLE）

非芳香类型，但在 n-1 个芳香原子的"强制"下可参与共轭。

| 类型 | 来源 | 化学含义 | 兼容原因 |
|------|------|---------|---------|
| C | Amber14sb | TYR CZ（芳香）或主链 C=O（非芳香） | 歧义类型，TYR CZ 实为芳香 |
| N | Amber14sb | 杂环氮（芳香）或主链酰胺氮（非芳香） | 歧义类型 |
| os | GAFF | 呋喃 O 或醚氧 | GAFF 不区分呋喃氧和醚氧 |
| ss | GAFF | 噻吩 S 或硫醚硫 | GAFF 不区分噻吩硫和硫醚硫 |
| cc | GAFF | 非纯芳香共轭环碳 | sp2 碳，可参与共轭 |
| cd | GAFF | 同 cc | 同上 |
| pc | GAFF | 共轭环内 sp2 磷 | sp2 磷，可参与共轭 |
| pd | GAFF | 同 pc | 同上 |

### 2.3 其他类型

不在 STRONG_AROMIC 且不在 COMPATIBLE 中的原子类型，不参与芳香判定。

---

## 3. 芳香性判定规则

对每个检测到的环（任意大小），同时满足以下三个条件即为芳香环：

### 3.1 条件 1：原子类型

环内至少 n-1 个原子的类型在 STRONG_AROMIC 中（n 为环大小）。

```
aromatic_count = 环内在 STRONG_AROMATIC 中的原子数
条件 1：aromatic_count >= n-1
```

### 3.2 条件 2：兼容原子

环内不在 STRONG_AROMIC 中的原子，必须全部在 COMPATIBLE 中。

```
non_aromatic = 环内不在 STRONG_AROMIC 中的原子
条件 2：所有 non_aromatic 原子的类型都在 COMPATIBLE 中
```

### 3.3 条件 3：平面性

环必须是平面的（所有环原子近似共面）。

**检测方法**：

```
对环内每个原子 a：
  找到 a 在环内的两个邻居 n1, n2
  计算两个向量：v1 = a→n1, v2 = a→n2
  计算法向量：normal = cross(v1, v2)

对所有法向量两两配对，计算夹角：
  if 任何一对夹角满足 5.0° < angle < 175.0° → 环不平面
  else → 环是平面的
```

**阈值**：5.0°（来源：PLIP `AROMATIC_PLANARITY`，经验值）

**物理意义**：
- 完美平面环：所有法向量平行，夹夹角 = 0° 或 180°
- 允许微小偏差：±5°（MD 热振动导致的偏离）
- 明显非平面：夹角在 5°~175° 之间，p 轨道无法有效重叠

### 3.4 完整判定

```
对每个检测到的环（任意大小）：
  n = 环大小
  aromatic_count = 环内在 STRONG_AROMATIC 中的原子数
  non_aromatic = 环内不在 STRONG_AROMIC 中的原子

  if aromatic_count >= n-1
     AND 所有 non_aromatic 原子类型 ∈ COMPATIBLE
     AND 环是平面的
     → 芳香环
  else
     → 非芳香环
```

---

## 4. 芳香环去重规则

环检测可能产生冗余环（如萘的外环由两个内环组成）。去重规则：

1. 将检测到的芳香环按原子数从小到大排序
2. 维护已接受原子集合（初始为空）
3. 对每个环，检查其所有原子是否已在已接受集合中
   - 若是 → 排除（该环由更小的芳香环组成）
   - 若否 → 保留，将其原子加入已接受集合

---

## 5. 验证

### 5.1 标准芳香环

| 分子 | 环原子类型 | STRONG | COMPATIBLE | 其他 | n-1 | 平面 | 判定 |
|------|-----------|--------|------------|------|-----|------|------|
| 苯 | 6×ca | 6 | 0 | 0 | 5 | ✅ | ✅ 芳香 |
| 吡啶 | 5×ca + 1×nb | 6 | 0 | 0 | 5 | ✅ | ✅ 芳香 |
| 吡咯 | 4×ca + 1×na | 5 | 0 | 0 | 4 | ✅ | ✅ 芳香 |
| PHE | 6×CA | 6 | 0 | 0 | 5 | ✅ | ✅ 芳香 |
| TYR | 5×CA + 1×C | 5 | 1 | 0 | 5 | ✅ | ✅ 芳香 |

### 5.2 GAFF 兼容原子

| 分子 | 环原子类型 | STRONG | COMPATIBLE | 其他 | n-1 | 平面 | 判定 |
|------|-----------|--------|------------|------|-----|------|------|
| 呋喃 | 4×ca + 1×os | 4 | 1 | 0 | 4 | ✅ | ✅ 芳香 |
| 噻吩 | 4×ca + 1×ss | 4 | 1 | 0 | 4 | ✅ | ✅ 芳香 |

### 5.3 非芳香环

| 环 | 环原子类型 | STRONG | COMPATIBLE | 其他 | n-1 | 判定 |
|----|-----------|--------|------------|------|-----|------|
| D927 含硫环 | 2×ca + 2×c2 + 1×ss | 2 | 1 | 2 | 4 | ❌ 有"其他" |
| 环己烷 | 6×c3 | 0 | 0 | 6 | 5 | ❌ 无芳香原子 |

### 5.4 稠合环去重

| 分子 | 检测到的环 | 去重后 |
|------|-----------|--------|
| 萘 | 环A(6) + 环B(6) + 外环(10) | 环A + 环B |

---

## 6. 参考文献

1. Hückel E. Quantentheoretische Beiträge zum Benzolproblem. Z Phys. 1931;70(3-4):204-286.
2. PLIP 源码: https://github.com/pharmai/plip — `basic/supplemental.py` (`ring_is_planar`)
3. PLIP 配置: `AROMATIC_PLANARITY = 5.0`
4. GAFF 原子类型: http://ambermd.org/antechamber/gaff.html
5. Wang J, et al. Development and testing of a general amber force field. J Comput Chem. 2004;25(9):1157-1174.

---

*文档结束*
