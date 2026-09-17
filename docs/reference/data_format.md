# HDF5 序列化格式

`dii run` 将检测结果保存为 HDF5 文件，`dii export` 读取该文件导出为可视化/表格。本文档说明 h5 的存储格式，帮助理解数据组织、跨工具复用或调试。

## 顶层结构

h5 文件顶层是一个 HDF5 容器，含属性与多个 interaction 子组。

```
<file.h5>/
├── attrs
│   ├── format_version = "1.0"      # 格式版本（不匹配则拒读）
│   └── n_interactions              # Interaction 对象数量
├── interaction_0/                  # 第一个 Interaction
├── interaction_1/                  # 第二个（若有）
└── ...
```

## Interaction 子组结构

每个 `interaction_<i>` 子组：

```
interaction_<i>/
├── attrs
│   ├── interaction_type    # "hydrogen_bond", "pi_stacking", ...（字符串）
│   ├── n_pairs             # 基团对数量
│   └── n_frames            # 帧数
├── existence/              # 存在性矩阵，(n_pairs, n_frames) bool
├── times/                  # 每帧时间，(n_frames,) float，单位 ps
├── metrics/                # 指标字典，每个指标一个数据集
│   ├── distance/           # 数值指标，(n_pairs, n_frames) float
│   └── pistacking_type/    # 字符串指标（如 π 堆积类型 P/T/N）保存为字符串数组
└── groups/                 # 基团对数据
    ├── group_id/           # (n_groups,) 全局基团 ID
    ├── pair_index/         # (n_groups,) 每个基团所属的 pair 索引
    ├── group_index_in_pair # (n_groups,) 基团在 pair 内的位置（0,1,...）
    ├── group_type/         # 基团类型（"aromatic_ring", "H_donor", ...）
    ├── molecule/           # 所属分子名
    ├── residue_name/       # 残基名
    ├── residue_id/         # 残基号
    ├── metadata_json/      # 基团附加信息的 JSON 字符串
    └── atoms/              # 原子数据（平铺存储）
        ├── atom_global_idx   # 全局原子索引
        ├── atom_idx_in_residue
        ├── atom_name
        ├── atom_type         # 力场类型（如 "ca"）
        ├── atom_element
        ├── atom_charge
        ├── atom_mass
        ├── pair_index        # 该原子所属 pair
        └── group_index_in_pair # 该原子所属基团在 pair 内的位置
```

## 关键设计

### Interaction 是矩阵式存储

一个 Interaction 对象 = **一种类型的全部结果**（所有 pair × 全部帧），而非逐条记录：

- `existence[i][j]`：第 i 个基团对在第 j 帧是否存在
- `metrics["distance"][i][j]`：第 i 个基团对在第 j 帧的距离
- `groups[i]`：第 i 个基团对（`existence[i]` 行与之对应）

### 原子/基团平铺 + 索引重建

`groups/atoms` 下所有原子**平铺**为一个长数组，通过 `pair_index` + `group_index_in_pair` 重建出每个 pair 内的基团及其原子。读取时据此重组为 `List[Tuple[Group, ...]]`。

### 数据类型

- 数值型指标存为 float32/float64
- 字符串指标（如 `pistacking_type`）存为 UTF-8 字符串数组，读取时还原为 `<U1` 数组
- 基团 `metadata` 序列化为 JSON 字符串（自动兜底 numpy 标量/数组，如 ndarray → list、np.int → int）

### 压缩

默认对数据集启用 gzip 压缩（`compress=True`），bool/float 数组压缩率高。小文件场景可关闭以换取速度。

## 无损性

`save_interactions` → `load_interactions` 往返保证无损：existence、times、全部 metrics（含字符串）、groups 及其原子、metadata 均可精确还原。项目单元测试 `test_io_h5.py` 覆盖此保证。

## 多 Interaction

`dii run` 每个类型存一个独立 h5 文件（`<类型>.h5`），因此单个 h5 通常只含一个 Interaction。但 `save_interactions` 支持列表，`dii export` 也支持遍历多个 Interaction（同类型重复自动加序号）。

## 与 Python 交互

```python
from DuIvyInteractions.io import load_interactions

its = load_interactions("out/salt_bridge.h5")  # List[Interaction]
it = its[0]
print(it.interaction_type, it.n_pairs, it.n_frames)
print(it.existence.shape, it.metrics.keys(), it.times)
```

## 兼容性

- 读取时校验 `format_version`，版本不匹配会拒绝并报错
- `path` 参数接受 `str` 或 `pathlib.Path`；传入 BytesIO/None 等非法类型会被拒绝（避免静默生成垃圾文件）