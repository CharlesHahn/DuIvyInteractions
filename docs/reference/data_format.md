# HDF5 序列化格式

`dii run` 将检测结果保存为 HDF5 文件，`dii export` 读取该文件导出为可视化/表格。本文档按 `DuIvyInteractions/io/h5.py` 的读写逻辑逐字段说明 h5 的存储格式，帮助理解数据组织、跨工具复用或调试。

## 顶层结构

h5 文件顶层是一个 HDF5 容器，含属性（attrs）与多个 interaction 子组。`save_interactions` 写入顶层属性 `format_version`（当前 `"1.0"`，见 `h5.py` 的 `FORMAT_VERSION`）与 `n_interactions`，随后每个 `Interaction` 写为一个 `interaction_<i>` 子组；读取时以 `n_interactions` 为循环上界逐个读回。

```
<file.h5>/
├── attrs
│   ├── format_version = "1.0"      # 格式版本（读取时校验，不匹配则拒读）
│   └── n_interactions              # Interaction 对象数量（读取循环上界）
├── interaction_0/                  # 第一个 Interaction
├── interaction_1/                  # 第二个（若有）
└── ...
```

## Interaction 子组结构

每个 `interaction_<i>` 子组由 `_write_interaction` 写出：

```
interaction_<i>/
├── attrs
│   ├── interaction_type    # "hydrogen_bond", "pi_stacking", ...（字符串）
│   ├── n_pairs             # 基团对数量
│   └── n_frames            # 帧数
├── existence/              # 存在性矩阵，(n_pairs, n_frames) bool（gzip）
├── times/                  # 每帧时间，(n_frames,) float（来自轨迹 ts.time，单位 ps）
├── metrics/                # 指标字典，每个指标一个数据集
│   ├── distance/           # 数值指标，(n_pairs, n_frames) float，原样存储
│   └── pistacking_type/    # 字符串指标（π 堆积类型 P/T/N）：展平为 1D 字符串列表存储，
│                           #   原始 shape 存于该数据集的 attrs['shape']，读取时 reshape 还原
└── groups/                 # 基团对数据（基团元数据 + 平铺原子数据）
    ├── group_id/           # (n_groups,) 全局基团 ID（int64）
    ├── pair_index/         # (n_groups,) 每个基团所属的 pair 索引（int64）
    ├── group_index_in_pair # (n_groups,) 基团在 pair 内的位置（0,1,...）
    ├── group_type/         # 基团类型（"aromatic_ring", "H_donor", ...，UTF-8 字符串）
    ├── molecule/           # 所属分子名（UTF-8 字符串）
    ├── residue_name/       # 残基名（UTF-8 字符串）
    ├── residue_id/         # 残基号（int64）
    ├── metadata_json/      # 基团 metadata 的 JSON 字符串（ensure_ascii=False + numpy 兜底）
    └── atoms/              # 原子数据（平铺存储，每条记录一个原子）
        ├── pair_index          # (n_atoms,) 该原子所属 pair
        ├── group_index_in_pair # (n_atoms,) 该原子所属基团在 pair 内的位置
        ├── atom_global_idx     # 全局原子索引（int64）
        ├── atom_idx_in_residue # 残基内索引（int64）
        ├── atom_name           # 原子名（UTF-8 字符串）
        ├── atom_type           # 力场类型，如 "ca"（UTF-8 字符串）
        ├── atom_element        # 元素符号（UTF-8 字符串）
        ├── atom_charge         # 电荷（float64）
        └── atom_mass           # 原子质量（float64）
```

> 注：`(gzip)` 标记表示该数据集在 `compress=True`（默认）时以 gzip 压缩写入。**属性（attrs）不做压缩**——`compression` 参数只作用于 `create_dataset` 创建的数据集。

## 关键设计

### Interaction 是矩阵式存储

一个 Interaction 对象 = **一种类型的全部结果**（所有 pair × 全部帧），而非逐条记录：

- `existence[i][j]`：第 i 个基团对在第 j 帧是否存在
- `metrics["distance"][i][j]`：第 i 个基团对在第 j 帧的距离
- `groups[i]`：第 i 个基团对（`existence[i]` 行与之对应）

### 原子/基团平铺 + 索引重建

`groups/atoms` 下所有原子**平铺**为一个长数组；基团元数据同样平铺（每条记录一个基团）。读取时以 `(pair_index, group_index_in_pair)` 为键重建：对每个 pair，按 `group_index_in_pair = 0, 1, ...` 连续取出基团组成 `Tuple[Group, ...]`，每个基团的原子从 `atoms` 平铺数据中按同一键取回，最终重组为 `List[Tuple[Group, ...]]`（`_read_groups`）。

### 指标的两种编码

`_write_interaction` 按 dtype 区分两种编码：

- **数值指标**（如 distance、angle）：`create_dataset(name, data=values)` 原样存储，shape 为 `(n_pairs, n_frames)`。
- **字符串指标**（dtype kind 为 `U`/`S`/`O`，如 π 堆积的 `pistacking_type`）：先 `flatten()` 为 1D 字符串列表，以 UTF-8 字符串 dtype（`h5py.string_dtype(encoding='utf-8')`）存储；**原始 shape 单独存于该数据集的 `attrs['shape']`**。读取时先把字符串解码为 Python str 列表（`_decode_strings`），再按 `shape` `reshape` 还原为同 shape 的 NumPy 字符串数组（如 `pistacking_type` 为单字符 `'P'`/`'T'`/`'N'` 数组）。

### 基团 metadata 序列化

`Group.metadata` 经 `json.dumps(g.metadata, ensure_ascii=False, default=_json_default)` 存为 `metadata_json` 字符串。`_json_default` 兜底 numpy 类型：ndarray → list、np.integer → int、np.floating → float、np.bool_ → bool，其余转 `str`。读取时 `json.loads` 还原。

### 压缩

默认对数据集启用 gzip 压缩（`compress=True`），bool/float 数组压缩率高。小文件场景可传 `compress=False` 关闭以换取速度；压缩参数只影响写出，不影响读取，且不影响 h5 的字段结构。

## 无损性

`save_interactions` → `load_interactions` 往返保证无损：existence、times、全部 metrics（含字符串指标及 `shape` 还原）、groups 及其原子、metadata（含 numpy 标量/数组兜底）均可精确还原。项目单元测试 `test_io_h5.py` 覆盖此保证。

## 多 Interaction

`dii run` 每类型存一个独立 h5 文件（`<output>/<类型>.h5`），因此单个 h5 通常只含一个 Interaction。但 `save_interactions` 接受 `List[Interaction]`，会写出多个 `interaction_<i>` 子组；`dii export` 也支持遍历多个 Interaction（同类型重复时输出文件名自动加序号，如 `salt_bridge_2_*`），并跳过 0 对或 0 帧的空 Interaction。

## 与 Python 交互

```python
from DuIvyInteractions.io import load_interactions

its = load_interactions("out/salt_bridge.h5")  # List[Interaction]
it = its[0]
print(it.interaction_type, it.n_pairs, it.n_frames)
print(it.existence.shape, it.metrics.keys(), it.times)
```

## 兼容性

- 读取时校验顶层 `format_version`，与 `FORMAT_VERSION = "1.0"` 不匹配则抛 `ValueError`（拒绝不同版本的旧/新文件）。
- `path` 参数接受 `str` 或 `pathlib.Path`；空路径抛 `ValueError`，传入 BytesIO/None 等非法类型抛 `TypeError`（避免静默生成垃圾文件）。
- 字符串读取按编码兜底：h5py 读回的 `bytes` 元素统一 `decode('utf-8')` 为 str，保证跨环境可读。