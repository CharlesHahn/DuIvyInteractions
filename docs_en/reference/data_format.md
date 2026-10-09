# HDF5 Storage Format

`dii run` saves the detection results as an HDF5 file, and `dii export` reads this file and exports it to visualizations/tables. This document describes the h5 storage format field by field, following the read/write logic of `DuIvyInteractions/io/h5.py`, to help understand the data organization, reuse across tools, or debugging.

## Top-Level Structure

The top level of an h5 file is an HDF5 container holding attributes and multiple interaction subgroups. `save_interactions` writes the top-level attributes `format_version` (currently `"1.0"`, see `FORMAT_VERSION` in `h5.py`) and `n_interactions`, then writes each `Interaction` as an `interaction_<i>` subgroup; on reading, `n_interactions` serves as the loop bound for reading them back.

```
<file.h5>/
├── attrs
│   ├── format_version = "1.0"      # format version (validated on reading; rejected if not matching)
│   └── n_interactions              # number of Interaction objects (read loop bound)
├── interaction_0/                  # first Interaction
├── interaction_1/                  # second (if any)
└── ...
```

## Interaction Subgroup Structure

Each `interaction_<i>` subgroup is written by `_write_interaction`:

```
interaction_<i>/
├── attrs
│   ├── interaction_type    # "hydrogen_bond", "pi_stacking", ... (string)
│   ├── n_pairs             # number of pairs
│   └── n_frames            # number of frames
├── existence/              # existence matrix, (n_pairs, n_frames) bool (gzip)
├── times/                  # time per frame, (n_frames,) float (from trajectory ts.time, unit ps)
├── metrics/                # metric dict, one dataset per metric
│   ├── distance/           # numeric metric, (n_pairs, n_frames) float, stored as-is
│   └── pistacking_type/    # string metric (π stacking type P/T/N): flattened to a 1D string list,
│                           #   the original shape is stored in this dataset's attrs['shape'],
│                           #   and restored by reshape on reading
└── groups/                 # pair data (group metadata + flattened atom data)
    ├── group_id/           # (n_groups,) global group ID (int64)
    ├── pair_index/         # (n_groups,) index of the pair each group belongs to (int64)
    ├── group_index_in_pair # (n_groups,) position of the group within the pair (0,1,...)
    ├── group_type/         # group type ("aromatic_ring", "H_donor", ..., UTF-8 string)
    ├── molecule/           # name of the molecule it belongs to (UTF-8 string)
    ├── residue_name/       # residue name (UTF-8 string)
    ├── residue_id/         # residue ID (int64)
    ├── metadata_json/      # JSON string of group metadata (ensure_ascii=False + numpy fallback)
    └── atoms/              # atom data (flattened storage, one record per atom)
        ├── pair_index          # (n_atoms,) the pair this atom belongs to
        ├── group_index_in_pair # (n_atoms,) the position of the group this atom belongs to within the pair
        ├── atom_global_idx     # global atom index (int64)
        ├── atom_idx_in_residue # in-residue atom index (int64)
        ├── atom_name           # atom name (UTF-8 string)
        ├── atom_type           # force field type, e.g., "ca" (UTF-8 string)
        ├── atom_element        # element symbol (UTF-8 string)
        ├── atom_charge         # charge (float64)
        └── atom_mass           # atomic mass (float64)
```

> Note: the `(gzip)` marker indicates that the dataset is written with gzip compression when `compress=True` (the default). **Attributes are never compressed** — the `compression` parameter applies only to datasets created by `create_dataset`.

## Key Design

### Interaction Is Stored as Matrices

An Interaction object = **all results of one type** (all pairs × all frames), rather than record-by-record entries:

- `existence[i][j]`: whether the i-th pair exists in frame j
- `metrics["distance"][i][j]`: the distance of the i-th pair in frame j
- `groups[i]`: the i-th pair (the row `existence[i]` corresponds to it)

### Flattened Atoms/Groups + Index Reconstruction

All atoms under `groups/atoms` are **flattened** into one long array; group metadata is likewise flattened (one record per group). On reading, groups are rebuilt using `(pair_index, group_index_in_pair)` as the key (`_read_groups`): for each pair, groups are taken consecutively with `group_index_in_pair = 0, 1, ...` to form a `Tuple[Group, ...]`, and each group's atoms are retrieved from the flattened atom data under the same key, finally reassembled into `List[Tuple[Group, ...]]`.

### Two Encodings of Metrics

`_write_interaction` distinguishes two encodings by dtype:

- **Numeric metrics** (e.g., distance, angle): stored as-is via `create_dataset(name, data=values)`, with shape `(n_pairs, n_frames)`.
- **String metrics** (dtype kind `U`/`S`/`O`, e.g., the `pistacking_type` of π stacking): first `flatten()`ed to a 1D list of strings and stored with the UTF-8 string dtype (`h5py.string_dtype(encoding='utf-8')`); **the original shape is stored separately in the dataset's `attrs['shape']`**. On reading, the strings are first decoded into a Python str list (`_decode_strings`) and then `reshape`d to the stored shape, restoring a NumPy string array of the same shape (e.g., `pistacking_type` as a single-character `'P'`/`'T'`/`'N'` array).

### Group Metadata Serialization

`Group.metadata` is stored as the `metadata_json` string via `json.dumps(g.metadata, ensure_ascii=False, default=_json_default)`. `_json_default` converts numpy types as a fallback: ndarray → list, np.integer → int, np.floating → float, np.bool_ → bool, anything else → `str`. On reading, it is restored with `json.loads`.

### Compression

gzip compression is enabled by default for datasets (`compress=True`), and bool/float arrays compress well. For small-file scenarios, pass `compress=False` to disable compression in exchange for speed; the compression parameter affects writing only, not reading, and does not change the field structure of the h5.

## Losslessness

The `save_interactions` → `load_interactions` round trip guarantees losslessness: existence, times, all metrics (including string metrics and their `shape` restoration), groups and their atoms, and metadata (including the numpy scalar/array fallback) can all be reproduced exactly. The project unit test `test_io_h5.py` covers this guarantee.

## Multiple Interactions

`dii run` stores one separate h5 file per type (`<output>/<type>.h5`), so a single h5 usually contains only one Interaction. However, `save_interactions` accepts a `List[Interaction]` and writes multiple `interaction_<i>` subgroups; `dii export` also supports iterating over multiple Interactions (duplicates of the same type automatically get a sequence number in the output file name, e.g., `salt_bridge_2_*`), and skips empty Interactions with 0 pairs or 0 frames.

## Interacting with Python

```python
from DuIvyInteractions.io import load_interactions

its = load_interactions("out/salt_bridge.h5")  # List[Interaction]
it = its[0]
print(it.interaction_type, it.n_pairs, it.n_frames)
print(it.existence.shape, it.metrics.keys(), it.times)
```

## Compatibility

- The top-level `format_version` is validated on reading; a mismatch with `FORMAT_VERSION = "1.0"` raises `ValueError` (older/newer files of a different version are rejected).
- The `path` parameter accepts `str` or `pathlib.Path`; an empty path raises `ValueError`, and illegal types such as BytesIO/None raise `TypeError` (to avoid silently producing junk files).
- String reading is encoding-tolerant: `bytes` elements read back by h5py are uniformly `decode('utf-8')`d to str, ensuring cross-environment readability.