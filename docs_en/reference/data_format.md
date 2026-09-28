# HDF5 Storage Format

`dii run` saves the detection results as an HDF5 file, and `dii export` reads this file and exports it to visualizations/tables. This document describes the h5 storage format, to help understand the data organization, reuse across tools, or debugging.

## Top-Level Structure

The top level of an h5 file is an HDF5 container holding attributes and multiple interaction subgroups.

```
<file.h5>/
├── attrs
│   ├── format_version = "1.0"      # format version (rejected if not matching)
│   └── n_interactions              # number of Interaction objects
├── interaction_0/                  # first Interaction
├── interaction_1/                  # second (if any)
└── ...
```

## Interaction Subgroup Structure

Each `interaction_<i>` subgroup:

```
interaction_<i>/
├── attrs
│   ├── interaction_type    # "hydrogen_bond", "pi_stacking", ... (string)
│   ├── n_pairs             # number of pairs
│   └── n_frames            # number of frames
├── existence/              # existence matrix, (n_pairs, n_frames) bool
├── times/                  # time per frame, (n_frames,) float, unit ps
├── metrics/                # metric dict, one dataset per metric
│   ├── distance/           # numeric metric, (n_pairs, n_frames) float
│   └── pistacking_type/    # string metric (e.g., π stacking type P/T/N) flattened to 1D storage, shape stored in attrs
└── groups/                 # pair data
    ├── group_id/           # (n_groups,) global group ID
    ├── pair_index/         # (n_groups,) index of the pair each group belongs to
    ├── group_index_in_pair # (n_groups,) position of the group within the pair (0,1,...)
    ├── group_type/         # group type ("aromatic_ring", "H_donor", ...)
    ├── molecule/           # name of the molecule it belongs to
    ├── residue_name/       # residue name
    ├── residue_id/         # residue ID
    ├── metadata_json/      # JSON string of additional group information
    └── atoms/              # atom data (flattened storage)
        ├── atom_global_idx   # global atom index
        ├── atom_idx_in_residue
        ├── atom_name
        ├── atom_type         # force field type (e.g., "ca")
        ├── atom_element
        ├── atom_charge
        ├── atom_mass
        ├── pair_index        # pair this atom belongs to
        └── group_index_in_pair # position of the group this atom belongs to within the pair
```

## Key Design

### Interaction Is Stored as Matrices

An Interaction object = **all results of one type** (all pairs × all frames), rather than record-by-record entries:

- `existence[i][j]`: whether the i-th pair exists in frame j
- `metrics["distance"][i][j]`: the distance of the i-th pair in frame j
- `groups[i]`: the i-th pair (the row `existence[i]` corresponds to it)

### Flattened Atoms/Groups + Index Reconstruction

All atoms under `groups/atoms` are **flattened** into one long array, and the groups within each pair and their atoms are reconstructed via `pair_index` + `group_index_in_pair`. On reading, they are reassembled into `List[Tuple[Group, ...]]` accordingly.

### Data Types

- Numeric metrics are stored as float32/float64
- String metrics (e.g., `pistacking_type`) are **flattened to 1D lists** and stored as UTF-8 string arrays, with the original shape stored separately in the dataset's `attrs['shape']`; on reading, they are restored to `<U1` arrays via `reshape` according to the shape
- Group `metadata` is serialized to a JSON string (numpy scalars/arrays are automatically converted as fallback, e.g., ndarray → list, np.int → int)

### Compression

gzip compression is enabled by default for datasets (`compress=True`), and bool/float arrays compress well. For small-file scenarios, it can be disabled in exchange for speed.

## Losslessness

The `save_interactions` → `load_interactions` round trip guarantees losslessness: existence, times, all metrics (including strings), groups and their atoms, and metadata can all be reproduced exactly. The project unit test `test_io_h5.py` covers this guarantee.

## Multiple Interactions

`dii run` stores one separate h5 file per type (`<type>.h5`), so a single h5 usually contains only one Interaction. However, `save_interactions` supports lists, and `dii export` also supports iterating over multiple Interactions (duplicates of the same type automatically get a sequence number).

## Interacting with Python

```python
from DuIvyInteractions.io import load_interactions

its = load_interactions("out/salt_bridge.h5")  # List[Interaction]
it = its[0]
print(it.interaction_type, it.n_pairs, it.n_frames)
print(it.existence.shape, it.metrics.keys(), it.times)
```

## Compatibility

- `format_version` is validated on reading; a version mismatch is rejected with an error
- The `path` parameter accepts `str` or `pathlib.Path`; illegal types such as BytesIO/None are rejected (to avoid silently producing junk files)
