# Changelog

## v0.1.0 (In Development)

### 2026-09-16

- **XPM bug fix**: bypass DuIvyTools `refresh_by_value_matrix` remapping, add `_build_discrete_xpm` to manually build Discrete XPM (value_matrix as color index, colors/notes aligned by index), fixing color misalignment in heatmaps when the value set is incomplete
- **dii export enhancement**: iterate over all Interactions (append sequence number for duplicate types), skip empty data (0 pairs / 0 frames), friendly error for corrupted h5, reject bool + dedupe pair_indices
- **h5 serialization hardening**: metadata numpy-type fallback serialization; path parameter type validation (reject BytesIO/None garbage files)
- **Dependencies**: add `h5py` / `scipy` / `DuIvyTools` to pyproject
- **Adversarial tests**: add XPM manual build (9), DII export (4), path validation (4)

### 2026-09-14

- **CLI tool**: add `dii run` (interaction detection + h5 saving) and `dii export` (export xvg/xpm/csv + overview)
- **Pipeline orchestration**: chain Reader → Identifier → Detector → h5 saving

### 2026-09-05 ~ 09-13

- **Result serialization**: lossless HDF5 store/load (full Interaction/Group/Atom roundtrip)
- **Exporters**: InteractionExporter base + 8 subclasses (xvg/xpm/CSV); CSV summary, π-stacking 3-value XPM
- **Result time series**: Interaction gains `times` field

### 2026-08-12 ~ 09-03

- **Interaction detection**: all 8 types × 3 strategies (PerTuple / PerFrame / TwoPass)
- **TwoPass performance**: water bridge reduced from 65h to ~5s via KDTree pre-filtering
- **Directory refactor**: `input_readers` → `system_readers`, add `io/` directory

### 2026-08-11 (Group Identification Complete)

- Group identification module complete (Amber force fields)
- End-to-end pipeline from tpr to functional groups validated (D927 system)
- Type mapping zero-conflict across full Amber family (amber03/94/96/99/99sb/99sb-ildn/GS/14sb + GAFF)