# Changelog

## v0.1.0 (In Development)

### 2026-09-17

- **Adversarial documentation revision**: Created the Chinese user documentation site (guide + reference); optimized the main index directory structure; fixed MyST internal cross-references; corrected the label formats for hydrophobic interactions/halogen bonds in Interpreting Results and the description of the h5 string metric storage mechanism

### 2026-09-16

- **XPM bug fix**: Bypassed the DuIvyTools `refresh_by_value_matrix` remapping and added `_build_discrete_xpm` to manually build the Discrete XPM (value_matrix is the color index; colors/notes are aligned by index), fixing heatmap color misalignment when the value set is incomplete
- **dii export enhancement**: Iterates over all Interactions (same-type duplicates get serial numbers), skips export for empty data/0 frames, reports corrupted h5 files gracefully, rejects bool for pair_indices plus order-preserving deduplication
- **h5 serialization hardening**: metadata supports fallback serialization for numpy types; type validation for the path parameter (rejects garbage files such as BytesIO/None)
- **Dependency completion**: Added `h5py` / `scipy` / `DuIvyTools` to pyproject
- **Adversarial testing**: Added tests for XPM manual construction (9 cases), DII export (4 cases), and path validation (4 cases)

### 2026-09-14

- **Command-line tools**: Added `dii run` (interaction detection + h5 saving) and `dii export` (export xvg/xpm/csv + overview)
- **Pipeline orchestration**: Chains Reader → Identifier → Detector → h5 saving

### 2026-09-05 ~ 09-13

- **Result serialization**: Lossless HDF5 storage/loading (full roundtrip of Interaction/Group/Atom)
- **Exporter**: InteractionExporter base class + 8 subclasses (xvg/xpm/CSV); CSV summary, three-value XPM for π stacking types
- **Result time series**: Added the times field to Interaction

### 2026-08-12 ~ 09-03

- **Interaction detection**: All 8 types × 3 strategies (PerTuple / PerFrame / TwoPass) implemented
- **TwoPass performance optimization**: After pre-filtering water bridges with KDTree, the runtime dropped from 65 h to ~5 s
- **Directory refactoring**: `input_readers` → `system_readers`; added the `io/` directory

### 2026-08-11 (Group identification completed)

- Group identification module completed (Amber force field)
- Full pipeline from tpr to functional groups validated (D927 system verification)
- The type mapping table shows zero conflicts across the entire Amber family (amber03/94/96/99/99sb/99sb-ildn/GS/14sb + GAFF)
