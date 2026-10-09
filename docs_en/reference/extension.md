# Extension Guide

This document explains how to add new force fields, new interaction types, or new criteria to DuIvyInteractions. It is intended for integrators/developers. All interfaces follow the registries in `DuIvyInteractions/core/interfaces.py` and `pipeline.py`.

## Architecture Overview

The tool adopts the **strategy pattern**: four layers — Reader → GroupIdentifier → InteractionDetector → io (serialization/export) — each layer is pluggable.

```
Input (tpr+xtc) → Reader → SystemData → GroupIdentifier → List[Group]
               → Detector (per_tuple/per_frame/two_pass) → List[Interaction]
               → io/h5 serialization → io/exporter export xvg/xpm/csv
```

The three strategy detectors share a consistent interface (`detect()` returns a unified matrix-style `List[Interaction]`), switched by the `STRATEGY_INDEX` of `Pipeline` (`{"two_pass": 0, "per_frame": 1, "per_tuple": 2}`); the detector registry is `pipeline.DETECTOR_CLASSES`.

## Adding a New Force Field

**Principle**: group determination is based on a feature-space mapping (type → {hybridization, aromaticity, polarity, has H, lone pair}) plus structural criteria, rather than hardcoded type names. The existing `GromosFFGroupIdentifier` / `CharmmFFGroupIdentifier` / `OplsFFGroupIdentifier` all **subclass `AmberFFGroupIdentifier`** to reuse force-field-independent logic (donors, halogen bonds, metals, water, metal binding, charged functional groups) and replace only the force-field-specific parts.

### Steps

1. **Subclass the base class**: the new identifier inherits from `GroupIdentifier` (full rewrite) or `AmberFFGroupIdentifier` (recommended, reuses the class-A logic) and implements the `name` property plus any necessary overrides.
2. **Define the type feature tables at the top of the new identifier module** (module-level `frozenset`/`dict` constants; **do not** modify the tables inside `amber_ff_identifier.py`), following the existing implementations:
   - `XXX_ACCEPTOR_TYPES`: H-bond acceptor types (with a lone pair; each force field's table has already been judged type by type against chemical facts, excluding ammonium/amide/H-bearing pyrrole/guanidinium N)
   - `XXX_STRONG_AROMATIC` / `XXX_COMPATIBLE_TYPES`: aromatic types / compatible types that participate in conjugation when forced by n-1 aromatic atoms
   - aromatic-ring strong signals/residue whitelists (e.g., GROMOS `GROMOS_RING_TYPES` / `GROMOS_AROMATIC_RESIDUES` tiered determination)
   - hydrophobic types (e.g., GROMOS `GROMOS_HYDROPHOBIC_TYPES` whitelist + `GROMOS_HYDROPHOBIC_EXCLUDED` polar-neighbor exclusion)
   - `METAL_IONS` (if a new metal is added) and protein positive/negative residue dictionaries (residue name → atom name list)
3. **Declare the `WATER_RESIDUES` class attribute**: it must be overridden with a non-empty set (the base class defaults to an empty `frozenset()`); the pipeline and the internal `_find_water`/`_find_metal_binding` methods all use `self.WATER_RESIDUES` for water exclusion. Reference implementations: amber=`{SOL,HOH,WAT}`, gromos=`{SOL}`, charmm=`{TIP3,HOH,SOL,WAT}`, opls=`{HOH,HO4,HO5,SOL,WAT}`.
4. **Override structural criteria when types are coarse**: if type names cannot distinguish chemical environments (e.g., GROMOS `N` covers backbone amide/Pro/aromatic N), override `_find_acceptors` with structural criteria — the GROMOS criterion is: ≥4 bonded neighbors (ammonium) or an H-bearing neighbor (ordinary amide/H-bearing pyrrole/guanidinium) → non-acceptor; no H (Pro N / His pyridine-type N) → acceptor (design basis: `doc/force_field_compatibility_survey.md` and `doc/gromos_identifier_design.md`). CHARMM/OPLS acceptor tables have been judged type by type (no ambiguous N types), so only the table needs to be replaced, without structural criteria. Donor (`_classify_dh_pair`: D=N/O/S/F and q(H)>0) and hydrophobic determinations depend on explicit H: **united-atom force fields** (no explicit aliphatic H) must be handled by type whitelist/blacklist (as GROMOS does).
5. **Register in the registry**: add `"force_field_name": identifier_class` to `IDENTIFIER_CLASSES` in `group_identifiers/__init__.py`.
6. **CLI takes effect automatically**: the `--ff` choices of `dii run` come directly from `IDENTIFIER_CLASSES` (`DII.py`); no CLI change is needed.

### Verified Force Fields

- **Amber family**: amber03/94/96/99/99SB/99SB-ildn/GS/14SB + GAFF/GAFF2 (cross-version type-name differences handled; the mapping has been verified to have zero conflicts)
- **GROMOS**: 53A6/54A7 (united atom, explicit polar H; tiered aromatic determination + hydrophobic blacklist, see `doc/gromos_identifier_design.md`). **Boundary**: GROMOS ligands are not supported (they require ATB automatic topology parameterization, as declared); SPC/SPC-E water model (residue name `SOL`)
- **CHARMM**: CHARMM36/C36m + CGenFF ligands (verbatim consistent with the CHARMM-GUI official rtf, see `doc/charmm_identifier_design.md`); TIP3/HOH water model
- **OPLS**: OPLS-AA/L (2001, see `doc/opls_identifier_design.md`); HOH/SPC, HO4/TIP4P, HO5/TIP5P water models

Before shipping a new force field, run an end-to-end validation on real data of that force field (see the "Testing and Validation" section below).

## Adding a New Interaction Type

### Steps

1. **Implement the detector**: inherit from `InteractionDetectorTwoPass` (recommended, best performance) or `PerFrame`/`PerTuple`, and implement the required methods:
   - `name`: the type name (e.g. `"pi_cation"`)
   - `required_group_types`: the required group types (Pipeline filters by these)
   - `metric_names`: list of metric names
   - `initialize_candidates` / `compute_pair_metrics` / `apply_threshold` (TwoPass); for PerTuple/PerFrame, `get_candidate_tuples` / `compute_metrics` (`compute_metrics_for_frame` for PerFrame) / `apply_threshold`
2. **Register with Pipeline**: add a `(TwoPass, PerFrame, PerTuple)` triple to `DETECTOR_CLASSES` in `pipeline.py`, and add the type name to the `ALL_INTERACTIONS` tuple (the source of valid `dii run --interactions` values)
3. **Implement the exporter**: inherit from `InteractionExporter` in `io/interaction_exporter.py`, implement the abstract properties `name`/`metric_labels` (`get_pair_label` has a default implementation that can be overridden), and export it in `io/__init__.py`
4. **Register with DII**: add `"type_name": exporter_class` to the `exporter_classes` dict in `_run_export` of `DII.py`
5. **Add tests**: create `Tests/unittests/test_<type>_per_frame.py` and `test_<type>_two_pass.py` (per the project test convention, the per_tuple-strategy `test_<type>.py` is in a "may be deprecated" state and is not run; see below)

#### Detector Code Skeleton (TwoPass)

```python
from ..core.interfaces import InteractionDetectorTwoPass

# Criterion threshold: module-level constant (consistent with PLIP)
MY_DIST_MAX = 4.0  # Å


class MyInteractionDetectorTwoPass(InteractionDetectorTwoPass):
    """New interaction detector (TwoPass strategy)."""

    @property
    def name(self) -> str:
        return "my_interaction"

    @property
    def required_group_types(self) -> list[str]:
        return ["H_donor", "H_acceptor"]

    @property
    def metric_names(self) -> list[str]:
        return ["distance"]

    def initialize_candidates(self, groups, trajectory, tuple_filter=None):
        # Before Pass1: generate candidate pairs (pre-filtering or overriding
        # run_pass1, e.g., with a KDTree, can be done here)
        return super().initialize_candidates(groups, trajectory, tuple_filter)

    def compute_pair_metrics(self, group_tuples, all_positions):
        # Compute per-frame metrics for the candidate pairs -> {"distance": (n_groups,)}
        ...

    def apply_threshold(self, metrics):
        # Judge per-frame existence by the criterion -> (n_groups,) bool
        return metrics["distance"] <= MY_DIST_MAX
```

#### Exporter Code Skeleton

```python
from typing import Dict
from .interaction_exporter import InteractionExporter


class MyInteractionExporter(InteractionExporter):
    """New interaction exporter."""

    @property
    def name(self) -> str:
        return "My Interaction"

    @property
    def metric_labels(self) -> Dict[str, str]:
        return {"distance": "Distance (Å)"}

    def get_pair_label(self, interaction, pair_idx: int) -> str:
        # Generate the pair label (optional override; the base default is
        # "ResidueNameResidueID-ResidueNameResidueID")
        g1, g2 = interaction.groups[pair_idx]
        return f"{g1.residue_name}{g1.residue_id}···{g2.residue_name}{g2.residue_id}"
```

## Adding a New Criterion (Multiple Criteria for the Same Type)

A single type can have multiple Detectors (e.g. `HBondStrict`, `HBondLoose`); just inherit the base class and override `apply_threshold`. Criterion thresholds are module-level constants, so tuning parameters only requires changing the constants; the unit tests of the affected types must then be rerun.

## Extending the Data Format

### Adding a New Metric

- Add a new key to the detector's `metrics` dictionary (shape `(n_pairs, n_frames)`) and list it in `metric_names`
- Add the corresponding label to the exporter's `metric_labels`. Numeric metrics are automatically exported to xvg; string metrics (dtype kind `U`/`S`/`O`, e.g., `pistacking_type`) are automatically skipped in xvg and stored in the h5 as a "flattened 1D string list + attrs['shape']" (see [HDF5 Storage Format](data_format.md))
- For type-specific visualization, add a method to the exporter (following `PiStackingExporter.save_xpm_stacking_type`) and call it by type in the export flow of `DII.py`

### Group Metadata

`Group.metadata` is a key-value dictionary whose keys must be strings and whose values are limited to JSON-serializable types (numpy scalars/arrays are automatically handled as a fallback). Identifiers can use it to carry additional chemical information (e.g., the `{"source": "element"}` of `_find_metal_binding`).

## Testing and Validation

- Run unit tests: per the project test convention, run the **explicit file list** `Tests/unittests/test_<type>_per_frame.py` and `test_<type>_two_pass.py` (the per_tuple-strategy `test_<type>.py` is in a "may be deprecated" state and is not run; do not run `pytest Tests/unittests/` indiscriminately)
- Real-data validation (end-to-end, 3 force fields):
  - `Tests/test_MD_case_amber/` (Amber KRAS–RBD D927)
  - `Tests/test_MD_case_gromos/` (GROMOS 53A6 protein + 6 ligands)
  - `Tests/test_MD_case_charmm36/` (CHARMM36 SMO–BST)
- Result validation: lossless h5 round-trip (`test_io_h5.py`), XPM value semantics consistent (`test_exporters.py`)