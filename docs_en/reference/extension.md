# Extension Guide

This document explains how to add new force fields, new interaction types, or new criteria to DuIvyInteractions. It is intended for integrators/developers.

## Architecture Overview

The tool adopts the **strategy pattern**: four layers — Reader → GroupIdentifier → InteractionDetector → io (serialization/export) — each layer is pluggable.

```
Input (tpr+xtc) → Reader → SystemData → GroupIdentifier → List[Group]
               → Detector (per_tuple/per_frame/two_pass) → List[Interaction]
               → io/h5 serialization → io/exporter export xvg/xpm/csv
```

The three strategy detectors share a consistent interface (`detect()` returns a unified matrix-style `List[Interaction]`), switched by the `STRATEGY_INDEX` of `Pipeline`, see `core/interfaces.py`.

## Adding a New Force Field

**Principle**: group determination is based on feature-space mapping (type → {hybridization, aromaticity, polarity, has H, lone pair}) rather than hardcoded type names. A new force field only needs to fill in the feature table, without rewriting logic.

### Steps

1. Add the types of the new force field to the feature sets at the top of `group_identifiers/amber_ff_identifier.py`:
   - `STRONG_AROMATIC`: aromatic types
   - `COMPATIBLE_TYPES`: compatible types (participate in conjugation when forced by n-1 aromatic atoms)
   - `ACCEPTOR_TYPES`: acceptor types (having a lone pair)
   - `METAL_IONS`: metal elements (if a new metal is added)
   - `WATER_RESIDUES`: water residue names (if different)
2. Register the new identifier class in `IDENTIFIER_CLASSES` in `group_identifiers/__init__.py`
3. The `--ff` argument of `dii run` can then select the new force field

**Verified compatible**: amber03/94/96/99/99sb/99sb-ildn/GS/14sb + GAFF/GAFF2. Cross-version differences in type names (e.g. CX vs CT) have been handled.

## Adding a New Interaction Type

### Steps

1. **Implement the detector**: inherit from `InteractionDetectorTwoPass` (recommended, best performance) or `PerFrame`/`PerTuple`, and implement the required abstract methods:
   - `name`: the type name (e.g. `"pi_cation"`)
   - `required_group_types`: the required group types (Pipeline filters by these)
   - `metric_names`: list of metric names
   - `initialize_candidates` / `compute_pair_metrics` / `apply_threshold` (TwoPass)
2. **Register with Pipeline**: add a `(TwoPass, PerFrame, PerTuple)` triple to `DETECTOR_CLASSES` in `pipeline.py`
3. **Implement the exporter**: inherit from `InteractionExporter` in `io/interaction_exporter.py`, implement `name`/`metric_labels`/`get_pair_label`, and register in `io/__init__.py`
4. **Register with DII**: add to `ALL_INTERACTIONS` and `exporter_classes` in `DII.py`
5. **Add tests**: create unit tests in `Tests/unittests/test_<type>_*.py`

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
        # Before Pass1: generate candidate pairs (pre-filtering can be done here)
        return super().initialize_candidates(groups, trajectory, tuple_filter)

    def compute_pair_metrics(self, group_tuples, all_positions):
        # Compute metrics over all frames for the candidate pairs -> {"distance": (n_pairs, n_frames)}
        ...

    def apply_threshold(self, metrics):
        # Judge per-frame existence by criterion -> (n_pairs, n_frames) bool
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
        # Generate the pair label
        g1, g2 = interaction.groups[pair_idx]
        return f"{g1.residue_name}{g1.residue_id}···{g2.residue_name}{g2.residue_id}"
```

## Adding a New Criterion (Multiple Criteria for the Same Type)

A single type can have multiple Detectors (e.g. `HBondStrict`, `HBondLoose`); just inherit the base class and override `apply_threshold`. Criterion thresholds are module-level constants, so tuning parameters only requires changing the constants.

## Extending the Data Format

### Adding a New Metric

- Add a new key to the detector's `metrics` dictionary
- Add the corresponding label to the exporter's `metric_labels` (numeric metrics are automatically exported to xvg; string metrics such as `pistacking_type` are automatically skipped in xvg, but can be extended with a dedicated xpm)

### Group Metadata

`Group.metadata` is a key-value dictionary whose keys must be strings and whose values are limited to JSON-serializable types (numpy scalars/arrays are automatically handled as a fallback). Identifiers can use it to carry additional chemical information.

## Testing and Validation

- Run unit tests: `pytest Tests/unittests/`
- Validate with real data: run `dii run` + `dii export` end-to-end on `Tests/test_MD_case/` (KRAS-RBD system, 1 ns trajectory)
- Validate results: lossless h5 round-trip (`test_io_h5.py`), XPM value semantics consistent (`test_exporters.py`)