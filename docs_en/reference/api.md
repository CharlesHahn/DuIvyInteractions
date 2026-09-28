# Python API

This document describes the Python interface of DuIvyInteractions for scripted use and secondary development. For command-line usage, see [Command Reference](../guide/command.md).

## Package Structure

```
DuIvyInteractions/
├── pipeline.py              # Pipeline main workflow
├── system_readers/          # Topology readers
├── group_identifiers/       # Group identifiers
├── interaction_detectors/   # Interaction detectors (8 types × 3 strategies)
├── io/                      # Serialization + export
└── core/                    # Data classes, interfaces, constants
```

## Top-Level Workflow: Pipeline

```python
from DuIvyInteractions.pipeline import Pipeline

# Construct: force field + strategy
pipeline = Pipeline(ff="amber", strategy="two_pass")

# Run: tpr + xtc → identify → detect → save h5
pipeline.run("md.tpr", "md.xtc", "out/", interactions=None)
# interactions=None means detecting all 8 types; a subset list can also be passed
```

Parameters:

| Parameter | Type | Description |
|:----------|:-----|:------------|
| `ff` | str | force field name, currently only `"amber"` |
| `strategy` | str | `"two_pass"` (default) / `"per_frame"` / `"per_tuple"` |
| `interactions` | list[str] \| None | subset of interaction types to detect, None = all |

## Data Reading (system_readers)

```python
from DuIvyInteractions.system_readers import GmxTprReader, GmxTprDumpReader

sd = GmxTprReader().read("md.tpr")        # read binary tpr (MDAnalysis)
# sd: SystemData (atoms/residues/bonds)
```

`GmxTprDumpReader` parses the `gmx dump` text output, serving as a cross-check/fallback for binary reading.

## Group Identification (group_identifiers)

```python
from DuIvyInteractions.group_identifiers import AmberFFGroupIdentifier

groups = AmberFFGroupIdentifier().identify(sd)   # -> List[Group]
```

`IDENTIFIER_CLASSES` is the identifier registry (the selection source of the `--ff` argument).

## Interaction Detection (interaction_detectors)

Detectors are named by type × strategy, e.g. `SaltBridgeDetectorTwoPass`, `HydrogenBondDetectorPerFrame`. Unified interface:

```python
from DuIvyInteractions.interaction_detectors import SaltBridgeDetectorTwoPass
import MDAnalysis as mda

# Filter by the group types required by the detector before detection (consistent with Pipeline internals)
det = SaltBridgeDetectorTwoPass()
filtered = [g for g in groups
            if g.group_type in det.required_group_types]

u = mda.Universe("md.tpr", "md.xtc")
results = det.detect(filtered, trajectory=u.trajectory)  # -> List[Interaction]
```

> Tip: for daily use, calling `Pipeline.run()` directly is recommended; internally it completes the entire "identify → filter by required_group_types → detect → save h5" workflow. Using detectors directly suits scenarios that require fine-grained control over detection parameters.

The three strategy base classes (`core/interfaces.py`): `InteractionDetectorPerTuple` / `PerFrame` / `TwoPass`, with a consistent `detect()` interface and results uniformly returned as `List[Interaction]`. `DETECTOR_CLASSES` in `pipeline.py` registers all 24 detectors (8 types × 3 strategies).

## Result Serialization (io.h5)

```python
from DuIvyInteractions.io import save_interactions, load_interactions

save_interactions(results, "out/salt_bridge.h5")   # List[Interaction] -> h5
its = load_interactions("out/salt_bridge.h5")      # h5 -> List[Interaction]
```

- `path` accepts `str` or `pathlib.Path`
- `compress=True` (default) enables gzip compression
- lossless round-trip (see [HDF5 Storage Format](data_format.md))

## Result Export (io.exporters)

```python
from DuIvyInteractions.io import SaltBridgeExporter

it = load_interactions("out/salt_bridge.h5")[0]
exp = SaltBridgeExporter()
exp.save_xvg_count(it, "sb_count.xvg")             # active pair count per frame
exp.save_xvg(it, "distance", "sb_distance.xvg")    # metric time series
exp.save_xpm(it, "sb_existence.xpm")               # existence heatmap
exp.to_csv_summary(it, "sb_summary.csv")           # per-pair summary
```

The 8 exporters (`HydrogenBondExporter` / `PiStackingExporter` / `SaltBridgeExporter` / `HydrophobicExporter` / `HalogenBondExporter` / `MetalCoordinationExporter` / `WaterBridgeExporter` / `PiCationExporter`) inherit from the `InteractionExporter` base class.

## Data Classes (core.datas)

| Class | Description |
|:------|:------------|
| `SystemData` | system data (residue list, atoms, inter-molecular bonds) |
| `ResidueData` | residue (residue name/number, molecule name, atom list) |
| `AtomData` | atom (global index, name, force field type, element, charge, mass) |
| `BondData` | bond (indices of the two bonded atoms, bond order) |
| `Group` | chemical group (group_id/group_type/molecule/residue/atoms/metadata) |
| `Interaction` | all detection results of one type (groups/existence/metrics/times) |
| `InteractionSparse` | sparse intermediate results of TwoPass Pass1 |

`Interaction` is stored in matrix form: `existence[i][j]` indicates whether the i-th pair exists in the j-th frame, and `metrics["distance"][i][j]` is the corresponding distance. See [HDF5 Storage Format](data_format.md).