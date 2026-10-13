# Python API

This document describes the Python interface of DuIvyInteractions for scripted use and secondary development. For command-line usage, see [Command Reference](../guide/command.md). All signatures follow `DuIvyInteractions/core/interfaces.py`, `core/datas.py`, `io/h5.py`, and `pipeline.py`.

## Package Structure

```
DuIvyInteractions/
├── pipeline.py              # Pipeline main workflow (Reader→Identifier→Detector→h5)
├── DII.py                   # Command-line entry (dii run / dii export)
├── selection.py             # User group selection & filtering (--select/--pairing)
├── system_readers/          # Topology readers (GmxTprReader / GmxTprDumpReader)
├── group_identifiers/       # Group identifiers (4 force fields: amber/gromos/charmm/opls; IDENTIFIER_CLASSES registry)
├── interaction_detectors/   # Interaction detectors (8 types × 3 strategies, 24 classes)
├── io/                      # h5 serialization + 8 exporters
└── core/                    # Data classes (datas), interfaces (interfaces), constants (constants)
```

## Top-Level Workflow: Pipeline

`Pipeline` orchestrates the complete "read tpr → identify groups → detect per type → write h5" workflow. The `ff` constructor argument supports 4 force fields:

```python
from DuIvyInteractions.pipeline import Pipeline

# ff: "amber" | "gromos" | "charmm" | "opls" (identical to dii run --ff)
pipeline = Pipeline(ff="amber", strategy="two_pass")

# Run: tpr + xtc → identify → detect → write <type>.h5 under output/
pipeline.run("md.tpr", "md.xtc", "out/", interactions=None, select=None, pairing="any")
# interactions=None means detecting all 8 types; a subset list can also be passed
# select=None means use the whole participation pool (no filtering);
# pairing="any" means no pairwise constraint between selects
```

| Parameter | Type | Description |
|:----------|:-----|:------------|
| `ff` | str | force field name: `"amber"` / `"gromos"` / `"charmm"` / `"opls"` (keys of `IDENTIFIER_CLASSES`; an unknown force field raises `ValueError`) |
| `strategy` | str | `"two_pass"` (default, best performance) / `"per_frame"` / `"per_tuple"` |
| `interactions` | list[str] \| None | subset of interaction types (names in `pipeline.ALL_INTERACTIONS`), None = all 8 types |
| `select` | list[str] \| None | optional MDAnalysis selection strings restricting the participation pool; each string is evaluated independently as its own atom set (never merged); a group is retained only if one select fully covers all its atoms; None = whole pool |
| `pairing` | str | `"any"` (default): all pool candidates, no pairing constraint. `"between"`: keep only candidate pairs whose two ends belong to different selects (requires ≥ 2 selects, otherwise `ValueError`) |

Behavior details of `run()`:

- Group identification is performed once and is frame-independent; the trajectory is loaded once (`mda.Universe(tpr, xtc)`).
- Water residues are excluded using the identifier's own `WATER_RESIDUES` (see below): `_filter_groups` drops groups whose residue name is in `WATER_RESIDUES` for every type except water bridges, which ensures correct cross-force-field exclusion (e.g., CHARMM TIP3, OPLS HO4/HO5). When `select` is given, a water-residue `H_donor`/`H_acceptor` group that is assigned to a concrete select index is exempted from this exclusion (explicit water participation, e.g. protein × water hydrogen bonds via `select=["segid seg_0_RBD_pro", "resname SOL"]`).
- When `select` is given, each expression is resolved once on the first frame (`u.trajectory[0]`) with MDAnalysis `select_atoms`; empty selects, pairwise-overlapping selects, and a pool left with zero non-water groups raise `ValueError`, which the CLI turns into a friendly `error: ...` message (`SystemExit`, no traceback). See the `tuple_filter` note under Interaction Detection below, and [Known Limitations](limitations.md) for the static-evaluation caveat.
- Each type is detected and saved as `<output>/<type>.h5`; a failure of one type only prints `[WARN]` and does not interrupt the remaining types.
- The full list of 8 type names is in `pipeline.ALL_INTERACTIONS`; the type→detector registry is `DETECTOR_CLASSES` (8 types × 3 strategies, 24 classes, keyed to `(TwoPass, PerFrame, PerTuple)` class triples), and strategy switching is handled by `STRATEGY_INDEX = {"two_pass": 0, "per_frame": 1, "per_tuple": 2}`.

### Selection Validation Errors (`select` / `pairing`)

All ambiguous or empty selection states are reported rather than silently degraded: `selection.py`/`pipeline.py` raise `ValueError`, and `dii run` converts it into a friendly `SystemExit` prefixed with `error: ` (no traceback). The message text is Chinese by design; the exact strings are:

| Case | Literal message (CLI output) |
|:-----|:------------------------------|
| select expression has a syntax error | `error: select 语法错误: '<expr>' → <MDAnalysis reason>` (select syntax error: `<expr>` → `<MDAnalysis reason>`) |
| a select matches 0 atoms | `error: select #N 未选中任何原子` (select #N matched no atoms) |
| overlap between two selects | `error: select #N 与 select #M 原子重叠 (K 个原子)，请用 and not 显式剔除后重试` (select #N overlaps select #M by K atoms; please exclude explicitly with `and not` and retry) |
| `between` with fewer than 2 selects | `error: --pairing between 需要至少 2 个 --select` (--pairing between requires at least 2 --select) |
| zero non-water groups left in the pool | `error: 参与池内无可用基团（非水基团 0 个）。若使用了 around/sphzone 几何选择，基团可能因部分原子在圈外被排除，请改用 byres around …` (no usable groups in the pool, 0 non-water groups; if you used an `around`/`sphzone` geometric selection, groups may have been excluded because some atoms fall outside the region — switch to `byres around …`) |

## User Group Selection (selection)

`DuIvyInteractions/selection.py` implements the core logic of `--select` / `--pairing` (design: `doc/group_selection_design.md`). It is called by `Pipeline.run(select=..., pairing=...)`, and the CLI layer (`DII.py`) passes the same parameters through.

```python
from DuIvyInteractions.selection import (
    WATER_SIDE, CROSS_MIN_SIDES, CONSTRAINT_POSITIONS,
    resolve_selection, validate_selections, assign_groups,
    filter_groups, build_tuple_filter,
)
```

Public functions:

| Function | Signature | Description |
|:---------|:----------|:------------|
| `resolve_selection` | `(universe, sel_str) -> frozenset[int]` | Resolve one selection expression to an atom-index set on the Universe (already seeked to the first frame); raises `ValueError` on a selection syntax error |
| `validate_selections` | `(sets) -> None` | Validate select sets: a select matching 0 atoms, or atomic overlap between two selects → raises `ValueError` |
| `assign_groups` | `(groups, sets) -> Assignment` | Determine each group's side: a non-water group → the select index that fully covers all of its atoms; a `group_type == "water"` group always belongs to the water side (`WATER_SIDE`); groups not fully covered by any select are omitted |
| `filter_groups` | `(groups, assignment) -> List[Group]` | Filter the group pool by side assignment: keep only groups with a side (including the water side) |
| `build_tuple_filter` | `(pairing, assignment, interaction_type=None) -> predicate \| None` | Assemble the tuple-level filter predicate: `pairing != "between"` → `None`; `between` → judge per `CONSTRAINT_POSITIONS` that the groups at the constrained positions belong to different selects |

Constants:

| Constant | Value | Meaning |
|:---------|:------|:--------|
| `WATER_SIDE` | `"water"` | The constant side marker for water (implicit select, cannot be disabled or selected by the user) |
| `CROSS_MIN_SIDES` | `2` | Minimum number of distinct sides required by `between` (cross-side = at least 2 different selects) |
| `CONSTRAINT_POSITIONS` | `{type: (positions, ...)}` | Fixed member positions participating in the pairing judgment per detector tuple (not user-configurable): the 7 pair-type detectors (`hydrogen_bond`/`pi_stacking`/`salt_bridge`/`hydrophobic`/`halogen_bond`/`metal_coordination`/`pi_cation`) → `(0, 1)`; `water_bridge` → `(0, 2)` — judging the (donor, acceptor) ends, water neutral |

Example — detect only the interactions between chains RBD_pro and KRAS_pro (across selects):

```python
from DuIvyInteractions.pipeline import Pipeline

Pipeline(ff="amber").run(
    "md.tpr", "md.xtc", "out/",
    interactions=["hydrogen_bond", "salt_bridge"],
    select=["segid seg_0_RBD_pro", "segid seg_2_KRAS_pro"],
    pairing="between",
)
```

## Data Reading (system_readers)

```python
from DuIvyInteractions.system_readers import GmxTprReader, GmxTprDumpReader

sd = GmxTprReader().read("md.tpr")        # read binary tpr (MDAnalysis)
# sd: SystemData (residues/atoms/bonds/inter-residue bonds)
```

- `GmxTprReader`: reads the binary tpr via MDAnalysis (main path).
- `GmxTprDumpReader`: parses the `gmx dump` text output, serving as a cross-check/fallback for binary reading.
- Both implement the same `Reader` interface (`core/interfaces.py`): a `name` property plus `read(source) -> SystemData`.

## Group Identification (group_identifiers)

Four force-field identifier classes (all under `DuIvyInteractions/group_identifiers/`; the latter three subclass `AmberFFGroupIdentifier` to reuse force-field-independent logic):

| Class | Force field | Water residue names (`WATER_RESIDUES` class attribute) |
|:------|:------------|:-------------|
| `AmberFFGroupIdentifier` | Amber family (amber03/94/96/99/99SB/99SB-ildn/GS/14SB) + GAFF/GAFF2 | `{"SOL", "HOH", "WAT"}` |
| `GromosFFGroupIdentifier` | GROMOS 53A6/54A7 (united atom, explicit polar H) | `{"SOL"}` |
| `CharmmFFGroupIdentifier` | CHARMM36/C36m + CGenFF | `{"TIP3", "HOH", "SOL", "WAT"}` |
| `OplsFFGroupIdentifier` | OPLS-AA/L | `{"HOH", "HO4", "HO5", "SOL", "WAT"}` |

```python
from DuIvyInteractions.group_identifiers import AmberFFGroupIdentifier

groups = AmberFFGroupIdentifier().identify(sd)   # -> List[Group]
```

- Interface: `GroupIdentifier` in `core/interfaces.py` — a `name` property plus `identify(SystemData) -> List[Group]`.
- **`WATER_RESIDUES` is a class attribute**: the base class defaults to an empty `frozenset()`, and subclasses must override it with a non-empty set. `_filter_groups` in the pipeline and the internal `_find_water`/`_find_metal_binding` methods all use `self.WATER_RESIDUES`, making it the single entry point for cross-force-field water exclusion (e.g., CHARMM TIP3, OPLS HO4/HO5).
- `IDENTIFIER_CLASSES` (`group_identifiers/__init__.py`) is the registry `{"amber", "gromos", "charmm", "opls"}`; it is the source of the `--ff` choices in `DII.py` and the selection source for `Pipeline(ff=...)`.

## Interaction Detection (interaction_detectors)

Detectors are named by type × strategy, e.g., `SaltBridgeDetectorTwoPass`, `HydrogenBondDetectorPerFrame`. Unified interface:

```python
from DuIvyInteractions.interaction_detectors import SaltBridgeDetectorTwoPass
import MDAnalysis as mda

det = SaltBridgeDetectorTwoPass()
# Filter by the group types required by the detector (consistent with Pipeline internals)
filtered = [g for g in groups
            if g.group_type in det.required_group_types]

u = mda.Universe("md.tpr", "md.xtc")
results = det.detect(filtered, trajectory=u.trajectory)  # -> List[Interaction]
```

Full `detect()` signature (identical across the three strategy base classes, see `core/interfaces.py`):

```python
detect(groups, trajectory=None, n_workers=1,
       topology_path=None, trajectory_path=None,
       tuple_filter=None) -> List[Interaction]
```

| Parameter | Description |
|:----------|:------------|
| `groups` | the group list already filtered by `required_group_types` |
| `trajectory` | an MDAnalysis trajectory object (required for serial execution; `PerFrame`/`TwoPass` ignore `n_workers` and require it) |
| `n_workers` | number of parallel workers (>1 requires `topology_path` and `trajectory_path`; parallel execution is implemented only in the `PerTuple` strategy) |
| `tuple_filter` | optional callback `(Tuple[Group, ...]) -> bool` applied to candidate group tuples for user-defined filtering (e.g., keep only pairs from different molecules); it is applied after candidate generation and before geometric evaluation |

> The same capability is exposed to end users: `Pipeline.run(..., select=..., pairing=...)` — and the CLI `--select` / `--pairing` of `dii run` — builds the `tuple_filter` internally (implementation in `selection.py`), so the detectors themselves are untouched. With `pairing="between"`, the predicate keeps only tuples whose two ends belong to *different* selects; the judged positions are fixed in `selection.CONSTRAINT_POSITIONS`: the 7 pair-type detectors (`hydrogen_bond`, `pi_stacking`, `salt_bridge`, `hydrophobic`, `halogen_bond`, `metal_coordination`, `pi_cation`) judge positions (0, 1), while `water_bridge` judges `(donor, water, acceptor)` → (0, 2) — the water molecule is neutral and never participates. `pairing="any"` attaches no filter (`None`, equivalent to the identity, zero overhead).

- The three strategy base classes — `InteractionDetectorPerTuple` / `InteractionDetectorPerFrame` / `InteractionDetectorTwoPass` (template method pattern) — share a consistent `detect()` interface and uniformly return `List[Interaction]`.
- Subclasses must implement the `name`, `required_group_types`, and `metric_names` properties plus the strategy-specific detection methods: `get_candidate_tuples`/`compute_metrics`/`apply_threshold` for PerTuple; `get_candidate_tuples`/`compute_metrics_for_frame`/`apply_threshold` for PerFrame; `initialize_candidates`/`compute_pair_metrics`/`apply_threshold` for TwoPass (the latter three have base-class default implementations returning empty values and must be overridden). All three may optionally override `filter_candidate_tuples` (first-frame pre-filtering) and `_post_process` (cross-pair post-processing hook).
- All 24 detectors (8 types × 3 strategies) are registered in `pipeline.DETECTOR_CLASSES`.
- Strategy notes: the `per_tuple` strategy is in a "may be deprecated" state and is not a current result (tests do not run it; see [Known Limitations](limitations.md)); the current strategies are `per_frame` and `two_pass`.

> Tip: for daily use, calling `Pipeline.run()` directly is recommended. Using detectors directly suits scenarios that require fine-grained control over candidate generation/filtering (e.g., `tuple_filter`).

## Result Serialization (io.h5)

```python
from DuIvyInteractions.io import save_interactions, load_interactions

save_interactions(results, "out/salt_bridge.h5")   # List[Interaction] -> h5
its = load_interactions("out/salt_bridge.h5")      # h5 -> List[Interaction]
```

- `save_interactions(interactions, path, compress=True)`: `compress` enables gzip compression (on by default).
- `path` accepts `str` or `pathlib.Path`; an empty path raises `ValueError`, a non-str/Path type raises `TypeError`.
- The top-level `format_version` is validated on reading; a mismatch with `FORMAT_VERSION = "1.0"` raises `ValueError`.
- Lossless round-trip (see [HDF5 Storage Format](data_format.md)); the project unit test `test_io_h5.py` covers this guarantee.

## Result Export (io.exporters)

```python
from DuIvyInteractions.io import SaltBridgeExporter

it = load_interactions("out/salt_bridge.h5")[0]
exp = SaltBridgeExporter()
exp.save_xvg_count(it, "sb_count.xvg")             # active pair count per frame
exp.save_xvg(it, "distance", "sb_distance.xvg")    # numeric metric time series
exp.save_xpm(it, "sb_existence.xpm")               # existence heatmap
exp.to_csv_summary(it, "sb_summary.csv")           # per-pair summary (occupancy + avg/std)
```

- The 8 exporters (`HydrogenBondExporter` / `PiStackingExporter` / `SaltBridgeExporter` / `HydrophobicExporter` / `HalogenBondExporter` / `MetalCoordinationExporter` / `WaterBridgeExporter` / `PiCationExporter`) inherit from `InteractionExporter` in `io/interaction_exporter.py` and are all exported in `io/__init__.py`.
- Subclasses must implement the two abstract properties `name` and `metric_labels`; `get_pair_label` has a default implementation (`"ResidueNameResidueID-ResidueNameResidueID"`, e.g., `ARG73-D927`) that can be overridden.
- `PiStackingExporter` additionally provides `save_xpm_stacking_type` (a P/T/N π-stacking type heatmap).
- `to_csv_summary` outputs numeric metrics only (string metrics such as `pistacking_type` are skipped); `dii export` looks up the exporter by the `interaction_type` stored in the h5 and exits with an error on unknown types.

## Data Classes (core.datas)

| Class | Description |
|:------|:------------|
| `SystemData` | system data: `system_name`, `residues` (residue list), `inter_residue_bonds` (inter-residue bonds such as peptide/disulfide bonds); validates global-index uniqueness of residues/atoms on construction |
| `ResidueData` | residue: residue name/global index/in-molecule number, molecule name, atom list, in-residue bond list |
| `AtomData` | atom: global index, in-residue index, atom name, force field type, element, charge, mass |
| `BondData` | in-residue bond: two in-residue endpoint indices + bond type (in `BOND_TYPES`) |
| `InterResidueBond` | inter-residue covalent bond (two residue global indices + in-residue atom indices + bond type) |
| `Group` | chemical group: `group_id`/`group_type`/`molecule`/`residue_name`/`residue_id`/`atoms`/`metadata`; properties `num_atoms`/`atom_indices`/`net_charge`; `group_type` must be a value of `GROUP_TYPES` and `atoms` must be non-empty |
| `Interaction` | all detection results of one type: `interaction_type`/`groups`/`existence`/`metrics`/`times`; properties `n_pairs`/`n_frames`, method `occupancy()` |
| `InteractionSparse` | sparse intermediate results of TwoPass Pass1 (keyed by `(group_id, ...)` tuples) |

- `Group.metadata`: a key-value dictionary whose keys must be strings and whose values are limited to JSON-serializable types (str/int/float/bool/None/list/dict; numpy scalars/arrays are automatically converted as a fallback).
- `Interaction` is stored in matrix form: `existence[i][j]` indicates whether the i-th pair exists in frame j, and `metrics["distance"][i][j]` is the corresponding geometric metric; `occupancy()` returns the per-pair existence ratio. See [HDF5 Storage Format](data_format.md).
