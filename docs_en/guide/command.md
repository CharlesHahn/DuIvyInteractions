# dii Command Reference

`dii` is the command-line entry point of DuIvyInteractions (registered in `[project.scripts]` in `pyproject.toml`, pointing to `DuIvyInteractions.DII:main`), providing two subcommands: `run` (detection) and `export` (export).

## Global Usage

```bash
dii [-h] {run,export} ...
```

`dii --help` prints the usage above; `dii run --help` and `dii export --help` show the full arguments of each subcommand.

## Subcommand: run — Run Interaction Detection

Detect interactions in the specified trajectory and save the results as HDF5 files (format version 1.0):

```bash
dii run -t TPR -f XTC -o OUTPUT --ff FF [--interactions LIST] [--strategy STRATEGY]
```

### Arguments

| Argument | Required | Default | Description |
|:-----|:----|:----|:-----|
| `-t, --tpr` | ✅ | — | GROMACS topology file (input to group identification) |
| `-f, --xtc` | ✅ | — | Trajectory file (any format readable by MDAnalysis, e.g. xtc) |
| `-o, --output` | ✅ | — | Output directory (created automatically) |
| `--ff` | ✅ | — | Force field: `amber` / `gromos` / `charmm` / `opls` (argparse choices generated from the `IDENTIFIER_CLASSES` registry) |
| `--interactions` | No | `all` | Interaction types to detect, comma-separated (case-insensitive), e.g. `hydrogen_bond,pi_stacking` |
| `--strategy` | No | `two_pass` | Detection strategy: `two_pass` / `per_frame` / `per_tuple` |

The four `--ff` values map to the identifiers (`DuIvyInteractions/group_identifiers/__init__.py`):

| Value | Identifier | Covered force fields |
|:-----|:-----|:-----|
| `amber` | `AmberFFGroupIdentifier` | Amber 03/94/96/99/99sb/99sb-ildn/GS/14sb protein + GAFF/GAFF2 ligands |
| `gromos` | `GromosFFGroupIdentifier` | GROMOS 53A6 / 54A7 (united-atom force field, SPC/SPC-E water) |
| `charmm` | `CharmmFFGroupIdentifier` | CHARMM36 / C36m (including CGenFF ligands) |
| `opls` | `OplsFFGroupIdentifier` | OPLS-AA / L |

An unknown `--ff` aborts and lists all available values; an unknown `--interactions` type aborts and lists all available types (`ALL_INTERACTIONS`).

### The 8 Supported Interaction Types

`hydrogen_bond` (hydrogen bond), `pi_stacking` (π-π stacking), `salt_bridge` (salt bridge), `hydrophobic` (hydrophobic interaction), `halogen_bond` (halogen bond), `metal_coordination` (metal coordination), `water_bridge` (water bridge), `pi_cation` (π-cation interaction).

- The `--interactions` value is lowercased, whitespace-stripped, and split on commas (`main` in `DII.py`).
- `all` is the default and is equivalent to "all 8 types"; **`all` cannot be mixed with other types** — mixing aborts with `'all' cannot be mixed with other types`.
- A per-type detection failure (e.g. the system contains no metals or halogens) only prints `[WARN] <type> detection failed: <reason>` without interrupting the other types.

### Internal Flow

1. The tpr is read with `GmxTprReader` and the identifier performs group identification once (frame-independent, done a single time);
2. The trajectory is loaded with `mda.Universe(tpr, xtc)`; the water-residue names are taken from the identifier's `WATER_RESIDUES` (water exclusion is correct across force fields);
3. For each type, the groups are filtered by the detector's `required_group_types` (water excluded except for the water bridge), detection runs frame by frame, and the result is saved as `<output>/<type>.h5`.

### Output

One h5 file is generated per type: `<output>/<interaction_type>.h5`. For example:

```bash
dii run -t md.tpr -f md.xtc -o out/ --ff amber
# generates out/hydrogen_bond.h5, out/pi_stacking.h5, out/salt_bridge.h5,
#           out/hydrophobic.h5, out/halogen_bond.h5, out/metal_coordination.h5,
#           out/water_bridge.h5, out/pi_cation.h5
```

The h5 format version is 1.0 (`FORMAT_VERSION` in `io/h5.py`), gzip-compressed; see {doc}`/reference/data_format`.

### Strategy Descriptions

| Strategy | Mechanism | Use Case |
|:-----|:-----|:---------|
| `two_pass` (default) | Pass1 discovers active group tuples frame by frame with exact criteria (KDTree pre-screening, sparse storage) → Pass2 fills full-frame metrics for the discovered tuples | Large systems and long trajectories; best performance (water bridge: 65 h → ~5 s after KDTree pre-screening) |
| `per_frame` | Vectorized frame-by-frame computation over all candidate tuples, with a pre-allocated dense matrix | Types with very large candidate sets (e.g. water bridge); note the dense-matrix memory usage on long trajectories |
| `per_tuple` | Iterates all frames per candidate tuple and vectorizes | Control implementation in the "possibly to be removed" state (not run by tests); extremely slow with many candidates |

## Subcommand: export — Export Results

Reads h5 files, exports them as xvg/xpm/csv, and prints an overview to the terminal:

```bash
dii export -i INPUT -o OUTPUT
```

### Arguments

| Argument | Required | Description |
|:-----|:----|:-----|
| `-i, --input` | ✅ | h5 file path |
| `-o, --output` | ✅ | Output directory (created automatically) |

### Behavior

- Supports h5 files containing multiple Interactions (`dii run` generates one file per type, which fits this scenario; `load_interactions` in `io/h5.py` returns a list of Interactions).
- If the h5 contains multiple Interactions of the same type, a sequence number is appended to the filenames to avoid overwrites (e.g. `salt_bridge_2_count.xvg`).
- Empty data (0 pairs or 0 frames) prints `[skip] <type>: nothing to export (pairs=0, frames=0)` and skips that Interaction.
- A missing or corrupted h5 aborts with `cannot read h5 file: <path> (file may be corrupted or missing)` (`SystemExit`).
- An h5 with no results at all aborts with `no interaction results in h5: <path>` (`SystemExit`).
- An h5 containing an unknown type aborts with `unknown interaction type: '<type>'. Available: ...` (`SystemExit`).
- Exported xvg files cover numeric metrics only (string metrics such as the π-stacking `pistacking_type` are not exported as xvg; they appear in the type.xpm instead).

### Overview Output

An overview is printed for each exported Interaction (`_print_overview` in `DII.py`):

```
===== Salt Bridge overview =====
Type:     salt_bridge
Pairs:    47
Frames:   101
Time range: 0.0 ~ 1000.0 ps

Top 5 occupancies:
  1. ARG210(3443-3451)···ASP211(3461-3463)  100.0%
  ...
```

`Pairs` is the number of group pairs of this type (`Interaction.n_pairs`), `Frames` is the number of frames, `Time range` is the first-to-last frame time (ps, one decimal place), and `Top 5` is sorted in descending occupancy (`OVERVIEW_TOP_N = 5`).

### Output Files

For each Interaction, the following files are generated in the output directory:

| File | Content |
|:-----|:-----|
| `<type>_count.xvg` | Number of active interactions per frame (X axis: time in ps) |
| `<type>_<metric>.xvg` | Time series of each numeric metric, one curve per column (corresponding to one group pair) |
| `<type>_existence.xpm` | Existence heatmap (rows = pairs, columns = frames; white = No, blue = Yes) |
| `<type>_type.xpm` | π-stacking only: stacking-type map (0=none / 1=T-shaped / 2=parallel) |
| `<type>_summary.csv` | Per-pair summary statistics (`pair_label, occupancy, avg_<metric>, std_<metric>`) |

Numeric metrics per type (they determine the `<type>_<metric>.xvg` filenames; data from the `metric_labels` of each exporter):

| Type | Metrics (corresponding xvg filenames) |
|:-----|:-----|
| `hydrogen_bond` | `distance` (D-A distance), `angle` (D-H···A angle) |
| `salt_bridge` | `distance` (charge-center distance) |
| `pi_stacking` | `distance` (ring-center distance), `angle` (normal-vector angle), `offset` |
| `hydrophobic` | `distance` |
| `halogen_bond` | `distance` (X···A distance), `don_angle` (C-X···A angle), `acc_angle` (X···A-R angle) |
| `metal_coordination` | `distance` (metal-ligand distance) |
| `water_bridge` | `dist_dw` (D-Ow), `dist_wa` (Ow-A), `theta` (O-H···Ow angle), `omega` (H-Ow···A angle) |
| `pi_cation` | `distance` (ring-charge distance), `offset` |

Filenames are based on the interaction type rather than the Interaction index in the h5 — the raw h5 filenames produced by `dii run` (`<type>.h5`) correspond one-to-one to the export filenames (`<type>_*`).

## Exit Codes

- `0`: success
- Non-zero: argument errors, missing/corrupted files, unknown types, empty results, etc., all raised via `SystemExit` (e.g. `unknown interaction type: 'xxx'. Available: ...`, `'all' cannot be mixed with other types`, `cannot read h5 file: ...`).

## Next Steps

- Format and meaning of each output file, see [Interpreting Results](result).
- A complete runnable example on real data, see [Quick Start](quickstart).