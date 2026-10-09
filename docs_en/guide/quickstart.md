# Quick Start

This guide runs the full `dii run → dii export` workflow with the real GROMACS test data bundled with the repository, and provides runnable commands for all 4 force fields. All example commands must be executed from the **repository root** (the directory containing `Tests/` and `pyproject.toml`); the sample outputs are aligned with the current code (default `two_pass` strategy, HDF5 v1.0 storage).

## 1. Prepare Data

The repository provides 3 sets of real GROMACS test data under `Tests/` (shared with the `Tests/unittests/` integration tests):

| Directory | Force field | System | Files |
|:-----|:-----|:-----|:-----|
| `test_MD_case_amber/` | Amber family (amber14sb protein + GAFF ligand) | KRAS-RBD D927 complex (contains the GNP/Mg²⁺ metal center), 116,383 atoms | `md.tpr` + `md1ns.xtc` (1 ns, 101 frames, 0–1000 ps) |
| `test_MD_case_gromos/` | GROMOS 53A6 | 130-residue protein + 6 ligands (ZIN1–6) | `gromos53a6_md.tpr` + `gromos53a6_md10ns.xtc` |
| `test_MD_case_charmm36/` | CHARMM36 | SMO-BST (Smoothened–β-sitosterol) complex | `SMO-BST_md_complex.tpr` + `SMO-BST_md_complex.xtc` |

> **OPLS-AA**: group identification is validated by `Tests/unittests/test_opls_identifier.py` (charged residues, aromatic rings, etc.), but no real trajectory test data is bundled yet; the usage of `dii run --ff opls` is identical to the other force fields.
>
> The CHARMM36 data also provides the all-atom pair `SMO-BST_md_fullatom.tpr` + `SMO-BST_md_fullatom_100ps.xtc` (used for water-bridge validation in `Tests/unittests/test_charmm_water_bridge_real.py`).

The remainder of this guide uses the Amber test data as the running example.

## 2. Run Interaction Detection

Use `dii run` to detect all 8 interaction types and save them as h5 files:

```bash
dii run -t Tests/test_MD_case_amber/md.tpr -f Tests/test_MD_case_amber/md1ns.xtc -o out_amber/ --ff amber
```

Argument notes:

- `-t` / `--tpr`: GROMACS topology file (required)
- `-f` / `--xtc`: trajectory file (required; any format readable by MDAnalysis)
- `-o` / `--output`: output directory (required; created automatically)
- `--ff`: force field, one of `amber` / `gromos` / `charmm` / `opls` (required; corresponds one-to-one to the keys of the `IDENTIFIER_CLASSES` registry in `DuIvyInteractions/group_identifiers/__init__.py`)
- By default, all 8 types are detected and saved as one h5 file per type: `<output>/<type>.h5`

After a successful run, the output directory contains:

```
$ ls out_amber/
hydrogen_bond.h5          metal_coordination.h5     salt_bridge.h5
hydrophobic.h5            pi_cation.h5              water_bridge.h5
halogen_bond.h5           pi_stacking.h5
```

> `dii run` prints no intermediate output on success; if a type fails, it prints `[WARN] <type> detection failed: <reason>` without interrupting the other types (per-type try/except in `Pipeline.run`, see `DuIvyInteractions/DII.py`).
>
> On this test set (101 frames), the `two_pass` strategy detects **47** salt-bridge pairs (this number is the assertion baseline of `Tests/unittests/test_saltbridge_two_pass.py`, with no missing distance values over all frames) and **10** π-stacking pairs (the pairs actually formed among the 38 aromatic rings; re-verified against the current code).

### Other Force Fields

GROMOS 53A6:

```bash
dii run -t Tests/test_MD_case_gromos/gromos53a6_md.tpr \
        -f Tests/test_MD_case_gromos/gromos53a6_md10ns.xtc \
        -o out_gromos/ --ff gromos
```

CHARMM36:

```bash
dii run -t Tests/test_MD_case_charmm36/SMO-BST_md_complex.tpr \
        -f Tests/test_MD_case_charmm36/SMO-BST_md_complex.xtc \
        -o out_charmm/ --ff charmm
```

> Water exclusion (so that water is not treated as a protein donor or coordination site in water-bridge / metal-coordination detection) is governed uniformly by the `WATER_RESIDUES` class attribute of each identifier (`DuIvyInteractions/core/interfaces.py`): Amber `{SOL, HOH, WAT}`, GROMOS `{SOL}`, CHARMM `{TIP3, HOH, SOL, WAT}`, OPLS `{HOH, HO4, HO5, SOL, WAT}`. No manual specification is needed — just make sure `--ff` matches the force field actually used for the system.

### Equivalent Python API

```python
from DuIvyInteractions.pipeline import Pipeline

Pipeline(ff="amber", strategy="two_pass").run(
    "Tests/test_MD_case_amber/md.tpr",
    "Tests/test_MD_case_amber/md1ns.xtc",
    "out_amber/",
    interactions=None,          # None = all 8 types; or pass a list of type names
)
```

## 3. Export Results

Use `dii export` to export the h5 results to xvg/xpm/csv and print an overview:

```bash
dii export -i out_amber/salt_bridge.h5 -o out_export/
```

Terminal overview (example; the concrete residues and occupancies depend on the system and version):

```
===== Salt Bridge overview =====
Type:     salt_bridge
Pairs:    47
Frames:   101
Time range: 0.0 ~ 1000.0 ps

Top 5 occupancies:
  1. ARG210(3443-3451)···ASP211(3461-3463)  100.0%
  2. LYS70(1157-1160)···ASP180(2983-2985)  100.0%
  3. ARG291(4737-4745)···ASP295(4795-4797)  100.0%
  4. ARG244(4009-4017)···ASP211(3461-3463)  100.0%
  ...
```

> Meaning of the overview: `Pairs` is the number of group pairs judged to exist in at least one frame (47 in this example, matching the unit-test baseline); `Frames` is the total number of trajectory frames (101); `Time range` is the first-to-last frame time; `Top 5 occupancies` is sorted in descending occupancy, where occupancy = number of frames in which the pair exists / total frames.

The following files are generated in the output directory:

```
$ ls out_export/
salt_bridge_count.xvg        # number of active salt bridges per frame
salt_bridge_distance.xvg     # distance time series per pair (charge-center distance, Å)
salt_bridge_existence.xpm    # existence heatmap (rows = pairs, columns = frames)
salt_bridge_summary.csv      # per-pair summary (pair_label, occupancy, avg/std)
```

> When exporting π-stacking (`pi_stacking.h5`), an additional `<type>_type.xpm` file (stacking-type map: none / T-shaped / parallel) is produced; if the h5 data is empty (0 pairs or 0 frames), `[skip] <type>: nothing to export` is printed and the export is skipped.

## 4. Common Arguments

Detect only a subset of interaction types (comma-separated; `all` cannot be mixed with other types — mixing aborts with an error):

```bash
dii run -t md.tpr -f md.xtc -o out/ --ff amber \
    --interactions hydrogen_bond,pi_stacking
```

Switch the detection strategy (default: `two_pass`, recommended for large systems; `per_frame` is the frame-by-frame vectorized implementation; `per_tuple` iterates the trajectory per candidate tuple and is for reference only):

```bash
dii run -t md.tpr -f md.xtc -o out/ --ff amber --strategy two_pass
```

## 5. Performance Reference

Reference timings for water-bridge detection on the Amber test system (KRAS-RBD, 116,383 atoms, 101 frames; numbers from the project design documents and `Tests/unittests/`, measured in the development environment; actual timings depend on the system size and machine):

| Detection Strategy | Water Bridge Time | Notes |
|:---------|:---------|:-----|
| `two_pass` (default) | ~5 s | Pass1 discovers active triples frame by frame with KDTree (pre-screening radius 8.2 Å) + Pass2 fills full-frame metrics; sparse storage |
| `per_frame` | ~5 s | First-frame KDTree pre-screening of candidate triples + frame-by-frame vectorized computation (~56 ms/frame); pre-allocates a dense matrix for all candidates |
| `per_tuple` | ~65 h | Iterates the trajectory per candidate triple; extremely slow with many candidates; reference only |

Water-bridge candidate triples can reach 249,000; `per_tuple` iterating the trajectory triple by triple is the root cause of its slowness. **For large systems / long trajectories / large candidate sets, use the default `two_pass`**.

## 6. Known Limitations and Further Reading

- The full list of known limitations (PBC not handled, missing Ow-A lower bound in the two-pass water bridge, hydrophobic-π de-duplication not implemented, PerFrame memory usage on long trajectories, etc.) is in {doc}`/reference/limitations`.
- Detection criteria and result interpretation for the 8 interaction types, see [Interpreting Results](result).
- All `dii` commands and arguments, see [Command Reference](command).
- Core principles (why tpr atom types are read directly), see {doc}`/reference/concepts`.