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
        [--select SELEXP ...] [--pairing any|between]
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
| `--select` | No | — | MDAnalysis selection string, `action="append", nargs="+"`: may be **repeated** or given as **space-separated multi-values** (e.g. `--select "segid A" "segid B"`); each value defines an **independent** atom set (never merged). A group must be **fully covered by one select** to enter the participation pool (see "User Group Selection & Filtering" below) |
| `--pairing` | No | `any` | Pairing rule: `any` = all pool candidates without a pairing constraint (default); `between` = keep only candidate pairs whose two ends belong to different selects (requires ≥ 2 selects, otherwise an error) |

The four `--ff` values map to the identifiers (`DuIvyInteractions/group_identifiers/__init__.py`):

| Value | Identifier | Covered force fields |
|:-----|:-----|:-----|
| `amber` | `AmberFFGroupIdentifier` | Amber 03/94/96/99/99sb/99sb-ildn/GS/14sb protein + GAFF/GAFF2 ligands |
| `gromos` | `GromosFFGroupIdentifier` | GROMOS 53A6 / 54A7 (united-atom force field, SPC/SPC-E water) |
| `charmm` | `CharmmFFGroupIdentifier` | CHARMM36 / C36m (including CGenFF ligands) |
| `opls` | `OplsFFGroupIdentifier` | OPLS-AA / L |

An unknown `--ff` aborts and lists all available values; an unknown `--interactions` type aborts and lists all available types (`ALL_INTERACTIONS`); `--select` / `--pairing` configuration errors (syntax errors, empty selects, overlapping selects, too few selects for `between`, an empty pool) abort with a friendly `error: ...` message (see "User Group Selection & Filtering" and "Exit Codes" below).

### The 8 Supported Interaction Types

`hydrogen_bond` (hydrogen bond), `pi_stacking` (π-π stacking), `salt_bridge` (salt bridge), `hydrophobic` (hydrophobic interaction), `halogen_bond` (halogen bond), `metal_coordination` (metal coordination), `water_bridge` (water bridge), `pi_cation` (π-cation interaction).

- The `--interactions` value is lowercased, whitespace-stripped, and split on commas (`main` in `DII.py`).
- `all` is the default and is equivalent to "all 8 types"; **`all` cannot be mixed with other types** — mixing aborts with `'all' cannot be mixed with other types`.
- A per-type detection failure (e.g. the system contains no metals or halogens) only prints `[WARN] <type> detection failed: <reason>` without interrupting the other types.

### User Group Selection & Filtering (--select / --pairing)

`--select` and `--pairing` restrict the scope of the analysis: `--select` defines the **participation pool** of groups that may participate in detection, and `--pairing` defines the pairing rule for candidate pairs. The selection syntax reuses the **MDAnalysis Atom selection language** (MDAnalysis is already a dependency of this project — no new syntax is introduced); selection and pairing act only on the group-pool and candidate-pair level, so the detectors' internal geometric criteria and thresholds are **left completely unchanged**. Core usage:

```bash
dii run -t md.tpr -f md.xtc -o out/ --ff amber \
       [--select "SELEXP" ["SELEXP" ...]] [--pairing any|between]
```

- `--select` is declared as `action="append", nargs="+"`: it may be **repeated**, or given as **space-separated multi-values** in a single occurrence; the values are flattened into one list of selection strings. Each value is an MDAnalysis selection string that defines an **independent** atom set, and multiple selects are **never merged**. The two forms below are equivalent:

```bash
dii run -t md.tpr -f md.xtc -o out/ --ff amber --select "segid A" --select "segid B"
dii run -t md.tpr -f md.xtc -o out/ --ff amber --select "segid A" "segid B"
```

- `--pairing` is `any` (default) or `between`: `any` = all pool candidates without a pairing constraint; `between` = keep only **cross-select** candidate pairs (requires ≥ 2 selects, otherwise an error).

#### Syntax Cheat-Sheet

Each `--select` value is an MDAnalysis selection-language string; the commonly used forms:

| Category | Syntax | Examples |
|:-----|:-----|:-----|
| Property selection | `segid` / `resname` / `resid` / `resnum` / `name` / `type` / `moltype` / `index` / `bynum` | `segid seg_0_RBD_pro`, `resname LIG`, `resid 1:30`, `resnum 1 to 5`, `name CA`, `type X`, `moltype M`, `index 0:100`, `bynum 1:100` |
| Standard sets | `protein` / `backbone` / `nucleic` / `nucleicbase` / `nucleicsugar` | `--select "protein"` |
| Boolean combinators | `and` / `or` / `not`, parentheses | `protein and not (resname ALA or resname LYS)` |
| Wildcards | `*` / `[ ]` / `?` | `GL*`, `GL[YN]`, `?` |
| Geometric keywords (evaluated statically at the first frame) | `around` / `byres` / `sphzone` / `point` / `prop` | `around 6 (resname LIG)`, `byres around 6 (resname LIG)`, `sphzone 6 (protein)`, `point X Y Z D`, `prop z > 5` |

> The authoritative reference for the full syntax and the semantics of each keyword is the MDAnalysis official Atom selection language documentation.

#### Semantics

**Group retention rule** (what makes up the participation pool):

- A group is retained ⇔ **all** of its atoms are fully covered by **one** select set (all atoms ∈ a single select);
- Neither "belongs only to the union of several selects" nor a group spanning across selects is accepted;
- Therefore **atom-level** selects (e.g. `name N` + `name H`) make groups span selects → not retained; this is a deterministic rule, not an error. **Residue/molecule-level** selections (`segid` / `resname` / `byres`) never trigger this branch.

**Water side (implicit and constant)**:

- Bridging-water groups with `group_type == "water"` always belong to the "water side" — an implicit select that the user cannot disable and does not need to select; water groups are always retained and skip select assignment;
- H_donor / H_acceptor groups identified on water residues (SOL etc.) do **not** belong to the water side; they are retained only if fully covered by a select (this prevents water-residue donors from polluting the participation pool);
- `between` pairing judgement **skips the water side**.

**Pairing rule (`between`)**:

- Keep only candidate pairs whose two ends belong to **different** selects; pairs with an end assigned to the water side (or not assigned at all) are excluded;
- The constraint-position table is fixed and not user-configurable: hydrogen_bond / pi_stacking / salt_bridge / hydrophobic / halogen_bond / metal_coordination / pi_cation → `(0, 1)` (both ends); water_bridge `(donor, water, acceptor)` → `(0, 2)`, with the water molecule neutral (never participating in pairing judgement);
- Intra-molecular analysis = a single `--select` + the default `any` (no pairing constraint); no extra argument is needed.

**Explicit water selection (scheme A)**:

- Without `--select` (`assignment=None`), water-residue groups are excluded via `WATER_RESIDUES` exactly as in previous versions;
- When water residues are explicitly selected (e.g. `--select "resname SOL"`), the water-residue H_donor / H_acceptor groups assigned to a **concrete select (not the water side)** are **exempted** from the `WATER_RESIDUES` exclusion for non-water-bridge detectors and can participate in hydrogen-bond-type detection (e.g. "protein × water hydrogen bonds").

#### Examples

```bash
# Protein–protein interface (Amber real data: RBD × KRAS, between)
dii run -t Tests/test_MD_case_amber/md.tpr -f Tests/test_MD_case_amber/md1ns.xtc \
        -o out_rbd/ --ff amber \
        --select "segid seg_0_RBD_pro" --select "segid seg_2_KRAS_pro" --pairing between

# Protein–ligand interface (CHARMM36 real data: SMO protein × BST ligand, between)
dii run -t Tests/test_MD_case_charmm36/SMO-BST_md_complex.tpr \
        -f Tests/test_MD_case_charmm36/SMO-BST_md_complex.xtc \
        -o out_smo/ --ff charmm \
        --select "protein" --select "resname BST" --pairing between

# Intra-molecular analysis (single select, pairing defaults to any)
dii run -t Tests/test_MD_case_amber/md.tpr -f Tests/test_MD_case_amber/md1ns.xtc \
        -o out_rbd_only/ --ff amber --select "segid seg_0_RBD_pro"
```

Generic forms (replace the placeholders): protein–ligand interface `--select "protein" --select "resname LIG" --pairing between`; trimer pairwise interfaces `--select "segid A" --select "segid B" --select "segid C" --pairing between` (all pairwise combinations of N selects are retained); ligand pocket residues `--select "byres around 6 (resname LIG)"`; explicit water (protein × water hydrogen bonds) `--select "segid seg_0_RBD_pro" --select "resname SOL" --pairing between`.

> **Geometric-selection note**: without `byres`, `around 6 (resname LIG)` may leave part of a contacting residue's atoms outside the region; such groups are excluded because they are not "fully covered", which can trigger the empty-pool error below. Adding `byres` guarantees that whole residues enter the selected set.

#### Error Messages

Selection configuration errors are always converted into **friendly errors** (`SystemExit`, printing `error: ...` without a traceback). The message text is Chinese by design (it is the literal CLI output); an English gloss is given in parentheses:

| Case | Printed message |
|:-----|:-----|
| select string has a syntax error | `error: select 语法错误: '<expression>' → <MDAnalysis reason>` (select syntax error: `<expression>` → `<MDAnalysis reason>`) |
| a select matches 0 atoms | `error: select #N 未选中任何原子` (select #N matched no atoms) |
| atomic overlap between selects | `error: select #N 与 select #M 原子重叠 (K 个原子)，请用 and not 显式剔除后重试` (select #N overlaps select #M by K atoms; please exclude them explicitly with `and not` and retry) |
| `between` with fewer than 2 selects | `error: --pairing between 需要至少 2 个 --select` (--pairing between requires at least 2 --select) |
| zero non-water groups in the participation pool | `error: 参与池内无可用基团（非水基团 0 个）。若使用了 around/sphzone 几何选择，基团可能因部分原子在圈外被排除，请改用 byres around …` (no usable groups in the pool, 0 non-water groups; if you used an `around`/`sphzone` geometric selection, groups may have been excluded because some atoms fell outside the region — switch to `byres around …`) |

The empty-pool check is based on **zero non-water groups**: bridging-water groups on the water side are always retained, so a wrong user selection (e.g. `--select "name CA"` matching bare atoms) cannot be masked by an all-water pool silently producing empty results. Atomic overlap between selects **must be fixed by the user in the selection string** (e.g. `protein and not ...`); the tool never auto-deduplicates or re-assigns.

#### Limitations

- **Geometric selections are a static approximation**: `around` / `sphzone` / `point` / `prop` are evaluated once against the **first-frame coordinates** and are not dynamic per frame; if the ligand drifts far from its first-frame position, the selected set is not updated. Selection evaluation is the only coordinate-dependent step and is fixed after the first frame — it does not conflict with the "group identification once, frame-independent" architecture;
- **Explicit water + water bridge**: in a water bridge the water molecule is the bridge, not an end; when water is explicitly selected, water still never becomes an end — if the donor and acceptor ends belong to the same select, there are no cross-side candidate pairs → an empty result (designed semantics, not a bug). For water-bridge analysis, keep the implicit water side (do not explicitly select water);
- Custom pairing rules beyond the "cross-select" semantics fall back to the Python API `tuple_filter` (see {doc}`/reference/api`).

### Internal Flow

1. The tpr is read with `GmxTprReader` and the identifier performs group identification once (frame-independent, done a single time);
2. The trajectory is loaded with `mda.Universe(tpr, xtc)`; the water-residue names are taken from the identifier's `WATER_RESIDUES` (water exclusion is correct across force fields);
3. (When `--select` is given) the trajectory is first advanced to the first frame (`u.trajectory[0]`), and each selection string is evaluated into an independent atom set (geometric selections are evaluated statically against the first-frame coordinates); after validating syntax / empty / pairwise overlap, each group is assigned to the select that fully covers it (`group_type == "water"` bridging-water groups always belong to the water side), and the participation pool is filtered accordingly (see "User Group Selection & Filtering" above);
4. For each type, the groups are filtered by the detector's `required_group_types` (water excluded except for the water bridge; water-residue groups explicitly assigned to a concrete select are exempted from the exclusion), with `between` keeping only cross-select candidate tuples via the tuple-level predicate; detection runs frame by frame, and the result is saved as `<output>/<type>.h5`.

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
- `--select` / `--pairing` errors also exit non-zero: a friendly message prefixed with `error:` is printed without a traceback (e.g. `error: select 语法错误: ...`, `error: select #1 未选中任何原子`, `error: --pairing between 需要至少 2 个 --select`, `error: 参与池内无可用基团（非水基团 0 个）...`).

## Next Steps

- Format and meaning of each output file, see [Interpreting Results](result).
- A complete runnable example on real data, see [Quick Start](quickstart).