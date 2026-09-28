# dii Command Reference

`dii` is the command-line entry point of DuIvyInteractions, providing two subcommands: `run` (detection) and `export` (export).

## Global Usage

```bash
dii [-h] {run,export} ...
```

## Subcommand: run — Run Interaction Detection

Detect interactions in the specified trajectory and save the results as h5 files.

```bash
dii run -t TPR -f XTC -o OUTPUT --ff FF [--interactions LIST] [--strategy STRATEGY]
```

### Arguments

| Argument | Required | Default | Description |
|:-----|:----|:----|:-----|
| `-t, --tpr` | ✅ | — | GROMACS topology file |
| `-f, --xtc` | ✅ | — | Trajectory file |
| `-o, --output` | ✅ | — | Output directory (created automatically) |
| `--ff` | ✅ | — | Force field type; currently only `amber` is supported |
| `--interactions` | No | `all` | Interaction types to detect, comma-separated (e.g. `hydrogen_bond,pi_stacking`) |
| `--strategy` | No | `two_pass` | Detection strategy: `two_pass` / `per_frame` / `per_tuple` |

### The 8 Supported Interaction Types

`hydrogen_bond` (hydrogen bond), `pi_stacking` (π-π stacking), `salt_bridge` (salt bridge), `hydrophobic` (hydrophobic interaction), `halogen_bond` (halogen bond), `metal_coordination` (metal coordination), `water_bridge` (water bridge), `pi_cation` (π-cation interaction)

### Output

One h5 file is generated per type: `<output>/<interaction_type>.h5`. For example:

```bash
dii run -t md.tpr -f md.xtc -o out/ --ff amber
# generates out/hydrogen_bond.h5, out/pi_stacking.h5, out/salt_bridge.h5, ...
```

### Strategy Descriptions

| Strategy | Mechanism | Use Case |
|:-----|:-----|:---------|
| `two_pass` (default) | Pass 1 discovers active pairs frame by frame (sparse) → Pass 2 completes all frames | Large systems, best performance (water bridge: 65h → ~5s after KDTree pre-screening) |
| `per_frame` | Processes all candidate pairs frame by frame, vectorized | Types with a very large number of candidate pairs (e.g. water bridge) |
| `per_tuple` | Loads all frames per candidate pair | Control implementation, for reference |

## Subcommand: export — Export Results

Reads h5 files, exports them as xvg/xpm/csv, and prints an overview.

```bash
dii export -i INPUT -o OUTPUT
```

### Arguments

| Argument | Required | Description |
|:-----|:----|:-----|
| `-i, --input` | ✅ | h5 file path |
| `-o, --output` | ✅ | Output directory (created automatically) |

### Behavior

- Supports h5 containing multiple Interactions (`dii run` generates one file per type, which fits this scenario)
- If the h5 contains multiple Interactions of the same type, a sequence number is appended to the filename to avoid overwrites (e.g. `salt_bridge_2_*`)
- Empty data (0 pairs or 0 frames) is automatically skipped with a message
- A corrupted or missing h5 produces a friendly error

### Output Files

For each Interaction, the following files are generated in the output directory:

| File | Content |
|:-----|:-----|
| `<type>_count.xvg` | Number of active interactions per frame |
| `<type>_<metric>.xvg` | Time series of each numeric metric (e.g. distance, angle) |
| `<type>_existence.xpm` | Existence heatmap (rows = pairs, columns = frames) |
| `<type>_type.xpm` | π-stacking only: stacking type map (0=none / 1=T-shaped / 2=P-type) |
| `<type>_summary.csv` | Summary statistics for each pair |

## Exit Codes

- `0`: success
- Non-zero: argument errors, missing/corrupted files, unknown types, etc. (raised via `SystemExit`)