# Quick Start

This guide runs the full workflow with the test data bundled with the repository (`Tests/test_MD_case/`), demonstrating the real usage and output of the `dii` command. The data is a 1 ns trajectory (101 frames) of a KRAS-RBD system.

## 1. Prepare Data

The test data is located under `Tests/test_MD_case/` in the repository and includes:

- `md.tpr` — GROMACS topology file (Amber-family force field: amber14sb protein + GAFF ligand)
- `md1ns.xtc` — 1 ns trajectory (101 frames)

## 2. Run Interaction Detection

Use `dii run` to detect all 8 interaction types and save them as h5 files:

```bash
dii run -t Tests/test_MD_case/md.tpr -f Tests/test_MD_case/md1ns.xtc -o out/ --ff amber
```

- `-t` / `-f`: topology and trajectory files
- `-o`: output directory (created automatically)
- `--ff`: force field type; currently only `amber` is supported
- By default, all 8 types are detected and saved as `out/<type>.h5`

`dii run` prints no intermediate output on success; after detection, one h5 file per type is generated in the output directory:

```
$ ls out/
hydrogen_bond.h5  pi_stacking.h5  salt_bridge.h5  hydrophobic.h5
halogen_bond.h5   metal_coordination.h5  water_bridge.h5  pi_cation.h5
```

> Note: if a type fails, `[WARN] <type> detection failed: <reason>` is printed without interrupting the other types. On the test data, 47 salt-bridge pairs and 10 π-stacking pairs were identified (101 frames).

## 3. Export Results

Use `dii export` to export the h5 results to xvg/xpm/csv and print an overview:

```bash
dii export -i out/salt_bridge.h5 -o out_export/
```

**Real output** (overview + generated files):

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
  5. ARG5(82-90)···PRO141(2352-2354)        100.0%

Output files:
salt_bridge_count.xvg         # number of active interactions per frame
salt_bridge_distance.xvg      # distance time series per pair
salt_bridge_existence.xpm     # existence heatmap (rows = pairs, columns = frames)
salt_bridge_summary.csv       # summary statistics per pair
```

> In the table above, Top 1 shows `ARG210···ASP211` with 100% occupancy, indicating that this salt bridge exists in all 101 frames — i.e., a salt bridge contact stably formed between RBD and the ligand/receptor.

## 4. More Parameters

Detect only a subset of interaction types:

```bash
dii run -t md.tpr -f md.xtc -o out/ --ff amber \
    --interactions hydrogen_bond,pi_stacking
```

Switch the detection strategy (default: `two_pass`, recommended for large systems; `per_tuple`/`per_frame` are reference implementations):

```bash
dii run -t md.tpr -f md.xtc -o out/ --ff amber --strategy two_pass
```

## Performance Reference

The following is reference performance on the test system (KRAS-RBD, 116,383 atoms, 101 frames):

| Detection Strategy | Water Bridge Time | Notes |
|:---------|:---------|:-----|
| `two_pass` (default) | ~5 s | Pass1 discovers active pairs frame by frame + Pass2 completes them; water bridge uses KDTree pre-screening |
| `per_frame` | ~5 s | Vectorized per-frame processing of all candidates |
| `per_tuple` | ~65 h | Iterates over all frames per candidate pair; extremely slow when there are many candidate tuples; for reference only |

Water bridge candidate tuples can reach 249,000, and the `per_tuple` strategy iterating the trajectory pair by pair is the root cause of its slowness. **For large systems/long trajectories, use the default `two_pass`**.

## Next Steps

- Learn about the detection criteria and result interpretation for the 8 interaction types, see [Interpreting Results](result)
- View all `dii` commands and parameters, see [Command Reference](command)
- Understand the core principles (why tpr types are read directly), see {doc}`/reference/concepts` (in the Reference)