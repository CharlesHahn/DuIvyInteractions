# Interpreting Results

`dii export` exports h5 results into three types of files: xvg / xpm / csv, and prints an overview to the terminal. This document describes the format and meaning of each file type and how to interpret them. The structure of the h5 file itself is described in {doc}`/reference/data_format`.

## Overview (Printed to the Terminal)

When running `dii export`, an overview is printed for each exported Interaction (`_print_overview` in `DII.py`):

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

- **Pairs**: the number of group pairs of this type (`Interaction.n_pairs`). The detectors keep only pairs judged to exist in at least one frame; candidate pairs that never occur do not enter the results or the exported files.
- **Frames / time range**: the total number of frames and the first-to-last frame time (ps, one decimal place).
- **Top N occupancy**: the N pairs with the highest occupancy (`OVERVIEW_TOP_N = 5`). Occupancy = number of frames in which the pair exists / total number of frames (`Interaction.occupancy()`), i.e. the fraction of the simulation during which the interaction is present, in the range 0–1.

## Output File Overview

Five file types are generated per Interaction (π-stacking additionally produces type.xpm):

| File | Content | Format |
|:-----|:-----|:-----|
| `<type>_count.xvg` | Number of active interactions per frame | GROMACS xvg time series |
| `<type>_<metric>.xvg` | Time series of each numeric metric (one curve per group pair) | GROMACS xvg time series |
| `<type>_existence.xpm` | Existence heatmap (rows = pairs, columns = frames) | GROMACS xpm matrix |
| `<type>_type.xpm` | π-stacking only: stacking-type map (none / T-shaped / parallel) | GROMACS xpm matrix |
| `<type>_summary.csv` | Per-pair summary statistics | CSV |

Filenames are based on the type rather than the index in the h5, corresponding one-to-one to the h5 filenames produced by `dii run` (`<type>.h5`).

## XVG (Time Series)

xvg is the standard GROMACS time series format and can be plotted with DuIvyTools or Xmgrace (the first column is time, in ps).

**count.xvg**: in `<type>_count.xvg` the Y axis is the number of active interactions per frame (frame-wise sum of `existence`), used to observe how the overall activity changes over time (e.g. whether the number of salt bridges decreases in a segment of the simulation).

**metric.xvg**: in `<type>_<metric>.xvg` the first column is time and each subsequent column is one curve corresponding to one group pair (legend = pair label), giving the value of the metric over **all frames**. Note:

- Metric values are written for all frames (`two_pass` fills full-frame metrics for the tuples discovered in Pass1; `per_frame` computes full-frame metrics for candidate pairs that exist in at least one frame). An individual geometric quantity may be NaN in a frame where the geometry is undefined (shown as a break in the plot).
- **The xvg itself does not mark in which frames a pair is "present"** — frame-level existence is judged by the `existence.xpm` and by the occupancy in `summary.csv`. To decide when a pair is established, consult the existence heatmap or the CSV summary.

## CSV Summary (Per-Pair Statistics)

`<type>_summary.csv` records the summary statistics of each pair in a row-column table, **including only pairs with occupancy > 0** (`to_csv_summary` skips rows with occupancy = 0):

| Column | Meaning |
|:---|:-----|
| `pair_label` | Pair label (format below) |
| `occupancy` | Occupancy (0–1, four decimal places) |
| `avg_<metric>` | Mean of the metric over **active frames** (frames where the pair exists) |
| `std_<metric>` | Standard deviation of the metric over active frames |

The statistics are computed over active frames only (the code takes `metrics[i][existence[i]]` and uses nanmean/nanstd); if all active-frame values of a metric are NaN, the corresponding cells are left empty. String metrics (e.g. the π-stacking `pistacking_type`) do not enter the CSV.

Example (illustrative values):

```
pair_label,occupancy,avg_distance,std_distance
ARG210(3443-3451)···ASP211(3461-3463),1.0000,3.6200,0.0548
```

Interpretation: the `ARG210···ASP211` salt bridge has 100% occupancy, an average distance of 3.62 Å, and a standard deviation of 0.055 Å — the salt bridge exists stably throughout the simulation with a compact conformation. Conversely, a pair with low occupancy and large std indicates an intermittent interaction or pronounced geometric fluctuations.

## XPM (Heatmaps)

xpm is the standard GROMACS matrix heatmap format, convenient for inspecting which pairs are present (or in which conformation) at which frames.

### Existence Heatmap

`<type>_existence.xpm`:

- **Rows**: pairs; **Columns**: frames (X axis is time)
- **Colors**: white = the interaction is absent in that frame (No), blue (#38A7D0) = present (Yes)
- **Use**: quickly inspect the temporal existence pattern of each pair — continuously stable (all-blue row), intermittent (alternating blue/white), or present only in part of the time.

### π-Stacking Type Heatmap

`<type>_type.xpm` (π-stacking only):

- Three-valued colors: white = no stacking (None), pink (#F67088) = T-shaped stacking, blue (#38A7D0) = parallel stacking
- **Use**: distinguish the geometric conformation of π stacking (T-shaped / parallel) over time; existence is judged by existence.xpm, and type.xpm is colored only on frames where the pair exists.

## Metrics per Type

The `<type>_<metric>.xvg` files and the metric columns of `summary.csv` are defined per type (units from the `metric_labels` of each exporter):

| Type | Metrics (meaning / unit) |
|:-----|:-----|
| `hydrogen_bond` | `distance`: D-A distance (Å); `angle`: D-H···A angle (°) |
| `salt_bridge` | `distance`: charge-center distance (Å) |
| `pi_stacking` | `distance`: ring-center distance (Å); `angle`: angle between ring normal vectors (°); `offset`: offset (Å) |
| `hydrophobic` | `distance`: atom-atom distance (Å) |
| `halogen_bond` | `distance`: X···A distance (Å); `don_angle`: C-X···A angle (°); `acc_angle`: X···A-R angle (°) |
| `metal_coordination` | `distance`: metal-ligand distance (Å) |
| `water_bridge` | `dist_dw`: D-Ow distance (Å); `dist_wa`: Ow-A distance (Å); `theta`: O-H···Ow angle (°); `omega`: H-Ow···A angle (°) |
| `pi_cation` | `distance`: ring-charge distance (Å); `offset`: offset (Å) |

## Pair Label Format

Pair labels are used to locate the atoms of an interaction on the structure (the `get_pair_label` implementation of each exporter). Atom numbers are **global atom indices** in the tpr:

| Type | Label format | Example |
|:-----|:---------|:-----|
| Salt bridge / π-cation / π-π stacking | `ResidueName ResidueNumber(startAtom-endAtom)···ResidueName ResidueNumber(startAtom-endAtom)` (the range is the minimum–maximum global index of the atoms in the group) | `ARG210(3443-3451)···ASP211(3461-3463)` |
| Hydrogen bond | `ResidueName:donorAtom(index)-hydrogen(index)···ResidueName:acceptorAtom(index)` | `ARG210:NE(3443)-HE(3444)···ASP211:OD1(3461)` |
| Water bridge | `donorResidue:donorAtom(index)-hydrogen(index)···waterResidue:OW(index)···acceptorResidue:atom(index)` | `ARG210:NE(3443)-HE(3444)···SOL:OW(5000)···ASP211:OD1(3461)` |
| Halogen bond | `donorResidue:carbon(index)-halogen(index)···acceptorResidue:atom(index)` (donor shown as C-X; the acceptor-side R is chosen dynamically per frame and does not appear in the label) | `D927:C(9)-CL(10)···ASP211:OD1(3461)` |
| Hydrophobic interaction | `ResidueName:atomName(index)···ResidueName:atomName(index)` (first atom of each group) | `D927:C7(100)···PHE232:CG(3700)` |
| Metal coordination | `metalResidue:metalAtom(index)···coordinatingResidue:coordinatingAtom(index)` | `MG:MG(100)···ASP211:OD1(3461)` |

## Interpretation Suggestions

- **Find stable interactions**: sort `summary.csv` by occupancy in descending order; pairs with high occupancy (e.g. > 0.8) and small std are stable candidate hotspots; for low-occupancy pairs, check existence.xpm to see whether they are intermittent or formed only late in the simulation.
- **Inspect temporal evolution**: use existence.xpm for the continuity/intermittency pattern of each pair over time, count.xvg for overall activity, and type.xpm (π-stacking) for conformational switching.
- **Check geometric plausibility**: compare the values in `<metric>.xvg` with the detection thresholds (e.g. hydrogen bond D-A ≤ 4.1 Å, D-H···A ≥ 100°; criteria for all types in {doc}`/reference/criteria`) to confirm that the detected pairs indeed fall within the criteria.

## Next Steps

- To understand the detection criteria and geometric definitions for the 8 interaction types, see {doc}`/reference/criteria`.
- To understand the result data format (h5), see {doc}`/reference/data_format`.