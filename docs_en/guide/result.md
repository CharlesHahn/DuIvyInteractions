# Interpreting Results

`dii export` exports h5 results into three types of files: xvg / xpm / csv. This document describes the format and meaning of each file type and how to interpret them.

## Overview (Printed to the Terminal)

When running `dii export`, an overview of the type is printed first:

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

- **Pair count**: the number of all candidate pairs under this type
- **Frames / time range**: the total number of frames and the time span of the trajectory
- **Top N occupancy**: the N pairs with the highest occupancy (default 5). Occupancy = number of frames in which the pair exists / total number of frames, i.e. the fraction of time this interaction appears throughout the simulation

## CSV Summary (Per-Pair Statistics)

`<type>_summary.csv` records the summary statistics of each pair in a row-column table. Only pairs with occupancy > 0 are included.

| Column | Meaning |
|:---|:-----|
| `pair_label` | Pair label (e.g. `ARG210(3443-3451)···ASP211(3461-3463)`) |
| `occupancy` | Occupancy (0~1) |
| `avg_<metric>` | Mean of the metric over active frames (e.g. average distance, in Å) |
| `std_<metric>` | Standard deviation of the metric over active frames |

Example:

```
pair_label,occupancy,avg_distance,std_distance
ARG210(3443-3443)···ASP211(3461-3463),1.0000,3.6200,0.0548
```

Interpretation: the `ARG210···ASP211` salt bridge has 100% occupancy, an average distance of 3.62 Å, and a standard deviation of 0.055 Å — indicating that this salt bridge exists stably throughout the simulation with a compact conformation.

## XVG (Time Series)

xvg is the standard GROMACS time series format and can be plotted with DuIvyTools or Xmgrace.

| File | Content |
|:-----|:-----|
| `<type>_count.xvg` | X-axis: time; Y-axis: the number of active interactions per frame |
| `<type>_<metric>.xvg` | One curve per column, the value of a metric (e.g. distance, angle) over all frames |

The metric of inactive frames is NaN (shown as breaks in the plot); it should be interpreted together with the existence heatmap.

## XPM (Heatmaps)

xpm is the standard GROMACS matrix heatmap format; rows and columns form a matrix, and colors represent values.

### Existence Heatmap

`<type>_existence.xpm`

- **Rows**: pairs; **Columns**: frames
- **Colors**: white = the interaction is absent in that frame (No), blue = present (Yes)
- **Use**: quickly inspect the temporal existence pattern of each pair (continuous / intermittent / present only in part of the time)

### π-Stacking Type Heatmap

`<type>_type.xpm` (π-stacking only)

- Three-valued colors: white = none, pink = T-shaped stacking, blue = parallel stacking
- **Use**: distinguish the geometric conformations of π stacking

## Pair Label Format

The meaning of the labels differs by type; they are used to locate the atoms of the interaction on the structure:

| Type | Label format | Example |
|:-----|:---------|:-----|
| Salt bridge / π-cation / π-π stacking | `ResidueName ResidueNumber(startAtom-endAtom)···ResidueName ResidueNumber(startAtom-endAtom)` | `ARG210(3443-3451)···ASP211(3461-3463)` |
| Hydrogen bond | `ResidueName:donorAtom(index)-hydrogen(index)···ResidueName:acceptorAtom(index)` | `ARG210:NE(3443)-HE(3444)···ASP211:OD1(3461)` |
| Water bridge | `donorResidue:donorAtom(index)-hydrogen(index)···waterResidue:OW(index)···acceptorResidue:atom(index)` | `ARG210:NE(3443)-HE(3444)···SOL:OW(5000)···ASP211:OD1(3461)` |
| Halogen bond | `donorResidue:carbon(index)-halogen(index)···acceptorResidue:atom(index)` (donor shown as C-X) | `D927:C(9)-CL(10)···ASP211:OD1(3461)` |
| Hydrophobic interaction | `ResidueName:atomName(index)···ResidueName:atomName(index)` | `D927:C7(100)···PHE232:CG(3700)` |
| Metal coordination | `metalResidue:metalAtom(index)···coordinatingResidue:coordinatingAtom(index)` | `MG:MG(100)···ASP211:OD1(3461)` |

## Next Steps

- To understand the detection criteria for the 8 interaction types, see {doc}`/reference/criteria` (Reference)
- To understand the result data format (h5), see {doc}`/reference/data_format` (Reference)