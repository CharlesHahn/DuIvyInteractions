# Interaction Criteria

Each interaction type is evaluated by its corresponding detector during the **geometric evaluation** stage (stage ② of the two-stage architecture, see [Core Concepts](concepts.md)), using distance/angle criteria to determine whether it exists in each frame. This document lists the criteria and thresholds of the 8 interaction types, for understanding the results and tuning parameters.

> All distances are in Å (Ångström) and all angles in degrees (°). The threshold constants are defined at the top of the detector files (`DuIvyInteractions/interaction_detectors/*_detector_{per_frame,two_pass}.py`).
>
> **Strategy-consistency note**: `per_frame` and `two_pass` use identical thresholds, with **one exception — the water bridge**: per_frame additionally applies a 2.5 Å lower bound on the distances, while two_pass has only the 4.1 Å upper bound (see "Water Bridge" below). The `per_tuple` strategy is no longer maintained (skipped in the test conventions); its water-bridge thresholds match per_frame.

## Geometric Quantities

The criteria involve the following geometric quantities (all computed per frame and per pair):

| Geometric quantity | Definition | Used for |
|:-------|:-----|:-----|
| `distance` | Distance between the reference points of two groups. The reference point varies by type: hydrogen bond = donor atom D and acceptor atom A; salt bridge = **partial-charge-weighted center** of the positive/negative group (centroid weighted by atomic partial charges); π stacking = geometric center of the aromatic ring; π-cation = ring center and cation charge center; halogen bond = halogen X and acceptor A; metal coordination = metal and coordinating atom | All types |
| `angle` (hydrogen bond) | Three-point angle of donor D, hydrogen H, acceptor A (vertex at H), i.e., the D-H···A angle | Hydrogen bond |
| `angle` (π stacking) | Angle between the normal vectors of two aromatic rings (the smaller of `raw` and `180°−raw`, handling orientation ambiguity; ring normal = average of cross products of adjacent ring atoms) | π stacking |
| `offset` (π stacking) | Lateral offset of one ring center relative to the other ring's plane: project the ring center onto the other ring's plane along its normal, and take the distance between the projected point and the other ring's center; each ring serves as reference for the other, and the smaller value is taken | π stacking, π-cation interaction |
| `planarity` | Maximum pairwise angle between local normal vectors of atoms in the ring (cross products of two adjacent bonds) | π stacking (optional) |
| `don_angle` | Three-point angle of carbon C, halogen X, acceptor A (vertex at X), i.e., the C-X···A angle | Halogen bond |
| `acc_angle` | Among the atoms adjacent to acceptor A, the X···A-R angle closest to 120° | Halogen bond |
| `theta` (water bridge) | Three-point angle of donor D, hydrogen H, water oxygen Ow (vertex at H), i.e., D-H···Ow | Water bridge |
| `omega` (water bridge) | Three-point angle of acceptor A, water oxygen Ow, hydrogen H (vertex at Ow), i.e., A-Ow···H | Water bridge |

## Hydrogen Bond (hydrogen_bond)

**Criterion**: the distance between donor D and acceptor A is ≤ 4.1 Å, and the D-H···A angle is ≥ 100°.

| Constant | Value | Meaning |
|:-----|:--|:-----|
| `HBOND_DIST_MAX` | 4.1 | Maximum D-A distance (inclusive) |
| `HBOND_DON_ANGLE_MIN` | 100 | Minimum D-H···A angle |

## π-π Stacking (pi_stacking)

**Criterion**: the distance between the centers of two aromatic rings ∈ (0.5, 5.5] Å, and the stacking type is not N (`pistacking_type != 'N'`, i.e., classified as T-type or P-type); planarity satisfaction is an optional condition (disabled by default, `check_planarity=False`).

| Constant | Value | Meaning |
|:-----|:--|:-----|
| `PISTACK_MIN_DIST` | 0.5 | Minimum distance (exclusive, excludes self) |
| `PISTACK_DIST_MAX` | 5.5 | Maximum ring-center distance (inclusive) |
| `PISTACK_OFFSET_MAX` | 2.0 | Upper limit of ring-center offset |
| `PISTACK_PLANARITY` | 5.0 | Upper limit of ring planarity deviation (°), disabled by default |
| `PISTACK_ANG_DEV` | 30 | Angle deviation for T/P classification |

**Stacking type classification** (`pistacking_type` metric; identical in both strategies), based on the angle `angle` between the two ring normal vectors and the offset `offset`:

- **P-type (parallel stacking)**: `angle ≤ 30°` and `offset < 2.0 Å`
- **T-type (edge-to-face stacking)**: `angle ≥ 60°` and `offset < 2.0 Å`
- **N (none)**: otherwise

In the output XPM: 0 = none (white `#FFFFFF`), 1 = T-type (pink `#F67088`), 2 = P-type (blue `#38A7D0`).

## Salt Bridge (salt_bridge)

**Criterion**: the distance between the **partial-charge-weighted centers** of a positive and a negative charged group is ≤ 5.5 Å (charge center = centroid of the group atoms weighted by partial charges).

| Constant | Value | Meaning |
|:-----|:--|:-----|
| `SALTBRIDGE_DIST_MAX` | 5.5 | Maximum charge-center distance (inclusive) |

## Hydrophobic Interaction (hydrophobic)

**Criterion**: the distance between two hydrophobic atoms ∈ (0.5, 4.0) Å (both endpoints exclusive).

| Constant | Value | Meaning |
|:-----|:--|:-----|
| `HYDROPH_MIN_DIST` | 0.5 | Minimum distance |
| `HYDROPH_DIST_MAX` | 4.0 | Maximum distance |

## Halogen Bond (halogen_bond)

**Criterion**: the distance between halogen X and acceptor A is ≤ 4.0 Å, the C-X···A angle ∈ [135°, 195°] (165°±30°), and the X···A-R angle ∈ [90°, 150°] (120°±30°, using the R neighbor of A closest to 120°).

| Constant | Value | Meaning |
|:-----|:--|:-----|
| `HALOGEN_DIST_MAX` | 4.0 | Maximum X···A distance (inclusive) |
| `HALOGEN_DON_ANGLE` | 165 | Optimal C-X···A angle |
| `HALOGEN_ACC_ANGLE` | 120 | Optimal X···A-R angle |
| `HALOGEN_ANGLE_DEV` | 30 | Upper limit of angle deviation |

## Metal Coordination (metal_coordination)

**Criterion**: the distance between the metal center (e.g., Mg²⁺) and coordinating atoms is < 3.0 Å. Coordinating atoms come from the `metal_binding` group (atoms of element ∈ {O, N, S}).

| Constant | Value | Meaning |
|:-----|:--|:-----|
| `METAL_DIST_MAX` | 3.0 | Maximum metal-coordinating atom distance (exclusive) |

> Note: coordination by pure water molecules is not reported — the `metal_binding` groups exclude water residues entirely at the identification stage (per-force-field `WATER_RESIDUES`); the first release applies only the distance criterion and does not match coordination geometries.

## Water Bridge (water_bridge)

**Criterion (two_pass, the default strategy)**: the distance from water oxygen Ow to polar atoms (donor D / acceptor A) is < 4.1 Å, with theta ≥ 100° and 71° ≤ omega ≤ 140°.

| Constant | Value | Meaning |
|:-----|:--|:-----|
| `WATER_BRIDGE_MAXDIST` | 4.1 | Maximum Ow-to-polar-atom distance |
| `WATER_BRIDGE_THETA_MIN` | 100 | Minimum D-H···Ow angle |
| `WATER_BRIDGE_OMEGA_MIN` | 71 | Minimum A-Ow···H angle |
| `WATER_BRIDGE_OMEGA_MAX` | 140 | Maximum A-Ow···H angle |

**Strategy difference (stated as-is)**:

- **per_frame**: `apply_threshold` additionally requires a **2.5 Å lower bound** on both distances — `2.5 < dist_dw < 4.1` and `2.5 < dist_wa < 4.1` (`WATER_BRIDGE_MINDIST = 2.5`, excluding overlapping/overly close atom pairs); candidate generation also requires the first-frame D···A distance ≥ 2.5 Å.
- **two_pass**: **lacks the Ow-D and Ow-A distance lower bound**, requiring only `dist_dw < 4.1` and `dist_wa < 4.1`. The angle thresholds (theta ≥ 100, 71 ≤ omega ≤ 140) are identical in both strategies.
- **per_tuple** (no longer maintained): same as per_frame (including the 2.5 Å lower bound).

Therefore, for the same system, per_frame and two_pass water-bridge results may differ in situations where a water molecule comes closer than 2.5 Å to a polar atom; the default two_pass strategy does not apply this lower bound.

## π-Cation Interaction (pi_cation)

**Criterion**: the distance from the ring center to the **charge center** of the cation (positively charged group) ∈ (0.5, 6.0) Å (both endpoints exclusive), with offset < 2.0 Å.

| Constant | Value | Meaning |
|:-----|:--|:-----|
| `PICATION_MIN_DIST` | 0.5 | Minimum distance |
| `PICATION_DIST_MAX` | 6.0 | Maximum ring-cation distance |
| `PICATION_OFFSET_MAX` | 2.0 | Upper limit of offset |

## Metrics by Type (metrics)

Metrics output by the detectors (stored in the h5 `metrics`, exported to xvg):

| Type | Metrics | Unit |
|:-----|:-----|:-----|
| hydrogen bond | distance, angle | Å, ° |
| π stacking | distance, angle, offset, pistacking_type (plus planarity_ring1/2 when planarity checking is enabled) | Å, °, Å, type |
| salt bridge | distance | Å |
| hydrophobic | distance | Å |
| halogen bond | distance, don_angle, acc_angle | Å, °, ° |
| metal coordination | distance | Å |
| water bridge | dist_dw, dist_wa, theta, omega | Å, Å, °, ° |
| π-cation interaction | distance, offset | Å, Å |

## Criterion Sources

The criteria and thresholds follow the interaction definition system of PLIP (Protein-Ligand Interaction Profiler; Salentin et al. 2015): the distance/angle cutoffs (e.g., hydrogen bond D-A 4.1 Å, salt bridge 5.5 Å, halogen bond 165°±30°) are all consistent with PLIP, making it easy to compare with existing PLIP results. The π stacking T/P-type classification standard is consistent with PLIP's `pi-stacking` classification (methodological reference: McGaughey et al. 1998); the water-bridge angle definitions follow Jiang et al. 2005; the halogen-bond geometry follows Auffinger et al. (as annotated in the detector comments). The concrete thresholds are governed by the constants at the top of the detector files.

## Interpreting Results: Key Points

- **existence**: whether a given pair exists in each frame (boolean matrix after criterion filtering)
- **metrics**: metric values per pair per frame (an (n_pairs, n_frames) array in `Interaction`; inactive frames are NaN-filled in two_pass)
- **occupancy**: number of frames in which a given pair exists / total number of frames (`Interaction.occupancy()`), measuring the persistence of that interaction