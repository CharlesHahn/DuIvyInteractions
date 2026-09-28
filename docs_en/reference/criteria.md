# Interaction Criteria

Each interaction type is evaluated by its corresponding detector during the **geometric evaluation** stage, using distance/angle criteria to determine whether it exists in each frame. This document lists the criteria and thresholds of the 8 interaction types, for understanding the results and tuning parameters.

> All distances are in Å (Ångström) and all angles in degrees (°). The threshold constants are defined at the top of the detector files, and the three detection strategies (two_pass / per_frame / per_tuple) use identical thresholds.

## Geometric Quantities

The criteria involve the following geometric quantities (all computed per frame and per pair):

| Geometric quantity | Definition | Used for |
|:-------|:-----|:-----|
| `distance` | Distance between the reference points of two groups. The reference point varies by type: hydrogen bond = donor D and acceptor A; salt bridge = positive/negative charge centers; π-related = aromatic ring center; metal coordination = metal and coordinating atom | All types |
| `angle` (hydrogen bond) | Three-point angle of donor D, hydrogen H, acceptor A (vertex at H), i.e., the D-H···A angle | Hydrogen bond |
| `angle` (π stacking) | Angle between the normal vectors of two aromatic rings | π stacking |
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
| `HBOND_DIST_MAX` | 4.1 | Maximum D-A distance |
| `HBOND_DON_ANGLE_MIN` | 100 | Minimum D-H···A angle |

## π-π Stacking (pi_stacking)

**Criterion**: the distance between the centers of two aromatic rings ∈ (0.5, 5.5] Å, and the stacking type is not N (`pistacking_type != 'N'`, i.e., classified as T-type or P-type); planarity satisfaction is an optional condition (disabled by default).

| Constant | Value | Meaning |
|:-----|:--|:-----|
| `PISTACK_MIN_DIST` | 0.5 | Minimum distance |
| `PISTACK_DIST_MAX` | 5.5 | Maximum ring-center distance |
| `PISTACK_OFFSET_MAX` | 2.0 | Upper limit of ring-center offset |
| `PISTACK_PLANARITY` | 5.0 | Upper limit of ring planarity deviation (°), disabled by default |
| `PISTACK_ANG_DEV` | 30 | Angle deviation for T/P classification |

**Stacking type classification** (`pistacking_type` metric), based on the angle `angle` between the two ring normal vectors and the offset `offset`:
- **P-type (parallel stacking)**: `angle ≤ 30°` and `offset < 2.0 Å`
- **T-type (edge-to-face stacking)**: `angle ≥ 60°` and `offset < 2.0 Å`
- **N (none)**: otherwise

In the output XPM: white = none, pink = T-type, blue = P-type.

## Salt Bridge (salt_bridge)

**Criterion**: the distance between a positive charge center (LYS/ARG, etc.) and a negative charge center (ASP/GLU, etc.) is ≤ 5.5 Å.

| Constant | Value | Meaning |
|:-----|:--|:-----|
| `SALTBRIDGE_DIST_MAX` | 5.5 | Maximum charge-center distance |

## Hydrophobic Interaction (hydrophobic)

**Criterion**: the distance between two hydrophobic atoms ∈ (0.5, 4.0) Å.

| Constant | Value | Meaning |
|:-----|:--|:-----|
| `HYDROPH_MIN_DIST` | 0.5 | Minimum distance |
| `HYDROPH_DIST_MAX` | 4.0 | Maximum distance |

## Halogen Bond (halogen_bond)

**Criterion**: the distance between halogen X and acceptor A is ≤ 4.0 Å, the C-X···A angle = 165°±30°, and the X···A-R angle = 120°±30°.

| Constant | Value | Meaning |
|:-----|:--|:-----|
| `HALOGEN_DIST_MAX` | 4.0 | Maximum X···A distance |
| `HALOGEN_DON_ANGLE` | 165 | Optimal C-X···A angle |
| `HALOGEN_ACC_ANGLE` | 120 | Optimal X···A-R angle |
| `HALOGEN_ANGLE_DEV` | 30 | Upper limit of angle deviation |

## Metal Coordination (metal_coordination)

**Criterion**: the distance between the metal center (e.g., Mg²⁺) and coordinating atoms is < 3.0 Å. Coordinating atoms come from the `metal_binding` group (atoms around the metal that can coordinate).

| Constant | Value | Meaning |
|:-----|:--|:-----|
| `METAL_DIST_MAX` | 3.0 | Maximum metal-coordinating atom distance |

> Note: coordination by pure water molecules is not reported (ignored when all coordinating atoms come from water); the first release only applies the distance criterion and does not match coordination geometries.

## Water Bridge (water_bridge)

**Criterion**: the distance from water oxygen Ow to polar atoms (donor D/acceptor A) is < 4.1 Å, and both angles are within range (theta ≥ 100°, 71° ≤ omega ≤ 140°).

| Constant | Value | Meaning |
|:-----|:--|:-----|
| `WATER_BRIDGE_MAXDIST` | 4.1 | Maximum Ow-to-polar-atom distance |
| `WATER_BRIDGE_THETA_MIN` | 100 | Minimum water O-donor H-donor D angle |
| `WATER_BRIDGE_OMEGA_MIN` | 71 | Minimum acceptor-water O-donor H angle |
| `WATER_BRIDGE_OMEGA_MAX` | 140 | Maximum acceptor-water O-donor H angle |

## π-Cation Interaction (pi_cation)

**Criterion**: the distance from the ring center to the cation (positively charged group) ∈ (0.5, 6.0) Å, with offset < 2.0 Å.

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
| π stacking | distance, angle, offset, pistacking_type | Å, °, Å, type |
| salt bridge | distance | Å |
| hydrophobic | distance | Å |
| halogen bond | distance, don_angle, acc_angle | Å, °, ° |
| metal coordination | distance | Å |
| water bridge | dist_dw, dist_wa, theta, omega | Å, Å, °, ° |
| π-cation interaction | distance, offset | Å, Å |

## Criterion Sources

The criteria and thresholds follow the interaction definition system of PLIP (Protein-Ligand Interaction Profiler): the distance/angle cutoffs (e.g., hydrogen bond D-A 4.1 Å, salt bridge 5.5 Å, halogen bond 165°±30°) are all consistent with PLIP, making it easy to compare with existing PLIP results. The π stacking T/P-type classification standard is consistent with PLIP's `pi-stacking` classification.

## Interpreting Results: Key Points

- **existence**: whether a given pair exists in each frame (boolean matrix after criterion filtering)
- **metrics**: metric values per pair per frame (NaN for inactive frames)
- **occupancy**: number of frames in which a given pair exists / total number of frames, measuring the persistence of that interaction
