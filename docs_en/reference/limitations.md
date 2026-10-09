# Known Limitations and Version Information

This document describes the tool's version information, supported scope, and known limitations to help users judge its applicability.

## Version Information

| Item | Value |
|:---|:---|
| Software version | v0.0.1 |
| HDF5 result format version | 1.0 |
| Supported Python | >= 3.9 |
| Package name (PyPI) | `duivyinteractions` |
| Command-line entry point | `dii` |

## Supported Scope

- **Force fields (4 families)**:
  - **Amber family**: amber03/94/96/99/99SB/99SB-ildn/GS/14SB proteins + GAFF/GAFF2 ligands (all-atom with explicit hydrogens)
  - **GROMOS 53A6/54A7**: united-atom force field with explicit polar H; paired with the SPC/SPC-E water model (`SOL`); coarse types compensated by structural criteria
  - **CHARMM36/C36m**: proteins + CGenFF ligands; default TIP3/HOH water model
  - **OPLS-AA/L**: includes `opls_XXX` number mapping; HOH/SPC, HO4/TIP4P, HO5/TIP5P water models
- **Topology**: GROMACS tpr (read via MDAnalysis)
- **Trajectory**: formats readable by MDAnalysis (xtc, etc.)
- **Structure**: donor determination relies on explicit H (a D–H bond with q(H)>0); in united-atom force fields (GROMOS) there is no explicit aliphatic H, so donor and hydrophobic-neighborhood determinations are correspondingly restricted (hydrophobicity is compensated by a type whitelist + polar-neighbor exclusion)

## Boundaries

- **GROMOS ligands are not supported**: GROMOS ligands require ATB (Automatic Topology Builder) topology parameterization, which this project does not cover; the identifier is validated for GROMOS protein residues.
- **The per_tuple strategy (strategy one) is in a "may be deprecated" state**: it is a per-candidate-tuple reference implementation whose tests are not run and whose results are not current; the current strategies are `per_frame` and `two_pass`.

## Known Limitations

### Detection Criteria Related

- **PBC not handled**: interactions under periodic boundary conditions (especially water bridges and hydrogen bonds) are not processed periodically; long trajectories or cross-boundary interactions may be inaccurate
- **Metal coordination geometry not matched**: only a distance criterion is used, with no matching of coordination geometries such as linear/trigonal/tetrahedral/octahedral; purely water-coordinated sites are excluded (`metal_binding` groups exclude `WATER_RESIDUES`)
- **Water bridge water deduplication not implemented**: when one water molecule participates in multiple hydrogen bonds, the two best H-O-H angles are not retained
- **Water bridge TwoPass strategy lacks a distance lower bound**: the TwoPass strategy does not apply an Ow-A distance lower bound (2.5 Å), which may cause small differences from the PerTuple/PerFrame strategies
- **Hydrophobic-aromatic not deduplicated**: π stacking and hydrophobic interactions may double-count the same physical contact
- **Cross-strategy differences**: the "same-residue deduplication" of hydrophobic interactions (keeping only the pair with the closest average distance per `(group_id, residue_id)`) is implemented only in the per_frame strategy and is missing in the two_pass strategy (a TODO remains in the detector code); the missing water-bridge lower bound in TwoPass is another such difference
- **Rough nitrogen group classification**: positively charged nitrogen groups are all classified as "tertiary amine" (reproducing the PLIP definition, which actually covers all sp3 N), without distinguishing primary/secondary/tertiary amines; the tertiary amine angle check in π-cation detection may trigger false positives
- **Group identification results not fully reviewed manually**: unit tests are based on code output; manual verification against "ground truth" is still needed (a pending item in `doc/TODO.md`)

### Performance Related

- **Long trajectory memory usage**: the PerFrame strategy pre-allocates a `(n_pairs, n_frames)` matrix for all candidate pairs; 1μs-scale long trajectories (hundreds of thousands of frames) have high memory usage (e.g., a distance matrix for 42,603 candidates × 500,000 frames can reach hundreds of GB); the TwoPass strategy does not have this problem (sparse storage + full metrics computed only for active pairs)

### Functional Scope

- **Visualization**: results are exported as xvg/xpm (can be plotted with DuIvyTools/Xmgrace); built-in plotting in the tool is not implemented
- **H-bond acceptor determination (fix A1, 2026-09-30)**: the cross-force-field issue of "ammonium/guanidinium/H-bearing pyrrole and other non-acceptor N being misjudged as acceptors" is fixed — all four force fields now exclude ordinary amide/ammonium/H-bearing pyrrole/guanidinium N type by type (based on force field types + chemical facts + literature), while retaining Pro N / H-free pyridine-type His N / neutral amines / nucleic-acid amino groups; measured effect (GROMOS system): acceptors 336→206, H-bond pairs 108→86, water bridges 2541→1698. The item-by-item justification is in `doc/acceptor_identification_evidence.md`

## Notes

- Thresholds are consistent with PLIP for easy result comparison, but specific systems may require adjustment as needed (constants are at the top of the detector files; rerun the unit tests of the affected types after changing them)
- Group identification is performed once and is frame-independent; geometric evaluation is performed per frame, and the computational cost grows linearly with the number of frames
- Water residue exclusion is force-field specific (the `WATER_RESIDUES` class attribute: Amber SOL/HOH/WAT, GROMOS SOL, CHARMM TIP3/HOH/SOL/WAT, OPLS HOH/HO4/HO5/SOL/WAT); always choose the `--ff` matching the system when running