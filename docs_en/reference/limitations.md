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

- **Force field**: Amber family (amber03/94/96/99/99sb/99sb-ildn/GS/14sb) proteins + GAFF/GAFF2 ligands
- **Topology**: GROMACS tpr (read via MDAnalysis)
- **Trajectory**: formats readable by MDAnalysis (xtc, etc.)
- **Structure**: all-atom force fields with explicit hydrogens (donor identification relies on explicit H)

## Known Limitations

### Detection Criteria Related

- **PBC not handled**: interactions under periodic boundary conditions (especially water bridges and hydrogen bonds) are not processed periodically; long trajectories or cross-boundary interactions may be inaccurate
- **Metal coordination geometry not matched**: only a distance criterion is used, with no matching of coordination geometries such as linear/trigonal/octahedral; purely water-coordinated sites are excluded
- **Water bridge water deduplication not implemented**: when one water molecule participates in multiple hydrogen bonds, the two best H-O-H angles are not retained
- **Water bridge TwoPass strategy lacks a distance lower bound**: the TwoPass strategy does not apply an Ow-A distance lower bound (2.5 Å), which may cause small differences from the PerTuple/PerFrame strategies
- **Hydrophobic-aromatic not deduplicated**: π stacking and hydrophobic interactions may double-count the same physical contact
- **Rough nitrogen group classification**: positively charged nitrogen groups are all classified as "tertiary amine", without distinguishing primary/secondary/tertiary amines; the tertiary amine angle check in π-cation detection may trigger false positives
- **Group identification results not fully reviewed manually**: unit tests are based on code output; manual verification against "ground truth" is still needed

### Performance Related

- **Long trajectory memory usage**: the PerFrame strategy pre-allocates a `(n_pairs, n_frames)` matrix for all candidate pairs; 1μs-scale long trajectories (hundreds of thousands of frames) have high memory usage; the TwoPass strategy does not have this problem

### Functional Scope

- **Visualization**: results are exported as xvg/xpm (can be plotted with DuIvyTools/Xmgrace); built-in plotting in the tool is not implemented
- **United-atom force fields**: force fields without explicit hydrogens such as GROMOS are not supported (donor identification relies on explicit H)

## Notes

- Thresholds are consistent with PLIP for easy result comparison, but specific systems may require adjustment as needed (constants are at the top of the detector files)
- Group identification is performed once and is frame-independent; geometric evaluation is performed per frame, and the computational cost grows linearly with the number of frames
