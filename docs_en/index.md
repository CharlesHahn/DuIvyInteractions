# DuIvyInteractions Documentation

An intermolecular interaction determination tool based on MD topology force field parameters.

**Core idea**: Read force-field atom types directly from the GROMACS tpr topology to identify chemical groups — group identification is based on the force-field parameters themselves, rather than rebuilding chemical structures from coordinates. This approach reuses the force-field type semantics encoded in the MD topology (e.g. GAFF `ca`, CHARMM `CG2R61`, OPLS `CA`): the outcome is **deterministic** (the same input always yields the same group set) and **self-consistent with the force field** (chemical information comes from the force field used in the simulation, not from third-party inference). Compared with the coordinate-based reconstruction of bond orders, aromaticity and hydrogens used by PLIP (OpenBabel) and ProLIF (RDKit), no MD topology chemistry is discarded and no reconstruction uncertainty is introduced.

## Feature Overview

**Multi-force-field support (4 force fields)**: Group identification covers the Amber family (amber03/94/96/99/99sb/99sb-ildn/GS/14sb + GAFF), GROMOS 53A6/54A7, CHARMM36/C36m (including CGenFF), and OPLS-AA/L. Identifiers are managed centrally by the `IDENTIFIER_CLASSES` registry (`DuIvyInteractions/group_identifiers/__init__.py`); `dii run --ff` accepts all four force fields. Water residue names (`WATER_RESIDUES`) are force-field specific (Amber SOL/HOH/WAT, GROMOS SOL, CHARMM TIP3, OPLS HO4/HO5) so that water bridges and metal-coordination detection exclude water correctly.

**8 interaction types**: hydrogen bonds, π–π stacking, salt bridges, hydrophobic contacts, halogen bonds, metal coordination, water bridges, and π–cation interactions. Each type provides per_frame and two_pass as the two current detection strategies (the per_tuple strategy is in a pending state); results are uniformly returned as `List[Interaction]` and stored in HDF5.

**Two-stage architecture**:

1. **Group identification** (frame-independent, performed once): tpr atom types + bonding graph + explicit hydrogens + charges → deterministic identification of groups (donors/acceptors/aromatic rings/charged/hydrophobic/halogen/metal/water). H-bond acceptor determination was revised in fix A1 (2026-09-30), which removed H-bearing non-acceptor N types; the evidence list is in `doc/acceptor_identification_evidence.md`.
2. **Geometric determination** (per frame): distance/angle/planarity criteria → per-frame interaction lists → temporal statistics.

**Real-data test cases (3 force fields)**: Amber KRAS–RBD D927 (with GNP/Mg metal centre), GROMOS 53A6 (130-residue protein + 6 ligands ZIN1–6), and CHARMM36 SMO–BST (Smoothened–β-sitosterol, Mendeley v94vzbwzf3, two tpr variants — a water-stripped complex and a full-atom system with 63,055 SOL waters — and 64 tests) are covered by integration tests (`Tests/test_MD_case_*` and `Tests/unittests/`); all eight interaction types are validated on real data (halogen bonds and metal coordination are asserted as zero-result in the SMO–BST system, which contains no halogen or metal centre).

## Quick Start

```bash
dii run -t md.tpr -f md.xtc -o out/ --ff amber      # detection → h5 (--ff: amber/gromos/charmm/opls)
dii export -i out/hydrogen_bond.h5 -o out_export/  # export xvg/xpm/csv + overview
```

```{toctree}
:maxdepth: 1
:caption: User Guide

guide/index
```

```{toctree}
:maxdepth: 1
:caption: Reference

reference/index
```

```{toctree}
:maxdepth: 1
:caption: Other

changelog
```