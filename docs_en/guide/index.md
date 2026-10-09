# User Guide

This section is for **end users** (computational chemistry and structural biology researchers): how to install, get started quickly, use the `dii` command, and interpret results.

## Overview

DuIvyInteractions is a tool for determining intermolecular interactions based on MD topology force-field parameters. Its core idea is to **read force-field atom types directly from the GROMACS tpr topology to identify chemical groups**, rather than rebuilding chemical structures (bond orders, aromaticity, protonation) from coordinates as PLIP/ProLIF do, so the chemical information already present in the MD topology is never discarded; the results are **deterministic** and **consistent with the force field used**.

Key features of the current version (package version v0.0.1, `pyproject.toml`):

- **Two-stage architecture**: ① group identification — performed once, frame-independent, deterministically identifying groups (donor / acceptor / aromatic ring / charged / hydrophobic / halogen / metal / water) from the tpr atom types + bonding graph + explicit hydrogens + charges; ② geometric determination — frame-by-frame distance/angle/planarity criteria, followed by temporal statistics.
- **4 force-field identifiers** (the `IDENTIFIER_CLASSES` registry in `DuIvyInteractions/group_identifiers/__init__.py`): Amber family (amber03/94/96/99/99sb/99sb-ildn/GS/14sb proteins + GAFF/GAFF2 ligands), GROMOS 53A6/54A7, CHARMM36/C36m (including CGenFF), and OPLS-AA/L. The choices of `dii run --ff` (`amber`/`gromos`/`charmm`/`opls`) correspond one-to-one.
- **8 interaction types**: hydrogen bond, π-π stacking, salt bridge, hydrophobic interaction, halogen bond, metal coordination, water bridge, π-cation (`ALL_INTERACTIONS` in `pipeline.py`).
- **Detection strategies**: `two_pass` (default; two-pass traversal with KDTree pre-screening and sparse storage; best performance) and `per_frame` (frame-by-frame vectorized) are the current strategies; `per_tuple` (per-candidate-tuple traversal over all frames) is in the "possibly to be removed" state and is provided for reference only.
- **Result storage and export**: HDF5 serialization (format version 1.0, `io/h5.py`) + xvg/xpm/CSV export (`dii export`, compatible with the GROMACS/DuIvyTools toolchain).
- **Real test cases (3 force fields)**: the repository bundles 3 sets of real GROMACS data — Amber KRAS-RBD D927 (116,383 atoms, contains GNP/Mg²⁺), GROMOS 53A6 (130-residue protein + 6 ligands), CHARMM36 SMO-BST — with `Tests/unittests/` integration tests covering all 8 interaction types; OPLS-AA group identification is validated by unit tests (no real trajectory data bundled yet).

Reliability notes (code state as of 2026-10-09):

- **H-bond acceptor determination (fix A1)**: across force fields, non-acceptor N types bearing H (ordinary amides, ammonium, H-bearing pyrrole, guanidinium) are excluded, while Pro N, His H-free pyridine N, neutral amines, and nucleic-acid amino N are retained; per-item evidence is in `doc/acceptor_identification_evidence.md`.
- **Cross-force-field water exclusion**: water-residue names are governed uniformly by the `GroupIdentifier.WATER_RESIDUES` class attribute (Amber `{SOL,HOH,WAT}`, GROMOS `{SOL}`, CHARMM `{TIP3,HOH,SOL,WAT}`, OPLS `{HOH,HO4,HO5,SOL,WAT}`), so water-bridge and metal-coordination detection exclude water correctly.
- Known limitations (PBC not handled, GROMOS ligands not supported, etc.) are listed in {doc}`/reference/limitations`.

## Reading Path

- **Installation**: see [Installation and Environment Setup](install).
- **Quick Start**: run the full workflow with real test data, see [Quick Start](quickstart).
- **Command Reference**: all arguments and behaviors of `dii run` / `dii export`, see [Command Reference](command).
- **Interpreting Results**: format and meaning of each xvg/xpm/csv file, see [Interpreting Results](result).

```{toctree}
:maxdepth: 1
:caption: Quick Start

install
quickstart
```

```{toctree}
:maxdepth: 1
:caption: Usage Reference

command
result
```
