# Glossary

This document summarizes the terms and symbols used in the tool for quick reference. The group-type and bond-type enumerations are consistent with `GROUP_TYPES` / `BOND_TYPES` in `DuIvyInteractions/core/constants.py`.

## Force Fields and General Concepts

| Term | Description |
|:-----|:-----|
| **tpr** | GROMACS binary topology file containing atom types, bonding graph, charges, and force field parameters |
| **xtc** | GROMACS compressed trajectory file containing per-frame atom coordinates |
| **MD** | Molecular dynamics (MD) simulation |
| **Force field** | Parameter set describing atomic interactions. This project supports **4 force field families**: the Amber family, GROMOS 53A6/54A7, CHARMM36/C36m, and OPLS-AA/L, managed centrally by the `IDENTIFIER_CLASSES` registry |
| **Amber family** | amber03/94/96/99/99SB/99SB-ildn/GS/14SB protein force fields + GAFF/GAFF2 ligand force fields; all-atom with complete explicit hydrogens |
| **GAFF / GAFF2** | General Amber Force Field, the general-purpose organic-molecule force field of the Amber family (commonly used for ligands) |
| **GROMOS 53A6/54A7** | GROMOS united-atom force field: hydrogens on aliphatic C are merged into heavy atoms, but **polar hydrogens (H on N/O/S) are explicit**. Paired with the SPC/SPC-E water model (residue name `SOL`). Ligands are not supported (they require ATB parameterization). Types are coarse (`N`/`C` cover many chemical environments) and are compensated by structural criteria |
| **CHARMM36 / C36m** | CHARMM protein force field; ligands use CGenFF (CHARMM General Force Field). Default TIP3 water model (residue names `TIP3`/`HOH`) |
| **OPLS-AA/L** | OPLS all-atom force field (2001 version, GROMACS `oplsaa.ff`). Water models HOH/SPC, HO4/TIP4P, HO5/TIP5P (residue names `HOH`/`HO4`/`HO5`) |
| **United atom** | A force-field style that merges hydrogens on aliphatic C into heavy atoms (e.g., GROMOS). Such force fields have no explicit aliphatic H, so determinations depending on explicit H (donors, hydrophobic neighborhood) are restricted |
| **Explicit H** | Hydrogen atoms modeled explicitly in the force field. H-bond donor determination relies on the D–H bond and q(H)>0, hence explicit H is required |
| **Water model** | Force-field parameters for water: SPC/SPC-E (GROMOS), TIP3 (common in CHARMM/Amber), TIP4P/TIP5P (OPLS HO4/HO5). Residue names: `SOL`/`HOH`/`WAT`/`TIP3`/`HO4`/`HO5`, collected in each force field's `WATER_RESIDUES` class attribute |

## Data Structures

| Term | Description |
|:-----|:-----|
| **SystemData** | System data parsed from the tpr (atoms, residues, bonds, inter-residue bonds) |
| **Group** | A chemical group that can participate in interactions (e.g., an aromatic ring, a hydrogen bond donor) |
| **pair** | A pair (or more, e.g., a water bridge triplet) of groups participating in an interaction |
| **Interaction** | All detection results of one interaction type (stored in matrix form) |
| **InteractionSparse** | Sparse intermediate results of the TwoPass strategy's Pass1 (keyed by `(group_id, ...)` tuples) |
| **existence** | Boolean matrix `(n_pairs, n_frames)` indicating whether a pair exists in a given frame |
| **metrics** | Geometric metric dictionary, e.g., distance, angle, with shape `(n_pairs, n_frames)` |
| **occupancy** | Occupancy = number of frames in which a pair exists / total number of frames |
| **h5** | HDF5 result file, the output of `dii run` (format version 1.0) |

## Group Types (GROUP_TYPES)

| Group Type | Meaning | Participating Interactions |
|:---------|:-----|:------------|
| `H_donor` | Hydrogen bond donor (D-H) | Hydrogen bond, water bridge |
| `H_acceptor` | Hydrogen bond acceptor (has lone pair) | Hydrogen bond, water bridge |
| `aromatic_ring` | Aromatic ring | π stacking, π-cation interaction, halogen bond |
| `charged_positive` | Positively charged group | Salt bridge, π-cation interaction |
| `charged_negative` | Negatively charged group | Salt bridge |
| `halogen_donor` | Halogen bond donor (C-X) | Halogen bond |
| `halogen_acceptor` | Halogen bond acceptor | Halogen bond |
| `metal` | Metal center | Metal coordination |
| `metal_binding` | Metal coordination atom | Metal coordination |
| `water` | Water molecule | Water bridge |
| `hydrophobic` | Hydrophobic atom | Hydrophobic interaction |

## Bond Types (BOND_TYPES)

| Bond Type | Meaning |
|:----------|:--------|
| `bond` | Unknown bond order (default) |
| `single` / `double` / `triple` / `aromatic` | Chemical bonds |
| `constrained` / `settle` / `virtual` | Constrained bonds (e.g., SHAKE/LINCS), SETTLE water bonds, virtual-atom bonds |

> Note: N–H donor bonds appear in the `Constraint:` section of the tpr, not in the `Bond:` section; donor identification must merge both. The `constrained` entry of `BOND_TYPES` corresponds to such bonds.

## Interaction Types

| Type | English | Description |
|:-----|:-----|:-----|
| Hydrogen bond | hydrogen_bond | D-H···A |
| π-π stacking | pi_stacking | Two aromatic rings stacked (T-shaped/parallel); the metric `pistacking_type` (P/T/N) classifies the geometry |
| Salt bridge | salt_bridge | Positively and negatively charged pair |
| Hydrophobic interaction | hydrophobic | Contact between nonpolar atoms |
| Halogen bond | halogen_bond | C-X···A |
| Metal coordination | metal_coordination | Metal-coordination atom |
| Water bridge | water_bridge | D-H···Ow···A |
| π-cation interaction | pi_cation | Aromatic ring-cation |

## Detection Strategies

| Strategy | Description |
|:---------|:------------|
| `two_pass` | Two passes: Pass1 discovers active pairs frame by frame (KDTree pre-screening, sparse storage) → Pass2 completes full-frame metrics. Default strategy, best performance |
| `per_frame` | Per-frame vectorized processing of all candidate pairs. Suitable when the number of candidate pairs is very large (e.g., water bridges); long trajectories have high memory usage due to pre-allocated matrices |
| `per_tuple` | Per-candidate-tuple traversal over all frames, vectorized. A reference implementation, currently in a "may be deprecated" state (tests do not run it) |

## Symbols

| Symbol | Meaning |
|:-----|:-----|
| Å | Angstrom, unit of length (1 Å = 0.1 nm) |
| ° | Degree, unit of angle |
| D | Hydrogen bond donor atom (N/O/S/F) |
| A | Hydrogen bond acceptor atom |
| H | Hydrogen atom |
| Ow | Oxygen atom of a water molecule |
| X | Halogen atom (F/Cl/Br/I) |
| q(H) | Partial charge of the hydrogen atom |
| P-type / T-type | Parallel / edge-to-face (T-shaped) π stacking |