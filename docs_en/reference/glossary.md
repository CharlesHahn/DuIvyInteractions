# Glossary

This document summarizes the terms and symbols used in the tool for quick reference.

## General Concepts

| Term | Description |
|:-----|:-----|
| **tpr** | GROMACS binary topology file containing atom types, bonding graph, charges, and force field parameters |
| **xtc** | GROMACS compressed trajectory file containing per-frame atom coordinates |
| **MD** | Molecular dynamics (MD) simulation |
| **Force field** | Parameter set describing atomic interactions (this project supports the Amber family) |
| **GAFF** | General Amber Force Field, a general-purpose force field for organic molecules (commonly used for ligands) |
| **ff14SB** | Amber protein force field (this project's test system uses amber14sb) |

## Data Structures

| Term | Description |
|:-----|:-----|
| **SystemData** | System data parsed from the tpr (atoms, residues, bonds) |
| **Group** | A chemical group that can participate in interactions (e.g., an aromatic ring, a hydrogen bond donor) |
| **pair** | A pair (or more, e.g., a water bridge triplet) of groups participating in an interaction |
| **Interaction** | All detection results of one interaction type (stored in matrix form) |
| **existence** | Boolean matrix `(n_pairs, n_frames)` indicating whether a pair exists in a given frame |
| **metrics** | Geometric metric dictionary, e.g., distance, angle, with shape `(n_pairs, n_frames)` |
| **occupancy** | Occupancy = number of frames in which a pair exists / total number of frames |
| **h5** | HDF5 result file, the output of `dii run` |

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

## Interaction Types

| Type | English | Description |
|:-----|:-----|:-----|
| Hydrogen bond | hydrogen_bond | D-H···A |
| π-π stacking | pi_stacking | Two aromatic rings stacked (T-shaped/parallel) |
| Salt bridge | salt_bridge | Positively and negatively charged pair |
| Hydrophobic interaction | hydrophobic | Contact between nonpolar atoms |
| Halogen bond | halogen_bond | C-X···A |
| Metal coordination | metal_coordination | Metal-coordination atom |
| Water bridge | water_bridge | D-H···Ow···A |
| π-cation interaction | pi_cation | Aromatic ring-cation |

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
