# Core Concepts

This document explains how the tool works: how force field atom types are read from the GROMACS tpr topology to identify chemical groups, and how intermolecular interactions are detected based on them. It is intended for computational chemistry and structural biology researchers.

## Method Overview

Intermolecular interaction analysis requires first determining the chemical structure of the molecules — which atoms form aromatic rings, which atoms are hydrogen bond donors/acceptors, which groups carry charges, and so on. Based on these chemical groups, the tool then performs geometric evaluation using trajectory coordinates.

There are two approaches to determining the chemical structure:

1. **Rebuild from coordinates**: Starting from atomic coordinates, infer bond orders, aromaticity, and hydrogen atom positions. Commonly used tools (e.g., PLIP with OpenBabel, ProLIF with RDKit) follow this approach.
2. **Read from force field parameters** (this tool): Read force field atom types directly from the GROMACS topology (tpr). These types were determined by antechamber/sobtop according to the molecular electronic structure during force field parameterization.

This tool adopts the second approach. The atom types in the tpr are "records of chemical judgment": force field developers have encoded electronic structure and chemical features into the type names (and residue definitions). Reading the types therefore yields chemical information directly, without rebuilding chemistry or discarding the chemical semantics present in the MD topology.

## Chemical Semantics of Force Field Atom Types

Force field atom types are not arbitrary labels; their naming encodes chemical features. Taking GAFF (General Amber Force Field) as an example, the suffixes of type names carry systematic meaning:

| Type | Chemical meaning | Naming basis |
|:-----|:---------|:---------|
| `c3` | sp3 carbon | c + number 3 |
| `c2` | sp2 olefinic carbon | c + number 2 |
| `ca` | aromatic carbon | c + a (aromatic) |
| `c` | sp2 carbonyl carbon | c without suffix |
| `na` | pyrrole-type aromatic nitrogen | n + a |
| `nb` | pyridine-type aromatic nitrogen | n + b |
| `os` | ether oxygen (single bond) | o + s (single) |

Amber protein types are similar: the type names used for ring atoms in residue definitions (rtp files) (e.g., `CA`) can be cross-referenced against known residue chemistry to deduce their meaning. The same holds for the CHARMM/CGenFF, OPLS-AA, and GROMOS families — for the former two, type names (`CG2R61`, `opls_145`, etc.) originate from the official rtf/parameter-table comments, so their semantics can be verified type by type (see [Force Field Type Mapping](force_field.md)).

Type names encode features such as hybridization state, aromaticity, polarity, and whether they carry hydrogens. Therefore, reading the types is sufficient to determine chemical groups, with no need to infer from coordinates. This guarantees: group identification depends only on topology, not on coordinates, and is consistent across frames and reproducible.

**Note**: type names are **not equivalent across force fields** (e.g., `CA` is an aromatic carbon in Amber but the α carbon in GROMOS); type names cannot be reused across force fields. The tool unifies semantics through each force field's **independent feature tables** (type → {hybridization, aromaticity, polarity, has H, lone pair}).

## Two-Stage Architecture

The analysis is divided into two stages:

```
① Group identification
   tpr force field parameters (atom types + bonding graph + explicit H + charges)
   → identify chemical groups (aromatic rings, H bond donors/acceptors,
     charged groups, hydrophobic atoms, metals, etc.)
   —— performed once, frame-independent

② Geometric evaluation
   with groups labeled, evaluate per frame with distance/angle/planarity criteria
   → per-frame interaction list (8 types) → temporal statistics
```

- **Group identification** depends only on topology, not on coordinates, so the results are frame-independent and reproducible
- **Geometric evaluation** computes geometric metrics (distance, angle, etc.) frame by frame for each pair, and determines whether an interaction exists according to the criteria (see [Interaction Criteria](criteria.md))

Group identification is the difficult stage: it must map tpr atom types to chemical groups correctly. The tool maintains per-force-field feature tables and the necessary structural criteria, and identifies groups once, to be reused for all frames.

## Explicit Hydrogens and Chemical Information

The following chemical information is explicitly present in the tpr topology and read directly from it, with no coordinate inference:

- **Hydrogen bond donor**: the bond entry between the donor atom D and hydrogen H exists (including Constraint-section bonds), and H carries a positive charge
- **Water molecules**: identified by residue names (`WATER_RESIDUES`, per-force-field sets), with oxygen atom OW and hydrogen atoms HW explicit
- **Metal center**: identified by element symbol (e.g., Mg); the reader infers the element from the atomic number

## Supported Scope

- **Force fields (4 families)**:
  - Amber family (amber03/94/96/99/99sb/99sb-ildn/GS/14sb) proteins + GAFF/GAFF2 ligands, with type mapping verified to have zero conflicts
  - GROMOS 53A6 / 54A7 proteins: supported, but as a **united-atom force field** — polar H (N-H, O-H) and aromatic C-H are explicit, while aliphatic nonpolar H are merged into CH1/CH2/CH3; the coarse-type parts are compensated by **structural criteria / tiered determination / anti-enumeration** (acceptors, aromatic rings, hydrophobicity); determinations relying on explicit aliphatic H (e.g., donor-neighborhood details) are restricted; GROMOS ligands (ATB parameterized) are not within the first-release support commitment
  - CHARMM36 / C36m proteins + CGenFF ligands
  - OPLS-AA/L (2001) proteins + ligands within the type tables
- **Missing bond orders**: even if the tpr does not preserve bond orders (all func=1), aromaticity determination is still based on type names; bond orders serve only as cross-validation, so aromatic determination is reliable
- **No bond-order/chemistry rebuild**: chemistry is not rebuilt with OpenBabel/RDKit, and no hydrogen addition or bond-order inference is performed

## Detailed Rules

- Complete mapping of type → chemical features: see [Force Field Type Mapping](force_field.md)
- Identification rules for each group: see [Group Identification Rules](group_rules.md)
- Geometric criteria for each interaction: see [Interaction Criteria](criteria.md)

## References

- GAFF force field: Wang J, Wolf RM, Caldwell JW, Kollman PA, Case DA. Development and testing of a general amber force field. *J Comput Chem*. 2004;25(9):1157-1174.
- Amber force field (ff14SB): Maier JA, Martinez C, Kasavajhala K, Wickstrom L, Hauser KE, Simmerling C. ff14SB: Improving the accuracy of protein side chain and backbone parameters from ff99SB. *J Chem Theory Comput*. 2015;11(8):3696-3713.
- PLIP interaction definitions: Salentin S, Schreiber S, Haupt VJ, Adasme MF, Schroeder M. PLIP: fully automated protein-ligand interaction profiler. *Nucleic Acids Res*. 2015;43(W1):W443-W447.
- MDAnalysis: Michaud-Agrawal N, Denning EJ, Woolf TB, Beckstein O. MDAnalysis: A toolkit for the analysis of molecular dynamics simulations. *J Comput Chem*. 2011;32(10):2319-2327.