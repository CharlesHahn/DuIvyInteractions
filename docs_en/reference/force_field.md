# Force Field Type Mapping

This project determines chemical groups directly from the force field atom types in the GROMACS tpr. This document lists the mapping rules from Amber family (Amber protein + GAFF ligand) types to chemical features, to help research users understand the basis of the determination.

## Mapping Principles

The type → chemical feature mapping comes from two independent, cross-validatable lines of evidence:

1. **GAFF types**: from the official antechamber naming rules (Wang et al., J Comput Chem 2004) — the suffix letters/digits of a type name encode chemical features (e.g., the `a` in `ca` = aromatic)
2. **Amber protein types**: inferred backwards from rtp residue definitions — the known chemical facts of residues (e.g., the TYR benzene ring is an aromatic ring) are compared against the type names (`CA`) used for ring atoms in the rtp file, yielding this type = aromatic carbon

The two lines are mutually independent, their conclusions agree, and they have been verified to have **zero conflicts** across the entire Amber family.

## Aromatic Types (STRONG_AROMATIC)

Atoms of these types are **directly determined to be aromatic by their type names**, and serve as strong evidence for ring detection.

| Category | Types | Notes |
|:---------|:------|:------|
| GAFF aromatic carbon | `ca, cg, ch, cm, cn, cp, cq, c1` | `c` + aromatic suffix |
| GAFF aromatic nitrogen | `na, nb, nh, ni, nj, n1, n2` | `na` = pyrrole-type, `nb` = pyridine-type |
| GAFF aromatic phosphorus | `pb` | — |
| Amber protein aromatic carbon | `CA, CB, CC, CK, CM, C5, C6, C7, C*, CW, CR, CN, CV, CQ` | rtp ring atom types |
| Amber protein aromatic nitrogen | `NA, NB, NC, N*` | — |

## Compatible Types (COMPATIBLE_TYPES)

Non-aromatic types that can participate in conjugation under the "coercion" of **n−1 atoms in the ring being aromatic types**. Used to handle ambiguous types and heterocycles.

| Type | Notes |
|:-----|:------|
| `C, N` | Amber protein ambiguous types (main-chain carbonyl carbon/aromatic ring carbon C; amide nitrogen/aromatic nitrogen N) |
| `os, ss` | GAFF furan oxygen / thiophene sulfur |
| `cc, cd` | GAFF conjugated ring carbons that are not purely aromatic |
| `pc, pd` | GAFF sp2 phosphorus in conjugated rings |

## Acceptor Types (ACCEPTOR_TYPES)

Candidate types for hydrogen bond acceptors (those with lone pairs).

| Category | Types |
|:---------|:------|
| GAFF oxygen | `o, o2, oh, os, oe, o1, ow` |
| GAFF nitrogen | `n, n2, n3, nb, ni, nj, nc, ne, nf, nk` (excluding `na, nh` — pyrrole-type can act as donors) |
| GAFF sulfur | `s, ss, sh, sx, s2` |
| Halogens | `f, cl, br, i` |
| Amber protein oxygen | `O, OH, O2, OS, OW` |
| Amber protein nitrogen | `N, N2, N3, NA, NB, N*, NC` |
| Amber sulfur | `S, SH` |

## Donor Determination

The **bond entries between a donor atom D (N/O/S/F) and its hydrogens H must exist explicitly** (Bond + Constraint merged, because N–H bonds are in the Constraint section), and H must carry a positive charge. This determination does not involve coordinate inference.

## Metal Ions (METAL_IONS)

From PLIP config.py: `Ca, Co, Mg, Mn, Fe, Cu, Zn, Li, Na, K, Rb, Sr, Cs, Ba, Cr, Ni, Ru, Rh, Pd, Ag, Cd, La, W, Os, Ir, Pt, Au, Hg, Ce, Pr, Sm, Eu, Gd, Tb, Yb, Lu, Al, Ga, In, Sb, Tl, Pb`.

## Water Molecules (WATER_RESIDUES)

Residue names `SOL, HOH, WAT`, with oxygen atoms OW and hydrogen atoms HW identified by name.

## Type Name Ambiguities and Handling

- **`C` dual identity**: serves both as the main-chain carbonyl carbon and as the Tyr CZ aromatic ring carbon. The determination requires promotion by "≥4 strong aromatic neighbors within the ring".
- **`N3` vs `n3`**: uppercase `N3` = protein positively charged amino group, lowercase `n3` = GAFF neutral amino group, meanings differ.
- **`CA` has different meanings across force fields**: amber = aromatic carbon, GROMOS = α carbon. Must be distinguished via the feature mapping table.
- **N–H bonds are in the Constraint section**: donor identification must merge both Bond + Constraint bond types.

## Compatibility Scope

Verified force fields: amber03, amber94, amber96, amber99, amber99sb, amber99sb-ildn, amberGS, amber14sb + GAFF/GAFF2 ligands. Type name differences across versions (e.g., CX vs CT) are all handled. If a new version introduces new types, simply add one entry to feature tables such as `STRONG_AROMATIC`/`COMPATIBLE_TYPES`/`ACCEPTOR_TYPES` (the design goal of the feature-space mapping).