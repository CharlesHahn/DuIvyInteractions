# Group Identification Rules

In the group identification stage, chemical groups capable of participating in interactions are identified from the atom types, bonding graph, charges, and explicit hydrogen atoms of the tpr. This document lists the determination rules and feasibility basis for each group, for understanding tool behavior and validating results.

## Complete Set of Group Types (GROUP_TYPES)

Group types the tool can identify (core/constants.py):

| Group Type | Participating Interactions |
|:-----------|:---------------------------|
| `H_donor` | hydrogen bond, water bridge |
| `H_acceptor` | hydrogen bond, water bridge |
| `aromatic_ring` | π-π stacking, π-cation interaction, halogen bond (π acceptor) |
| `charged_positive` | salt bridge, π-cation interaction |
| `charged_negative` | salt bridge |
| `halogen_donor` | halogen bond |
| `halogen_acceptor` | halogen bond |
| `metal` | metal coordination |
| `metal_binding` | metal coordination |
| `water` | water bridge |
| `hydrophobic` | hydrophobic interaction |

## Per-Group Feasibility Matrix

| Group | Required Chemical Features | Evidence in tpr | Certainty |
|:------|:---------------------------|:----------------|:----------|
| **Aromatic ring** (π-π/π-cation/halogen-π) | aromatic ring | types `ca/na/nb/cp/cg` → graph-theoretic ring detection (edge deletion + BFS + fused ring deduplication) → all atoms in ring aromatic; optional cross-validation with bond order func=4/5 and planarity | high |
| **H bond donor** | D–H (D=N/O/S/F) | D–H bond entries exist explicitly (all-atom, Bond + Constraint merged) + q(H)>0 | no coordinate inference |
| **H bond acceptor** | lone pair available | type rule table (`o/oh/os/nb/n/f`…) + H count (`nb` with no H is necessarily an acceptor) + negative charge validation | high |
| **Salt bridge** | formal charge pair | protein residue name dictionary (LYS/ARG/ASP/GLU/HIP encode protonation state); ligands by type + H count | high for protein, medium for ligand |
| **Hydrophobic** | nonpolar C/S/X | type set (`c3/c2`…) + no polar substituents | high |
| **Halogen bond** | σ-hole halogen | halogen types (`f/cl/br/i`) + neighboring carbon type distinguishes aromatic/alkyl halides | high |
| **Metal coordination** | metal center | element/type + charge (Mg²⁺ in this project); uses atomic number (not dependent on type names) | no coordinate inference |
| **Water bridge** | water molecules | residue names SOL/HOH + OW/HW atom names | no coordinate inference |

## Key Points of Identification Methods

### Aromatic Ring Detection

1. Determine strong aromatic atoms from type names (`STRONG_AROMATIC`)
2. Graph-theoretic ring detection (edge deletion + BFS) to find all rings
3. Filter: at least n−1 atoms in the ring are aromatic (the rest may be `COMPATIBLE_TYPES` compatible types, participating in conjugation under the "coercion" of strong aromatic neighbors)
4. Fused ring redundancy removal (`_deduplicate_aromatic_rings`)
5. Optional: cross-validation with bond order func=4/5 and planarity

Ring atoms are stored in ring order (BFS path order) for subsequent normal vector computation.

On the D927 validation system: 3 aromatic rings + 1 non-aromatic sulfur-containing ring were identified (a 2,3-dihydrothiophene-type fused heterocycle, with the C22=C23 double bond being a non-aromatic thiophene); the determination of RBD's 27 rings (Pro×9 + Tyr×11 + His×3 + Trp×2 + Phe×1) agrees with chemical facts.

### Donor Determination

Criterion = D–H bond (Bond + Constraint merged) + q(H)>0. Because N–H bonds are in the Constraint section of the tpr (not the Bond section), **the two must be merged**. On the D927 validation system: the charges of RBD's 263 donor H atoms are all +0.19~+0.45, with no q(H)≤0.

### Protein Charged Groups

Protonation states are encoded by the residue name dictionary: LYS (+), ARG (+), ASP (−), GLU (−), HIP (positively charged histidine state). Ligand charged groups are judged by type + H count. After identification, validation is performed by net charge: positively charged groups must have net charge > +0.1, negatively charged groups must have net charge < −0.1 (`CHARGE_THRESHOLD`).

### Hydrophobic Atom Determination

A hydrophobic atom = a carbon atom whose all bonded neighbors are C or H (i.e., no polar substituents).

### Metal Coordination

Metal centers are identified by element (the `METAL_IONS` set), not dependent on type names (using atomic number). Coordination atoms are the atoms around the metal that can coordinate (`metal_binding`).

## Type Name Ambiguities and Handling

- **`C` dual identity**: main-chain carbonyl carbon vs Tyr CZ aromatic ring carbon. Requires promotion to aromatic by "≥4 strong aromatic neighbors within the ring".
- **`N3` vs `n3`**: protein positively charged amino group vs GAFF neutral amino group; the case differs in meaning.
- **`CA` has different meanings across force fields**: amber = aromatic carbon, GROMOS = α carbon. Must go through the feature mapping table.
- **N–H bonds are in the Constraint section**: donor identification must merge Bond + Constraint.

## Supported Boundaries

- All-atom force fields (Amber family) rely on explicit hydrogen atoms; for united-atom force fields (GROMOS has no explicit H) donor identification fails and is declared unsupported
- Robust to missing bond orders: even if the tpr does not retain bond orders, type names already encode aromaticity, so the determination is reliable