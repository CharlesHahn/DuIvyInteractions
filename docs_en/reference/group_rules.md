# Group Identification Rules

In the group identification stage (stage ① of the two-stage architecture, see [Core Concepts](concepts.md)), chemical groups capable of participating in interactions are identified from the atom types, bonding graph, charges, and explicit hydrogen atoms of the tpr, **once and frame-independently**. This document lists the determination rules for each group and the implementation differences across the four force field identifiers (Amber / GROMOS / CHARMM / OPLS), for understanding tool behavior and validating results. The complete per-force-field type feature tables are given in [Force Field Type Mapping](force_field.md).

## Complete Set of Group Types (GROUP_TYPES)

Group types the tool can identify (`GROUP_TYPES` in `core/constants.py`):

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

## Determination Basis at a Glance

| Group | Evidence in tpr | Certainty |
|:------|:----------------|:----------|
| Aromatic ring | types ∈ `STRONG_AROMATIC`/compatible → graph-theoretic ring detection (BFS) → two-condition filter (GROMOS: tiered determination) | high |
| H bond donor | D–H bond entries exist explicitly (Bond/Constraint/Settle merged) + q(H)>0 | no coordinate inference |
| H bond acceptor | type table + structural criterion + q<0 (positive-charge members / in-ring H-bearing N excluded) | high |
| Charged group | residue-name dictionary / functional-group pattern + net-charge validation (threshold ±0.1) | high for proteins |
| Hydrophobic | C with all neighbors ∈ {C,H} (GROMOS: type whitelist + no polar neighbors) | high |
| Halogen bond | halogen element/type + neighboring carbon (distinguishes aromatic/alkyl halides) | high |
| Metal | element ∈ `METAL_IONS` (inferred by reader) | 100% |
| Water | residue name ∈ `WATER_RESIDUES` + OW/HW atom names | 100% |

## Identification Flow

`identify()` identifies residues one by one in a fixed order (base class `AmberFFGroupIdentifier._identify_residue`; GROMOS/CHARMM/OPLS are subclasses that override the force-field-dependent parts), and D–H pairs across residue boundaries are checked separately:

1. Aromatic rings (`_find_aromatic_rings`) → 2. Donors (`_find_donors`) → 3. Charged groups (`_find_charged`, identified early so that positively charged members can be excluded from acceptors) → 4. Acceptors (`_find_acceptors`, with an exclusion set) → 5. Halogen donors → 6. Halogen acceptors → 7. Metals → 8. Water → 9. Hydrophobic atoms → 10. Metal-binding atoms

**Acceptor exclusion set** (built in `_identify_residue` for acceptor determination): ① positively charged group members (no lone pair, cannot be acceptors); ② in-ring H-bearing nitrogens (pyrrole-type NH, lone pair in aromatic conjugation, cannot be acceptors).

## Per-Group Determination Rules

### Aromatic Ring (aromatic_ring)

**Purpose**: provide ring geometry for π-π / π-cation / halogen bond (π acceptor).

- **Ring detection**: graph-theoretic; BFS from each edge finds all rings (no size limit, `_detect_rings`), then rings are sorted by size and deduplicated — larger rings fully covered by an accepted smaller ring are removed (`_deduplicate_aromatic_rings`)
- **Aromatic filtering (Amber / CHARMM / OPLS, two conditions)**: **≥ n−1 atoms** in the ring belong to `STRONG_AROMATIC`, and all remaining atoms must belong to `COMPATIBLE_TYPES`
- **GROMOS tiered determination** (coarse types, `_filter_aromatic_rings` override): ring composition must be ⊆ `{C, CR1, NR}`; containing `NR` (aromatic nitrogen) is directly aromatic; all-C rings fall back to the residue-name whitelist (PHE/TYR/TRP/His variants, see [Force Field Type Mapping](force_field.md))
- Ring atoms are stored in ring order (BFS path order); the π-π detector computes ring normals from cross products of adjacent atoms

**Note**: aromaticity uses **only type names and the bonding graph, not bond orders** (reliable even if the tpr does not retain bond orders); planarity is not a group-identification criterion but an optional geometric check of the π-π detector ([Interaction Criteria](criteria.md)).

### H Bond Donor (H_donor)

**Criterion** = D–H bond entries exist (D ∈ {N, O, S, F}) + q(H) > 0 (`_classify_dh_pair`). Because N–H bonds reside in the **Constraint section** of the tpr, the reader merges Bond / Constraint / Settle into a unified bond list, and **the merged bond graph must be used** (scanning only the Bond section misses N–H). D–H pairs across residue boundaries are checked separately (`_find_inter_residue_donors`). No coordinate inference is involved.

### H Bond Acceptor (H_acceptor)

**Criterion** = type ∈ acceptor type table + **partial charge q < 0** + not in the exclusion set (positive-charge members, in-ring H-bearing N).

**No "H counting" is used** (it was never implemented); the exclusion of H-bearing environments is done by each force field's type table or structural criteria:

| Force field | Acceptor determination implementation |
|:------------|:--------------------------------------|
| Amber | type table `ACCEPTOR_TYPES` + structural criterion for the ambiguous type `N` (H-bearing ordinary amide → excluded; H-free Pro N → kept; `N2`/`NB`/`NC` kept directly by the type table) + exclusion set |
| GROMOS | type table `GROMOS_ACCEPTOR_TYPES` + N structural criterion (`_find_acceptors` override, for element N: **≥ 4 bonds (ammonium RNH₃⁺) → excluded**; **has an H neighbor (the union of ordinary amide / side-chain amide / H-bearing pyrrole / guanidinium) → excluded**; **H-free (Pro N / H-free pyridine-type His N) → kept**) |
| CHARMM | type table `CHARMM_ACCEPTOR_TYPES` already carries a per-type judgment (protein: NH1/NH2/NH3/NC2/NY/NR1/NR3 excluded, N/NR2 kept; CGenFF: 10 excluded / 6 kept), no ambiguous N types |
| OPLS | type table `OPLS_ACCEPTOR_TYPES` already carries a per-type judgment (nitrogen keeps only opls_239/511/900, 10 excluded) |

Acceptor qualification (nine N environments + literature) is detailed in the "Acceptor Qualification" section of [Force Field Type Mapping](force_field.md) and in `doc/acceptor_identification_evidence.md`.

### Charged Groups (charged_positive / charged_negative)

Three layers (base class `_find_charged`; GROMOS/CHARMM/OPLS replace dictionaries or override validation):

1. **Layer 1: residue-name dictionary** — produces charged groups directly from residue name + atom-name lists (`_identify_protein_charged`)

   | Force field | Positive dictionary | Negative dictionary |
   |:------------|:--------------------|:--------------------|
   | Amber | ARG, LYS, HIP, ORN, DAB, M3L, MLY | ASP, GLU, CYM, KCX, PCA, SEP, TPO, PTR |
   | GROMOS | ARG, LYSH, HISH | ASP, GLU |
   | CHARMM | LYS, ARG, HSP (includes imidazole carbons CG/CE1/CD2) | ASP, GLU, CYM |
   | OPLS | ARG, LYSH, HISH (includes imidazole carbons CG/CE1/CD2) | ASP, GLU |

2. **Layer 2: functional-group pattern matching** (following PLIP `is_functional_group`, force-field independent): positive — quaternary ammonium (N with 4 bonds, no H), tertiary amine (N with ≥3 bonds), guanidine (C bonded to 3 N, at least one N bonded only to that C), sulfonium (S with 3 bonds, no H); negative — phosphate (P with all O neighbors), sulfonic acid (S with 3 O neighbors), sulfate (S with 4 O neighbors), carboxylate (C with 2 O + exactly 1 C neighbor)

3. **Layer 3: net-charge validation** — a positive group must have net charge > +0.1 and a negative group < −0.1 (`CHARGE_THRESHOLD = 0.1`)

**CHARMM / OPLS specifics**: ammonium nitrogens (N-terminal NH₃⁺, LYSH side chain) carry negative partial charges, so single-atom validation fails; both force fields override the functional-group layer and, only for nitrogens satisfying the **terminal-ammonium structural criterion** (`_is_terminal_ammonium`: N neighbors ⊆ {C,H} and heavy neighbors not bonded to ≥2 N), validate with charge atoms expanded to N + bonded H.

**GROMOS specifics**: the positively charged N atoms of its protonated residues themselves carry positive partial charges, so the functional-group layer would additionally emit a per-N positive group that is a subset of the dictionary-layer complete group; `_deduplicate_charged` is overridden to remove groups strictly contained in a more complete group of the same type (Amber has no such issue).

### Halogen Donor / Acceptor (halogen_donor / halogen_acceptor)

- **Donor**: a halogen (F/Cl/Br/I) bonded to a carbon atom → group atoms = [C, X] (first the carbon, then the halogen)
- **Acceptor**: a central atom A ∈ {C, P, S} with at least one neighbor R ∈ {O, P, N, S} → group atoms = [A, R₁, R₂, …] (geometric evaluation uses the X···A-R angle)

### Hydrophobic Atoms (hydrophobic)

| Force field | Criterion |
|:------------|:----------|
| Amber / CHARMM / OPLS | carbon atom whose **all bonded neighbors ∈ {C, H}** (no polar substituent) |
| GROMOS (anti-enumeration) | type ∈ `GROMOS_HYDROPHOBIC_TYPES` (`C, CH0, CH1, CH2, CH3, CH4, CH2r`, CH3p/CR1 excluded) and **no neighbor of element O/N/S** — because united-atom aliphatic C has no H neighbors |

### Metal Coordination (metal / metal_binding)

- **Metal center**: element symbol ∈ `METAL_IONS` (42 elements, inferred by the reader, independent of type names)
- **Coordinating atoms**: atoms of element ∈ {O, N, S} (`_find_metal_binding`); **water residues are excluded entirely** (returns empty for residue names in `self.WATER_RESIDUES`) — so "coordination by pure water molecules" is not reported, and exclusion is correct per force field via `WATER_RESIDUES` (CHARMM TIP3 / OPLS HO4/HO5)

### Water (water)

Residue name ∈ `WATER_RESIDUES` (base-class class attribute; per-force-field sets in [Force Field Type Mapping](force_field.md)); group atoms = all atoms of the water residue (OW/HW identified by atom names).

## Type Name Ambiguities and Handling

- **`C` dual identity** (Amber): main-chain carbonyl carbon vs Tyr CZ aromatic ring carbon; `C` participates in conjugation as a compatible type when the ring has ≥ n−1 strong aromatic atoms (5 of 6 for a six-membered ring).
- **`N3` vs `n3`**: protein positively charged amino group (ammonium) vs GAFF neutral amine — the case differs in meaning.
- **`CA` has different meanings across force fields**: Amber = aromatic carbon, GROMOS = α carbon; must go through each force field's independent feature tables.
- **N–H bonds are in the Constraint section**: donor identification must rely on the merged Bond + Constraint (+ Settle) bond graph.

## Validation Baseline

On the Amber real test system (KRAS-RBD D927, `Tests/test_MD_case_amber`):

- **Aromatic rings**: the ligand yields 3 aromatic rings + 1 non-aromatic sulfur-containing ring (a 2,3-dihydrothiophene-type fused heterocycle whose C22=C23 double bond is a non-aromatic thiophene); the determination of RBD's 27 rings (Pro×9 + Tyr×11 + His×3 + Trp×2 + Phe×1) agrees with chemical facts
- **Donors**: the charges of RBD's 263 donor H atoms are all +0.19~+0.45, with no q(H) ≤ 0

Real integration tests are in `Tests/unittests/test_*_real_identifier.py` and `test_<ff>_interactions_{per_frame,two_pass}.py` (three force fields: Amber / GROMOS / CHARMM).

## Supported Boundaries

- **All-atom force fields** (Amber / CHARMM / OPLS) rely on explicit hydrogen atoms for donor and hydrophobic-neighborhood determination; fully usable
- **GROMOS** is supported (53A6/54A7), but as a united-atom force field: polar H (N-H, O-H) and aromatic C-H are explicit, while aliphatic nonpolar H are merged into CH1/CH2/CH3 — coarse-type parts are compensated by structural criteria / tiered determination / anti-enumeration (acceptors, aromatic rings, hydrophobicity), and determinations depending on aliphatic H are restricted; GROMOS ligands (ATB parameterized) are not within the first-release support commitment
- **Robust to missing bond orders**: even if the tpr does not retain bond orders (all func=1), type names already encode aromaticity, so the determination is reliable; bond order/planarity are not group-identification criteria