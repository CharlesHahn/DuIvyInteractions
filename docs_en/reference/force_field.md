# Force Field Type Mapping

This project determines chemical groups directly from the force field atom types in the GROMACS tpr (see [Core Concepts](concepts.md)). This document lists the type-to-chemical-feature mapping rules for the **four registered force field families** — the Amber family (Amber proteins + GAFF ligands), GROMOS 53A6/54A7, CHARMM36/C36m (proteins + CGenFF ligands), and OPLS-AA/L — so that research users can understand the basis of the determination and assess applicability.

All feature-table constants are defined at the top of the identifier modules: `DuIvyInteractions/group_identifiers/{amber,gromos,charmm,opls}_ff_identifier.py`. The full multi-force-field compatibility survey is archived in `doc/force_field_compatibility_survey.md`.

## Mapping Principles and Feature-Space Design

The type → chemical feature mapping comes from two independent, cross-validatable lines of evidence:

1. **GAFF types**: from the official antechamber naming rules (Wang et al., *J Comput Chem* 2004) — the suffix letters/digits of a type name encode chemical features (e.g., the `a` in `ca` = aromatic)
2. **Protein types**: inferred backwards from rtp/rtf residue definitions — the known chemical facts of residues (e.g., the TYR benzene ring is aromatic) are compared against the type names used for ring atoms in the rtp file (e.g., Amber's `CA`), yielding that type = aromatic carbon

The design is a **feature-space mapping**: each type maps to a feature vector {hybridization, aromaticity, polarity, has H, lone-pair availability}, and groups are determined from feature combinations (e.g., "aromatic ring = ring with ≥ n−1 strong aromatic atoms"). Registering a new force field = filling in one feature table, without changing the determination logic (see [Extension Guide](extension.md)).

**Registered force fields and boundaries**:

| Force field family | Proteins | Ligands | Boundary |
|:-------------------|:---------|:--------|:---------|
| Amber (amber03/94/96/99/99sb/99sb-ildn/GS/14sb) | ✓ | GAFF/GAFF2, mapping verified zero conflicts | — |
| GROMOS 53A6/54A7 | ✓ (coarse types, compensated by structural criteria / tiered ring determination / anti-enumeration) | ✗ not supported (ATB parameterized, declared for the first release) | coarse types (`N`/`C` cover many chemical environments); united-atom, so determinations that depend on explicit aliphatic H (donors, hydrophobic neighborhood) are restricted |
| CHARMM36/C36m | ✓ | CGenFF (`NG*` acceptor qualification judged per type) | — |
| OPLS-AA/L (2001) | ✓ | types within tables | requires `opls_XXX` number mapping |

## Amber Family Feature Tables

Constants are defined at the top of `DuIvyInteractions/group_identifiers/amber_ff_identifier.py`.

### Aromatic Types (STRONG_AROMATIC)

Atoms of these types are **directly determined to be aromatic by their type names** and serve as strong evidence for ring detection (a ring is aromatic when ≥ n−1 of its atoms belong to this set; see [Group Identification Rules](group_rules.md)):

| Category | Types | Notes |
|:---------|:------|:------|
| GAFF aromatic carbon | `ca, cg, ch, cm, cn, cp, cq, c1` | `c` + aromatic suffix |
| GAFF aromatic nitrogen | `na, nb, nh, ni, nj, n1, n2` | `na` = pyrrole-type, `nb` = pyridine-type |
| GAFF aromatic phosphorus | `pb` | — |
| Amber protein aromatic carbon | `CA, CB, CC, CK, CM, C5, C6, C7, C*, CW, CR, CN, CV, CQ` | rtp ring atom types |
| Amber protein aromatic nitrogen | `NA, NB, NC, N*` | — |

### Compatible Types (COMPATIBLE_TYPES)

Non-aromatic types that can participate in conjugation under the "coercion" of **n−1 atoms in the ring being aromatic types**. Used to handle ambiguous types and heterocycles:

| Type | Notes |
|:-----|:------|
| `C, N` | Amber protein ambiguous types (main-chain carbonyl carbon/aromatic ring carbon `C`; amide nitrogen/aromatic nitrogen `N`) |
| `os, ss` | GAFF furan oxygen / thiophene sulfur |
| `cc, cd` | GAFF conjugated ring carbons that are not purely aromatic |
| `pc, pd` | GAFF sp2 phosphorus in conjugated rings |

### Acceptor Types (ACCEPTOR_TYPES)

Candidate types for hydrogen bond acceptors (those with an available lone pair). **Note**: this table was revised by the 2026-09-30 acceptor-qualification fix (A1), which removed the H-bearing non-acceptor nitrogen types (see "[Acceptor Qualification](#acceptor-qualification)" below):

| Category | Types |
|:---------|:------|
| GAFF oxygen | `o, o2, oh, os, oe, o1, ow` |
| GAFF nitrogen | `n2, n3, nb, ni, nj, nc, ne, nf, nk` — `n2` imine-type (sp2 N, disubstituted, lone pair available, acceptor; GAFF official type-table definition), `n3` neutral amine; `n` (ligand amide) and `na/nh` (pyrrole-type, lone pair in aromatic conjugation) removed |
| GAFF sulfur | `s, ss, sh, sx, s2` |
| Halogens | `f, cl, br, i` |
| Amber protein oxygen | `O, OH, O2, OS, OW` |
| Amber protein nitrogen | `N, N2, NB, NC` — `N` serves both ordinary amide (H-bearing, non-acceptor) and Pro N (H-free, acceptor); the same type is ambiguous and is distinguished by H-neighbor count in `_find_acceptors`; `N2` nucleic-acid amino (adenine N6, 2H, acceptor); `NB/NC` H-free pyridine/in-ring nitrogen; `N3` (ammonium), `NA` (H-bearing pyrrole), `N*` (nucleoside glycosidic N, 3-bond pyrrole-type) removed |
| Amber sulfur | `S, SH` |

### Other Constants

- **METAL_IONS**: set of 42 metal elements (from PLIP config.py), see "Metal Ions" below
- **WATER_RESIDUES**: `{SOL, HOH, WAT}` (class attribute, see "Water Molecules" below)
- **CHARGE_THRESHOLD**: 0.1 (net-charge validation threshold for charged groups)
- **Positive residue dictionary** (residue name → atom-name lists): ARG, LYS, HIP, ORN, DAB, M3L, MLY
- **Negative residue dictionary**: ASP, GLU, CYM, KCX, PCA, SEP, TPO, PTR

## GROMOS Feature Tables (53A6 / 54A7)

Constants are defined in `DuIvyInteractions/group_identifiers/gromos_ff_identifier.py` (measured against the local gromos54a7.ff; common to 53A6/54A7). GROMOS types are **coarse** (`N`/`C` cover many chemical environments) and cannot be distinguished by type alone; acceptors, aromatic rings, and hydrophobic atoms are therefore determined with compensating structural criteria / tiered ring determination / anti-enumeration.

### Acceptor Types (GROMOS_ACCEPTOR_TYPES) and the N Structural Criterion

```text
GROMOS_ACCEPTOR_TYPES = {O, OM, OA, OE, OW,          # oxygen: carbonyl/carboxyl/hydroxyl/ether/water
                         N, NT, NL, NR, NZ, NE,     # nitrogen: peptide/terminal/aromatic/guanidinium
                         S,                          # sulfur
                         F, CL, BR}                  # halogens
```

`_find_acceptors` applies a structural criterion to atoms of element N (chemistry + literature; see `doc/acceptor_identification_evidence.md`):

- **≥ 4 neighbors (ammonium RNH₃⁺) → excluded** (no lone pair, non-acceptor)
- **has H neighbor → excluded** (ordinary amide / side-chain amide / H-bearing pyrrole / guanidinium, non-acceptor)
- **H-free → kept** (Pro N / H-free pyridine-type His N, acceptors)

### Tiered Aromatic Ring Determination

GROMOS has no per-type aromatic table; determination is tiered (`_filter_aromatic_rings` override):

| Set | Contents | Role |
|:----|:---------|:-----|
| `GROMOS_RING_TYPES` | `{C, CR1, NR}` | full set of ring carbon/nitrogen types; ring composition filter |
| `GROMOS_AROMATIC_STRONG` | `{NR}` | strong signal: aromatic nitrogen; ring containing NR is directly aromatic |
| `GROMOS_AROMATIC_RESIDUES` | PHE, TYR, TRP, HISA, HISB, HISH, HIS1, HIS2 | weak signal: all-C rings only for standard aromatic protein residues |

Determination order: ring composition ⊆ `{C, CR1, NR}` → containing `NR` is directly aromatic → all-C rings fall back to the residue-name whitelist.

### Hydrophobic Anti-Enumeration

GROMOS is a united-atom force field (no H neighbors on aliphatic C), so hydrophobicity is not determined by the Amber rule "C with all neighbors C/H"; instead a **type whitelist + polar-neighbor exclusion** is used (`_find_hydrophobic` override):

- `GROMOS_HYDROPHOBIC_TYPES` = `{C, CH0, CH1, CH2, CH3, CH4, CH2r}` (**excluding** CH3p polar choline N⁺ and CR1 aromatic/alkene sp2)
- `GROMOS_HYDROPHOBIC_EXCLUDED` = `{O, N, S}` — any neighbor in this set excludes the atom

### Water Residues and Charged Residues

- **WATER_RESIDUES**: `{SOL}` (GROMOS is paired with SPC/SPC-E)
- **Positive residue dictionary**: ARG, LYSH (protonated ε-ammonium NH₃⁺), HISH (doubly protonated His) — GROMOS splits LYS into neutral LYS / protonated LYSH; only the latter carries +1
- **Negative residue dictionary**: ASP, GLU
- Charged-group deduplication override `_deduplicate_charged`: removes groups strictly contained in a more complete group of the same type (the positively charged N atoms of GROMOS protonated residues themselves carry positive partial charges, so the functional-group layer would additionally emit a per-N charged group; see [Group Identification Rules](group_rules.md))

## CHARMM Feature Tables (CHARMM36 / C36m + CGenFF)

Constants are defined in `DuIvyInteractions/group_identifiers/charmm_ff_identifier.py` (verbatim consistent with CHARMM-GUI official `top_all36_prot.rtf` / `top_all36_cgenff.rtf`).

### Acceptor Types (CHARMM_ACCEPTOR_TYPES)

| Category | Kept (acceptors) | Excluded (non-acceptors) |
|:---------|:-----------------|:-------------------------|
| Protein oxygen | `O, OB, OC, OH1, OS` | — |
| CGenFF oxygen | `OG2D1–5, OG2P1, OG2R50, OG301–304, OG311, OG312, OG3C51, OG3C61, OG3R60` | — |
| Protein nitrogen | `N` (Pro N, H-free, acceptor), `NR2` (H-free pyridine-type His N, acceptor) | `NH1` (peptide), `NH2` (amide), `NH3` (ammonium), `NC2` (guanidinium), `NY` (pyrrole), `NR1/NR3` (protonated His) |
| CGenFF nitrogen | `NG2D1` (neutral imine/Schiff base), `NG2R50` (purine N7), `NG2R60/NG2R62` (6-membered pyridine-type), `NG2S3` (exocyclic amine/aniline-type), `NG3N1` (hydrazine N, sp3 amine) — **6 kept** | `NG2O1` (nitrobenzene N), `NG2P1` (protonated imine), `NG2R51` (5-membered singly bonded sp2 N = His/Trp pyrrole with H), `NG2R52` (protonated Schiff base/amidine/guanidine), `NG2R61` (6-membered singly bonded imine N), `NG2RC0` (bridgehead N), `NG2S0` (N,N-disubstituted amide), `NG2S1` (peptide N), `NG2S2` (terminal amide N), `NG2S4` (hydroxamic-acid N) — **10 excluded** |
| Sulfur | `S, SM, SS` | — |
| Halogens | `FGA1–3, FGR1, CLGA1, CLGA3, CLGR1, BRGA1–3, BRGR1, IGR1` (CGenFF names) | — |

> Each CGenFF `NG*` type was judged item by item against the official MASS-section comments (6 kept / 10 excluded); `NG2S0` (N,N-disubstituted amide) is registered as a known exception: it would include ligand Pro-type weak acceptors, but they are inseparable from ordinary tertiary amides in the same type (see `doc/acceptor_identification_evidence.md` §3).

### Aromatic Ring Types

- **CHARMM_STRONG_AROMATIC**: protein aromatic carbon `CA, CAI, CPH1, CPH2, CPT, CY`; protein aromatic nitrogen `NR1, NR2, NR3, NY`; CGenFF aromatic carbon `CG2R51–53, CG2R57, CG2R61–64, CG2R66, CG2R67, CG2R71, CG2RC0`; CGenFF aromatic nitrogen `NG2R50–52, NG2R57, NG2R60–62, NG2R67, NG2RC0`; heterocyclic O/S `OG2R50, SG2R50`
- **CHARMM_COMPATIBLE_TYPES**: `C, CC, CD` (in-ring main-chain carbonyl/carboxyl/amide carbon), `CG2D1, CG2D2, CG2DC1–3` (conjugated alkenes), `CG2O1–6` (in-ring carbonyl)

### Water Residues and Charged Residues

- **WATER_RESIDUES**: `{TIP3, HOH, SOL, WAT}` (CHARMM defaults TIP3/HOH; also accepts the GROMACS-conventional SOL)
- **Positive residue dictionary**: LYS, ARG, HSP (doubly protonated His — the HSP imidazole nitrogens carry negative partial charges, so the imidazole carbons CG/CE1/CD2 (types CPH1/CPH2) must be included to make the charge center positive)
- **Negative residue dictionary**: ASP, GLU, CYM (thiolate Cys⁻)
- CHARMM protein N atoms carry negative partial charges; N-terminal NH₃⁺ is identified via a **terminal-ammonium structural criterion** (`_is_terminal_ammonium`: N neighbors ⊆ {C,H} and heavy neighbors not bonded to ≥2 N; charge atoms = N + bonded H); see [Group Identification Rules](group_rules.md)

## OPLS Feature Tables (OPLS-AA/L 2001)

Constants are defined in `DuIvyInteractions/group_identifiers/opls_ff_identifier.py` (from protein rtp measurements + ffnonbonded symbol mapping; version = OPLS-AA/L 2001).

### Acceptor Types (OPLS_ACCEPTOR_TYPES)

| Category | Kept (acceptors) | Excluded (non-acceptors) |
|:---------|:-----------------|:-------------------------|
| Oxygen | `opls_236, opls_272, opls_154, opls_167, opls_268, opls_269` (O/O2/OH/O_3) | — |
| Nitrogen | `opls_239` (Pro N, H-free), `opls_511` (H-free pyridine-type His N), `opls_900` (LYS neutral amine NZ) — **3 kept** | `opls_238` (main-chain N), `opls_237` (side-chain amide N), `opls_287` (LYSH ammonium), `opls_300/opls_303` (Arg guanidinium), `opls_503` (H-bearing pyrrole: HISD ND1, TRP NE1), `opls_512` (doubly protonated His N), `opls_749/750/751` (ARGN neutral guanidine) — **10 excluded** |
| Sulfur | `opls_202, opls_200` (S/SH) | — |
| Halogens | ionic halides `opls_400–403`; organic halogens `opls_123, 151, 164, 226, 264, 709, 719, 721, 722, 726, 728, 730, 732, 786, 956, 965` | — |

### Aromatic Ring Types

- **OPLS_STRONG_AROMATIC**: CA family `opls_145, 166, 302, 752`; 5-membered rings `opls_500–502, 506–510, 514`; aromatic nitrogen `opls_503, 511, 512` (NA/NB)
- **OPLS_COMPATIBLE_TYPES**: `opls_235, opls_267` (in-ring carbonyl carbon), `opls_271` (in-ring carboxylic-acid carbon)

### Water Residues and Charged Residues

- **WATER_RESIDUES**: `{HOH, HO4, HO5, SOL, WAT}` (HOH/SPC, HO4/TIP4P, HO5/TIP5P)
- **Positive residue dictionary**: ARG, LYSH (protonated ε-ammonium), HISH (doubly protonated His, includes imidazole carbons CG/CE1/CD2) — OPLS LYS is neutral, LYSH protonated
- **Negative residue dictionary**: ASP, GLU
- OPLS ammonium nitrogen carries a negative partial charge; the same **terminal-ammonium structural criterion** as CHARMM is reused (`_is_terminal_ammonium`, identical logic, independently maintained) to identify the LYSH side-chain ammonium and the N-terminal NH₃⁺

(acceptor-qualification)=
## Acceptor Qualification (A1 fix, 2026-09-30)

Hydrogen bond acceptor qualification depends on whether the N in a given environment has an **available lone pair**. The fix removed the H-bearing non-acceptor N types (ordinary amide / ammonium / H-bearing pyrrole / guanidinium) and kept Pro N, H-free pyridine-type His N, neutral amines, and nucleic-acid amino groups. Determination for the nine N environments (the full evidence list is archived in `doc/acceptor_identification_evidence.md`):

| N environment | Lone-pair state | Acceptor? | Basis |
|:--------------|:----------------|:----------|:------|
| Ordinary amide N (H-bearing, bonded to C=O) | in amide resonance | no | Eildal 2013, *JACS* 135:12998; InterMap analysis (critique of ProLIF peptide-N misjudgment, `doc/intermap_paper_analysis.md`) |
| Ammonium N (4 bonds, RNH₃⁺) | none | no | textbook (four-bond ammonium) |
| H-bearing pyrrole N (in aromatic ring) | in aromatic conjugation | no | textbook (heterocyclic chemistry) |
| Guanidinium N (bonded to C(N)(N)) | in guanidinium resonance | no | chemical facts (only "guanidine as donor" literature found; no "guanidine as acceptor" literature) |
| Pyridine-type N (H-free, in aromatic ring) | lone pair in sp2 orbital | yes | textbook (pyridine can be protonated/coordinated) |
| Neutral amine N (H-bearing, sp3) | free lone pair | yes | Luisi 1998, *J Mol Biol* |
| Nucleic-acid amino N (adenine N6, guanine N2, 2H) | conjugated with ring but still available | yes | Luisi 1998; Baik 2003, *JACS* |
| Pro N (H-free, bonded to carbonyl C) | more basic than ordinary amide | yes | Deepak 2016, *Biophys J* |
| Purine glycosidic N (N9/N1, 3-bond in-ring sp2 N) | in purine aromaticity (pyrrole-type) | no | textbook heterocyclic chemistry; base-pairing acceptor sites are N1/N3/N7 (revised by the 2026-09-30 adversarial review, evidence list §2/§3) |

**Implementation per force field**:

- **Amber**: `ACCEPTOR_TYPES` table + type `N` distinguished by H-neighbor count (H-bearing ordinary amide skipped, H-free Pro N kept) + exclusion set (positively charged group members, in-ring H-bearing N)
- **GROMOS**: `GROMOS_ACCEPTOR_TYPES` table + N structural criterion (≥4 bonds or H-bearing excluded, H-free kept)
- **CHARMM / OPLS**: the type tables already encode per-type judgment (no ambiguous N types), applied directly

In all force fields, acceptors additionally require **partial charge q < 0**.

**Known boundaries** (registered in evidence list §3): GROMOS's simplified structural criterion "H-bearing → excluded" also excludes nucleic-acid amino groups (which Luisi judges as acceptors and Amber keeps via `N2`); GROMOS ligand/nucleic-acid scenarios are therefore restricted. CHARMM `NG2S0` (N,N-disubstituted amide) is inseparable from ligand Pro-type weak acceptors, so the table exclusion loses that weak-acceptor class.

## Donor Determination

The **bond entries between a donor atom D (N/O/S/F) and its hydrogens H must exist explicitly**, and H must carry a positive charge (q(H) > 0). In the tpr, N–H bonds reside in the **Constraint section** (not the Bond section); the reader **merges the Bond, Constraint (and Settle) sections into a unified per-residue bond list** (`_parse_molblock` in `gmx_tpr_dump_reader.py` writes all three ilist sections into `bonds`; the MDAnalysis reader merges via `u.bonds` as well). Donor identification depends on this merged bond graph — scanning only the Bond section is insufficient; D–H pairs across residue boundaries are additionally checked by `_find_inter_residue_donors`. This determination involves no coordinate inference.

## Metal Ions (METAL_IONS)

From PLIP config.py, 42 elements (initial capital, matching element symbols):

`Ca, Co, Mg, Mn, Fe, Cu, Zn, Li, Na, K, Rb, Sr, Cs, Ba, Cr, Ni, Ru, Rh, Pd, Ag, Cd, La, W, Os, Ir, Pt, Au, Hg, Ce, Pr, Sm, Eu, Gd, Tb, Yb, Lu, Al, Ga, In, Sb, Tl, Pb`

Metal centers are identified by **element symbol** (the reader infers the element via `ATOMIC_NUMBER_TO_ELEMENT` in `core/constants.py`), independent of type names.

## Water Molecules (WATER_RESIDUES)

The `GroupIdentifier.WATER_RESIDUES` class attribute (base class default is an empty `frozenset()`, overridden by each identifier subclass) is the **unified cross-force-field entry point**: the pipeline excludes water using this set per force field, and the identifier-internal `_find_water` / `_find_metal_binding` also use `self.WATER_RESIDUES`:

| Force field | WATER_RESIDUES (class attribute) | Notes |
|:------------|:---------------|:------|
| Amber (`amber_ff_identifier.py`) | `{SOL, HOH, WAT}` | SOL is the main name |
| GROMOS (`gromos_ff_identifier.py`) | `{SOL}` | paired with SPC/SPC-E |
| CHARMM (`charmm_ff_identifier.py`) | `{TIP3, HOH, SOL, WAT}` | CHARMM defaults TIP3/HOH |
| OPLS (`opls_ff_identifier.py`) | `{HOH, HO4, HO5, SOL, WAT}` | HOH/SPC, HO4/TIP4P, HO5/TIP5P |

Non-`SOL` water residue names such as CHARMM TIP3 and OPLS HO4/HO5 are therefore identified and excluded correctly (fix of 2026-10-08).

## Type Name Ambiguities and Handling

- **`C` dual identity**: in Amber, `C` serves both as the main-chain carbonyl carbon and as the Tyr CZ aromatic ring carbon. When the ring has ≥ n−1 aromatic atoms (5 of 6 for a six-membered ring), `C` is "promoted" as a compatible type to participate in conjugation.
- **`N3` vs `n3`**: uppercase `N3` = Amber protein positively charged amino group (ammonium, non-acceptor); lowercase `n3` = GAFF neutral amine (acceptor) — the case differs in meaning.
- **`CA` has different meanings across force fields**: in the Amber type system `CA` is an aromatic carbon; in the GROMOS system `CA` denotes the α carbon. They must be distinguished through each force field's independent feature tables; type names cannot be reused across force fields.
- **N–H bonds are in the Constraint section**: donor identification must merge Bond + Constraint (+ Settle) bond types.

## Compatibility Scope

- **Amber family**: amber03, amber94, amber96, amber99, amber99sb, amber99sb-ildn, amberGS, amber14sb + GAFF/GAFF2 ligands; all cross-version type-name differences (e.g., CX vs CT) are handled, and the mapping has been verified to have zero conflicts
- **GROMOS**: 53A6 / 54A7 family (measured against local gromos54a7.ff)
- **CHARMM**: CHARMM36 / C36m proteins + CGenFF ligands (verbatim consistent with CHARMM-GUI official rtf)
- **OPLS**: OPLS-AA/L 2001 (GROMACS oplsaa.ff)

If a new version introduces new types, simply add one entry to feature tables such as `STRONG_AROMATIC`/`COMPATIBLE_TYPES`/`ACCEPTOR_TYPES` (the design goal of the feature-space mapping; full procedure in [Extension Guide](extension.md)).
