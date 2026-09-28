# Core Concepts

This document explains how the tool works: how force field atom types are read from the GROMACS tpr topology to identify chemical groups, and how intermolecular interactions are detected based on them.

## Method Overview

Intermolecular interaction analysis requires first determining the chemical structure of the molecules — which atoms form aromatic rings, which atoms are hydrogen bond donors/acceptors, which groups carry charges, and so on. Based on these chemical groups, the tool then performs geometric evaluation using trajectory coordinates.

There are two approaches to determining the chemical structure:

1. **Rebuild from coordinates**: Starting from atomic coordinates, infer bond orders, aromaticity, and hydrogen atom positions. Commonly used tools (e.g., PLIP with OpenBabel, ProLIF with RDKit) follow this approach.
2. **Read from force field parameters** (this tool): Read force field atom types directly from the GROMACS topology (tpr). These types were determined by antechamber/sobtop according to the molecular electronic structure during force field parameterization.

This tool adopts the second approach.

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

Amber protein types are similar: the type names used for ring atoms in residue definitions (rtp files) (e.g., `CA`) can be cross-referenced against known residue chemistry to deduce their meaning.

Type names encode features such as hybridization state, aromaticity, polarity, and whether they carry hydrogens. Therefore, reading the types is sufficient to determine chemical groups, with no need to infer from coordinates. This guarantees: group identification depends only on topology, not on coordinates, and is consistent across frames and reproducible.

## Two-Stage Architecture

The analysis is divided into two stages:

```
① 基团鉴定
   tpr 力场参数（原子类型 + 键合图 + 显式 H + 电荷）
   → 识别化学基团（芳香环、H 键供受体、带电基团等）
   —— 只做一次，与帧无关

② 几何判定
   基团带标签后，逐帧用距离/角度/平面临近判据
   → 每帧相互作用列表 → 时间统计
```

- **Group identification** depends only on topology, not on coordinates, so the results are frame-independent and reproducible
- **Geometric evaluation** computes geometric metrics (distance, angle, etc.) frame by frame for each pair, and determines whether an interaction exists according to the criteria

## Explicit Hydrogens and Chemical Information

In all-atom force fields (Amber family), the following chemical information is explicitly present in the topology:

- **Hydrogen bond donor**: the bond entry between the donor atom D and hydrogen H exists, and H carries a positive charge
- **Water molecules**: identified by residue names (SOL/HOH), with oxygen atom OW and hydrogen atoms HW explicit
- **Metal center**: identified by element (e.g., Mg) and charge

This information is read directly from the topology and involves no coordinate inference.

## Supported Scope

- **Force fields**: Amber family (amber03/94/96/99/99sb/99sb-ildn/GS/14sb) proteins + GAFF/GAFF2 ligands, with type mapping verified to have zero conflicts
- **Requires explicit hydrogens**: the tool relies on the explicit hydrogen atoms of all-atom force fields; donor identification does not apply to united-atom force fields (e.g., GROMOS, no explicit H)
- **Missing bond orders**: even if the tpr does not preserve bond orders (all func=1), aromaticity determination is still based on type names; bond orders serve only as cross-validation

## Detailed Rules

- Complete mapping of type → chemical features: see [Force Field Type Mapping](force_field.md)
- Identification rules for each group: see [Group Identification Rules](group_rules.md)
- Geometric criteria for each interaction: see [Interaction Criteria](criteria.md)

## References

- GAFF force field: Wang J, Wolf RM, Caldwell JW, Kollman PA, Case DA. Development and testing of a general amber force field. *J Comput Chem*. 2004;25(9):1157-1174.
- Amber force field (ff14SB): Maier JA, Martinez C, Kasavajhala K, Wickstrom L, Hauser KE, Simmerling C. ff14SB: Improving the accuracy of protein side chain and backbone parameters from ff99SB. *J Chem Theory Comput*. 2015;11(8):3696-3713.
- PLIP interaction definitions: Salentin S, Schreiber S, Haupt VJ, Adasme MF, Schroeder M. PLIP: fully automated protein-ligand interaction profiler. *Nucleic Acids Res*. 2015;43(W1):W443-W447.
- MDAnalysis: Michaud-Agrawal N, Denning EJ, Woolf TB, Beckstein O. MDAnalysis: A toolkit for the analysis of molecular dynamics simulations. *J Comput Chem*. 2011;32(10):2319-2327.
