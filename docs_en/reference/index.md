# Reference

This section is intended for **researchers who want to understand the principles and customize the tool**: core concepts, force-field type mapping, group identification rules, interaction criteria, data format, the Python API, and extension methods.

## Current-State Overview

The reference implementation state of DuIvyInteractions v0.0.1 is as follows (each item corresponds to the code in `DuIvyInteractions/`):

- **Two-stage architecture**: ① group identification — performed once and frame-independent, deterministically identifying groups from tpr atom types + bonding graph + explicit hydrogens + charges; ② geometric determination — per-frame distance/angle/planarity criteria for 8 interaction types. See [Core Concepts](concepts), [Group Identification Rules](group_rules), and [Interaction Criteria](criteria).
- **4 force-field identifiers**: the Amber family (proteins + GAFF/GAFF2), GROMOS 53A6/54A7, CHARMM36/C36m (including CGenFF), and OPLS-AA/L, managed by the `IDENTIFIER_CLASSES` registry; water exclusion uses each identifier's own `WATER_RESIDUES` class attribute (correctly excluding CHARMM TIP3, OPLS HO4/HO5, etc.). See [Force-Field Type Mapping](force_field) and [Python API](api).
- **8 interaction types × 3 strategies**: `two_pass` (default, best performance) / `per_frame` (current) / `per_tuple` (strategy one, in a "may be deprecated" state); H-bond acceptor determination has been revised by fix A1 (removing non-acceptor N such as ammonium/guanidinium/H-bearing pyrrole; evidence in `doc/acceptor_identification_evidence.md`).
- **Storage and export**: lossless HDF5 serialization (format version 1.0) + 8 exporters producing xvg/xpm/csv; command line `dii run` / `dii export`. See [HDF5 Storage Format](data_format) and [Command Reference](../guide/command.md).
- **Real-data test cases**: three real MD datasets — Amber KRAS–RBD D927, GROMOS 53A6 (protein + 6 ligands), and CHARMM36 SMO–BST (`Tests/test_MD_case_*`), covering all 8 interaction types.
- **Scope and limitations**: PBC, metal coordination geometry, water-bridge deduplication, hydrophobic-aromatic deduplication, and visualization are not implemented; see [Known Limitations](limitations).

```{toctree}
:maxdepth: 1
:caption: Core Concepts

concepts
glossary
```

```{toctree}
:maxdepth: 1
:caption: Scientific Details

force_field
group_rules
criteria
```

```{toctree}
:maxdepth: 1
:caption: Data and Extension

data_format
api
extension
```

```{toctree}
:maxdepth: 1
:caption: Scope and Limitations

limitations
```