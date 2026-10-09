# Installation and Environment Setup

This document explains how to install DuIvyInteractions and verify the environment. The released package version is v0.0.1 (`pyproject.toml`), and result files use HDF5 format version 1.0.

## System Requirements

- **Python** >= 3.9 (`pyproject.toml` declares `requires-python = ">=3.9"`; classifiers cover 3.9–3.12; 3.10+ recommended)
- **Operating system**: Linux / macOS / Windows (Linux or macOS recommended; tpr/xtc parsing is pure Python + MDAnalysis and is OS-independent)
- **Input data**: one GROMACS topology file (tpr) and one trajectory file (xtc, or any format readable by MDAnalysis). `dii run` requires both: the tpr provides the force-field atom types and bonding graph (the sole basis of group identification), and the trajectory provides per-frame coordinates (geometric determination)

> **GROMACS is not required**: the default reading path, `GmxTprReader`, parses the tpr binary directly with MDAnalysis and never invokes `gmx`. GROMACS is needed only to use `GmxTprDumpReader` (a separate `Reader` implementation that parses `gmx dump` text output), and the current `Pipeline` always uses the direct binary path — there is no automatic fallback.

## Dependencies

Installing DuIvyInteractions automatically installs the following Python dependencies (exactly matching the `dependencies` in `pyproject.toml`):

| Dependency | Version Requirement | Purpose |
|:-----|:---------|:-----|
| numpy | >= 1.20 | Numerical computation (matrix-based detection, statistics) |
| MDAnalysis | >= 2.0 | tpr binary parsing (`gmx_tpr_reader`), trajectory loading (`mda.Universe(tpr, xtc)` in `Pipeline.run`) |
| h5py | >= 3.0 | HDF5 serialization of results (format version 1.0, gzip compression) |
| scipy | >= 1.7 | Spatial indexing KDTree (pre-screening in the hydrogen-bond / water-bridge TwoPass and PerFrame detectors) |
| DuIvyTools | >= 0.6.0 | Parsing, building, and later plotting of xvg/xpm result files |

Optional dependencies:

- **Development/testing**: `pip install -e ".[dev]"` additionally installs pytest (the `[project.optional-dependencies].dev` entry).

### MDAnalysis Version and tpr Compatibility

Binary tpr parsing depends on MDAnalysis support for the corresponding GROMACS version (the tpx version number). tpr files produced by newer GROMACS versions require newer MDAnalysis; if parsing reports `Your tpx version is XXX, which this parser does not support`, upgrade MDAnalysis (the development environment uses 2.7–2.10; 2.10 has been verified to parse the GROMACS 2018.1 tpr files bundled with the repository).

## Installation

### Method 1: Install from PyPI (Recommended)

```bash
pip install duivyinteractions
```

The package name is `duivyinteractions` (v0.0.1 published on PyPI). All dependencies are installed automatically, and the `dii` command is registered on PATH (entry: `[project.scripts]` in `pyproject.toml`, pointing to `DuIvyInteractions.DII:main`).

### Method 2: Install from Source (Development)

Run this in the project root directory (the directory containing `pyproject.toml`):

```bash
pip install -e .
```

The `-e` flag performs an editable install so that source changes take effect immediately, which is suitable for development and debugging. For a release build, run `pip install .` without `-e`.

## Verify Installation

Check whether the command-line tool is available:

```bash
dii --help
```

It should output (with the two subcommands `run` and `export` registered):

```
usage: dii [-h] {run,export} ...
```

Inspect the subcommand arguments (note the choices of `--ff` and `--strategy`):

```bash
dii run --help
dii export --help
```

Verify that the module can be imported in Python:

```python
from DuIvyInteractions.pipeline import Pipeline
from DuIvyInteractions.io import load_interactions

print("DuIvyInteractions OK")
```

With an installed package you can also inspect the version and dependencies:

```bash
pip show duivyinteractions
```

## Next Steps

After installation, proceed to [Quick Start](quickstart) and run the full workflow with the real test data bundled with the repository (Amber / GROMOS / CHARMM36).
