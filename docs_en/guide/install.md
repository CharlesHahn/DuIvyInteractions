# Installation and Environment Setup

DuIvyInteraction is an intermolecular interaction analysis tool based on MD topology force field parameters. This document explains how to install and verify the environment.

## System Requirements

- **Python** >= 3.9 (3.10+ recommended)
- **Operating system**: Linux / macOS / Windows (Linux or macOS recommended; full support for GROMACS tpr parsing)

## Dependencies

Installing DuIvyInteractions automatically installs the following Python dependencies:

| Dependency | Version Requirement | Purpose |
|:-----|:---------|:-----|
| numpy | >= 1.20 | Numerical computation |
| MDAnalysis | >= 2.0 | tpr binary reading, trajectory loading |
| h5py | >= 3.0 | HDF5 serialization of results |
| scipy | >= 1.7 | Spatial indexing (KDTree, TwoPass detector pre-screening) |
| DuIvyTools | >= 0.6.0 | Parsing and visualization of xvg/xpm result files |

In addition, to parse the tpr in text format (`gmx_tpr_dump_reader`), **GROMACS** must be installed and the `gmx dump` command must be available. For daily use, it is recommended to use `gmx_tpr_reader` directly (MDAnalysis reads the binary tpr), which does not require GROMACS.

## Installation

### Method 1: Install from PyPI (Recommended)

```bash
pip install duivyinteractions
```

Release versions have been packaged and uploaded to PyPI (current version v0.0.1), and all dependencies are installed automatically.

### Method 2: Install from Source (Development)

Run this in the project root directory:

```bash
pip install -e .
```

`-e` is an editable install, suitable for development and debugging; in a release environment you can omit `-e` and run `pip install .` directly.

## Verify Installation

Check whether the command-line tool is available:

```bash
dii --help
```

It should output `usage: dii [-h] {run,export} ...`, indicating that the `dii` command has been installed.

Verify that the module can be imported in Python:

```python
from DuIvyInteractions.pipeline import Pipeline
from DuIvyInteractions.io import load_interactions
print("DuIvyInteractions OK")
```

## Next Steps

After installation, proceed to [Quick Start](quickstart) and run the full workflow with real data.