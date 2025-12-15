# ChemVCS Python Layer

Python domain layer for ChemVCS, providing chemistry-specific abstractions on top of the Go-based VCS core.

## Features

- **Repository Interface**: Interact with ChemVCS repositories from Python
- **Domain Objects**: Chemistry-specific objects (Structure, Run, Workflow)
- **File Parsers**: Read/write common chemistry file formats (XYZ, etc.)
- **Type-safe API**: Dataclass-based models with validation

## Installation

```bash
# From source
cd python
pip install -e .

# With development dependencies
pip install -e ".[dev]"
```

## Quick Start

```python
from chemvcs_py import Repo, Structure
from chemvcs_py.io import read_xyz

# Open repository
repo = Repo()  # Finds .chemvcs in current or parent directories

# Read structure from XYZ file
structure = read_xyz("molecule.xyz")

# Convert to core object and inspect
obj = structure.to_core_object()
print(obj.type)  # "structure"
print(obj.meta)  # {"formula": "H2O", "positions": [...], ...}

# List all structures in repository
objects = repo.list_objects(type_filter="structure")
for obj_info in objects:
    print(f"{obj_info['Hash'][:8]}: {obj_info['Type']}")
```

## Package Structure

```
chemvcs_py/
├── __init__.py
├── core/
│   ├── repo.py         # Repository interface
│   └── objects.py      # Core object representations
├── domain/
│   ├── structure.py    # Structure domain object
│   ├── run.py          # Run domain object (planned)
│   └── workflow.py     # Workflow domain object (planned)
├── io/
│   └── xyz.py          # XYZ file format support
└── util/
    └── errors.py       # Custom exceptions
```

## Development

```bash
# Install in development mode
pip install -e ".[dev]"

# Run tests
pytest

# Format code
black chemvcs_py/

# Type checking
mypy chemvcs_py/
```

## Requirements

- Python >=3.8
- numpy >=1.20.0
- ChemVCS CLI (Go binary) in PATH or ../go/

## Status

**Version**: 0.1.0 (Alpha)

Currently implemented:
- ✅ Repository interface (Repo)
- ✅ Core object models (CoreObject, Reference)
- ✅ Structure domain object
- ✅ XYZ file parser/writer
- ✅ List and inspect objects via CLI

Planned for M5:
- 🔲 Run domain object
- 🔲 Workflow domain object
- 🔲 POSCAR parser (VASP)
- 🔲 CIF parser (crystallography)
- 🔲 Integration tests

## License

[To be determined]
