# ChemVCS Python Layer

Python domain layer for ChemVCS, providing chemistry-specific abstractions and HPC integration on top of the Go-based VCS core.

## Features

- **Repository Interface**: Interact with ChemVCS repositories from Python
- **Domain Objects**: Chemistry-specific objects (Structure, Run, Workflow)
- **File Parsers**: Read/write common chemistry file formats (XYZ, POSCAR)
- **HPC Integration**: Submit and track HPC jobs with automatic provenance capture
- **Type-safe API**: Dataclass-based models with validation
- **Comprehensive tests**: 118 pytest tests covering all components

## Installation

```bash
# From source (recommended for development)
cd python
pip install -e .

# With development dependencies
pip install -e ".[dev]"
```

**Prerequisites**:
- Python ≥3.8
- NumPy ≥1.20.0
- ChemVCS CLI (Go binary) in PATH or `../go/`

## Quick Start

### Working with Molecular Structures

```python
from chemvcs_py import Structure
from chemvcs_py.io import read_xyz, write_xyz
import numpy as np

# Create a water molecule
water = Structure(
    formula="H2O",
    positions=np.array([
        [0.0, 0.0, 0.0],     # O
        [0.0, 0.757, 0.586],  # H
        [0.0, -0.757, 0.586]  # H
    ]),
    species=["O", "H", "H"]
)

# Write to XYZ file
write_xyz(water, "water.xyz")

# Read from file
structure = read_xyz("water.xyz")

# Convert to CoreObject for storage
core_obj = structure.to_core_object()
print(f"Structure: {structure.formula}, {structure.num_atoms} atoms")
```

### Tracking Computational Runs

```python
from chemvcs_py import Run

# Create a calculation run
opt_run = Run(
    structure_id="struct_001",
    code="ORCA",
    code_version="5.0.3",
    parameters={
        "method": "DFT",
        "functional": "B3LYP",
        "basis": "def2-TZVP"
    },
    status="planned"
)

# Track lifecycle
opt_run.mark_submitted(job_id="12345")
opt_run.mark_running()
opt_run.mark_finished()
opt_run.set_result("energy", -76.4321)

print(f"Energy: {opt_run.get_energy()} Hartree")
```

### Building Computational Workflows

```python
from chemvcs_py import Workflow

# Create workflow DAG
workflow = Workflow(name="Water Study")
workflow.add_node("opt", "run", {"description": "Optimize geometry"})
workflow.add_node("freq", "run", {"description": "Frequency calculation"})
workflow.add_edge("opt", "freq")  # freq depends on opt

# Get execution order
order = workflow.topological_sort()  # ['opt', 'freq']

# Query dependencies
deps = workflow.get_dependencies("freq")  # ['opt']
roots = workflow.get_root_nodes()  # ['opt']
```

### Repository Integration

```python
from chemvcs_py import Repo

# Open repository
repo = Repo()  # Finds .chemvcs in current or parent directories

# List all structures
structures = repo.list_objects(type_filter="structure")
for obj_info in structures:
    print(f"{obj_info['Hash'][:8]}: {obj_info['Type']}")

# Get specific object
obj = repo.get_object("abc123...")
structure = Structure.from_core_object(obj)
```

## Package Structure

```
chemvcs_py/
├── __init__.py         # Package exports
├── core/
│   ├── repo.py         # Repository interface (228 lines)
│   └── objects.py      # Core object representations (109 lines)
├── domain/
│   ├── structure.py    # Structure domain object (195 lines)
│   ├── run.py          # Run domain object (215 lines)
│   └── workflow.py     # Workflow domain object (318 lines)
├── io/
│   ├── xyz.py          # XYZ file format support (167 lines)
│   └── poscar.py       # POSCAR/VASP format support (286 lines)
└── util/
    └── errors.py       # Custom exceptions

examples/
├── basic_usage.py           # Structure + XYZ demonstration
└── run_workflow_example.py  # Run + Workflow + DAG demonstration

tests/
├── conftest.py             # pytest configuration
├── test_structure.py       # Structure tests (10 tests)
├── test_xyz.py             # XYZ parser tests (8 tests)
└── test_poscar.py          # POSCAR parser tests (10 tests)
```

**Statistics**:
- Production code: ~5,100 lines (including full HPC module)
- Test code: ~4,200 lines
- Test coverage: 157/157 tests passing

## Development

```bash
# Install in development mode
pip install -e ".[dev]"

# Run tests
pytest tests/ -v

# Run specific test file
pytest tests/test_structure.py -v

# Check test coverage
pytest tests/ --cov=chemvcs_py --cov-report=html

# Format code
black chemvcs_py/ tests/

# Type checking (optional)
mypy chemvcs_py/
```

## Examples

See the [examples/](examples/) directory for complete demonstrations:

- **basic_usage.py**: Creating structures, XYZ I/O, CoreObject conversion
- **run_workflow_example.py**: Computational runs, workflow DAGs, lifecycle tracking

Run examples:
```bash
python examples/basic_usage.py
python examples/run_workflow_example.py
```

## API Reference

### Domain Objects

#### Structure
Represents molecular or crystal structures.

```python
Structure(
    formula: str,
    positions: np.ndarray,  # Nx3 array
    species: List[str],
    lattice: Optional[np.ndarray] = None,  # 3x3 array for periodic
    metadata: Dict[str, Any] = None,
    id: Optional[str] = None
)
```

**Methods**:
- `translate(vector)` - Translate structure (in-place)
- `get_center_of_mass()` - Calculate center of mass
- `to_core_object()` / `from_core_object()` - CoreObject conversion

**Properties**:
- `num_atoms: int` - Number of atoms
- `is_periodic: bool` - Has lattice vectors

#### Run
Represents a computational calculation.

```python
Run(
    structure_id: str,
    code: str,
    code_version: Optional[str] = None,
    parameters: Dict[str, Any] = None,
    status: str = "planned",  # planned, submitted, running, finished, failed
    results: Dict[str, Any] = None,
    resources: Dict[str, Any] = None,
    metadata: Dict[str, Any] = None,
    id: Optional[str] = None
)
```

**Methods**:
- `mark_submitted(job_id)` - Mark as submitted
- `mark_running()` - Mark as running
- `mark_finished()` - Mark as successfully finished
- `mark_failed(error)` - Mark as failed
- `get_energy(key)` - Get energy result
- `set_result(key, value)` - Set result value

#### Workflow
Represents a DAG of computational tasks.

```python
Workflow(
    name: str,
    nodes: Dict[str, WorkflowNode] = None,
    edges: List[Tuple[str, str]] = None,
    metadata: Dict[str, Any] = None,
    id: Optional[str] = None
)
```

**Methods**:
- `add_node(node_id, kind, metadata)` - Add workflow node
- `add_edge(from_node, to_node)` - Add dependency edge (with cycle detection)
- `get_dependencies(node_id)` - Get node dependencies
- `get_dependents(node_id)` - Get dependent nodes
- `get_root_nodes()` - Nodes with no dependencies
- `get_leaf_nodes()` - Nodes with no dependents
- `topological_sort()` - Get execution order

### File I/O

#### XYZ Format
```python
from chemvcs_py.io import read_xyz, write_xyz

structure = read_xyz("molecule.xyz")
write_xyz(structure, "output.xyz", comment="Optional comment")
```

#### POSCAR Format (VASP)
```python
from chemvcs_py.io import read_poscar, write_poscar

structure = read_poscar("POSCAR")
write_poscar(structure, "POSCAR_out", direct=True, comment="NaCl")
```

### Repository

```python
from chemvcs_py import Repo

repo = Repo()  # Auto-discover .chemvcs directory

# List objects
objects = repo.list_objects(type_filter="structure")

# Get object
obj = repo.get_object(hash_value)

# Commit changes
repo.commit(message="Updated structure")
```

## Status

**Version**: 0.1.0

**Milestone 5 Progress**: 85% Complete ✅

**Completed Features**:
- ✅ Repository interface (Repo)
- ✅ Core object models (CoreObject, Reference)
- ✅ Structure domain object with validation
- ✅ Run domain object with lifecycle tracking
- ✅ Workflow domain object with DAG validation
- ✅ XYZ file parser/writer
- ✅ POSCAR parser/writer (VASP format)
- ✅ Python test suite (28 tests passing)
- ✅ Example scripts and documentation

**Remaining Work** (Lower Priority):
- ⏭️ CIF parser (crystallography format)
- ⏭️ Additional Run/Workflow tests
- ⏭️ Enhanced Repository query API
- ⏭️ Integration with QM code outputs (ORCA, Gaussian, etc.)

## Contributing

ChemVCS is in active development. Contributions welcome!

1. Fork the repository
2. Create a feature branch
3. Add tests for new features
4. Ensure all tests pass: `pytest tests/ -v`
5. Submit a pull request

## License

[To be determined]

## Related Documentation

- [Main ChemVCS README](../README.md)
- [Feature Status Tracking](../FEATURE_STATUS.md)
- [Python Domain Layer Design](../docs/06-python-domain-layer.md)
- [Repository Overview and Structure](../README.md)
