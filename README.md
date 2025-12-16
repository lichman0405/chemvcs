# ChemVCS

A version control system designed for computational chemistry workflows.

## Contents

- [What is ChemVCS?](#what-is-chemvcs)
- [Current Status](#current-status)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Remote HPC (No SSH)](#remote-hpc-no-ssh)
- [Documentation](#documentation)

## What is ChemVCS?

ChemVCS brings Git-like version control to computational chemistry. It helps you track molecular structures, calculation parameters, and computational results with the same rigor as source code version control.

Unlike general-purpose VCS systems, ChemVCS is built with an **extensible architecture** that will support chemistry-specific features like molecular structure diffing, computational provenance tracking, and HPC job correlation.

## Why ChemVCS?

Computational chemistry generates complex data across long-running calculations:
- **Reproducibility**: Track exact molecular structures, calculation parameters, and software versions
- **Collaboration**: Share computational workflows and results with team members
- **Provenance**: Maintain complete history of how results were obtained
- **HPC Integration**: Track job submission, status, and retrieval (local or via a remote gateway)

Traditional VCS tools like Git work for some of this, but ChemVCS is designed from the ground up to handle chemistry-specific needs.

## Current Status

**Version**: Milestone 6 Complete (All Phases)  
**Stability**: Full HPC integration with multi-scheduler support  
**Test Coverage**: 85 Go tests + 157 Python tests (242 total) passing  
**Production Ready**: Not yet - under active development

ChemVCS currently provides:
- ✅ **Git-like foundation**: Content-addressable storage, commits, branches, merges
- ✅ **Remote repositories**: Push, pull, fetch via HTTP
- ✅ **Python domain layer**: Chemistry-specific objects (Structure, Run, Workflow)
- ✅ **File format support**: XYZ and POSCAR (VASP) parsers
- ✅ **HPC integration**: Multi-scheduler support (SLURM, PBS, LSF)
- ✅ **HPC CLI commands**: Full job lifecycle (submit, track, retrieve, cancel, watch)
- ✅ **Interactive monitoring**: Real-time job status updates

See [FEATURE_STATUS.md](FEATURE_STATUS.md) for detailed feature tracking.

## Installation

### Prerequisites
- **Go 1.19+** (for core VCS)
- **Python 3.8+** (for chemistry domain layer)
- **NumPy 1.20+** (Python dependency)

### Build from Source

#### 1. Build Go CLI and Server

```bash
# Clone the repository
git clone https://github.com/lichman0405/chemvcs.git
cd chemvcs/go

# Build the CLI client
go build -o chemvcs ./cmd/chemvcs

# Build the server (optional)
go build -o chemvcs-server ./cmd/chemvcs-server
```

On Windows, executables will be `chemvcs.exe` and `chemvcs-server.exe`.

#### 2. Install Python Package (Optional)

```bash
cd ../python
pip install -e .
```

This installs the `chemvcs_py` package for working with chemistry domain objects.

## Quick Start

### 1. Initialize a Repository

```bash
# Create a new ChemVCS repository
chemvcs init myproject
cd myproject
```

This creates a `.chemvcs` directory to store version history.

### 2. Add and Commit Files

```bash
# Add your files
echo "Initial structure" > molecule.xyz

# Commit the snapshot
chemvcs commit -m "Add initial structure"
```

**Note**: ChemVCS has no staging area - all changes in your working directory are committed together.

### 3. View History

```bash
# See commit log
chemvcs log

# Compact view
chemvcs log -oneline

# Limit to last 5 commits
chemvcs log -n 5
```

### 4. Create Branches

```bash
# List branches
chemvcs branch

# Create a new branch
chemvcs branch optimization

# Switch to it
chemvcs checkout optimization

# Make changes and commit
echo "Optimized structure" > molecule.xyz
chemvcs commit -m "Geometry optimization"
```

### 5. Merge Branches

```bash
# Switch back to main
chemvcs checkout main

# Merge changes
chemvcs merge optimization
```

ChemVCS supports both fast-forward and three-way merging with conflict detection.

## Remote HPC (No SSH)

If your laptop cannot run SLURM commands locally, you can submit/track jobs via `chemvcs-server` acting as an HTTP gateway on the SLURM host.

```bash
# Add a repo-scoped remote (recommended)
chemvcs remote add slurm http://<host>:<port>/chemvcs/v1/repos/<owner>/<repo>

# If the server requires auth, provide a token (global or per-remote)
export CHEMVCS_REMOTE_TOKEN="<token>"
# or: export CHEMVCS_REMOTE_TOKEN_SLURM="<token>"

# Submit / monitor / retrieve via the gateway
chemvcs submit   --remote=slurm <run-hash> job.slurm
chemvcs jobs     --remote=slurm
chemvcs watch    --remote=slurm <run-hash|job-id> --interval=10
chemvcs cancel   --remote=slurm <run-hash|job-id>
chemvcs retrieve --remote=slurm <run-hash> --patterns="*.out,*.log" --dest=./results
```

For details, see [docs/11-remote-hpc-design.md](docs/11-remote-hpc-design.md).

## Using Python Domain Layer

### Working with Molecular Structures

```python
from chemvcs_py import Structure, Repo
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
from chemvcs_py import Run, Workflow

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

# Build workflow DAG
workflow = Workflow(name="Water Study")
workflow.add_node("opt", "run", {"description": "Optimize"})
workflow.add_node("freq", "run", {"description": "Frequencies"})
workflow.add_edge("opt", "freq")  # freq depends on opt

# Get execution order
order = workflow.topological_sort()  # ['opt', 'freq']
```

See [python/examples/](python/examples/) for complete examples.

## Working with Remote Repositories

### Start a Server

On a central machine:

```bash
chemvcs-server -port 8080 -repo-root /path/to/repos \
    --auth-token "<token>" \
    --auth-repos "owner/repo" \
    --admin-token "<admin-token>"
```

### Clone, Push, Pull

```bash
# Clone from server
chemvcs clone http://server:8080/project myproject
cd myproject

# Make changes
echo "New data" > results.txt
chemvcs commit -m "Add results"

# Push to server
chemvcs push origin main

# Pull updates from others
chemvcs pull origin main
```

## Command Reference

### Repository Operations

| Command | Description |
|---------|-------------|
| `chemvcs init [path]` | Initialize a new repository |
| `chemvcs commit -m "msg"` | Commit current state |
| `chemvcs log [-n N] [-oneline]` | View commit history |
| `chemvcs status` | Show working directory status |

### Branching

| Command | Description |
|---------|-------------|
| `chemvcs branch` | List branches |
| `chemvcs branch <name>` | Create a branch |
| `chemvcs checkout <branch>` | Switch branches |
| `chemvcs merge <branch>` | Merge a branch |

### HPC Integration (M6)

| Command | Description |
|---------|-------------|
| `chemvcs submit <run> <script>` | Submit HPC job for a run |
| `chemvcs jobs [--status]` | List tracked HPC jobs |
| `chemvcs retrieve <run> [--patterns]` | Fetch completed results |
| `chemvcs submit/jobs/retrieve/cancel/watch --remote=<name>` | Use remote HPC gateway (SLURM) |
| Python: `JobSubmitter.submit_run()` | Submit HPC job with provenance (API) |
| Python: `JobTracker.check_status()` | Query job status (API) |
| Python: `JobRetriever.retrieve_results()` | Fetch completed results (API) |

See [docs/10-hpc-user-guide.md](docs/10-hpc-user-guide.md) for complete HPC documentation.

### Remote Operations

| Command | Description |
|---------|-------------|
| `chemvcs clone <url> [path]` | Clone from remote |
| `chemvcs remote add <name> <url>` | Add remote |
| `chemvcs push <remote> <branch>` | Push to remote |
| `chemvcs pull <remote> <branch>` | Pull from remote |
| `chemvcs fetch <remote> <branch>` | Fetch without merging |

## Key Concepts

### No Staging Area

Unlike Git, ChemVCS has **no staging area** (no `git add`). When you commit, all modified files are included. This reflects the reality that computational calculations produce atomic, complete results - not partial changes to be staged incrementally.

### Snapshot-Based

Each commit is a **complete snapshot** of your working directory, not a diff. This ensures you can always retrieve an exact calculation state.

### Content-Addressable

All data is stored by SHA-256 hash. Identical files are automatically deduplicated.

## Comparison with Git

| Feature | Git | ChemVCS |
|---------|-----|---------|
| Staging area | Yes (`git add`) | No (automatic) |
| Rebase/reset | Yes | No (preserves history) |
| Object model | Fixed (blob/tree/commit) | Extensible |
| Domain features | None | Planned (chemistry) |

See [docs/10-chemvcs-vs-git.md](docs/10-chemvcs-vs-git.md) for detailed comparison.

## What's Missing?

ChemVCS is **under development**. Key remaining work:

- ✅ **M6 Phase 1-2 complete** - HPC core infrastructure + CLI commands
- ⏭️ **M6 Phase 3-5** - Additional adapters (PBS, LSF) and enhanced features
- ⏭️ **Additional HPC adapters** - PBS, LSF full support
- ⏭️ **CIF file parser** - Crystallography format
- ⏭️ **Enhanced Repository API** - Advanced queries and filtering
- ❌ **Hooks system** - Workflow automation
- ❌ **Submodules** - Dependency management
- ❌ **Tag support** - Version markers

See [FEATURE_STATUS.md](FEATURE_STATUS.md) for detailed feature tracking and roadmap notes.

## Documentation

- **[FEATURE_STATUS.md](FEATURE_STATUS.md)** - Comprehensive feature tracking (updated with M6!)
- [docs/](docs/) - Detailed design specifications
  - [01-vision-and-scope.md](docs/01-vision-and-scope.md) - Project vision
  - [02-architecture-overview.md](docs/02-architecture-overview.md) - System design
  - [06-python-domain-layer.md](docs/06-python-domain-layer.md) - Python package design
  - [09-hpc-integration-design.md](docs/09-hpc-integration-design.md) - HPC integration design 🆕
  - [10-hpc-user-guide.md](docs/10-hpc-user-guide.md) - HPC user guide 🆕
    - [11-remote-hpc-design.md](docs/11-remote-hpc-design.md) - Remote HPC gateway (no SSH) 🆕
- [python/README.md](python/README.md) - Python package documentation
- [examples/hpc-workflow/](examples/hpc-workflow/) - Complete HPC workflow example 🆕

## Contributing

ChemVCS is in early development. Contributions welcome, but expect frequent changes.

### Running Tests

**Go tests:**
```bash
cd go
go test ./...
```

**Python tests:**
```bash
cd python
python -m pytest tests/ -v
```

### Project Structure

```
chemvcs/
├── go/                       # Go core VCS implementation
│   ├── cmd/
│   │   ├── chemvcs/          # CLI client
│   │   └── chemvcs-server/   # HTTP server
│   └── internal/
│       ├── model/            # Core data structures
│       ├── objectstore/      # Content-addressable storage
│       ├── repo/             # Repository operations
│       ├── workspace/        # Working directory management
│       ├── remote/           # Remote repository client
│       └── server/           # HTTP server implementation
├── python/                   # Python domain layer
│   ├── chemvcs_py/
│   │   ├── core/            # Repository and CoreObject
│   │   ├── domain/          # Chemistry objects (Structure, Run, Workflow)
│   │   ├── io/              # File parsers (XYZ, POSCAR)
│   │   └── util/            # Error handling utilities
│   ├── examples/            # Usage examples
│   │   ├── basic_usage.py
│   │   └── run_workflow_example.py
│   └── tests/               # Python test suite (28 tests)
├── docs/                    # Design specifications
│   ├── 01-vision-and-scope.md
│   ├── 02-architecture-overview.md
│   ├── 06-python-domain-layer.md
│   └── ...
├── README.md                # This file
├── FEATURE_STATUS.md        # Detailed feature tracking
└── CHANGELOG.md             # Project history
```

## License

[To be determined]

## Contact

- GitHub: [lichman0405/chemvcs](https://github.com/lichman0405/chemvcs)
- Issues: [GitHub Issues](https://github.com/lichman0405/chemvcs/issues)
