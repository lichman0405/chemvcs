# ChemVCS - Computational Chemistry Version Control System

ChemVCS is a domain-specific version control system designed for computational chemistry workflows. It provides Git-like version control optimized for managing complex molecular structures, calculation outputs, and computational workflows with **integrated HPC job management**.

## Features

### Core VCS Features (Milestone 1-4) ✅
- **Content-addressable storage**: SHA-256 hashed objects
- **Commit history**: Full snapshot-based versioning with parent tracking
- **Branching**: Create, switch, and manage multiple branches
- **Merging**: Fast-forward and three-way merge with conflict detection
- **Remote repositories**: Push, pull, fetch operations via HTTP
- **Working directory management**: Clean checkout with automatic file cleanup

### Python Domain Layer (Milestone 5) ✅ Complete
- **Chemistry domain objects**: Structure, Run (calculations), Workflow (DAGs)
- **File format parsers**: XYZ (molecular), POSCAR (VASP)
- **Repository integration**: Python API for ChemVCS operations
- **Comprehensive tests**: 73 Python tests covering all components
- **Example scripts**: Demonstrating molecular structures and computational workflows

### HPC Integration (Milestone 6 - All Phases) ✅ Complete
- **Multi-scheduler support**: SLURM, PBS/Torque, LSF adapters
- **Full job lifecycle**: Submit, track, retrieve, cancel, monitor
- **Provenance capture**: Automatic tracking of modules, environment, resources
- **CLI commands**: `chemvcs submit/jobs/retrieve/cancel/watch` for HPC operations
- **Python API**: JobSubmitter, JobTracker, JobRetriever, JobWatcher classes
- **Interactive monitoring**: Real-time job status updates with `watch` command
- **Complete testing**: 84 Python HPC tests + 5 Go CLI tests
- **Documentation**: Comprehensive user guide with examples

## Status

**Current Version**: Milestone 6 Complete (All Phases)

**Test Coverage**: 
- Go: 85 tests passing across 8 packages
- Python: 157 tests passing across 10 test files
- **Total: 242 automated tests**

### Completed Milestones
- ✅ **M1**: Local VCS core (objects, snapshots, refs) - 38 tests
- ✅ **M2**: Working directory and status tracking - 51 tests
- ✅ **M3**: Fast-forward merge - 57 tests
- ✅ **M4**: Remote repositories (client + HTTP server) - 72 tests
- ✅ **M5**: Python domain layer (Structure, Run, Workflow, parsers) - 73 tests
- ✅ **M6 Phase 1**: HPC core infrastructure (Python) - 45 tests
- ✅ **M6 Phase 2**: HPC CLI commands (Go) - 5 tests
- ✅ **M6 Phase 3**: Additional adapters (PBS, LSF) - 39 tests
- ✅ **M6 Phase 4**: Enhanced CLI (cancel, watch) - 0 new tests (integrated)

### Recent Additions (M6 Phase 1-4)

**M6 Phase 1: Python HPC Core**
1. ✅ **Extended Run object** - 7 HPC fields (job_id, modules, resources, etc.)
2. ✅ **JobAdapter interface** - Abstract base for scheduler adapters
3. ✅ **SlurmAdapter** - Complete SLURM integration (220 lines)
4. ✅ **Provenance capture** - Modules, env vars, script parsing (180 lines)
5. ✅ **JobSubmitter/Tracker/Retriever** - High-level HPC API (350 lines)
6. ✅ **HPC test suite** - 45 tests with mock-based testing
7. ✅ **Design document** - 53KB comprehensive specification
8. ✅ **Example workflow** - Complete VASP workflow demonstration

**M6 Phase 2: Go CLI Commands**
1. ✅ **Go-Python integration** - Python executor with PYTHONPATH management
2. ✅ **chemvcs submit** - Submit HPC jobs from CLI
3. ✅ **chemvcs jobs** - List and filter tracked jobs
4. ✅ **chemvcs retrieve** - Fetch completed results
5. ✅ **CLI tests** - 5 test suites for CLI integration
6. ✅ **Updated user guide** - CLI command reference and examples

**M6 Phase 3: Additional Adapters**
1. ✅ **PbsAdapter** - PBS/Torque integration (209 lines)
2. ✅ **LsfAdapter** - LSF integration (212 lines)
3. ✅ **Comprehensive tests** - 18 PBS tests + 21 LSF tests (39 total)
4. ✅ **Extended JobInfo** - Added cores field for LSF
5. ✅ **Updated documentation** - PBS/LSF examples in user guide

**M6 Phase 4: Enhanced CLI Features**
1. ✅ **Job cancellation** - `chemvcs cancel` command
   - Cancel by run hash or job ID
   - Automatic adapter detection
   - Updates Run status to cancelled
2. ✅ **Interactive monitoring** - `chemvcs watch` command
   - Real-time status updates with timestamps
   - Configurable polling interval
   - Optional timeout support
   - Clean emoji indicators (⏳⏸️✅❌🚫)
3. ✅ **JobWatcher class** - Python monitoring API (195 lines)
4. ✅ **Updated help messages** - Complete CLI reference

### Optional Future Extensions
- ⏭️ Job arrays and dependencies
- ⏭️ Batch workflow submission
- ⏭️ Web dashboard for job monitoring
- ⏭️ CIF parser (crystallographic format)
- ⏭️ Enhanced Repository API (advanced queries)

## Architecture

### Core Components (Go)
- **model**: Object, Snapshot, Reference data structures
- **objectstore**: Content-addressable storage (blobs and objects)
- **repo**: Repository operations (commit, branch, merge)
- **workspace**: Working directory scanning and restoration
- **remote**: Client for remote repositories
- **server**: HTTP server implementing ChemVCS protocol

### Python Domain Layer
- **core**: Repository API and CoreObject models
- **domain**: Chemistry objects (Structure, Run, Workflow)
- **io**: File parsers (XYZ, POSCAR, future CIF)
- **util**: Error handling and validation utilities

### Binaries
- **chemvcs**: Main CLI tool for local and remote operations
- **chemvcs-server**: HTTP server for hosting remote repositories

### Python Package
- **chemvcs_py**: Python library for chemistry domain objects and file I/O

## Installation

### Prerequisites
- Go 1.19 or higher
- Git

### Building

Build the client:
```bash
cd gLocal Repository Operations

#### Initialize Repository
```bash
chemvcs init [path]
```

#### Commit Changes
```bash
chemvcs commit -m "Your commit message"
```

#### View History
```bash
chemvcs log [-n <count>] [-oneline]
```

#### Check Status
```bash
chemvcs status
```

### Branching and Merging

```bash
# List branches
chemvcs branch

# Create branch
chemvcs branch feature-x

# Switch branch
chemvcs checkout feature-x

# Merge branch (supports fast-forward and three-way merge)
chemvcs merge feature-x
```

### Remote Operations

```bash
# Add remote
chemvcs remote add origin http://server:8080/repo

# Push changes
chemvcs push origin main

# Pull changes
chemvcs pull origin main

# Fetch objects without updating local ref
chemvcs fetch origin main
```

### HPC Operations (M6)

```bash
# Submit HPC job
chemvcs submit <run-hash> script.slurm

# List tracked jobs
chemvcs jobs
chemvcs jobs --status=RUNNING

# Retrieve completed results
chemvcs retrieve <run-hash>
chemvcs retrieve <run-hash> --patterns="*.out,*.log"
```

### Start Server

```bash
# Default port 8080
chemvcs-server

# Custom configuration
chemvcs-server -port 9000 -repo-root /path/to/repos
```

On Windows:
```powershell
.\chemvcs-server.exe -port 9000 -repo-root C:\chemvcs\repos
```

## Testing

### Go Tests

Run tests with coverage:
```bash
go test -cover ./...
```

### Go Tests

Run all tests:
```bash
cd go
go test ./...
# 85 tests passing
```

Run with verbose output:
```bash
go test ./... -v
```

Run tests for a specific package:
```bash
go test ./internal/hpc -v
```

### Python Tests

Run all tests:
```bash
cd python
pytest tests/ -v
# 118 tests passing
```

Run with coverage:
```bash
pytest tests/ --cov=chemvcs_py --cov-report=html
```

Run specific test file:
```bash
pytest tests/test_run.py -v
```

## Project Structure

```
chemvcs/
├── go/
│   ├── cmd/
│   │   ├── chemvcs/          # CLI entry point
│   │   └── chemvcs-server/   # HTTP server entry point
│   ├── internal/
│   │   ├── model/            # Core data structures (Object, Snapshot, Reference)
│   │   ├── objectstore/      # Content-addressable storage with sharding
│   │   ├── repo/             # Repository operations (commit, branch, merge)
│   │   ├── workspace/        # Working directory scanning and restoration
│   │   ├── remote/           # Remote client (push/pull/fetch)
│   │   ├── server/           # HTTP server implementation
│   │   └── hpc/              # HPC integration (Go-Python bridge)
│   └── go.mod
├── python/
│   ├── chemvcs_py/
│   │   ├── core/             # Repository and CoreObject interfaces
│   │   ├── domain/           # Structure, Run, Workflow objects
│   │   ├── io/               # File parsers (XYZ, POSCAR)
│   │   ├── util/             # Utilities and error handling
│   │   └── hpc/              # HPC integration (adapters, submission, tracking)
│   ├── tests/                # Python test suite
│   └── examples/             # Usage examples
├── examples/
│   └── hpc-workflow/         # Complete HPC workflow example
├── docs/                      # Design documentation
│   ├── 01-vision-and-scope.md
│   ├── 02-architecture-overview.md
│   ├── 03-object-model-and-storage.md
│   ├── 04-cli-spec-mvp.md
│   ├── 05-remote-protocol-and-server.md
│   ├── 06-python-domain-layer.md
│   ├── 07-hpc-adapter-design.md
│   ├── 08-testing-and-quality.md
│   ├── 09-hpc-integration-design.md
│   └── 10-hpc-user-guide.md
├── README.md
├── FEATURE_STATUS.md
├── PROJECT_STATUS.md
├── TODO.md
└── CHANGELOG.md
```

## Development

### Running Tests
```bash
cd go
go test ./...
```

### Test Coverage by Package
- **model**: Object and Snapshot core types
- **objectstore**: Storage and retrieval
- **repo**: Repository operations and merge logic
- **workspace**: Directory scanning and restoration
- **remote**: Client protocol implementation
- **server**: HTTP server endpoints

## ChemVCS vs Git

While ChemVCS shares core concepts with Git (content-addressable storage, Merkle DAG, snapshots), it diverges in important ways to better serve computational science:

### Key Differences

**1. Extensible Object Model** ⭐
- **Git**: Fixed types (blob, tree, commit, tag)
- **ChemVCS**: Flexible `Object` with extensible `Type` field
  ```go
  type Object struct {
      Type string                 // "file", "folder", or future: "structure", "calculation"
      Meta map[string]interface{} // Domain-specific metadata (energy, method, etc.)
      Refs []Reference            // Flexible references
  }
  ```
- **Benefit**: Can represent chemistry-specific entities natively in future milestones

**2. No Staging Area (Index)**
- **Git**: Three-stage workflow (working → staging → repository)
- **ChemVCS**: Two-stage workflow (working → repository)
- **Rationale**: Scientific calculations are atomic operations; partial staging contradicts "complete calculation state" model
- **UX**: Simpler mental model for researchers

**3. Terminology Aligned with Science**
- **Git**: "commit" (suggests code changes)
- **ChemVCS**: "snapshot" (emphasizes capturing computational state)
- **Philosophy**: Reflects that each version is a complete scientific result, not just incremental changes

**4. Deliberately Excluded Features**
- ❌ **Rebase/Cherry-pick**: Violates scientific provenance principles
  - Science requires authentic history, not beautified commits
  - Real calculation order matters for reproducibility
- ❌ **Interactive staging**: No staging area in ChemVCS
- ❌ **Reset modes**: Unnecessary complexity without staging

**5. Planned Chemistry-Specific Features** (M5+)
- Domain-aware diff for molecular structures (RMSD, graph comparison)
- Automatic provenance capture (compute environment, resources)
- HPC job tracking and result correlation
- Chemistry file format support (XYZ, CIF, quantum chemistry outputs)

### Current State

**Today (M1-M4)**: ChemVCS is a "cleaner Git" with extensible architecture
- ✅ Simpler workflow (no staging)
- ✅ Architecture ready for domain extensions
- ⚠️ Chemistry features not yet implemented

**Tomorrow (M5-M6)**: True differentiation emerges
- Chemistry-specific object types
- Computational provenance tracking
- HPC workflow integration

See [TODO.md](TODO.md) for detailed roadmap and design decisions.

---

## Design Principles

1. **Content-addressable**: All objects identified by SHA-256 hash
2. **Immutable objects**: Once stored, objects never change
3. **Snapshot-based**: Each commit captures full state
4. **Domain-agnostic core**: Chemistry semantics in domain layer
5. **HTTP protocol**: Simple REST-like API for remotes
6. **Scientific integrity**: Authentic history over convenient rewrites

## Documentation

Detailed specifications in `docs/`:
- `01-vision-and-scope.md` - Project goals and scope
- `02-architecture-overview.md` - System architecture
- `03-object-model-and-storage.md` - Core data model
- `04-cli-spec-mvp.md` - CLI command specifications
- `05-remote-protocol-and-server.md` - HTTP protocol spec
- `06-python-domain-layer.md` - Python integration plan
- `07-hpc-adapter-design.md` - HPC integration plan
- `08-testing-and-quality.md` - Testing strategy

## License

[Your license here]

## Contributing

[Your contribution guidelines here]
