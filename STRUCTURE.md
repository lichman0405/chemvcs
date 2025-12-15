# ChemVCS Project Structure

Complete overview of the ChemVCS repository organization.

## Directory Layout

```
chemvcs/
├── go/                           # Go implementation (core VCS)
│   ├── cmd/
│   │   ├── chemvcs/              # Main CLI client
│   │   │   └── main.go           # Entry point + command handlers
│   │   └── chemvcs-server/       # HTTP server
│   │       └── main.go           # Server implementation
│   └── internal/
│       ├── model/                # Core data structures
│       │   ├── object.go         # Object, Snapshot, Reference
│       │   ├── object_test.go
│       │   └── ...
│       ├── objectstore/          # Content-addressable storage
│       │   ├── store.go          # Storage operations
│       │   ├── store_test.go
│       │   └── ...
│       ├── repo/                 # Repository operations
│       │   ├── repo.go           # Init, commit, branch, merge
│       │   ├── repo_test.go
│       │   └── ...
│       ├── workspace/            # Working directory management
│       │   ├── workspace.go      # Scan, restore, status
│       │   └── workspace_test.go
│       ├── remote/               # Remote repository client
│       │   ├── client.go         # Push, pull, fetch operations
│       │   └── client_test.go
│       └── server/               # HTTP server
│           ├── server.go         # Request handlers
│           └── server_test.go
│
├── python/                       # Python domain layer
│   ├── chemvcs_py/               # Main package
│   │   ├── __init__.py           # Package exports
│   │   ├── core/                 # Core VCS integration
│   │   │   ├── __init__.py
│   │   │   ├── repo.py           # Repository API (CLI integration)
│   │   │   └── objects.py        # CoreObject, Reference models
│   │   ├── domain/               # Chemistry domain objects
│   │   │   ├── __init__.py
│   │   │   ├── structure.py      # Molecular/crystal structures (195 lines)
│   │   │   ├── run.py            # Computational calculations (215 lines)
│   │   │   └── workflow.py       # Workflow DAGs (318 lines)
│   │   ├── io/                   # File format parsers
│   │   │   ├── __init__.py
│   │   │   ├── xyz.py            # XYZ molecular format (167 lines)
│   │   │   └── poscar.py         # VASP POSCAR format (286 lines)
│   │   └── util/                 # Utilities
│   │       ├── __init__.py
│   │       └── errors.py         # Custom exceptions
│   ├── examples/                 # Usage examples
│   │   ├── basic_usage.py        # Structure + XYZ demo (126 lines)
│   │   └── run_workflow_example.py # Run + Workflow demo (197 lines)
│   ├── tests/                    # Python test suite
│   │   ├── __init__.py
│   │   ├── conftest.py           # pytest configuration
│   │   ├── test_structure.py     # Structure tests (195 lines, 10 tests)
│   │   ├── test_xyz.py           # XYZ parser tests (186 lines, 8 tests)
│   │   └── test_poscar.py        # POSCAR parser tests (201 lines, 10 tests)
│   ├── setup.py                  # Package configuration
│   ├── requirements.txt          # Python dependencies
│   └── README.md                 # Python package documentation
│
├── docs/                         # Design specifications
│   ├── 01-vision-and-scope.md    # Project vision and goals
│   ├── 02-architecture-overview.md # System architecture
│   ├── 03-object-model-and-storage.md # Object storage design
│   ├── 04-cli-spec-mvp.md        # CLI command specifications
│   ├── 05-remote-protocol-and-server.md # Remote protocol design
│   ├── 06-python-domain-layer.md # Python package design
│   ├── 07-hpc-adapter-design.md  # HPC integration (planned)
│   ├── 08-testing-and-quality.md # Testing strategy
│   ├── 10-chemvcs-vs-git.md      # Comparison with Git
│   └── chemvcs-dev-plan.md       # Development roadmap
│
├── README.md                     # Main project README
├── FEATURE_STATUS.md             # Detailed feature tracking
├── PROJECT_STATUS.md             # Milestone progress
├── TODO.md                       # Development roadmap
└── .gitignore                    # Git ignore rules

```

## Statistics

### Go Codebase
- **Production code**: ~3,512 lines
- **Test code**: ~2,837 lines
- **Test coverage**: 80 tests passing
- **Packages**: 7 (model, objectstore, repo, workspace, remote, server, + cmd)

### Python Codebase
- **Production code**: ~1,800 lines
- **Test code**: ~1,100 lines
- **Test coverage**: 28 tests passing
- **Packages**: 4 (core, domain, io, util)

### Documentation
- **Design docs**: 9 markdown files
- **Examples**: 2 Python scripts
- **Total docs**: ~15,000 words

## Key Files

### Go Entry Points
- `go/cmd/chemvcs/main.go` - CLI client (862 lines)
- `go/cmd/chemvcs-server/main.go` - HTTP server (138 lines)

### Core Go Packages
- `go/internal/model/object.go` - Data structures (237 lines)
- `go/internal/objectstore/store.go` - Storage layer (452 lines)
- `go/internal/repo/repo.go` - Repository operations (588 lines)
- `go/internal/workspace/workspace.go` - Working directory (334 lines)

### Python Core
- `python/chemvcs_py/core/repo.py` - Repository API (228 lines)
- `python/chemvcs_py/domain/structure.py` - Structures (195 lines)
- `python/chemvcs_py/domain/run.py` - Calculations (215 lines)
- `python/chemvcs_py/domain/workflow.py` - Workflows (318 lines)

### Python I/O
- `python/chemvcs_py/io/xyz.py` - XYZ parser (167 lines)
- `python/chemvcs_py/io/poscar.py` - POSCAR parser (286 lines)

## Build Artifacts

### Executables (not in repo)
- `go/chemvcs` (or `chemvcs.exe` on Windows)
- `go/chemvcs-server` (or `chemvcs-server.exe`)

### Python Installation
- Installed via `pip install -e python/`
- Creates `chemvcs_py.egg-info/` directory

### Repository Data (not in repo)
- `.chemvcs/` directory in user projects
  - `objects/` - Content-addressable storage
  - `refs/` - Branch and HEAD references
  - `config` - Repository configuration

## Development Tools

### Testing
```bash
# Go tests
cd go && go test ./...

# Python tests
cd python && python -m pytest tests/ -v
```

### Building
```bash
# Build CLI
cd go && go build -o chemvcs ./cmd/chemvcs

# Build server
cd go && go build -o chemvcs-server ./cmd/chemvcs-server

# Install Python package
cd python && pip install -e .
```

### Code Quality
```bash
# Go formatting
cd go && go fmt ./...

# Go linting (requires golangci-lint)
cd go && golangci-lint run

# Python formatting (requires black)
cd python && black chemvcs_py/ tests/

# Python linting (requires pylint)
cd python && pylint chemvcs_py/
```

## Dependencies

### Go Dependencies
- Standard library only
- No external dependencies

### Python Dependencies
- **numpy** (>=1.20.0) - Array operations for molecular structures
- **pytest** (dev) - Testing framework

## Git Workflow

### Branches
- `main` - Stable development branch
- Feature branches as needed

### Commit History
- Initial commit: Go core implementation
- Recent: Python domain layer (M5)
- Latest: f8e341f - Complete M5 with Run, Workflow, tests

## Future Additions (Planned)

### M6: HPC Integration
```
go/internal/hpc/          # HPC adapter interface
python/chemvcs_py/hpc/    # Python HPC client
```

### Additional File Formats
```
python/chemvcs_py/io/
├── cif.py                # CIF (crystallography)
├── gaussian.py           # Gaussian output
└── vasp_output.py        # VASP calculation results
```

### Enhanced Features
```
go/internal/
├── diff/                 # Molecular structure diff
└── hooks/                # Hook system

python/chemvcs_py/
├── analysis/             # Data analysis tools
└── visualization/        # Plotting utilities
```
