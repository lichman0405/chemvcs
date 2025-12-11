# ChemVCS - Computational Chemistry Version Control System

ChemVCS is a domain-specific version control system designed for computational chemistry workflows. It provides Git-like version control optimized for managing complex molecular structures, calculation outputs, and computational workflows.

## Features

### Core VCS Features (Milestone 1-4) ✅
- **Content-addressable storage**: SHA-256 hashed objects
- **Commit history**: Full snapshot-based versioning with parent tracking
- **Branching**: Create, switch, and manage multiple branches
- **Merging**: Fast-forward and three-way merge with conflict detection
- **Remote repositories**: Push, pull, fetch operations via HTTP
- **Working directory management**: Clean checkout with automatic file cleanup

### Recent Improvements
- ✅ **Checkout cleanup**: Automatically removes files not in target snapshot
- ✅ **Three-way merge**: Full support for diverged branches with intelligent merging
- ✅ **Conflict detection**: Identifies conflicting changes with clear reporting

## Status

**Current Version**: Milestone 4 Complete + Quality Improvements

**Test Coverage**: 80 tests passing across 7 packages

### Completed Milestones
- ✅ **M1**: Local VCS core (objects, snapshots, refs) - 38 tests
- ✅ **M2**: Working directory and status tracking - 51 tests
- ✅ **M3**: Fast-forward merge - 57 tests
- ✅ **M4**: Remote repositories (client + HTTP server) - 72 tests
- ✅ **Quality improvements**: Checkout cleanup, three-way merge, conflict detection - 80 tests

### Recent Quality Improvements (Post-M4)
1. ✅ **Checkout cleanup** - Files from previous branches are now automatically removed
2. ✅ **Three-way merge** - Diverged branches can be merged automatically with intelligent merging
3. ✅ **Conflict detection** - Conflicting changes are detected and reported clearly

### Upcoming Milestones
- **M5**: Python domain layer (chemistry-specific semantics)
- **M6**: HPC adapter (SLURM integration, job tracking)

## Architecture

### Core Components
- **model**: Object, Snapshot, Reference data structures
- **objectstore**: Content-addressable storage (blobs and objects)
- **repo**: Repository operations (commit, branch, merge)
- **workspace**: Working directory scanning and restoration
- **remote**: Client for remote repositories
- **server**: HTTP server implementing ChemVCS protocol

### Binaries
- **chemvcs**: Main CLI tool for local and remote operations
- **chemvcs-server**: HTTP server for hosting remote repositories

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

### Start Server

```bash
# Default port 8080
chemvcs-server

# Custom configuration
chemvcs-server -port 9000 -repo-root /path/to/repos
```

On Windows:
```powershell
.\chemvcs-server.exe -port 9000 -repo-root C:\chemvcs\
chemvcs commit -m "Your commit message"
```

### Branching
```bash

Run all tests:
```bash
cd go
go test ./...
```

Run tests with coverage:
```bash
go test -cover ./...
```

Run tests for a specific package:
```bash
go test ./internal/model -v
```

Run tests matching a pattern:
```bash
go test -run TestMerge ./internal/repo
```

### Project Structure

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
│   │   └── server/           # HTTP server implementation
│   └── go.mod
├── docs/                      # Design documentation
│   ├── 01-vision-and-scope.md
│   ├── 02-architecture-overview.md
│   ├── 03-object-model-and-storage.md
│   ├── 04-cli-spec-mvp.md
│   ├── 05-remote-protocol-and-server.md
│   ├── 06-python-domain-layer.md
│   ├── 07-hpc-adapter-design.md
│   └── 08-testing-and-quality.md
└── README.md
```
# Add remote
chemvcs remote add origin http://server:8080/repo

# Push changes
chemvcs push origin main

# Pull changes
chemvcs pull origin main
```

### Start Server
```bash
chemvcs-server -port 8080 -repo-root /path/to/repos
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

## Design Principles

1. **Content-addressable**: All objects identified by SHA-256 hash
2. **Immutable objects**: Once stored, objects never change
3. **Snapshot-based**: Each commit captures full state
4. **Domain-agnostic core**: Chemistry semantics in domain layer
5. **HTTP protocol**: Simple REST-like API for remotes

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
