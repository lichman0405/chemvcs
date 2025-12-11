# ChemVCS

A domain-aware, distributed version control system for computational chemistry and materials science projects.

## Project Status

**Current Phase:** Milestone 2 - Working Directory and Status ✅ **COMPLETE**

### Completed
- ✅ Phase 0: Project foundation
  - Go module initialised
  - Directory structure created
  - CI/CD ready
- ✅ Milestone 1: Local VCS Core MVP
  - Core data structures (Blob, Object, Snapshot, Reference)
  - Content-addressable storage with sharding
  - Repository operations (init, commit, log, branch, checkout)
  - CLI with 6 commands
  - 38 tests passing, 65-77% coverage
- ✅ Milestone 2: Working directory and status
  - Workspace scanning for recursive directory trees
  - Folder and file object creation with blob references
  - Diff engine for comparing object trees
  - Status command showing added/modified/deleted files
  - Working directory restoration on checkout
  - Real file tracking in commits
  - **Total: 51 tests passing, 65-76% coverage**

### In Progress
- 🚧 Milestone 3: Branches and merge

### Planned
- Milestone 3: Branch merging and conflict resolution
- Milestone 4: Remote repositories
- Milestone 5: Python domain layer
- Milestone 6: HPC integration

## Documentation

Comprehensive design documentation is available in the `docs/` directory:

- [01-vision-and-scope.md](../docs/01-vision-and-scope.md) - Project vision and goals
- [02-architecture-overview.md](../docs/02-architecture-overview.md) - System architecture
- [03-object-model-and-storage.md](../docs/03-object-model-and-storage.md) - Data model specification
- [04-cli-spec-mvp.md](../docs/04-cli-spec-mvp.md) - CLI specification
- [05-remote-protocol-and-server.md](../docs/05-remote-protocol-and-server.md) - Remote sync protocol
- [06-python-domain-layer.md](../docs/06-python-domain-layer.md) - Python domain layer design
- [07-hpc-adapter-design.md](../docs/07-hpc-adapter-design.md) - HPC scheduler integration
- [08-testing-and-quality.md](../docs/08-testing-and-quality.md) - Testing strategy
- [09-execution-plan.md](../docs/09-execution-plan.md) - Detailed implementation plan

## Development

### Prerequisites

- Go 1.19 or higher
- Git

### Building

```bash
cd go
go build -o chemvcs ./cmd/chemvcs
```

On Windows:
```bash
go build -o chemvcs.exe ./cmd/chemvcs
```

### Running

The `chemvcs` executable provides the following commands:

```bash
# Initialise a new repository
chemvcs init [path]

# Create a snapshot (commit)
chemvcs commit -m "Your message"

# View history
chemvcs log [-n <count>] [-oneline]

# Manage branches
chemvcs branch              # List branches
chemvcs branch <name>       # Create branch

# Switch branches or snapshots
chemvcs checkout <target>

# Show working directory status
chemvcs status

# Show version
chemvcs version
```

Set author information via environment variables:
```bash
export CHEMVCS_AUTHOR_NAME="Your Name"
export CHEMVCS_AUTHOR_EMAIL="your.email@example.com"
```

Or on Windows:
```powershell
$env:CHEMVCS_AUTHOR_NAME="Your Name"
$env:CHEMVCS_AUTHOR_EMAIL="your.email@example.com"
```

### Testing

Run all tests:
```bash
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

### Project Structure

```
go/
├── cmd/
│   └── chemvcs/          # CLI entry point
├── internal/
│   ├── model/            # Core data structures
│   ├── objectstore/      # Content-addressable storage
│   ├── repo/             # Repository operations
│   └── workspace/        # Working directory mapping
└── go.mod
```

## Design Principles

1. **Content-addressable and immutable** - All entities are addressed by cryptographic hashes
2. **Domain-neutral core** - Generic VCS with chemistry-specific extensions
3. **Provenance first** - Capture complete computational lineage
4. **HPC-friendly** - Designed for scientific computing environments
5. **Incremental adoption** - Start local, scale to distributed

## Licence

[To be determined]

## Authors

- Development team (see CONTRIBUTORS.md)

## Acknowledgements

This project draws inspiration from Git's elegant design whilst addressing the specific needs of computational science workflows.
