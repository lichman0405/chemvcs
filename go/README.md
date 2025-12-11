# ChemVCS

A domain-aware, distributed version control system for computational chemistry and materials science projects.

## Project Status

**Current Phase:** Milestone 1 - Local VCS Core MVP

### Completed
- ✅ Phase 0: Project foundation
  - Go module initialised
  - Directory structure created
  - CI/CD ready
- ✅ Milestone 1, Step 1: Model layer
  - Core data structures (Blob, Object, Snapshot, Reference)
  - Deterministic hashing with canonical JSON serialisation
  - Comprehensive unit tests (9/9 passing)

### In Progress
- 🚧 Milestone 1, Step 2: ObjectStore layer

### Planned
- Milestone 1, Step 3: Repository layer
- Milestone 1, Step 4: CLI implementation
- Milestone 2: Working directory and status
- Milestone 3: Branches and merge

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
go build ./cmd/chemvcs
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
