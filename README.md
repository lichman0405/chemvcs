# ChemVCS

A version control system designed for computational chemistry workflows.

## What is ChemVCS?

ChemVCS brings Git-like version control to computational chemistry. It helps you track molecular structures, calculation parameters, and computational results with the same rigor as source code version control.

Unlike general-purpose VCS systems, ChemVCS is built with an **extensible architecture** that will support chemistry-specific features like molecular structure diffing, computational provenance tracking, and HPC job correlation.

## Why ChemVCS?

Computational chemistry generates complex data across long-running calculations:
- **Reproducibility**: Track exact molecular structures, calculation parameters, and software versions
- **Collaboration**: Share computational workflows and results with team members
- **Provenance**: Maintain complete history of how results were obtained
- **HPC Integration**: (Planned) Correlate version control with batch job systems

Traditional VCS tools like Git work for some of this, but ChemVCS is designed from the ground up to handle chemistry-specific needs.

## Current Status

**Version**: Milestone 4 Complete  
**Stability**: Core features working, 80 tests passing  
**Production Ready**: Not yet - under active development

ChemVCS currently provides a solid Git-like foundation with:
- Content-addressable storage
- Full commit history and branching
- Three-way merge with conflict detection
- Remote repositories via HTTP

**Chemistry-specific features** (M5-M6) are planned but not yet implemented. See [PROJECT_STATUS.md](PROJECT_STATUS.md) for detailed progress.

## Installation

### Prerequisites
- Go 1.19 or higher

### Build from Source

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

## Working with Remote Repositories

### Start a Server

On a central machine:

```bash
chemvcs-server -port 8080 -repo-root /path/to/repos
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

ChemVCS is **under development**. Currently missing:

- ❌ Chemistry-specific features (molecular diffs, structure tracking)
- ❌ HPC integration (job tracking, SLURM correlation)
- ❌ Hooks system
- ❌ Submodules
- ❌ Tag support (planned)

See [TODO.md](TODO.md) for the roadmap.

## Documentation

- [PROJECT_STATUS.md](PROJECT_STATUS.md) - Development progress and milestones
- [TODO.md](TODO.md) - Roadmap and planned features
- [docs/](docs/) - Detailed design specifications
  - [01-vision-and-scope.md](docs/01-vision-and-scope.md) - Project vision
  - [02-architecture-overview.md](docs/02-architecture-overview.md) - System design
  - [10-chemvcs-vs-git.md](docs/10-chemvcs-vs-git.md) - ChemVCS vs Git comparison

## Contributing

ChemVCS is in early development. Contributions welcome, but expect frequent changes.

### Running Tests

```bash
cd go
go test ./...
```

### Project Structure

```
chemvcs/
├── go/
│   ├── cmd/
│   │   ├── chemvcs/          # CLI client
│   │   └── chemvcs-server/   # HTTP server
│   └── internal/
│       ├── model/            # Core data structures
│       ├── objectstore/      # Content-addressable storage
│       ├── repo/             # Repository operations
│       ├── workspace/        # Working directory management
│       ├── remote/           # Remote client
│       └── server/           # HTTP server
└── docs/                      # Design documentation
```

## License

[To be determined]

## Contact

- GitHub: [lichman0405/chemvcs](https://github.com/lichman0405/chemvcs)
- Issues: [GitHub Issues](https://github.com/lichman0405/chemvcs/issues)
