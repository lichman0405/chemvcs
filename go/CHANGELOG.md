# Changelog

All notable changes to this project will be documented in this file.

## [M5 - In Progress] - 2025-12-15

### Added

- **JSON API for Python Integration**
  - `chemvcs inspect-object <hash> [--format=json]` - Inspect object details with optional JSON output
  - `chemvcs list-objects [--type=<type>] [--format=json]` - List all objects with optional type filtering and JSON output
  - `Store.ListObjects(typeFilter string)` method for programmatic object listing
  - Foundation for Python layer to query ChemVCS objects via CLI

### Example Usage

```bash
# List all objects
chemvcs list-objects

# List only file objects
chemvcs list-objects --type=file

# Get JSON output for Python parsing
chemvcs list-objects --format=json

# Inspect a specific object
chemvcs inspect-object <hash>

# Get object as JSON
chemvcs inspect-object --format=json <hash>
```

## [M4 Complete] - 2025-12-11

### Added

- Remote repository support via HTTP
- `chemvcs-server` binary for hosting repositories
- Push, pull, fetch, clone operations
- Three-way merge algorithm with conflict detection
- Checkout cleanup (removes files not in target snapshot)

### Changed

- Documentation restructure: README.md is now user-facing, PROJECT_STATUS.md tracks development

## [M3 Complete] - 2025-12-10

### Added

- Fast-forward merge support
- Branch management commands

## [M2 Complete] - 2025-12-09

### Added

- Working directory management
- Status tracking
- Checkout command

## [M1 Complete] - 2025-12-08

### Added

- Local VCS core with content-addressable storage
- Basic commit, log, branch commands
- SHA-256 hashed objects
