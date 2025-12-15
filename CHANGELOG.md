# ChemVCS Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Planned
- CIF file parser (crystallography format)
- HPC integration (M6: SLURM adapter, job tracking)
- Enhanced Repository query API
- Hooks system
- Tag support

## [0.5.0] - 2025-12-15

### Added - Milestone 5: Python Domain Layer (85% Complete)

#### Python Package (`chemvcs_py`)
- **Structure domain object** (195 lines)
  - Molecular and crystal structure representation
  - NumPy array storage for atomic positions
  - Support for periodic (lattice) and non-periodic systems
  - Operations: translate(), get_center_of_mass()
  - Validation for positions (Nx3), species, and lattice (3x3)
  - CoreObject conversion for VCS storage

- **Run domain object** (215 lines)
  - Computational calculation representation
  - Lifecycle tracking: planned → submitted → running → finished/failed/cancelled
  - Status transition methods with timestamps
  - Result storage and retrieval (get_energy, set_result)
  - Link to Structure via structure_id reference

- **Workflow domain object** (318 lines)
  - Directed Acyclic Graph (DAG) for computational workflows
  - WorkflowNode dataclass (id, kind, metadata)
  - Cycle detection via depth-first search
  - Graph operations: add_node(), add_edge() with validation
  - Dependency queries: get_dependencies(), get_dependents()
  - Execution order: topological_sort() using Kahn's algorithm
  - Root/leaf node identification

- **XYZ file parser** (167 lines)
  - read_xyz(): Parse molecular XYZ files
  - write_xyz(): Write XYZ with optional custom comment
  - Automatic formula generation from species
  - Comprehensive error handling and validation

- **POSCAR file parser** (286 lines)
  - read_poscar(): VASP 5.x format support
  - Support for Direct (fractional) and Cartesian coordinates
  - Scaling factor and lattice vector handling
  - write_poscar(): Write in both coordinate modes
  - Element symbol and count parsing

- **Repository API** (228 lines)
  - Repo class for repository discovery
  - CLI integration via subprocess
  - Methods: list_objects(), get_object(), commit()
  - Automatic chemvcs executable discovery

- **Core object models** (109 lines)
  - CoreObject dataclass: version, type, meta, refs, hash
  - Reference dataclass: kind, id
  - JSON serialization/deserialization

#### Testing
- **Python test suite** (582 lines, 28 tests)
  - test_structure.py: 10 tests (creation, validation, operations, conversion)
  - test_xyz.py: 8 tests (read, write, errors, round-trip)
  - test_poscar.py: 10 tests (Direct/Cartesian, validation, round-trip)
  - pytest configuration with fixtures and temp files
  - 100% test pass rate

#### Examples
- **basic_usage.py** (126 lines)
  - Structure creation and manipulation
  - XYZ file I/O demonstration
  - CoreObject conversion
  - Repository integration

- **run_workflow_example.py** (197 lines)
  - Run object lifecycle demonstration
  - Workflow DAG construction (simple and complex)
  - Dependency analysis and topological sorting
  - Cycle detection showcase
  - Multi-conformer workflow example (10 nodes)

#### Go Core Enhancements
- **JSON API support**
  - `chemvcs inspect-object <hash> --format=json`
  - `chemvcs list-objects --type=<type> --format=json`
  - Store.ListObjects() method with type filtering
  - ObjectInfo struct (Hash, Type)
  - Fixed: empty results return [] instead of null

#### Documentation
- **FEATURE_STATUS.md** (251 lines)
  - Comprehensive feature tracking document
  - M1-M4: 100% complete
  - M5: 85% complete with detailed breakdown
  - Clear remaining work items

- **STRUCTURE.md** (250+ lines)
  - Complete directory layout with line counts
  - Statistics (Go: 3,512 + 2,837 lines, Python: 1,800 + 1,100 lines)
  - Key files reference
  - Development tools and commands
  - Dependencies listing

- **Updated README.md**
  - Python installation instructions
  - Python usage examples (Structure, Run, Workflow)
  - Current status updated to M5 85%
  - Links to comprehensive documentation

- **Updated python/README.md**
  - Expanded from 107 to 250+ lines
  - Complete API reference
  - Multiple quick start examples
  - Development guide

- **Updated PROJECT_STATUS.md**
  - M5 features and progress
  - Architecture with Python components
  - Test coverage statistics

- **Updated TODO.md**
  - M5 tasks marked complete (15 items)
  - Remaining lower-priority tasks listed

### Changed
- **Structure.translate()**: Changed from returning new Structure to in-place modification
- **Structure.to_core_object()**: Now includes num_atoms and is_periodic in meta
- **write_xyz()**: Added optional comment parameter

### Fixed
- **list-objects JSON**: Fixed null return for empty results (now returns [])
- **Structure.from_core_object()**: Properly excludes num_atoms and is_periodic from metadata

### Statistics
- **Go codebase**: 3,512 production + 2,837 test lines (80 tests passing)
- **Python codebase**: 1,800 production + 1,100 test lines (28 tests passing)
- **Documentation**: 15,000+ words across multiple documents
- **Total commits**: 100+ (f8e341f, c459c18, 8318858 for M5)

## [0.4.0] - 2024-12-11

### Added - Milestone 4: Remote Repositories
- Remote repository support via HTTP
- Client commands: push, pull, fetch, clone
- Server implementation: chemvcs-server
- Remote management: remote add, remote list
- Complete remote synchronization protocol
- 72 tests passing

### Added - Quality Improvements
- Checkout cleanup: Automatic removal of files not in target snapshot
- Three-way merge: Intelligent merging for diverged branches
- Conflict detection: Clear reporting of conflicting changes
- 80 tests passing

## [0.3.0] - 2024-12-08

### Added - Milestone 3: Merging
- Fast-forward merge support
- Merge conflict detection
- Branch divergence handling
- 57 tests passing

## [0.2.0] - 2024-12-05

### Added - Milestone 2: Working Directory
- Working directory management
- File scanning and status tracking
- Restore operations
- 51 tests passing

## [0.1.0] - 2024-12-01

### Added - Milestone 1: Core VCS
- Content-addressable storage (SHA-256)
- Object model: blobs, objects, snapshots
- Commit history with parent tracking
- Branch management (create, list, switch)
- Reference management (HEAD, branches)
- 38 tests passing

### Initial Features
- `chemvcs init` - Initialize repository
- `chemvcs commit` - Create commits
- `chemvcs log` - View history
- `chemvcs branch` - Manage branches
- `chemvcs checkout` - Switch branches

---

## Release Notes

### Version 0.5.0 (Current)
This release completes 85% of Milestone 5, adding a comprehensive Python domain layer for chemistry-specific operations. The package provides high-level abstractions for molecular structures, computational calculations, and workflow management, with full integration to the Go-based VCS core.

Key highlights:
- Three domain objects: Structure, Run, Workflow
- Two file format parsers: XYZ, POSCAR
- 28 comprehensive tests (100% passing)
- Complete documentation and examples
- Ready for chemistry workflow integration

### Version 0.4.0
Added complete remote repository support with push/pull/fetch/clone operations. Includes quality improvements for checkout cleanup and three-way merging.

### Version 0.3.0
Introduced merge capabilities with fast-forward support and conflict detection.

### Version 0.2.0
Added working directory management and status tracking.

### Version 0.1.0
Initial release with core VCS functionality: commits, branches, and history tracking.
