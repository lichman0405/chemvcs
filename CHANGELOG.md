# ChemVCS Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Remote HPC gateway over HTTP (no SSH): server-side SLURM execution via repo-scoped `/hpc/*` endpoints
- Remote client + CLI support for HPC commands via `--remote=<name>`
- `chemvcs retrieve --commit` (optional snapshot after pulling outputs)
- `chemvcs jobs [<run-hash|job-id>]` filtering
- Bearer token authentication for `chemvcs-server` (optional) with repo scoping and an optional admin token
- Per-remote client token support via `CHEMVCS_REMOTE_TOKEN[_<REMOTE>]`

### Security
- Hardened HPC result retrieval to prevent path traversal and escaping the repo root
- Safe ZIP extraction on the client side to prevent zip-slip style writes outside destination
- Enforced optimistic concurrency for ref updates (updates must provide expected `old_target`)

### Changed
- Remote HPC retrieve now enforces file-count/size limits and uses command timeouts for SLURM calls
- Remote HPC retrieve streams to disk on the client (avoids buffering the entire zip in memory)
- Server `/hpc/jobs` list includes `updated_at` and `working_dir`

### Documentation
- Updated docs and README to describe Remote HPC (no SSH) usage and remove stale references


### Planned
- Job arrays and dependencies
- Batch workflow submission
- Web dashboard for job monitoring
- CIF file parser (crystallography format)
- Enhanced Repository query API
- Hooks system
- Tag support

## [0.6.0] - 2025-12 (M6: HPC Integration - All Phases Complete)

### Added - Milestone 6 Phase 1-2: Full HPC Integration

#### Phase 1: HPC Core Infrastructure (Python)

#### HPC Module (`chemvcs_py.hpc`)

**JobAdapter Interface** (`adapter.py` - 110 lines)
- Abstract base class for HPC scheduler adapters
- `JobStatus` enum: PENDING, RUNNING, COMPLETED, FAILED, CANCELLED, TIMEOUT, UNKNOWN
- `JobInfo` dataclass: comprehensive job information (nodes, state, runtime, etc.)
- Methods: `submit()`, `get_status()`, `get_info()`, `cancel()`

**SlurmAdapter** (`slurm_adapter.py` - 220 lines)
- Complete SLURM workload manager integration
- `submit()`: Parse sbatch output, extract job ID
- `get_status()`: Query via squeue, fallback to sacct for completed jobs
- `get_info()`: Extract detailed job information (nodes, walltime, memory, etc.)
- `cancel()`: Job cancellation via scancel
- Comprehensive status mapping from SLURM states
- Timeout and error handling

**Exception Hierarchy** (`exceptions.py`)
- `HpcError`: Base exception for all HPC operations
- `JobSubmissionError`: Job submission failures
- `JobNotFoundError`: Missing or purged jobs
- `JobCancellationError`: Cancellation failures
- `JobInfoError`: Info query failures

**Provenance Capture** (`provenance.py` - 180 lines)
- `capture_modules()`: Parse `module list -t` output for loaded modules
- `capture_env_vars()`: Extract relevant environment variables (OMP_NUM_THREADS, SLURM_*, etc.)
- `parse_slurm_script()`: Extract #SBATCH directives (nodes, walltime, memory, etc.)
- `parse_pbs_script()`: Extract #PBS directives
- `detect_script_type()`: Auto-detect scheduler from script content
- Resource extraction: nodes, ntasks, walltime, memory
- Graceful handling of missing commands (FileNotFoundError, TimeoutExpired)

**Job Submission** (`submission.py` - 80 lines)
- `JobSubmitter` class for run submission
- `submit_run()`: Submit with automatic provenance capture
- Environment capture on submission (modules, env vars, resources)
- Run status update after submission
- Integration with Repository for persistence

**Job Tracking** (`tracking.py` - 130 lines)
- `JobTracker` class for status monitoring
- `check_status()`: Query current job state
- `wait_for_completion()`: Blocking wait with configurable polling interval
- `update_run_status()`: Sync Run object with job state
- Automatic Run persistence after status updates

**Result Retrieval** (`retrieval.py` - 140 lines)
- `JobRetriever` class for output fetching
- `retrieve_results()`: Fetch files from working directory
- Pattern-based file filtering (glob patterns)
- Automatic Run status update to RETRIEVED
- Support for custom destination paths
- Optional stdout/stderr capture

#### Extended Run Domain Object

**New HPC Fields** (`domain/run.py`)
- `job_id`: Cluster job identifier (e.g., "12345")
- `job_system`: Scheduler type ("slurm", "pbs", "lsf", etc.)
- `queue_name`: Submission queue/partition name
- `submit_script`: Full submission script snapshot
- `modules_loaded`: List of environment modules with versions
- `environment_vars`: Dict of relevant environment variables
- `job_resources`: Dict of allocated resources (nodes, walltime, memory)

**New Methods**
- `mark_queued(job_id, job_system, queue_name)`: Update status after job queuing
- `mark_retrieved()`: Mark results as retrieved from cluster

**Enhanced Serialization**
- `to_core_object()`: Serialize all HPC fields to VCS storage
- `from_core_object()`: Deserialize with backward compatibility

#### Comprehensive Testing (+45 tests, 118 total)

**Adapter Tests** (`test_slurm_adapter.py` - 21 tests)
- Submit success/failure scenarios
- Status queries for all job states (PENDING/RUNNING/COMPLETED/FAILED/etc.)
- Job info parsing and validation
- Cancel operations
- Edge cases: timeouts, purged jobs, invalid output parsing
- Mock-based using `unittest.mock.patch('subprocess.run')`

**Provenance Tests** (`test_provenance.py` - 16 tests)
- Module capture and parsing (empty, single, multiple modules)
- Environment variable extraction (OMP_NUM_THREADS, SLURM_CPUS_PER_TASK, etc.)
- SLURM script parsing: nodes, ntasks, walltime, memory
- PBS script parsing: equivalent directives
- Script type detection (SLURM vs PBS)
- Missing command handling (FileNotFoundError)

**Run HPC Tests** (`test_run.py` - 8 new tests in TestRunHPC class)
- HPC field validation (job_id, job_system, modules, env_vars, resources)
- `mark_queued()` method with all fields
- `mark_retrieved()` method
- Serialization round-trip with HPC data
- Complete HPC workflow: planned → queued → retrieved

**Testing Strategy**
- Mock-based: No real SLURM cluster required
- Subprocess mocking for all CLI commands (sbatch, squeue, sacct, scancel)
- Comprehensive edge case coverage
- All 118 Python tests passing

#### Documentation and Examples

**Design Document** (`docs/09-hpc-integration-design.md` - 53KB)
- Motivation and use cases
- Architecture overview with component diagrams
- Data model extensions
- JobAdapter interface specification
- CLI commands design (for Phase 2)
- Provenance tracking strategy
- Security considerations
- Testing strategy
- Complete implementation checklist

**User Guide** (`docs/10-hpc-user-guide.md` - 18KB)
- Installation and prerequisites
- Quick start guide
- Complete command reference (Python API)
- Provenance tracking examples
- Advanced usage: custom adapters, polling strategies
- Troubleshooting common issues
- Best practices for HPC workflows
- MockAdapter for testing without SLURM
- FAQ section

**Example Workflow** (`examples/hpc-workflow/`)
- `README.md`: Step-by-step tutorial (9KB)
- `water.xyz`: H2O molecular structure
- `vasp_relax.slurm`: VASP geometry optimization script
- `vasp_static.slurm`: VASP single-point calculation script
- `workflow.py`: Complete Python workflow demonstration (210 lines)
- Demonstrates: initialization, submission, tracking, retrieval, provenance viewing
- Multi-step workflow: relaxation → static calculation

#### Phase 2: Go CLI Commands

**Go-Python Integration Layer** (`go/internal/hpc/hpc.go` - 440 lines)
- `findPythonExecutable()`: Locate Python interpreter (python3 or python)
- `findChemVCSPyModule()`: Find chemvcs_py module in repository structure
- `runPythonScript()`: Execute Python scripts with proper PYTHONPATH
- `SubmitJob()`: Submit HPC job by calling Python JobSubmitter
- `ListJobs()`: List tracked jobs from repository
- `CheckJobStatus()`: Query job status via Python JobTracker
- `RetrieveResults()`: Fetch completed results via Python JobRetriever
- Data structures: JobStatus, JobInfo, SubmitOptions, RetrieveOptions

**CLI Command Handlers** (`go/cmd/chemvcs/main.go`)
- `chemvcs submit <run-hash> <script> [--capture-env]`: Submit HPC job
- `chemvcs jobs [--status=<status>] [-v]`: List tracked jobs
- `chemvcs retrieve <run-hash> [--patterns] [--dest]`: Retrieve results
- Updated help message with HPC Commands section
- User-friendly output formatting

**CLI Integration Tests** (`go/internal/hpc/hpc_test.go` - 330 lines)
- TestHPCPackageStructure: Data structure validation
- TestPythonAvailability: Python interpreter detection
- TestCLIIntegration: End-to-end CLI command testing
- TestMockWorkflow: Mock workflow demonstration
- All 85 Go tests passing (80 pre-M6 + 5 M6)

**Documentation Updates**
- Updated `docs/10-hpc-user-guide.md` with CLI commands section
- Added Quick Start with CLI and Python API options
- Complete CLI command reference with examples and output
- Side-by-side comparison of CLI vs Python API usage

### Changed

**Run Domain Object Enhancements**
- Added 7 HPC-related fields
- Added 2 new lifecycle methods (mark_queued, mark_retrieved)
- Updated serialization for backward compatibility
- Enhanced RunStatus enum with RETRIEVED state

**Python Package Structure**
- New `hpc/` subpackage (7 modules)
- Updated `__init__.py` exports
- Enhanced test suite organization

**CLI Help Message**
- Added "HPC Commands:" section
- Submit, jobs, and retrieve command descriptions

**User Guide Organization**
- Quick Start now shows CLI commands first (recommended approach)
- Python API section remains for advanced usage

### Code Statistics (M6 Phase 1-2 Combined)

**New Python Code** (Phase 1): ~2,800 lines
- Production: ~1,100 lines (7 HPC modules)
- Tests: ~1,700 lines (3 test files, 45 tests)

**New Go Code** (Phase 2): ~770 lines
- Production: ~440 lines (1 HPC module)
- Tests: ~330 lines (1 test file, 5 test suites)

**Total Codebase**:
- Python: ~8,400 lines (4,600 production + 3,800 tests)
- Go: ~7,120 lines (3,952 production + 3,168 tests)
- Combined: ~15,520 lines

**Documentation**: +85KB
- Design document: 53KB
- User guide: 32KB (updated)

**Total Tests**: 203 passing (Phase 1-2)
- Go: 85 tests (M1-M4: 80 + M6: 5)
- Python: 118 tests (M5: 73 + M6: 45)

### Git Commits (Phase 1-2)
- Phase 1: `979a40e` - "feat(M6): Phase 1 - HPC core infrastructure complete"
- Phase 1 docs: `24ea0b1` - "docs: Add M6 Phase 1 documentation"
- Phase 2: `2cfd282` - "feat(M6): Phase 2 - Go CLI commands for HPC integration"
- Phase 2 docs: `8cbb376` - "docs: Update user guide with CLI commands"

#### Phase 3: Additional HPC Scheduler Adapters

**Features**:
- PBS/Torque adapter (`pbs_adapter.py`)
  - Command integration: qsub, qstat, qdel
  - Status mappings: Q/W/H→PENDING, R/E→RUNNING, C→COMPLETED
- LSF adapter (`lsf_adapter.py`)
  - Command integration: bsub, bjobs, bkill
  - LSF-specific: slots/cores model, shell redirection for bsub
  - Status mappings: PEND/WAIT→PENDING, RUN→RUNNING, DONE→COMPLETED
- Extended JobInfo dataclass with `cores` field
- Updated exports in hpc/__init__.py
- Comprehensive test coverage (18 PBS + 21 LSF tests)

**Code Statistics** (Phase 3): +421 lines Python
- PbsAdapter: 209 lines
- LsfAdapter: 212 lines
- Tests: 39 new tests (all passing)

#### Phase 4: Enhanced CLI Features (Job Management)

**Features**:
- **Job Cancellation**
  - JobTracker.cancel_job() and cancel_run() methods
  - JobCancellationError exception
  - Go CancelJob() function (~70 lines)
  - CLI command: `chemvcs cancel <identifier>`
  - Updates Run status to "cancelled"
  
- **Interactive Job Monitoring**
  - JobWatcher class (195 lines, monitoring.py)
  - watch_job() for single job, watch_multiple() for batch monitoring
  - Real-time status updates with timestamp logging
  - Emoji status indicators: ⏳ PENDING, ⏸️ SUSPENDED, ✅ COMPLETED, ❌ FAILED, 🚫 CANCELLED
  - Configurable polling interval (default 30s) and optional timeout
  - Go WatchJob() function (~50 lines)
  - CLI command: `chemvcs watch <identifier> [--interval] [--timeout]`

**Code Statistics** (Phase 4): ~405 lines total
- Python: ~265 lines (monitoring.py, tracking.py, exceptions.py updates)
- Go: ~140 lines (CancelJob and WatchJob functions, CLI handlers)
- Tests: All 242 tests passing (85 Go + 157 Python)

### Git Commits (Phase 3-4)
- Phase 3: `7aee1d8` - "feat(M6): Phase 3 - PBS and LSF adapter implementations"
- Phase 4a: `fb4c3bc` - "feat(M6): Phase 4a - Job cancellation functionality"
- Phase 4b: `61f6323` - "feat(M6): Phase 4b - Interactive job monitoring"
- Documentation: `17c1c32` - "docs(M6): Update PROJECT_STATUS with Phase 3-5 completion"

**M6 Complete Statistics**:
- Python HPC: ~3,700 lines (9 modules: adapter, slurm_adapter, pbs_adapter, lsf_adapter, tracking, monitoring, validation, template, exceptions)
- Go HPC: ~510 lines (hpc.go with 5 operations)
- Total Tests: 242 passing (85 Go + 157 Python)
- CLI Commands: 5 HPC operations (submit, jobs, retrieve, cancel, watch)
- HPC Adapters: 3 schedulers (SLURM, PBS, LSF)

**Total Python Codebase**: ~8,400 lines
- Production: ~4,600 lines (pre-M6: 1,800 + M6: 2,800)
- Tests: ~3,800 lines (pre-M6: 1,100 + M6: 2,700)

**Documentation**: +71KB
- Design document: 53KB
- User guide: 18KB

**Total Tests**: 198 passing
- Go: 80 tests (M1-M4)
- Python: 118 tests (M5: 73 + M6: 45)

### Assessment: Essential Innovation Achieved

**Before M6**: ChemVCS was a well-implemented Git clone with chemistry data structures. No essential differentiation.

**After M6 Phase 1**: ChemVCS achieves essential innovation through:
1. **HPC Job Lifecycle Integration**: Links version-controlled structures to cluster computations
2. **Computational Provenance**: Captures environment modules, scheduler directives, resources
3. **Reproducibility Beyond Code**: Tracks not just *what* was calculated but *how* it was computed
4. **Chemistry-Specific Workflow**: End-to-end tracking from molecule → calculation → cluster job → results

**Key Differentiation**: Generic VCS (Git) tracks file changes. ChemVCS tracks the full computational lifecycle of chemistry calculations with complete provenance for reproducibility.

---

## [0.5.0] - 2025-12-15

### Added - Milestone 5: Python Domain Layer (100% Core Complete)

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
- **Python test suite** (1,100+ lines, 73 tests)
  - test_structure.py: 10 tests (creation, validation, operations, conversion)
  - test_run.py: 18 tests (lifecycle, status transitions, result management, conversion)
  - test_workflow.py: 27 tests (node/edge, cycle detection, topological sort, conversion)
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
- **FEATURE_STATUS.md** (280 lines)
  - Comprehensive feature tracking document
  - M1-M4: 100% complete
  - M5: 100% core complete (optional extensions remaining)
  - Clear remaining work items and priorities
  - Fixed duplicate content issues

- **STRUCTURE.md** (250+ lines)
  - Complete directory layout with line counts
  - Statistics (Go: 3,512 + 2,837 lines, Python: 1,800 + 1,100 lines)
  - Key files reference
  - Development tools and commands
  - Dependencies listing

- **Updated README.md**
  - Python installation instructions
  - Python usage examples (Structure, Run, Workflow)
  - Current status updated to M5 100% core complete
  - Test coverage: 80 Go + 73 Python tests
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
