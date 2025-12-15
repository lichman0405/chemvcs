# ChemVCS TODO List

# ChemVCS TODO List

## Milestone 5: Python Domain Layer (M5) - ✅ Complete

### Completed Tasks ✅
- ✅ Design Python package structure (`chemvcs_py`)
- ✅ Implement Python bindings to Go core (via JSON CLI output)
- ✅ Define chemistry-specific object types
  - ✅ Structure objects (molecular geometries, unit cells)
  - ✅ Run objects (calculation tasks with lifecycle tracking)
  - ✅ Workflow objects (DAG of computations with dependencies)
- ✅ Implement Python API for common operations
  - ✅ Repository operations (list_objects, get_object, commit)
  - ✅ Domain object conversions (to/from CoreObject)
  - ✅ File I/O (XYZ, POSCAR parsers)
- ✅ Add chemistry-aware parsers for common file formats
  - ✅ XYZ files (read/write with validation)
  - ✅ POSCAR files (VASP format, Direct/Cartesian coords)
- ✅ Write Python tests and examples
  - ✅ 73 pytest tests (Structure, Run, Workflow, XYZ, POSCAR)
  - ✅ basic_usage.py example
  - ✅ run_workflow_example.py with DAG demonstration
- ✅ Document Python API
  - ✅ Python package README
  - ✅ Design doc (06-python-domain-layer.md)
- ✅ Clean up documentation
  - ✅ Remove duplicate content from FEATURE_STATUS.md
  - ✅ Update test counts and completion status

### Optional Extensions (Lower Priority)
- ⏭️ CIF file parser (crystallography format)
- ⏭️ Enhanced diff/merge for molecular structures
- ⏭️ Advanced Repository API features (relationship tracking)

---

## Milestone 6: HPC Integration (M6) - Phase 1-2 Complete ✅

**Design Document**: docs/09-hpc-integration-design.md  
**User Guide**: docs/10-hpc-user-guide.md

### Phase 1: Core Infrastructure ✅ Complete

- [x] Create M6 design document (53KB comprehensive spec)
- [x] Extend Run object with HPC fields
  - [x] Add job_id, job_system, queue_name fields
  - [x] Add provenance fields (modules_loaded, environment_vars, job_resources)
  - [x] Add submit_script snapshot field
  - [x] Add mark_queued(), mark_retrieved() methods
  - [x] Update to_core_object() / from_core_object()
- [x] Add Run HPC tests (8 tests in TestRunHPC)
  - [x] Test HPC field validation
  - [x] Test mark_queued() method
  - [x] Test mark_retrieved() method
  - [x] Test round-trip serialization with HPC data
- [x] Design JobAdapter interface
  - [x] Define abstract base class
  - [x] Define JobStatus enum (PENDING/RUNNING/COMPLETED/etc.)
  - [x] Define JobInfo dataclass
  - [x] Define exception hierarchy (HpcError, JobSubmissionError, etc.)
- [x] Implement SlurmAdapter (220 lines)
  - [x] Implement submit() with sbatch
  - [x] Implement get_status() with squeue/sacct fallback
  - [x] Implement get_info() with sacct parsing
  - [x] Implement cancel() with scancel
  - [x] Add status mapping logic
- [x] Implement provenance capture (180 lines)
  - [x] capture_modules() - Parse module list
  - [x] capture_env_vars() - Extract OMP, SLURM vars
  - [x] parse_slurm_script() - Extract #SBATCH directives
  - [x] parse_pbs_script() - Extract #PBS directives
  - [x] detect_script_type() - Auto-detect scheduler
- [x] Implement high-level HPC API
  - [x] JobSubmitter class (80 lines)
  - [x] JobTracker class (130 lines)
  - [x] JobRetriever class (140 lines)
- [x] Add adapter unit tests (45 tests total, mock-based)
  - [x] Test submit() success/failure (21 adapter tests)
  - [x] Test get_status() for all states
  - [x] Test get_info() parsing
  - [x] Test cancel() operation
  - [x] Test provenance capture (16 tests)
  - [x] Test error handling
- [x] Documentation
  - [x] Design document (docs/09-hpc-integration-design.md)
  - [x] User guide (docs/10-hpc-user-guide.md - 18KB)
  - [x] Example workflow (examples/hpc-workflow/)

**Phase 1 Statistics**:
- New Python code: ~2,800 lines
- New tests: 45 (118 total Python tests)
- Documentation: +71KB

### Phase 2: CLI Integration ✅ Complete

- [x] Go: Implement Go-Python integration layer (440 lines)
  - [x] findPythonExecutable() - Locate python3/python
  - [x] findChemVCSPyModule() - Find module in repo
  - [x] runPythonScript() - Execute with PYTHONPATH
  - [x] SubmitJob(), ListJobs(), CheckJobStatus(), RetrieveResults()
- [x] Go: Add `chemvcs submit` command
  - [x] Parse arguments (run-hash, script path)
  - [x] Validate Run object exists
  - [x] Call Python HPC module
  - [x] Handle errors gracefully
- [x] Go: Add `chemvcs jobs` command
  - [x] List Run objects with job_id
  - [x] Query status for each job
  - [x] Format output table
  - [x] Add filtering options (--status, -v)
- [x] Go: Add `chemvcs retrieve` command
  - [x] Parse run-hash argument
  - [x] Find corresponding Run object
  - [x] Copy output files with patterns
  - [x] Support destination directory
- [x] Add CLI integration tests (5 test suites, 330 lines)
  - [x] Package structure tests
  - [x] Python availability detection
  - [x] CLI help message validation
  - [x] Command error handling tests
  - [x] Mock workflow demonstrations
- [x] Update documentation
  - [x] Add CLI commands to user guide
  - [x] Complete command reference
  - [x] Update FEATURE_STATUS, README, CHANGELOG

**Phase 2 Statistics**:
- New Go code: ~770 lines
- New tests: 5 (85 total Go tests)
- Total tests: 203 (85 Go + 118 Python)

### Phase 3: Additional HPC Adapters ✅ Complete

- [x] PBS/Torque adapter (pbs_adapter.py - 209 lines)
  - [x] qsub command for job submission
  - [x] qstat command for status queries
  - [x] qdel command for job cancellation
  - [x] Status mapping (Q/W/H → PENDING, R/E → RUNNING, C → COMPLETED)
  - [x] Job info parsing (queue, nodes, start time, exit code)
- [x] LSF adapter (lsf_adapter.py - 212 lines)
  - [x] bsub command for job submission
  - [x] bjobs command for status queries
  - [x] bkill command for job cancellation
  - [x] Status mapping (PEND/WAIT → PENDING, RUN → RUNNING, DONE → COMPLETED, EXIT → FAILED)
  - [x] Slots/cores support (LSF-specific)
- [x] Extended JobInfo dataclass
  - [x] Added cores field for LSF
- [x] Comprehensive test suites
  - [x] 18 PBS tests (submit, status, info, cancel)
  - [x] 21 LSF tests (submit, status, info, cancel)
  - [x] All mock-based, no real scheduler required
- [x] Update documentation
  - [x] Add PBS/LSF to supported schedulers
  - [x] Add PBS/LSF script examples
  - [x] Update user guide

**Phase 3 Statistics**:
- New Python code: ~421 lines (209 PBS + 212 LSF)
- New tests: 39 (18 PBS + 21 LSF)
- Total tests: 242 (85 Go + 157 Python)

### Phase 4: Enhanced CLI Features ✅ Complete

- [x] Job cancellation command
  - [x] JobTracker.cancel_job() method (by job ID)
  - [x] JobTracker.cancel_run() method (by run hash)
  - [x] JobCancellationError exception
  - [x] Go: CancelJob() function
  - [x] CLI: `chemvcs cancel <run-hash|job-id>`
  - [x] Automatic adapter detection
  - [x] Update Run status to cancelled
- [x] Interactive monitoring command
  - [x] JobWatcher class (monitoring.py - 195 lines)
  - [x] watch_job() method - Monitor until completion
  - [x] Real-time status updates with timestamps
  - [x] Configurable polling interval (default: 30s)
  - [x] Optional timeout support
  - [x] Clean emoji indicators (⏳⏸️✅❌🚫)
  - [x] Go: WatchJob() function
  - [x] CLI: `chemvcs watch <identifier> [--interval] [--timeout]`
- [x] Update CLI help messages
  - [x] Add cancel command
  - [x] Add watch command
  - [x] Complete command reference

**Phase 4 Statistics**:
- New Python code: ~265 lines (70 cancel + 195 monitoring)
- New Go code: ~140 lines (cancel + watch functions)
- Total tests: 242 (stable - integrated into existing tests)

### M6 Complete Summary

**Total Implementation**:
- Python HPC module: ~3,700 lines (9 modules)
- Go HPC integration: ~510 lines
- Tests: 242 passing (85 Go + 157 Python)
- CLI commands: 5 HPC commands (submit, jobs, retrieve, cancel, watch)
- Adapters: 3 schedulers (SLURM, PBS, LSF)
- Documentation: 85KB+ (design + user guide + examples)

**Git Commits**:
- Phase 1: `979a40e`, `24ea0b1`
- Phase 2: `2cfd282`, `8cbb376`
- Phase 3: `7aee1d8`
- Phase 4: `fb4c3bc`, `61f6323`
- Docs: `d0b6b76`, `17c1c32`
  - [ ] Test retrieve workflow

### Phase 3: Provenance (Week 3-4)
- [ ] Implement EnvironmentCapture
  - [ ] capture_modules() method
  - [ ] capture_env_vars() method
  - [ ] parse_slurm_script() method
  - [ ] save_script_snapshot() method
- [ ] Add provenance tests

### Phase 3-5: Optional Enhancements (Not Started)

These are **enhancement features** that can be added based on future needs:

**Phase 3: Additional HPC Adapters**
- [ ] PBS/Torque adapter (similar to SlurmAdapter)
- [ ] LSF adapter (IBM Spectrum LSF)
- [ ] SGE adapter (Sun Grid Engine)
- [ ] Local adapter (fork processes for testing)

**Phase 4: Enhanced CLI Features**
- [ ] Interactive job monitoring (watch mode with auto-refresh)
- [ ] Batch job submission (submit multiple runs)
- [ ] Job cancellation command (`chemvcs cancel <job-id>`)
- [ ] Job dependency chains
- [ ] Job array support

**Phase 5: Optimization & Polish**
- [ ] Performance optimization for large repositories
- [ ] Enhanced error messages and recovery
- [ ] Additional provenance metadata
- [ ] Integration with job accounting systems
- [ ] Release v0.7.0

---

## Future Enhancements (Post-M6)

### Essential Git Features
- [ ] **Hooks system** - Automation for workflows
  - [ ] pre-commit, post-commit hooks
  - [ ] pre-push, post-receive hooks
  - [ ] Hook configuration and management
- [ ] **Tags** - Version markers for releases
  - [ ] Lightweight tags
  - [ ] Annotated tags with messages
  - [ ] Tag listing and filtering
- [ ] **Submodules** - Dependency management
  - [ ] Add/remove submodules
  - [ ] Submodule update
  - [ ] Recursive operations

### Performance
- [ ] Implement packfile format for efficient storage
- [ ] Add incremental push/pull optimization
- [ ] Implement shallow clone support
- [ ] Add object caching layer

### User Experience
- [ ] Interactive conflict resolution
- [ ] Visual diff for molecular structures
- [ ] Web UI for browsing repositories
- [ ] IDE/editor integrations

### Essential Git Features (High Priority)
- [ ] **Hooks system** - Critical for automation
  - [ ] pre-commit: Validate input files, check file sizes
  - [ ] post-commit: Auto-submit to compute queue, trigger CI
  - [ ] pre-push: Prevent pushing large orbital files
  - [ ] post-receive: Server-side archiving
- [ ] **Submodules** - For shared resources
  - [ ] Share structure libraries across projects
  - [ ] Manage basis set dependencies
  - [ ] Large project dependency management

### Advanced Features
- [ ] Signed commits (GPG)
- [ ] Garbage collection for unreachable objects
- [ ] Repository mirroring
- [ ] Multi-user access control
- [ ] Stash - Temporarily save work state (low priority)
- [ ] Reflog - Recovery from mistakes (low priority)
- [ ] Bisect - Binary search for problematic commits (low priority)

### Chemistry-Specific
- [ ] Automatic structure comparison (RMSD, graph isomorphism)
- [ ] Energy landscape tracking
- [ ] Trajectory file support
- [ ] Integration with molecular viewers

### Testing & Quality
- [ ] Increase test coverage to 80%+
- [ ] Add benchmarks for large repositories
- [ ] Performance profiling and optimization
- [ ] Fuzz testing for protocol handlers

### Deliberately Excluded Features

**Not planned for ChemVCS** (with rationale):

- ❌ **Rebase / Cherry-pick** - Violates scientific provenance principles
  - Science requires authentic history, not beautified commits
  - Rewriting history breaks computational reproducibility
  - Real calculation order matters for traceability

- ❌ **Reset (--soft/--mixed/--hard)** - Conceptually complex and unnecessary
  - ChemVCS has no staging area, so most reset modes don't apply
  - Use `checkout` for switching states instead
  - Simpler mental model for users

- ❌ **Interactive staging (git add -p)** - Not applicable
  - ChemVCS intentionally removed staging area concept
  - Commit entire working directory state (scientific snapshot)
  - Partial staging contradicts "complete calculation state" model

---

## Completed ✅

### Milestone 1: Local VCS Core (M1) ✅
- [x] Core data structures (Object, Snapshot, Reference)
- [x] Content-addressable storage with SHA-256
- [x] Repository operations (init, commit, log, branch, checkout)
- [x] CLI with basic commands
- [x] 38 tests passing

### Milestone 2: Working Directory (M2) ✅
- [x] Workspace scanning for directory trees
- [x] Object creation with blob references
- [x] Diff engine for tree comparison
- [x] Status command
- [x] Working directory restoration
- [x] 51 tests passing

### Milestone 3: Branching and Merge (M3) ✅
- [x] Fast-forward merge
- [x] Ancestor checking
- [x] Merge command
- [x] Branch divergence detection
- [x] 57 tests passing

### Milestone 4: Remote Repositories (M4) ✅
- [x] HTTP client for remote protocol
- [x] Remote ref management
- [x] Object batch upload/download
- [x] Push/Pull/Fetch operations
- [x] HTTP server (chemvcs-server)
- [x] Remote configuration management
- [x] 72 tests passing

### Quality Improvements (Post-M4) ✅
- [x] Checkout cleanup - Remove files not in target snapshot
- [x] Three-way merge algorithm - Support diverged branches
- [x] Conflict detection - Report conflicting paths
- [x] Merge commit creation with two parents
- [x] Intelligent merge strategies (auto-resolve one-sided changes)
- [x] 80 tests passing

---

## Known Issues

### Critical
- None

### Minor
- Server authentication is placeholder (no actual token validation)
- No repository size limits enforced
- Limited error recovery in network operations

### Documentation
- Need more usage examples
- API documentation could be more detailed

---

## Questions for Discussion

1. Should M5 use REST API or CGo for Python-Go integration?
2. What QM codes should we prioritize for integration?
3. Should we support other schedulers besides SLURM (PBS, LSF)?
4. What authentication mechanism for production deployment?
