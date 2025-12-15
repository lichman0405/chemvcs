# ChemVCS TODO List

## Milestone 5: Python Domain Layer (M5) - 100% Core Complete ✅

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

### Documentation
- ✅ Python API usage examples
- ⏭️ Comprehensive API reference documentation
- ⏭️ Integration guide with QM codes (ORCA, VASP, Gaussian)

---

## Milestone 6: HPC Integration (M6) - In Progress

**Design Document**: docs/09-hpc-integration-design.md

### Phase 1: Core Infrastructure (Week 1-2)
- [x] Create M6 design document
- [ ] Extend Run object with HPC fields
  - [ ] Add job_id, job_system, queue_name fields
  - [ ] Add provenance fields (modules, resources, script)
  - [ ] Add mark_queued(), mark_retrieved() methods
  - [ ] Update to_core_object() / from_core_object()
- [ ] Add Run HPC tests
  - [ ] Test HPC field serialization
  - [ ] Test mark_queued() method
  - [ ] Test mark_retrieved() method
  - [ ] Test round-trip with HPC data
- [ ] Design JobAdapter interface
  - [ ] Define abstract base class
  - [ ] Define JobStatus enum
  - [ ] Define JobInfo dataclass
  - [ ] Define exception hierarchy
- [ ] Implement SlurmAdapter
  - [ ] Implement submit() with subprocess mock
  - [ ] Implement get_status() with squeue
  - [ ] Implement get_info() with sacct
  - [ ] Implement cancel() with scancel
  - [ ] Add status parsing logic
- [ ] Add adapter unit tests (mock-based)
  - [ ] Test submit() success/failure
  - [ ] Test get_status() for all states
  - [ ] Test get_info() parsing
  - [ ] Test cancel() operation
  - [ ] Test error handling

### Phase 2: CLI Integration (Week 2-3)
- [ ] Go: Add `chemvcs submit` command
  - [ ] Parse arguments (run-id, script path)
  - [ ] Validate Run object exists
  - [ ] Call Python HPC module
  - [ ] Handle errors gracefully
- [ ] Go: Add `chemvcs jobs` command
  - [ ] List Run objects with job_id
  - [ ] Query status for each job
  - [ ] Format output table
  - [ ] Add filtering options
- [ ] Go: Add `chemvcs retrieve` command
  - [ ] Parse job-id argument
  - [ ] Find corresponding Run object
  - [ ] Copy output files
  - [ ] Create commit
- [ ] Python: Implement JobSubmitter class
  - [ ] submit_run() method
  - [ ] Capture environment
  - [ ] Update Run object
  - [ ] Store to repository
- [ ] Python: Implement JobTracker class
  - [ ] list_jobs() method
  - [ ] get_job_status() method
  - [ ] find_run_by_job_id() method
- [ ] Python: Implement JobRetriever class
  - [ ] retrieve_results() method
  - [ ] Copy output files
  - [ ] Update Run object
  - [ ] Create commit
- [ ] Add CLI integration tests
  - [ ] Test submit workflow
  - [ ] Test jobs listing
  - [ ] Test retrieve workflow

### Phase 3: Provenance (Week 3-4)
- [ ] Implement EnvironmentCapture
  - [ ] capture_modules() method
  - [ ] capture_env_vars() method
  - [ ] parse_slurm_script() method
  - [ ] save_script_snapshot() method
- [ ] Add provenance tests
  - [ ] Test module capture (mocked)
  - [ ] Test script parsing
  - [ ] Test resource extraction
- [ ] Integrate into submit workflow
  - [ ] Capture before submission
  - [ ] Store in Run object
  - [ ] Verify in tests

### Phase 4: Testing & Documentation (Week 4-5)
- [ ] End-to-end workflow tests
  - [ ] Test complete submit → query → retrieve
  - [ ] Test error scenarios
  - [ ] Test concurrent jobs
- [ ] Create user guide (docs/10-hpc-user-guide.md)
  - [ ] Quick start guide
  - [ ] Command reference
  - [ ] Example workflows
  - [ ] Troubleshooting
- [ ] Create example repository
  - [ ] examples/hpc-workflow/
  - [ ] Sample structure files
  - [ ] Sample job scripts
  - [ ] Tutorial walkthrough
- [ ] Update project documentation
  - [ ] README.md with M6 features
  - [ ] FEATURE_STATUS.md mark M6 complete
  - [ ] CHANGELOG.md add v0.6.0

### Phase 5: Polish (Week 5-6)
- [ ] Code review and refactoring
- [ ] Performance optimization
- [ ] Additional error handling
- [ ] Documentation improvements
- [ ] Release v0.6.0

### Optional Extensions
- [ ] PBS/Torque adapter
- [ ] SGE adapter
- [ ] Local adapter (fork processes)
- [ ] Job array support
- [ ] Job dependency chains

---

## Future Enhancements (Post-MVP)

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
