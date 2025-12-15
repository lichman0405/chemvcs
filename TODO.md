# ChemVCS TODO List

## Milestone 5: Python Domain Layer (M5) - 85% Complete ✅

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
  - ✅ 28 pytest tests (Structure, XYZ, POSCAR)
  - ✅ basic_usage.py example
  - ✅ run_workflow_example.py with DAG demonstration
- ✅ Document Python API
  - ✅ Python package README
  - ✅ Design doc (06-python-domain-layer.md)

### Remaining Tasks (Lower Priority)
- ⏭️ CIF file parser (crystallography format)
- ⏭️ Additional domain object tests (Run, Workflow specific)
- ⏭️ Enhanced diff/merge for molecular structures
- ⏭️ Advanced Repository API features

### Documentation
- ✅ Python API usage examples
- ⏭️ Comprehensive API reference documentation
- ⏭️ Integration guide with QM codes (ORCA, VASP, Gaussian)

---

## Milestone 6: HPC Integration (M6)

### Core Tasks
- [ ] Design SLURM adapter interface
- [ ] Implement job submission tracking
  - [ ] Capture job ID in commit metadata
  - [ ] Link outputs back to commits
- [ ] Add HPC-specific commands
  - [ ] `chemvcs submit` - Submit job and commit
  - [ ] `chemvcs jobs` - List tracked jobs
  - [ ] `chemvcs retrieve` - Fetch completed job outputs
- [ ] Implement job status monitoring
- [ ] Add provenance metadata for computational jobs
  - [ ] Compute environment (modules, versions)
  - [ ] Job script
  - [ ] Resource usage
- [ ] Write integration tests with SLURM
- [ ] Document HPC workflow patterns

### Documentation
- [ ] HPC setup guide
- [ ] SLURM integration examples
- [ ] Best practices for computational provenance

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
