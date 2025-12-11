# ChemVCS – Testing and Quality Strategy (v0.1)

## 1. Introduction

This document describes the testing and quality assurance (QA) strategy for
ChemVCS, covering:

- Testing goals and principles;
- Test taxonomy (unit, integration, end-to-end, performance);
- Strategies for the Go core, remote server, Python domain layer, and HPC layer;
- Continuous integration (CI) and quality gates;
- Data integrity, compatibility, and regression safeguards.

The intent is to ensure that ChemVCS remains reliable and maintainable as it
evolves, and that changes can be made confidently without breaking existing
repositories or workflows.

---

## 2. Testing Goals and Principles

### 2.1 Goals

The testing strategy aims to:

1. **Protect data integrity**
   - Prevent corruption of `.chemvcs` repositories;
   - Detect hash mismatches, invalid refs, and inconsistent object graphs early.

2. **Ensure functional correctness**
   - Confirm that core commands (`init`, `commit`, `log`, `status`, `branch`,
     `checkout`) and remote operations behave as specified;
   - Validate that domain and HPC layers correctly reflect the state of runs,
     structures, and workflows.

3. **Support refactoring and evolution**
   - Allow internal changes (e.g. refactors, optimisations) to be made with
     confidence that observable behaviour remains correct.

4. **Enable reproducibility**
   - Ensure that given the same inputs and repository state, operations produce
     the same results across platforms and versions.

### 2.2 Principles

Key principles guiding test design:

- **Fast feedback** – unit tests and core integration tests should run quickly
  in CI to provide feedback on each change.
- **Determinism** – tests should avoid flaky behaviour (e.g. dependence on wall
  clock time or external services without proper control).
- **Isolation** – tests should not depend on global state; they should create
  and clean up their own repositories and temporary directories.
- **Representative scenarios** – end-to-end tests should mirror realistic usage
  patterns (e.g. “researcher runs a DFT calculation and collects results”).

---

## 3. Test Taxonomy

Tests are categorised as follows:

1. **Unit tests**
   - Validate individual functions or small modules in isolation.
   - Use mocks or in-memory structures where appropriate.

2. **Integration tests**
   - Exercise interactions between modules within a language layer.
   - For example, Go CLI + repo + objectstore, or Python domain + IO parsers.

3. **End-to-end (E2E) tests**
   - Treat ChemVCS as a black box from a user perspective.
   - Use temporary directories, real CLI commands, and (optionally) remote server
     instances to simulate complete workflows.

4. **Performance and stress tests (later)**
   - Measure performance for operations such as committing large trees, pushing
     large repositories, or scanning many runs.
   - Identify and prevent regressions in critical paths.

5. **Compatibility and migration tests**
   - Validate reading/writing of repositories created by older versions;
   - Test migration routines if/when on-disk formats evolve.

---

## 4. Go Core Testing Strategy

The Go core comprises the `model`, `objectstore`, `repo`, `workspace`, `remote`
(client), and CLI packages.

### 4.1 Unit Tests

Each internal package should have unit tests in `*_test.go` files.

#### 4.1.1 `model`

- Verify JSON encoding/decoding of `Blob` metadata (if any), `Object`, and
  `Snapshot` structs.
- Confirm canonicalisation where required (e.g. stable ordering of fields).
- Validate round-trip: encode → decode yields an equivalent in-memory
  representation.

#### 4.1.2 `objectstore`

- Test storage and retrieval of blobs and objects using temporary directories.
- Verify sharding layout (`objects/aa/bb...`).
- Confirm that the stored hash matches content for both writes and reads.
- Test behaviour on invalid input (e.g. malformed hashes, missing files).

#### 4.1.3 `repo`

- Test repository initialisation and opening.
- Validate creation of snapshots:
  - Correct handling of parent lists;
  - Proper ref updates.
- Verify `Log` traversal:
  - Simple linear chains;
  - Multi-parent (merge) scenarios when implemented.

#### 4.1.4 `workspace`

- For the generic file tree mapping:
  - Test conversion between file system trees and object graphs;
  - Verify that round-trip (object graph → filesystem → object graph) yields
    equivalent graphs.
- Test diff logic used by `status` for simple and nested directory scenarios.

### 4.2 Integration Tests (Go)

Integration tests exercise CLI commands and multiple packages together.

Typical patterns:

- Use `os.MkdirTemp` to create isolated project directories;
- Run CLI entrypoints via test harness or `exec.Command` in a separate process;
- Inspect `.chemvcs` contents and CLI output.

Test scenarios:

- `init` followed by `commit` and `log`:
  - Create files, run `chemvcs commit -m`, check that a snapshot is created and
    logged as expected.
- `status` and `checkout`:
  - Modify files between commits and verify that `status` reports correct changes;
  - Use `checkout` to revert to prior snapshots and confirm file contents.
- Error paths:
  - Running `commit` without a repository;
  - Using invalid options or arguments.

### 4.3 Remote Client and Server Tests

Remote-related tests include:

- Unit tests for client-side `remote` package:
  - Mock HTTP server to simulate responses;
  - Verify object negotiation logic and ref update handling.
- Unit tests for server handlers:
  - Use an in-memory or temporary filesystem-backed repository root;
  - Test endpoints for object existence, upload/download, and ref updates.

End-to-end tests (Go + server):

- Start a `chemvcs-server` instance bound to a random local port;
- Use CLI commands (or direct client library calls) to:
  - Push a local repository to the server;
  - Pull from the server into a fresh repository;
  - Verify that resulting repositories have identical content and history.

---

## 5. Python Domain Layer Testing Strategy

The Python domain layer (`chemvcs_py`) includes domain models, mapping logic,
IO parsers, and query utilities.

### 5.1 Unit Tests (Python)

Use `pytest` as the primary test framework.

#### 5.1.1 Domain Models

- Test dataclass constructors and simple invariants for `Structure`, `Run`,
  `Workflow`, `Dataset`.
- Validate that minimal required fields are present and consistent (e.g.
  positions and species length match, etc.).

#### 5.1.2 Mapping Functions

- Test `structure_to_object` / `object_to_structure`, `run_to_object` /
  `object_to_run`, etc., for round-trip equivalence.
- Use controlled metadata and references to ensure deterministic hashes and
  serialisation.

#### 5.1.3 IO Parsers

- For each supported code (e.g. VASP), include small test fixtures:
  - Example POSCAR/INCAR/OUTCAR files with known contents;
  - Test parsing functions (`parse_inputs`, `parse_outputs`) to ensure expected
    metadata and results.
- Validate error handling for malformed or incomplete files.

### 5.2 Integration Tests (Python)

Integration tests in Python exercise the Python layer against real `.chemvcs`
repositories.

Scenarios:

- Create a temporary repository using the Go CLI or Python helper;
- Use Python domain APIs to:
  - Import a structure from a POSCAR/ CIF file;
  - Create a `Run` object linked to that structure;
  - Convert and store these domain objects as core `Object`s via a `Repo` helper;
  - Verify that objects can be reloaded from disk and reconstructed as domain
    objects with equivalent data.

- Query tests:
  - Populate a repository with multiple runs (different codes, parameters);
  - Use query utilities to filter by code, parameter ranges, status;
  - Validate results against expectations.

### 5.3 Type Checking and Static Analysis

- Use `mypy` to type-check the Python codebase where type annotations are present.
- Enforce a baseline of type coverage for core modules (domain models, mapping,
  repo utilities).
- Optionally use `ruff` or `flake8` for linting and style checks.

---

## 6. HPC Layer Testing Strategy

The HPC layer interacts with external schedulers and the filesystem, which can
make pure unit testing challenging. The strategy is to combine:

- **Unit tests with mocks** for scheduler adapters;
- **Integration tests** that shell out to dummy scheduler binaries or test
  harnesses.

### 6.1 Unit Tests

- Implement mock `SchedulerAdapter` implementations that simulate submission,
  status transitions, and cancellation without calling real schedulers.
- Test run lifecycle functions (`create_run`, `submit_run`, `update_run_status`,
  `collect_run_outputs`) using these mocks and temporary run directories.
- Confirm that run metadata and statuses change as expected.

### 6.2 Integration Tests with Dummy Schedulers

To avoid relying on a real SLURM/PBS environment in CI, introduce dummy scheduler
scripts that mimic scheduler behaviour at a basic level:

- `dummy_sbatch`, `dummy_squeue`, etc., which:
  - Accept job scripts;
  - Write “job IDs” to a state file;
  - Simulate status transitions based on simple rules or triggers.

Environment variables and config can direct the HPC layer to use these dummy
commands during integration tests.

### 6.3 Real-Cluster Testing (Outside CI)

For real HPC environments, separate manual or scheduled tests can be run:

- Deploy ChemVCS on a test cluster;
- Submit real jobs (small test computations);
- Verify end-to-end behaviour in a realistic environment.

These tests are not part of regular CI but provide valuable validation before
releases.

---

## 7. Data Integrity and Compatibility Tests

### 7.1 Hash and Ref Integrity

Add dedicated tests to ensure that:

- All stored objects pass hash verification when scanned;
- Refs always point to existing snapshots after typical operations;
- Snapshots reference existing root objects and that object graphs are complete.

A “repository validator” tool or function can be used in tests:

- Traverse `objects/` and `refs/`;
- Validate hashes and JSON schemas;
- Report any inconsistencies as test failures.

### 7.2 Cross-Version Compatibility

As ChemVCS evolves, compatibility tests should ensure that new versions can:

- Open and operate on repositories created by earlier versions;
- Read and interpret previous `version` fields in object and snapshot JSON.

Test approach:

- Maintain a set of “golden” repositories produced by specific released versions:
  - Small but representative collections of snapshots, refs, and objects.
- For each new version:
  - Load and validate these repositories;
  - Run a subset of typical operations (log, status, commit) and verify outcomes.

If migration routines are introduced (e.g. format upgrades), tests should cover:

- Migration of older repositories to newer layouts;
- Safety of migration in the presence of partial or corrupted data.

---

## 8. Continuous Integration (CI) and Quality Gates

### 8.1 CI Pipeline Overview

A typical CI pipeline (e.g. GitHub Actions) should include:

1. **Go checks**:
   - `go test ./...` for unit and integration tests;
   - `go vet ./...` for static analysis;
   - `gofmt` or equivalent check for code formatting.

2. **Python checks**:
   - `pytest` for unit and integration tests;
   - `mypy` for type checking (where enabled);
   - `ruff` or `flake8` for linting.

3. **End-to-end smoke tests**:
   - A small set of E2E workflows that run the CLI and Python toolchain together
     on temporary repositories.

4. **Optional coverage reporting**:
   - Code coverage for Go and Python to monitor test completeness over time.

### 8.2 Supported Platforms

At minimum, CI should test on:

- Linux (primary deployment target for HPC and servers).

Where feasible, additional platforms can be added:

- macOS (developer machines);
- Windows (if the architecture aims to support it).

Tests should be written such that they do not depend on platform-specific paths
or shell features, or handle them explicitly.

### 8.3 Quality Gates

To maintain code quality, merges to main branches should require:

- Successful execution of all tests (Go and Python);
- No `go vet` errors;
- No critical linter errors;
- For significant changes, test coverage should not drop below an agreed threshold
  (e.g. 70–80% for core modules).

Review guidelines can also include:

- Ensuring that new features include appropriate tests;
- Ensuring that changes touching on-disk formats or network protocols are
  accompanied by updated documentation and migration notes.

---

## 9. Developer Guidelines for Testing

### 9.1 Writing New Tests

Developers should:

- Add unit tests for new functions and components, especially for public APIs;
- Extend integration tests where new interactions are introduced;
- Update or add E2E tests when user-facing behaviour changes.

### 9.2 Test Data Management

- Keep test fixtures small and focused (e.g. minimal POSCAR/OUTCAR files);
- Use synthetic or anonymised data where possible;
- Avoid embedding large binary files directly in the repository unless necessary
  (consider generating them during tests if feasible).

### 9.3 Debugging Failing Tests

When tests fail:

- Use log output and temporary directories to inspect intermediate states;
- For Go tests, leverage `t.Log`/`t.Helper` and verbose modes;
- For Python tests, use `pytest -vv` and fixtures for better diagnostics.

Where tests reveal actual design gaps or unclear behaviour, specifications
should be updated accordingly.

---

## 10. Summary

ChemVCS’s testing and quality strategy combines:

- Comprehensive unit tests for core components in Go and Python;
- Integration and end-to-end tests that simulate realistic workflows;
- Explicit data integrity checks for hashes, refs, and object graphs;
- CI pipelines that enforce consistent quality gates across the codebase.

This foundation is intended to keep the project robust as it grows in scope,
enabling contributors to evolve the architecture, extend domain and HPC features,
and optimise performance while maintaining confidence in correctness and
compatibility.
