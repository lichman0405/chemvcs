# ChemVCS Development Plan (Draft v0.1)

> Goal: Build a distributed version control system (ChemVCS) tailored for computational chemistry / materials projects.
> The foundation is a general VCS core, on top of which domain-specific chemistry and HPC integration layers can be added.

---

## 1. Project Goals and Scope

### 1.1 Goals

1. Provide a **content-addressable, Merkle DAG–based version control core** that supports:
   - Local repository initialization and management;
   - Storage and retrieval of core entities (Blob / Object / Snapshot);
   - Branch and history (refs and snapshots) management.

2. For computational chemistry / materials projects, progressively provide:
   - Domain-level entities for structures, calculation runs, and workflows (DAGs of runs);
   - Full lifecycle traceability of calculations (input, output, parameters, environment);
   - Integration with HPC schedulers (submit, monitor, collect).

3. Support basic collaboration:
   - Remote repositories (central repos) with push / pull / fetch capabilities;
   - Shared object storage and access control for multi-user / multi-project scenarios.

### 1.2 Out-of-Scope (for the initial phases)

For the v0.x phase (at least the first 6–12 months), the project will **not** aim to:

- Fully replace git for **source code management**;
- Provide advanced three-way automatic merge for arbitrary text files;
- Implement a complete public hosting platform comparable to GitHub;
- Deliver a full-featured web UI platform (early stages may include minimal APIs and simple visualisation only).

---

## 2. Milestones and Phased Delivery

The project will be developed in phases. Each milestone should produce a usable, self-contained increment.

### Milestone 1: Local VCS Core MVP

Scope:

- Repository initialization: `chemvcs init`;
- Local object storage for Blob / Object / Snapshot;
- Branch and HEAD / refs management;
- Creating snapshots: `chemvcs commit` (initially with a simplified root object model);
- Viewing history: `chemvcs log`;
- Basic `.chemvcs/` directory layout.

Out of scope for this milestone:

- Full working directory scanning and domain object construction;
- `diff` / `status`;
- Remote sync;
- HPC integration.

### Milestone 2: Working Directory View and Change Detection

On top of Milestone 1:

- Define mapping rules between the working directory and the object graph;
- Provide:
  - `chemvcs status`: detect added / modified / removed objects;
  - `chemvcs checkout`: restore a working directory view from a snapshot.
- Implement an initial object-level diff capability (list added / removed / modified Objects).

### Milestone 3: Branches and Basic Merge Strategy

On top of Milestone 2:

- Branch management commands: `chemvcs branch` / `chemvcs checkout <branch>`;
- Support fast-forward merges only:
  - `chemvcs merge <branch>` updates the target branch when it is an ancestor of the source;
- No complex three-way merge yet; only simple history advancement is supported.

### Milestone 4: Remote Repositories and Sync Protocol

Once the core is stable, implement remote sync:

- Remote repository server (HTTP API):
  - Object upload / download;
  - Ref read / update;
  - Basic authentication / token support.
- Client-side commands:
  - `chemvcs remote add` / `remote list`;
  - `chemvcs push` / `chemvcs pull` / `chemvcs fetch`.

### Milestone 5: Computational Chemistry Domain Layer (Python)

On top of the generic VCS core, add a domain layer for computational chemistry:

- Python package: `chemvcs_py`;
- Define domain-level entities:
  - Structure;
  - Run (single calculation);
  - Workflow (DAG of runs);
- Project layout conventions:
  - Agreed directory structure, e.g. `structures/`, `runs/`, etc.;
- Code parser plugins for VASP / QE (and others):
  - Input / output parsing;
  - Parameter extraction;
  - Result summaries.

### Milestone 6: HPC Scheduler and Run Integration

On top of Milestone 5, provide HPC integration:

- Scheduler plugins: SLURM / PBS / others (extensible);
- Run lifecycle management:
  - `chemvcs run create`: generate inputs and job scripts;
  - `chemvcs run submit`: submit jobs and record job IDs;
  - `chemvcs run status`: query job status;
  - `chemvcs run collect`: collect outputs and register results.

---

## 3. Technology Stack and Overall Architecture

### 3.1 Technology Choices

- Core language: **Go**
  - VCS core;
  - Command-line interface (CLI);
  - Remote repository server.

- Domain and integration layer: **Python**
  - Computational chemistry domain entities and parsing;
  - HPC scheduler integration;
  - Interfaces for analysis workflows (Jupyter / scripts).

- Storage and serialization:
  - Local objects: file system + SHA-256;
  - Metadata serialization: JSON (for the MVP), with the option to switch to Msgpack/CBOR later if needed.

### 3.2 Repository Layout (Suggested)

```text
chemvcs/
├─ go/
│  ├─ cmd/
│  │  └─ chemvcs/            # CLI entry point
│  └─ internal/
│     ├─ model/              # Blob/Object/Snapshot/Ref types
│     ├─ objectstore/        # Local object storage
│     ├─ repo/               # Repository operations: init/commit/log/branch/...
│     ├─ workspace/          # Working directory view and status/checkout
│     └─ remote/             # Remote sync client logic
├─ python/
│  └─ chemvcs_py/            # Python domain layer and HPC integration
├─ server/
│  └─ main.go                # Remote repository server entry (shares logic with internal packages)
├─ docs/
│  ├─ 01-vision-and-scope.md
│  ├─ 02-architecture-overview.md
│  ├─ 03-object-model-and-storage.md
│  ├─ 04-cli-spec-mvp.md
│  ├─ 05-remote-protocol-and-server.md
│  ├─ 06-python-domain-layer.md
│  ├─ 07-hpc-adapter-design.md
│  └─ 08-testing-and-quality.md
└─ .github/workflows/        # CI configuration
```

---

## 4. Object Model and Storage Design (Overview)

A detailed specification will be provided in `03-object-model-and-storage.md`. This section summarises the key ideas.

### 4.1 Core Types

#### Blob

- Represents arbitrary binary data (e.g. file contents).
- Fields:
  - `id`: string, SHA-256 hash (hex);
  - `size`: byte length.

#### Object

- Represents a typed node with metadata.
- Fields:
  - `id`: SHA-256 hash (hex);
  - `type`: string, e.g. `generic` / `structure` / `run` / `workflow` / `folder`;
  - `meta`: JSON object, key-value pairs defined by higher layers;
  - `refs`: array of strings, each a hash pointing to other Objects or Blobs.

#### Snapshot

- Represents a snapshot of the project state (similar to a git commit).
- Fields:
  - `id`: SHA-256 hash;
  - `root`: hash of the root Object for this snapshot;
  - `parents`: array of parent snapshot hashes;
  - `author`: string;
  - `timestamp`: UTC timestamp;
  - `message`: commit message.

#### Ref

- A named pointer to a Snapshot.
- Layout:
  - `refs/heads/<branch>` for branches;
  - `HEAD` to point to a ref (e.g. `ref: refs/heads/main`).

### 4.2 Local Repository Layout

In the project root directory:

```text
.chemvcs/
├─ objects/
│  ├─ ab/
│  │  └─ cdef1234...        # object or blob file, named by full hash
│  └─ ...
├─ refs/
│  └─ heads/
│     └─ main               # stores snapshot hash
├─ HEAD                     # current HEAD, e.g. "ref: refs/heads/main"
└─ config                   # repository configuration (author info, etc.)
```

Object file contents (MVP):

- Blob: raw binary data;
- Object / Snapshot: JSON (optionally compressed).

---

## 5. CLI Specification (MVP Overview)

Detailed command specifications will be documented in `04-cli-spec-mvp.md`. This section outlines the main commands for Milestones 1–3.

### 5.1 `chemvcs init`

- Purpose: initialise a new repository in the current or specified directory.
- Behaviour:
  - Create `.chemvcs/objects`, `.chemvcs/refs/heads`, etc.;
  - Create `HEAD` pointing to `refs/heads/main` by default;
  - Do not create any Snapshots yet.

### 5.2 `chemvcs commit`

- Purpose: create a new snapshot based on the current working directory state.
- MVP implementation:
  - Initially, use a minimal root Object representation (e.g. a simple file tree model);
  - Domain-aware object graph construction will be added later.
- Options:
  - `-m, --message`: commit message (required).
- Behaviour:
  - Use the current HEAD snapshot as the parent (if any);
  - Create the new Snapshot;
  - Update the current branch ref.

### 5.3 `chemvcs log`

- Purpose: display the snapshot history for the current branch.
- Behaviour:
  - Start from the Snapshot pointed to by HEAD and follow parent links backwards;
  - Print snapshot ID, author, timestamp, and message.
- Options:
  - `--oneline`, `--max-count` etc. can be added in later iterations.

### 5.4 `chemvcs branch`

- Purpose: list, create, or delete branches.
- Behaviour:
  - Without arguments: list all branches under `refs/heads/`;
  - `chemvcs branch <name>`: create a new branch at the current HEAD.

### 5.5 `chemvcs checkout`

- Purpose: switch the current branch or snapshot view.
- MVP behaviour:
  - Update what `HEAD` points to (branch or detached snapshot);
  - Working directory restoration is implemented in Milestone 2.

### 5.6 `chemvcs status` (Milestone 2)

- Purpose: show differences between working directory and the current snapshot.
- Behaviour:
  - List added / modified / removed objects;
  - When the domain layer is integrated, provide domain-aware diffs.

---

## 6. Detailed Tasks for Milestones 1–3

### 6.1 Milestone 1: Local VCS Core MVP

#### 6.1.1 Engineering Setup

- Initialise the Go module: `go mod init`;
- Create the base directory structure: `go/cmd/chemvcs`, `go/internal/...`;
- Set up basic CI:
  - `go test ./...`;
  - Code formatting with `gofmt`.

#### 6.1.2 Implement Object Model

- In `internal/model`:
  - Define Go structs for `Blob`, `Object`, `Snapshot`, `Ref`;
  - Implement JSON marshalling / unmarshalling;
  - Add unit tests for normal and error cases.

#### 6.1.3 Implement Object Store

- In `internal/objectstore`:
  - Use SHA-256 as the object ID, layout `objects/aa/bb...`;
  - APIs:
    - `PutBlob` / `GetBlob`;
    - `PutObject` / `GetObject`;
  - Ensure idempotent writes for existing objects;
  - Validate hash inputs and handle invalid hashes gracefully;
  - Use temporary directories for tests.

#### 6.1.4 Repository Structure and Refs

- In `internal/repo`:
  - `Init(path)`: create `.chemvcs` directory layout;
  - `Open(path)`: open an existing repo;
  - Refs API:
    - Read/write `refs/heads/<branch>`;
    - Parse and resolve `HEAD`.

#### 6.1.5 Snapshots and History

- In `repo`:
  - Implement `CreateSnapshot` and `GetSnapshot`;
  - Implement `Log(ref, limit)` to traverse the snapshot chain;
- Testing:
  - Construct a chain of snapshots, verify parent relationships and log order.

#### 6.1.6 CLI MVP

- In `cmd/chemvcs`:
  - `chemvcs init`;
  - `chemvcs commit -m` (initially using a simple root object representation);
  - `chemvcs log`.
- End-to-end tests:
  - Use temporary directories to run CLI commands;
  - Inspect `.chemvcs` contents and `log` output.

### 6.2 Milestone 2: Working Directory and Status

#### 6.2.1 Working Directory / Object Mapping

- In `internal/workspace` define:
  - How to construct an object tree from the working directory;
  - How to materialise a snapshot (object graph) into the working directory.
- For the MVP, use a generic file tree → folder+blob object mapping (domain-agnostic).

#### 6.2.2 Implement `status`

- Compare the current snapshot with working directory state to compute:
  - Added / modified / removed objects;
- Connect results to CLI output.

#### 6.2.3 Implement `checkout`

- Restore the working directory from a given snapshot;
- Handle cases where the working directory has uncommitted changes (warnings / protection).

### 6.3 Milestone 3: Branches and Basic Merge

#### 6.3.1 Branch Management

- Implement `chemvcs branch`:
  - List branches;
  - Create branches at current HEAD.
- Implement `chemvcs checkout <branch>`:
  - Switch the current branch ref and working directory view.

#### 6.3.2 Fast-Forward Merge

- Implement `chemvcs merge <branch>`:
  - Check branch relationship (fast-forward only);
  - Update the target branch ref to the new snapshot;
  - When not a fast-forward case, return an explicit error (no attempt to auto-merge).

---

## 7. Testing and Quality Assurance

### 7.1 Testing Strategy

- Unit tests:
  - Cover `model`, `objectstore`, `repo` and other core packages;
- Integration tests:
  - End-to-end tests for CLI commands executed in temporary directories;
- Invariant checks:
  - Verify object hash matches content after each write;
  - Validate data integrity when reading objects.

### 7.2 Static Analysis and Code Style

- Use Go tooling:
  - `go vet` for static analysis;
  - `gofmt` for formatting;
- Optionally introduce `golangci-lint` for stricter linting.

### 7.3 Data Safety and Failure Modes

- Separate object writes from ref updates:
  - First write objects, then update refs;
  - Ensure failures do not leave inconsistent state;
- Provide simple backup guidance:
  - Copying the `.chemvcs` directory is sufficient to back up a repository.

---

## 8. Future Work and Extensions (High Level)

After Milestones 1–3 are complete, the project can progress to:

1. Remote repository server (in `server/main.go`):
   - HTTP API and access control;
2. Python domain layer:
   - Computational chemistry object model;
   - Parsers and query interfaces;
3. HPC integration:
   - Scheduler plugins;
   - Run lifecycle management commands.

These topics will be detailed in subsequent design documents.
