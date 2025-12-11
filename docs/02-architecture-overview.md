# ChemVCS – Architecture Overview (v0.1)

## 1. Introduction

This document describes the overall software architecture of ChemVCS. It focuses on
the structure and responsibilities of the core components, their interactions, and
the main data flows for key operations such as `init`, `commit`, `log`, `status`,
`checkout`, and (in later milestones) remote synchronisation and domain-specific
extensions.

The goals of this document are to:

- Provide a clear mental model of how ChemVCS is organised;
- Define module boundaries and responsibilities in the Go and Python layers;
- Show how local repositories, working directories, and remote servers interact;
- Establish a foundation for more detailed design and implementation documents
  (object model, CLI spec, remote protocol, domain layer, HPC adapter).

This document describes the **target architecture** for Milestones 1–6. Some parts
will be implemented later, but the high-level structure is intended to remain stable.

---

## 2. High-Level Architecture

### 2.1 Layered View

ChemVCS is organised into four major layers:

1. **VCS Core (Go, local):**
   - Manages repositories under `.chemvcs/`;
   - Implements the content-addressable object store (blobs, objects, snapshots);
   - Manages refs, branches, and history;
   - Provides core CLI commands (`init`, `commit`, `log`, `branch`, `checkout`, `status`).

2. **Remote Services (Go, server-side):**
   - Hosts one or more repositories on a shared server;
   - Provides HTTP APIs for push / pull / fetch, object upload/download, and ref updates;
   - Handles basic authentication and access control;
   - Shares internal libraries with the local core where appropriate.

3. **Domain Layer (Python):**
   - Encodes the computational chemistry domain model: structures, runs, workflows, datasets;
   - Provides parsers and utilities for specific electronic-structure and atomistic codes;
   - Exposes a Python API for querying and manipulating domain objects backed by ChemVCS;
   - Interacts with the Go core via the CLI and/or well-defined file formats.

4. **Integration Layer (Python, HPC and tooling):**
   - Integrates with HPC schedulers (SLURM, PBS, etc.) through pluggable adapters;
   - Provides higher-level commands / scripts for run lifecycle management;
   - Enables researchers and tooling authors to embed ChemVCS into their workflows.

These layers are loosely coupled: the Go core does not depend on Python, and the
remote services do not require the domain or integration layers to function. This
allows the core to be used in non-chemistry contexts and supports incremental
adoption of domain-specific features.

### 2.2 Component Overview Diagram (Conceptual)

At a conceptual level, the main components can be visualised as:

- **Local Machine / HPC Login Node**
  - `chemvcs` CLI (Go)
  - VCS core libraries (Go `internal/...`)
  - `.chemvcs/` repository
  - Working directory (project files)
  - Python domain & integration tools

- **Remote Server**
  - `chemvcs-server` process (Go)
  - Repository storage (shared filesystem / object store)
  - Authentication / access control layer

Data flows between these components primarily via:

- On-disk repository structures (`.chemvcs/`);
- CLI commands and their JSON or text output;
- HTTP APIs for remote synchronisation.

---

## 3. Local VCS Core (Go)

The local core is implemented in Go and is responsible for all fundamental version
control operations.

### 3.1 Package Structure (Go)

The core Go codebase is organised as follows:

```text
go/
├─ cmd/
│  └─ chemvcs/           # CLI entry point (main)
└─ internal/
   ├─ model/             # Blob, Object, Snapshot, Ref types
   ├─ objectstore/       # Storage and retrieval of blobs/objects
   ├─ repo/              # Repository operations (init, open, commit, log, refs)
   ├─ workspace/         # Mapping between working directory and object graph
   └─ remote/            # Client-side logic for remote sync (later milestones)
```

Each package has a specific responsibility:

- `model`:
  - Defines core data types and their serialisation (JSON or similar);
  - Has no external dependencies on filesystem or CLI.

- `objectstore`:
  - Implements the content-addressable store under `.chemvcs/objects/`;
  - Provides APIs to store and fetch blobs and objects by hash;
  - Enforces hashing, layout, and basic integrity checks.

- `repo`:
  - Implements repository-level operations:
    - `Init`, `Open`;
    - Snapshot creation and retrieval;
    - Branch and ref management;
    - History traversal (log);
  - Uses `objectstore` and `model` but does not know about domain concepts.

- `workspace`:
  - Implements mapping from the working directory file tree to an object graph
    (folder + file objects in the MVP);
  - Computes differences between working directory and a snapshot (for `status`);
  - Applies a snapshot to the working directory (for `checkout`).

- `remote` (client):
  - Implements push / pull / fetch mechanics over HTTP;
  - Handles negotiation of missing objects and ref updates.

The `cmd/chemvcs` package wires these components together behind a consistent CLI.

### 3.2 Repository Responsibilities

A local repository, represented by a `.chemvcs/` directory, is responsible for:

- Storing all objects and snapshots needed to reconstruct its history;
- Storing refs (branch pointers) and the HEAD reference;
- Maintaining on-disk invariants such as:
  - Every ref either points to a valid snapshot or is empty;
  - Every snapshot references existing objects;
  - Every object references blobs/objects that exist or can be fetched.

The Go core ensures these invariants by design, using a strict sequence of
operations: write data first, update refs last.

### 3.3 CLI Responsibilities

The `chemvcs` CLI is the main user-facing entry point for the core. It:

- Translates user commands and flags into calls to `repo`, `workspace`, and `remote` packages;
- Handles error reporting, exit codes, and human-readable output;
- Optionally provides structured output (e.g. JSON) for programmatic clients.

The CLI itself does not embed domain logic; operations such as “what is a run?”
or “how to parse a VASP output?” live in the Python domain layer.

---

## 4. Remote Services (Go)

Remote services allow multiple machines to share repositories and synchronise
objects and refs.

### 4.1 Server Process

A remote server is a long-running Go process, e.g. `chemvcs-server`, which exposes
an HTTP API for repository operations:

- Listing repositories and creating new ones (administrative);
- Uploading / downloading objects and snapshots;
- Reading / updating refs for branches and tags;
- Negotiating missing objects during push / pull.

Internally, the server reuses the same `model` and `objectstore` concepts as the
local core, but operates over a storage root managed by the server (e.g. on a
shared filesystem or object storage backend).

### 4.2 Storage Layout on the Server

Each repository on the server has a layout analogous to the local `.chemvcs/`
directory structure, but rooted under a server-managed path, e.g.:

```text
/var/lib/chemvcs/repos/<repo-id>/
  ├─ objects/
  ├─ refs/
  ├─ HEAD
  └─ config
```

The server may optionally support:

- Additional metadata (e.g. owner, ACLs) stored in a server-side database;
- Background tasks such as garbage collection and packfile creation (later).

### 4.3 Client–Server Interaction

The client-side `remote` package in Go encapsulates the communication patterns
with the server. Typical flows include:

- **Push:**
  - Client determines local objects that are not present on the server;
  - Client uploads those objects in batches;
  - Client requests ref update on the server (e.g. move `refs/heads/main`).

- **Pull / Fetch:**
  - Client asks server which objects are needed to update certain refs;
  - Client downloads missing objects;
  - Client updates local refs accordingly.

The protocol is defined in the remote protocol document and aims to be simple
and robust rather than fully general-purpose.

---

## 5. Domain Layer (Python)

The domain layer provides computational chemistry–specific semantics on top of
the core VCS. It is implemented as a Python package, tentatively `chemvcs_py`.

### 5.1 Responsibilities

The domain layer is responsible for:

- Defining domain models:
  - `Structure`: atomic positions, cell, metadata;
  - `Run`: a single calculation, with input, output, parameters, environment;
  - `Workflow`: a DAG of runs and possibly sub-workflows;
  - `Dataset`: groupings of runs and derived quantities.

- Mapping domain models to and from core `Object`s and `Blob`s:
  - Serialising domain metadata into `Object.meta` as structured JSON;
  - Managing references (`Object.refs`) between domain entities and raw files.

- Parsing and summarising inputs/outputs for specific codes:
  - VASP (INCAR, POSCAR, OUTCAR, etc.);
  - Quantum ESPRESSO, CP2K, and others (as plugins);
  - Capturing key parameters and results in a structured, queryable way.

- Providing a Python API for higher-level operations:
  - Querying runs by structure, method, parameters, etc.;
  - Building workflows programmatically;
  - Integrating with notebooks, scripts, or other tools.

### 5.2 Interaction with the Go Core

The domain layer interacts with the Go core primarily through:

- The local filesystem (`.chemvcs/` directory);
- CLI commands (`chemvcs`) invoked as subprocesses, with text or JSON output;
- Potential direct object reading/writing via the documented file formats.

This approach keeps the Go core independent of Python, while still allowing
Python tools to drive and extend ChemVCS behaviour.

### 5.3 Plugin Mechanism

The domain layer is expected to support plugin-like extension points for:

- **Code parsers:**
  - A uniform interface for “given a run directory, extract metadata and results”;
  - Separate implementations for different quantum chemistry / materials codes.

- **Structure formats:**
  - Loaders and writers for CIF, POSCAR, XYZ, etc.;
  - Optional integration with external libraries if available.

These plugins are purely Python-level and do not affect the core Go architecture.

---

## 6. HPC Integration Layer (Python)

The integration layer connects ChemVCS to HPC schedulers and run management.

### 6.1 Scheduler Adapters

Scheduler adapters are Python components that encapsulate interactions with
specific schedulers, for example:

- SLURM (via `sbatch`, `squeue`, `scancel`);
- PBS / Torque;
- Other batch systems as needed.

Each adapter provides a minimal interface, such as:

- `submit(job_script_path) -> job_id`;
- `status(job_id) -> enum (pending, running, finished, failed, unknown)`;
- `cancel(job_id)`.

The integration layer uses these adapters to implement higher-level commands:

- `chemvcs run create` (prepare run directories and inputs);
- `chemvcs run submit` (submit and record run);
- `chemvcs run status` (monitor runs);
- `chemvcs run collect` (gather outputs and update the VCS).

### 6.2 Separation of Concerns

The HPC integration layer does **not** modify core repository semantics. Instead, it:

- Creates and updates domain-level objects in the Python layer;
- Triggers core VCS operations (commits) through the CLI or direct file writes;
- Records scheduler-specific IDs and statuses as domain metadata.

This separation allows ChemVCS to remain usable for non-HPC or non-chemistry
projects if desired, while still providing a rich integration story for HPC users.

---

## 7. Key Data Flows

This section gives a high-level view of how the main operations flow through the
architecture.

### 7.1 `chemvcs init`

1. User runs `chemvcs init` in a directory.
2. CLI resolves the target path and calls into `repo.Init(path)`.
3. `repo.Init` creates the `.chemvcs/` directory and subdirectories:
   - `objects/`;
   - `refs/heads/`;
   - `HEAD` (pointing to `refs/heads/main`);
   - `config` (optional initial configuration).
4. CLI prints a confirmation message.

No domain or remote components are involved.

### 7.2 `chemvcs commit` (Generic MVP)

1. User runs `chemvcs commit -m "message"`.
2. CLI:
   - Opens the repository (`repo.Open`);
   - Uses `workspace` to scan the working directory and produce a generic object graph representing the file tree (folders + file blobs);
   - Calls `repo.CreateSnapshot(root_object_hash, parents, message, author)`.
3. `repo.CreateSnapshot`:
   - Stores the root object and all required sub-objects via `objectstore`;
   - Creates a `Snapshot` object and stores it;
   - Updates the current branch ref (e.g. `refs/heads/main`);
   - Leaves `HEAD` pointing to that ref.
4. CLI prints the new snapshot hash and returns.

In later milestones, the domain layer may assist in constructing richer object
graphs (runs, workflows) before committing.

### 7.3 `chemvcs log`

1. User runs `chemvcs log`.
2. CLI:
   - Resolves `HEAD` to a branch and snapshot;
   - Calls `repo.Log(ref, limit)` to traverse the snapshot chain.
3. `repo.Log`:
   - Follows parent links from the given snapshot backwards;
   - Returns a sequence of `Snapshot` records.
4. CLI formats and prints the log (optionally in short or JSON form).

### 7.4 `chemvcs status` and `chemvcs checkout` (Later Milestone)

- `status`:
  1. CLI opens the repo and resolves the current snapshot.
  2. `workspace` constructs an object graph representing the current working directory.
  3. The snapshot’s root object is loaded and compared to the working directory graph.
  4. Differences (added/modified/removed) are reported.

- `checkout`:
  1. CLI resolves the target snapshot or branch.
  2. `workspace` materialises the snapshot’s root object into the working directory
     (subject to conflict checks and user confirmation).
  3. `HEAD` and refs are updated as appropriate.

### 7.5 Push / Pull (Later Milestone)

- **Push:**
  1. CLI identifies remote and target branch.
  2. `remote` client gathers local object IDs reachable from the local branch.
  3. Client contacts server and negotiates which objects are missing.
  4. Client uploads missing objects, then requests a ref update on the server.

- **Pull:**
  1. CLI identifies remote and branch to sync from.
  2. `remote` client asks server for the objects reachable from the remote ref.
  3. Client downloads missing objects.
  4. Client updates local refs and optionally working directory (if requested).

---

## 8. Cross-Cutting Concerns

### 8.1 Error Handling and Robustness

- The core must treat repository updates in an atomic manner at the logical level:
  - Objects are written first;
  - Refs are updated last;
  - Incomplete operations should not leave invalid refs.
- CLI should surface clear error messages and exit codes.
- The server should validate incoming data and reject malformed or inconsistent updates.

### 8.2 Performance and Scalability

Initial implementations will favour simplicity and correctness. Over time, the
architecture allows for:

- Packfile or chunk-based storage in `objectstore` to reduce filesystem overhead;
- Caching and lazy loading for large repositories;
- More efficient negotiation protocols between client and server.

These optimisations will reuse the same logical architecture; no fundamental
restructuring should be required.

### 8.3 Security and Access Control

- Local repositories rely on OS-level permissions.
- Remote server integrates authentication (e.g. token-based) and access control.
- Server-side validation prevents ref hijacking and unauthorised updates.

### 8.4 Observability

- Logging in both CLI and server supports debugging of repository operations.
- Optional metrics (e.g. number of objects, repository sizes, request rates)
  can be added on the server side without changing core logic.

---

## 9. Summary

The ChemVCS architecture is built around a small, composable VCS core written in
Go, with a domain-specific layer and HPC integration written in Python. Local
repositories hold all essential VCS data under `.chemvcs/`; the CLI orchestrates
repository operations by invoking core packages. Remote servers host repositories
and participate in object and ref synchronisation. The domain and integration
layers add computational chemistry semantics and HPC capabilities without
modifying the core.

This architecture is intended to be:

- **Modular:** clear boundaries between core, remote, domain, and integration layers;
- **Extensible:** new domain plugins and schedulers can be added without changes
  to the core;
- **Incremental:** local-only usage is possible, with remote and domain features
  added as needed;
- **Robust:** content-addressable storage, immutable snapshots, and conservative
  update semantics protect data integrity.

Subsequent documents (object model and storage, CLI specification, remote protocol,
domain layer design, HPC adapter design) will refine and detail specific parts of
this architecture.
