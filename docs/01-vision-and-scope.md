# ChemVCS – Vision and Scope (v0.1)

## 1. Introduction

ChemVCS is a domain-aware, distributed version control system designed for computational chemistry and computational materials projects. It aims to provide a robust, content-addressable data backbone for managing simulations, workflows, and derived data across local machines, HPC clusters, and shared servers.

Traditional source-code–oriented VCS tools (e.g. git) work well for text and small files, but they do not model the semantics and lifecycle of scientific calculations. ChemVCS focuses on:

- Tracking **scientific objects** rather than just files;
- Capturing the **provenance** of computational runs and workflows;
- Integrating with **HPC schedulers** and scientific codes;
- Enabling reproducible and shareable computational projects.

This document describes the vision, high-level goals, non-goals, stakeholders, primary use cases, scope for the initial versions, and key constraints for ChemVCS. It serves as the foundation for more detailed design documents and implementation work.

---

## 2. Vision

### 2.1 Vision Statement

ChemVCS should become the **standard backbone for managing computational chemistry and materials data**, in the same way that git became the standard for source code. It should make it natural and convenient to:

- Capture every meaningful calculation as an immutable, versioned object;
- Reconstruct and reproduce previous computational workflows with minimal effort;
- Share and publish selected results, with machine-readable provenance;
- Build higher-level tools (e.g. workflow managers, web UIs, data explorers) on top of a clean, well-defined core.

In other words, ChemVCS is envisioned as **a domain-aware data and history layer** that sits beneath scripts, workflow engines, notebooks, and GUIs, providing a single source of truth for what was computed, how, when, and where.

### 2.2 Core Principles

1. **Content-addressable and immutable core**
   - All fundamental entities (blobs, objects, snapshots) are addressed by cryptographic hashes.
   - Objects are immutable; new states are expressed as new objects and new snapshots.
   - This guarantees reproducibility and makes caching and deduplication natural.

2. **Domain-neutral core, domain-specific layers**
   - The core VCS engine remains domain-agnostic (no built-in assumptions about chemistry).
   - Computational chemistry logic (structures, runs, workflows) is implemented in a domain layer on top of the core.
   - This separation preserves generality and makes future reuse in other domains possible.

3. **Provenance first**
   - Every calculation is represented as a node in a graph with explicit dependencies and outputs.
   - Provenance (code version, parameters, environment, HPC resources) is captured as metadata, not an afterthought.
   - Queries over “what was computed” and “how it was computed” become first-class operations.

4. **HPC- and data-friendly**
   - Designed with many small files and occasional large files in mind.
   - Integrates with schedulers (SLURM, PBS, etc.) in a configurable, pluggable way.
   - Minimises assumptions about local filesystem layouts on HPC clusters.

5. **Incremental adoption**
   - Users can start by using ChemVCS as a local data/history manager for a single project.
   - Remote repositories, collaboration, and advanced integration can be added later without breaking early adopters.
   - The system should fail gracefully and preserve data integrity even when partially configured.

---

## 3. High-Level Goals

### 3.1 Functional Goals

1. **Local VCS Core**
   - Initialise and manage local repositories (`chemvcs init`).
   - Store and retrieve blobs, typed objects, and snapshots under a `.chemvcs` directory.
   - Manage branches and snapshot histories (HEAD, refs).
   - Provide basic CLI operations: `init`, `commit`, `log`, `branch`, `checkout`, `status` (later milestone).

2. **Working Directory / Object Graph Mapping**
   - Define and implement a generic mapping from file trees to object graphs and back.
   - Detect changes in the working directory relative to a snapshot (`chemvcs status`).
   - Restore the working directory from a snapshot (`chemvcs checkout`).

3. **Remote Synchronisation**
   - Provide a remote server component to host repositories.
   - Implement push / pull / fetch semantics with efficient object transfer.
   - Support basic authentication and access control.

4. **Domain Layer for Computational Chemistry (Python)**
   - Define domain entities: structures, runs, workflows, datasets.
   - Provide tools to parse and summarise inputs / outputs for common codes (e.g. VASP, QE).
   - Attach domain metadata to core objects in a structured, queryable manner.

5. **HPC Integration**
   - Abstract interfaces for schedulers (submit, query status, cancel).
   - Commands to create, submit, track, and collect runs in a way that updates the VCS state.
   - Support multiple HPC environments via plugins.

### 3.2 Non-Functional Goals

1. **Data integrity and robustness**
   - Hash-based addressing for all stored data.
   - Clear separation between object writes and ref updates to prevent inconsistent states.
   - Conservative behaviour in the presence of errors; never silently discard or corrupt data.

2. **Deployability in HPC environments**
   - Minimal runtime dependencies for the Go binaries.
   - No hard dependency on system-wide Python or external databases for the core VCS.
   - Work well under typical HPC constraints (limited home quota, different modules, etc.).

3. **Extensibility**
   - Clear internal package boundaries in the Go codebase.
   - Plugin-like mechanisms for code parsers and schedulers in the Python layer.
   - Well-documented CLI and file formats, enabling external tools to integrate with ChemVCS.

4. **Predictable and inspectable behaviour**
   - CLI commands should be explicit and unsurprising.
   - Repositories should remain readable and understandable on disk (`.chemvcs` contents are not opaque black boxes).
   - Diagnostics and logs should help users understand what the system is doing.

---

## 4. Non-Goals

In order to keep the project focused and feasible, ChemVCS will explicitly **not** target the following areas in the initial phases:

1. **Full-featured source code VCS replacement**
   - Although ChemVCS can store source files, it is not intended to replace git for day-to-day code management.
   - Advanced features like fine-grained text merge, rebase, cherry-pick, etc., are out of scope for the core roadmap.

2. **Generic big data lake or database system**
   - ChemVCS is not a general-purpose data lake or distributed database.
   - Its concern is **versioned scientific artefacts and provenance**, not arbitrary analytics workloads.

3. **Rich graphical web platform (initially)**
   - A full-blown web UI with dashboards, rich visualisations, and multi-tenant management is not part of the MVP.
   - Early web components, if any, will focus on debugging and basic browsing of repositories.

4. **Fully automated workflow manager**
   - ChemVCS is a foundation for tracking workflows, not an opinionated workflow engine.
   - Integration with existing workflow tools (e.g. custom Python scripts, schedulers, or external engines) is expected.

5. **Highly optimised performance for extreme-scale use cases (initially)**
   - The design should not preclude scaling, but the first versions will prioritise correctness, clarity, and robustness over micro-optimisations.

---

## 5. Stakeholders and Use Cases

### 5.1 Stakeholders

1. **Individual researchers**
   - Run calculations on local workstations and 1–2 HPC clusters.
   - Need reproducible records of their computational experiments.
   - Want to avoid losing track of which input produced which result.

2. **Research groups / labs**
   - Multiple students and postdocs share structures, scripts, and HPC resources.
   - Need a common, standardised way to organise calculations and provenance.
   - Require lightweight collaboration and internal result sharing.

3. **Infrastructure / data managers**
   - Maintain HPC clusters and internal data services.
   - Care about storage efficiency, backup, and enforceable policies.
   - Need a system that is inspectable and maintainable over time.

4. **Tooling and workflow developers**
   - Build higher-level tools (web dashboards, workflow managers, data explorers).
   - Need stable, documented APIs and data models to build on top of ChemVCS.

### 5.2 Representative Use Cases

1. **Track a project’s calculations over time**
   - A researcher initialises a ChemVCS repository for a project.
   - Each new calculation run is registered and committed as a new snapshot.
   - Months later, they can see exactly which inputs, parameters, and code version produced a published result.

2. **Compare different parameter choices**
   - For a given structure, the user runs several calculations with varying functionals or cutoffs.
   - ChemVCS records each run as an object with metadata describing parameters.
   - The user or a script can query and compare the runs via the domain layer.

3. **Move computations between clusters**
   - A project is partly computed on cluster A and then continued on cluster B.
   - ChemVCS synchronises the repository via a remote server.
   - Required objects (inputs, minimal outputs) are fetched on-demand to the new environment.

4. **Publish a reproducible dataset**
   - A subset of calculations are tagged for publication.
   - ChemVCS can export a self-contained snapshot of necessary inputs, metadata, and selected outputs.
   - This snapshot can be archived in an external repository or attached to a paper.

5. **Automated run tracking for a workflow**
   - A Python script, using the domain layer, creates runs, submits jobs to the scheduler, and collects outputs.
   - Each run is recorded in ChemVCS as a node in a provenance graph.
   - The workflow can be rerun or extended while still using the same underlying history.

---

## 6. Scope of Initial Releases

### 6.1 Scope for Milestone 1–3 (Core VCS)

Within Milestones 1–3, the following are in scope:

- Local VCS core (Go):
  - Repository layout (`.chemvcs/`);
  - Object model and storage;
  - Refs and HEAD management;
  - Snapshot creation and history traversal.
- CLI for core operations:
  - `init`, `commit`, `log`, `branch`, `checkout` (with limited functionality for working directory updates);
  - `status` implemented in a generic, file-tree–based way.
- Working directory mapping for generic file trees (no chemistry-specific assumptions).

The following are **not** in scope for Milestones 1–3:

- Remote server and push/pull;
- Domain-specific metadata for chemistry codes;
- HPC schedulers and run lifecycle management.

### 6.2 Scope for Milestone 4–6 (Extension Phases)

Later milestones gradually extend scope:

- Milestone 4:
  - Remote repository server and sync protocol;
  - Remote management commands in the CLI.
- Milestone 5:
  - Python domain layer (`chemvcs_py`);
  - Chemistry code plugins and metadata extraction.
- Milestone 6:
  - HPC scheduler integration;
  - Run lifecycle commands (`run create/submit/status/collect`).

Each step should preserve backward compatibility of the core repository format where feasible.

---

## 7. Constraints and Assumptions

### 7.1 Technical Constraints

- The core VCS must run on typical Linux HPC environments.
- The Go binaries should have minimal external dependencies.
- The core must not depend on a central database; all essential repository data must live within `.chemvcs/`.

### 7.2 Usage Assumptions

- Users are comfortable with basic CLI tools.
- Git remains the primary tool for source code; ChemVCS focuses on data and provenance.
- HPC administrators may impose restrictions on long-running daemons and background services; ChemVCS should work primarily as a batch or interactive CLI tool.

### 7.3 Interoperability Assumptions

- The `.chemvcs` layout and object formats will be documented and stable across minor versions.
- Other tools (Python scripts, web services) can interact via documented CLI commands and/or well-specified file formats.

---

## 8. Success Criteria and Risks

### 8.1 Success Criteria

ChemVCS can be considered successful in its initial stages if:

1. **A single researcher** can:
   - Initialise a repository;
   - Run calculations and record snapshots regularly;
   - Recover and inspect historical results without confusion.

2. **A small research group** can:
   - Share repositories across multiple machines and at least one remote server;
   - Use the same conventions for structures, runs, and workflows;
   - Avoid duplicating or losing important calculation data.

3. **Tool developers** can:
   - Integrate ChemVCS into their scripts and tooling using the CLI and Python layer;
   - Rely on stable object models and repository layouts across versions.

### 8.2 Risks

Key risks include:

1. **Scope creep**
   - Trying to replicate all git features or to cover too many domains too early.
   - Mitigation: keep the core VCS simple, focus on a limited set of operations and clear milestones.

2. **Usability barriers**
   - If the CLI and concepts feel too heavy or complex, users may fall back to ad hoc solutions.
   - Mitigation: provide simple defaults, clear documentation, and predictable behaviour.

3. **Integration complexity with HPC environments**
   - Different clusters and schedulers may impose various constraints.
   - Mitigation: design scheduler integration as a thin, pluggable layer with clear boundaries.

4. **Data model rigidity**
   - Overly rigid designs in the domain layer may not accommodate new codes or workflows.
   - Mitigation: keep the core object model generic and push domain-specific details to metadata and plugins.

---

This vision and scope document provides the overarching direction for ChemVCS. More detailed design and specification documents (architecture overview, object model, CLI specification, remote protocol, domain layer, and HPC adapter design) will build on this foundation.
