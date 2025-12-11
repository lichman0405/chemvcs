# ChemVCS – Python Domain Layer Design (v0.1)

## 1. Introduction

This document describes the design of the Python domain layer for ChemVCS,
tentatively implemented as a package named `chemvcs_py`.

The Python domain layer sits on top of the Go-based VCS core and is responsible
for:

- Expressing computational chemistry concepts (structures, runs, workflows,
  datasets) as domain objects;
- Mapping these domain objects to and from the core `Object` and `Blob` entities
  stored under `.chemvcs/`;
- Providing parsers and utilities for specific electronic-structure codes
  (e.g. VASP, Quantum ESPRESSO);
- Offering a Python API for querying, manipulating, and analysing ChemVCS-backed
  projects;
- Integrating with HPC adapters for run creation, submission, and collection
  (in coordination with the HPC design).

The domain layer is optional from the core’s perspective: the Go core and remote
server can operate without it. However, it is central to ChemVCS’s value for
computational chemistry workflows.

---

## 2. Goals and Non-Goals

### 2.1 Goals

The Python domain layer aims to:

1. Provide a **clear domain model** for computational chemistry projects,
   focused initially on DFT and related electronic-structure calculations.

2. Offer **convenient Python APIs** for:
   - Registering structures and runs in ChemVCS;
   - Inspecting and filtering runs by metadata (code, parameters, status);
   - Building and traversing workflows (DAGs of runs).

3. Enable **robust mapping** between domain entities and core VCS objects:
   - Well-defined conversion functions;
   - Stable metadata schemata in `Object.meta` where it matters;
   - Freedom to evolve domain-specific conventions within a versioned framework.

4. Support **plugin-like extensibility** for new codes and data formats:
   - VASP, QE, CP2K, etc.;
   - Parsers can be distributed separately while sharing the same core abstractions.

### 2.2 Non-Goals

The Python domain layer is **not** intended to:

- Replace the Go core as the source of truth for repository integrity;
- Provide a full workflow engine (scheduling, branching logic, error recovery);
- Dictate a rigid project structure for all users (but it will provide sensible,
  documented conventions).

Workflow orchestration and complex logic remain the responsibility of higher-level
tools or user scripts built on top of `chemvcs_py`.

---

## 3. Package Structure

The Python library is organised as follows (initial proposal):

```text
python/
  chemvcs_py/
    __init__.py
    config.py
    core/
      repo.py           # repository discovery and CLI integration
      objects.py        # domain object wrappers around core objects
      query.py          # query utilities
    domain/
      structures.py     # Structure model and utilities
      runs.py           # Run model and run-related utilities
      workflows.py      # Workflow model and DAG operations
      datasets.py       # Dataset/grouping utilities
    io/
      vasp.py           # VASP-specific parsers and serializers
      qe.py             # Quantum ESPRESSO parsers
      cp2k.py           # CP2K parsers (future)
      formats.py        # CIF/POSCAR/XYZ, etc.
    util/
      logging.py
      errors.py
      types.py
```

This structure separates concerns into:

- `core` – interactions with the Go core and core-level objects;
- `domain` – conceptual models (Structure, Run, Workflow, Dataset);
- `io` – input/output for specific codes and file formats;
- `util` – common utilities, type definitions, and error classes.

The exact structure can be refined during implementation, but the main goal is
to keep domain logic clearly separated from repository integrity mechanisms.

---

## 4. Domain Models

This section defines high-level domain entities and their relationships.

### 4.1 Structure

Represents an atomic structure (crystal, molecule, supercell).

Key attributes:

- Element types and counts;
- Atomic positions;
- Lattice vectors and cell parameters (for periodic systems);
- Symmetry-related information (optional);
- Provenance (source, transformations applied).

Python representation (illustrative):

```python
@dataclass
class Structure:
    id: Optional[str]  # hash of underlying core Object, if persisted
    formula: str
    lattice: Optional[np.ndarray]  # 3x3 matrix or None for molecules
    positions: np.ndarray          # N x 3
    species: List[str]             # length N
    metadata: Dict[str, Any] = field(default_factory=dict)
```

Mapping to core `Object`:

- `type = "structure"`;
- `meta` contains:
  - `formula`;
  - `lattice` (serialised as lists of lists);
  - `positions` (serialised as lists of lists);
  - `species`;
  - additional `metadata` fields;
- `refs` may include:
  - `kind="blob"` pointing to original structure file (CIF/POSCAR/etc.);
  - `kind="object"` linking to derived structures (e.g. relaxed versions).

### 4.2 Run

Represents a single computational run (e.g. a DFT calculation with a given code).

Key attributes:

- Which structure it uses (Structure reference);
- Which code and version;
- Input parameters (numerical, flags, k-point grid, etc.);
- HPC resource request (if applicable);
- Status (planned, submitted, running, finished, failed);
- Outputs (energies, forces, stress, band gap, etc.).

Illustrative Python representation:

```python
@dataclass
class Run:
    id: Optional[str]  # hash of core Object
    structure_id: str
    code: str
    code_version: Optional[str]
    parameters: Dict[str, Any]
    resources: Dict[str, Any]
    status: str                # e.g. "planned", "running", "finished", "failed"
    results: Dict[str, Any]    # energies, band gaps, etc.
    metadata: Dict[str, Any] = field(default_factory=dict)
```

Mapping to core `Object`:

- `type = "run"`;
- `meta` contains `code`, `code_version`, `parameters`, `resources`, `status`,
  `results`, plus `metadata`;
- `refs` include:
  - `structure` object;
  - `input` blobs (INCAR, POSCAR, KPOINTS, etc.);
  - `output` blobs (OUTCAR, DOSCAR, etc.) or folder objects representing outputs;
  - possibly other linked runs (e.g. ancestor runs in a workflow).

### 4.3 Workflow

Represents a workflow as a DAG of runs and sub-workflows (e.g. relax → SCF → band).

Key attributes:

- Set of nodes (runs or sub-workflows);
- Directed edges representing dependencies;
- Workflow type (e.g. “standard DFT workflow”, “phonon calculation”);
- Parameters or templates applied to the workflow as a whole.

Illustrative representation:

```python
@dataclass
class WorkflowNode:
    id: str                   # run or sub-workflow id (core Object hash)
    kind: str                 # "run" or "workflow"
    metadata: Dict[str, Any]

@dataclass
class Workflow:
    id: Optional[str]
    name: str
    nodes: Dict[str, WorkflowNode]
    edges: List[Tuple[str, str]]  # (from_node_id, to_node_id)
    metadata: Dict[str, Any] = field(default_factory=dict)
```

Mapping to core `Object`:

- `type = "workflow"`;
- `meta` contains `name`, `metadata`, and a serialised representation of nodes
  and edges;
- `refs` point to the underlying run or workflow objects for each node.

### 4.4 Dataset

Represents a logical grouping of runs (and/or structures) forming a dataset.

Examples:

- All converged calculations for a certain material system;
- All runs that contribute to a training dataset for ML.

Representation:

```python
@dataclass
class Dataset:
    id: Optional[str]
    name: str
    run_ids: List[str]
    structure_ids: List[str]
    metadata: Dict[str, Any] = field(default_factory=dict)
```

Mapping to core `Object`:

- `type = "dataset"`;
- `meta` includes `name` and `metadata`;
- `refs` list referenced run and structure objects by hash.

---

## 5. Mapping Between Domain Objects and Core Objects

The domain layer must provide stable, well-defined conversion between Python
domain models and core `Object` / `Blob` entities.

### 5.1 Design Principles

- Domain ↔ core mappings should be **pure** functions where possible:
  - Given a domain object, produce deterministic core representations;
  - Given a core object (with appropriate `type`), reconstruct the domain object.

- All non-trivial fields should be stored in `meta`, using JSON-serialisable types.

- References (`refs`) should be used to express:
  - Entity relationships (e.g. run → structure);
  - File attachments (e.g. run → output blobs);
  - Workflow topology (via workflow node refs).

### 5.2 Core Conversion APIs

Proposed functions in `chemvcs_py.core.objects`:

```python
def structure_to_object(struct: Structure) -> CoreObject: ...
def object_to_structure(obj: CoreObject) -> Structure: ...

def run_to_object(run: Run) -> CoreObject: ...
def object_to_run(obj: CoreObject) -> Run: ...

def workflow_to_object(wf: Workflow) -> CoreObject: ...
def object_to_workflow(obj: CoreObject) -> Workflow: ...

def dataset_to_object(ds: Dataset) -> CoreObject: ...
def object_to_dataset(obj: CoreObject) -> Dataset: ...
```

`CoreObject` here is a Python representation mirroring the Go `Object` struct
with fields `id`, `type`, `meta`, and `refs`.

### 5.3 Blob Management

The domain layer may need to create or read blobs for:

- Structure files (POSCAR, CIF);
- Input decks (INCAR, KPOINTS, QE input files);
- Output files (OUTCAR, logs, etc.).

Blob operations are mediated through repository-level APIs, e.g. in `core.repo`
or via CLI calls:

```python
class Repo:
    def put_blob(self, data: bytes) -> str: ...
    def get_blob(self, blob_id: str) -> bytes: ...

    def put_object(self, obj: CoreObject) -> str: ...
    def get_object(self, obj_id: str) -> CoreObject: ...
```

For each domain entity, the domain layer decides which parts are stored as meta
vs. as attached blobs.

---

## 6. Interaction with the Go Core

### 6.1 Repository Discovery

The Python layer needs to locate and operate on the associated ChemVCS repository.
It follows the same discovery rules as the CLI:

1. Use an explicit path if provided (`Repo(path=...)`);
2. Otherwise, search upwards from the current working directory for `.chemvcs`;
3. Optionally honour an environment variable (e.g. `CHEMVCS_REPO`).

`chemvcs_py.core.repo` may wrap these rules in a `Repo` class:

```python
class Repo:
    def __init__(self, root: Optional[Path] = None): ...
    @property
    def path(self) -> Path: ...
```

### 6.2 CLI Integration vs Direct File Access

There are two main interaction patterns:

1. **Via CLI (subprocess)**:

   - Python invokes `chemvcs` commands (`commit`, `log`, etc.) as subprocesses;
   - Captures output (e.g. JSON from future `--format` options) for integration.

   Pros:
   - Keeps Python code decoupled from low-level details of `.chemvcs` layout;
   - Reuses CLI’s validation and behaviour.

   Cons:
   - Overhead of subprocess calls;
   - Parsing text or JSON output.

2. **Direct file access**:

   - Python directly reads/writes `.chemvcs/objects/*`, refs, etc.;
   - Uses the object model specification and simple helper functions.

   Pros:
   - Lower overhead for bulk operations;
   - More control for advanced tooling.

   Cons:
   - Requires careful adherence to core invariants;
   - Risk of divergence if the Go core changes format.

In the MVP, the domain layer can adopt a **hybrid approach**:

- Use CLI for high-level operations (commits, branch management);
- Use direct access for read-only queries and object introspection;
- Provide a small, well-tested abstraction (`Repo`) that centralises format
  assumptions.

### 6.3 Committing Domain Changes

Typical pattern for committing changes from Python:

1. Construct or modify domain objects (`Structure`, `Run`, etc.).
2. Convert them to core `Object`s and `Blob`s using conversion APIs.
3. Store these objects/blobs in the repository (via `Repo.put_object`/`put_blob`).
4. Build or update a root object representing the project/file tree.
5. Instruct the Go core to create a new snapshot:
   - Either by calling `chemvcs commit` with an appropriate message; or
   - By writing the snapshot JSON and updating refs directly (advanced mode).

The recommended approach for early versions is to call the CLI for snapshot
creation, letting the Go core enforce integrity rules.

---

## 7. Code Parsers and IO Layer

### 7.1 Goals

The `io` subpackage provides parsers and serializers for:

- Electronic-structure codes (VASP, QE, CP2K, etc.);
- Structure formats (POSCAR, CIF, XYZ, etc.).

These parsers are used to:

- Extract domain metadata from raw input/output files;
- Populate `Run.parameters`, `Run.results`, `Structure` fields;
- Optionally generate input decks from domain models.

### 7.2 Parser Interface

A generic parser interface for runs might look like:

```python
class RunParser(Protocol):
    code: str

    def parse_inputs(self, run_dir: Path) -> Dict[str, Any]: ...
    def parse_outputs(self, run_dir: Path) -> Dict[str, Any]: ...
    def summarise_results(self, outputs: Dict[str, Any]) -> Dict[str, Any]: ...
```

Code-specific modules (e.g. `vasp.py`, `qe.py`) implement these methods.

The domain layer uses these parsers to construct `Run` objects from existing
run directories or from HPC job outputs.

### 7.3 Structure Format Utilities

For structures:

```python
class StructureFormatHandler(Protocol):
    def load(self, path: Path) -> Structure: ...
    def dump(self, structure: Structure, path: Path) -> None: ...
```

Handlers may support:

- POSCAR/CONTCAR;
- CIF;
- XYZ and others.

Internally, the domain layer uses these handlers when importing existing
structures into ChemVCS or exporting structures for further calculations.

### 7.4 Plugin Registration

To support new codes and formats without modifying the core package, a simple
registry mechanism can be used:

```python
RUN_PARSERS: Dict[str, RunParser] = {}
STRUCTURE_FORMATS: Dict[str, StructureFormatHandler] = {}

def register_run_parser(parser: RunParser) -> None: ...
def get_run_parser(code: str) -> Optional[RunParser]: ...

def register_structure_format(name: str, handler: StructureFormatHandler) -> None: ...
def get_structure_format(name: str) -> Optional[StructureFormatHandler]: ...
```

Third-party modules can register custom parsers at import time.

---

## 8. Query and Analysis APIs

### 8.1 Query Goals

The domain layer should provide convenient functions to:

- Retrieve all runs in a repository;
- Filter runs based on code, parameter ranges, status, or structure attributes;
- Traverse workflows and datasets.

These operations may rely on:

- Walking the object graph via core objects; and/or
- Maintaining a lightweight index (e.g. SQLite) for faster queries (future).

### 8.2 Example Query API

In `core.query` or a dedicated module:

```python
def list_runs(repo: Repo, code: Optional[str] = None) -> List[Run]: ...

def find_runs(
    repo: Repo,
    filters: Dict[str, Any]
) -> List[Run]: ...

def get_workflow_runs(repo: Repo, workflow_id: str) -> List[Run]: ...
```

Filters could support simple operators (equals, in-range) for fields like:

- `code`, `status`;
- `parameters.encut` or `parameters.kpoints`;
- `results.energy`, etc.

The MVP may implement queries by scanning all objects in the repository and
filtering in memory; later versions can add indexing.

### 8.3 Integration with Notebooks and Analysis

Because the domain layer is Python-based, it can integrate directly with:

- Jupyter notebooks for exploratory analysis;
- DataFrames (pandas) for tabular operations;
- Plotting libraries for visualising trends across runs.

The domain layer itself does not need to depend heavily on these libraries,
but it should return data in forms that are easy to consume (lists of dataclasses,
dicts, etc.).

---

## 9. Configuration and Environment

### 9.1 Repository Configuration

The domain layer may read and write configuration values in the repository’s
`config` file (or a sub-namespace of it), such as:

- Default code and code version for new runs;
- Default structure format for import/export;
- Paths or patterns for locating run directories.

Configuration management should be conservative and documented, avoiding
unexpected changes to core behaviour.

### 9.2 Global and User Configuration

In addition to repository-level config, user-level configuration files (e.g. in
the home directory) may store:

- Default HPC scheduler;
- Default remote for publishing;
- Default run templates.

The exact mechanism (e.g. `~/.chemvcs/config`) can be introduced incrementally.

---

## 10. Error Handling and Logging

### 10.1 Error Classes

The domain layer should define clear exception classes for common failure modes, e.g.:

- `ChemVCSDomainError` – base class;
- `RepositoryNotFoundError`;
- `ObjectDecodingError`;
- `ParserError`;
- `ConfigError`.

These exceptions distinguish domain-level issues from lower-level I/O errors,
which can be propagated as `OSError` or similar.

### 10.2 Logging

The domain layer may use Python’s `logging` module with a dedicated namespace
(e.g. `chemvcs_py`). Logging should be:

- Minimal by default (info and above);
- Configurable via environment or config file (debug for troubleshooting).

---

## 11. Versioning and Compatibility

### 11.1 Schema Evolution

Domain schemas (e.g. the contents of `Object.meta` for `run` or `structure`) will
evolve over time. To manage this:

- Metadata may include a `schema_version` field for each domain type;
- Conversion functions should be backward compatible:
  - Able to read older versions of metadata and upgrade to the current model;
  - Avoid breaking changes where possible.

### 11.2 Python Package Versioning

The `chemvcs_py` package versioning should be:

- Semantic where feasible (MAJOR.MINOR.PATCH);
- Coordinated with core ChemVCS versions when domain changes depend on core changes;
- Documented in a changelog, especially for schema changes.

---

## 12. Summary

The Python domain layer (`chemvcs_py`) provides the computational chemistry
semantics for ChemVCS. It defines domain models (Structure, Run, Workflow,
Dataset), maps them to the underlying VCS object model, integrates code-specific
parsers, and offers query and analysis APIs. It interacts with the Go core via
CLI and/or direct file access, ensuring that domain-level operations respect
repository integrity constraints.

This design aims to balance:

- A clear and structured domain model;
- Flexibility for different codes and workflows;
- Minimal coupling to the core so that both can evolve at an appropriate pace.

Subsequent implementation work will refine the exact APIs and structures, but
this document provides the initial blueprint for the Python layer.
