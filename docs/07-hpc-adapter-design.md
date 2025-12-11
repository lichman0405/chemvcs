# ChemVCS – HPC Adapter and Run Lifecycle Design (v0.1)

## 1. Introduction

This document describes the design of the HPC integration layer for ChemVCS,
with a focus on:

- Scheduler adapters (SLURM, PBS, etc.);
- The run lifecycle (create, submit, monitor, collect);
- Interaction patterns between the HPC layer, the Python domain layer, and the
  Go-based VCS core;
- Configuration, error handling, and extensibility.

The HPC layer is implemented in Python, on top of the `chemvcs_py` domain layer.
It does not modify the core repository semantics, but rather uses them to record
and track the state and provenance of HPC jobs.

---

## 2. Goals and Scope

### 2.1 Goals

The HPC integration layer aims to:

1. Provide a common abstraction for interacting with different batch schedulers
   (SLURM, PBS, etc.).

2. Connect the lifecycle of HPC jobs to ChemVCS run objects, so that:

   - Each run has a corresponding job submission record;
   - Status changes in the scheduler are reflected in run metadata;
   - Outputs collected from HPC are captured as blobs/objects in ChemVCS.

3. Offer a small set of high-level operations that can be used from both:

   - A future CLI command group (`chemvcs run ...`); and
   - Python scripts and notebooks.

### 2.2 Scope (MVP)

The MVP HPC layer focuses on:

- SLURM support (primary target);
- Basic PBS/Torque support where possible;
- Synchronous and polling-based status checks (no long-running daemons);
- Handling single-step runs (one job per run object).

Multi-step workflows, job arrays, and advanced scheduler features can be added
later, leveraging the same abstractions.

---

## 3. Architectural Overview

### 3.1 Layers and Responsibilities

The HPC integration layer sits above the domain layer and uses a scheduler
adapter interface to interact with different batch systems:

- **Domain layer (`chemvcs_py.domain.runs`)**
  - Defines `Run` objects and maps them to/from core objects.
  - Contains logic for representing run parameters and results.

- **HPC layer (`chemvcs_py.hpc`)**
  - Defines the abstract scheduler interface and concrete adapters.
  - Implements high-level run lifecycle operations.

- **Go core / CLI (`chemvcs`)**
  - Ensures repository integrity and snapshot creation.
  - Remains unaware of HPC-specific details.

### 3.2 Data Flow

1. User uses Python or CLI to define a `Run` (code, parameters, resources).
2. Domain layer creates a `Run` object and commits it into ChemVCS.
3. HPC layer uses a scheduler adapter to generate and submit a job script.
4. Scheduler returns a job ID; HPC layer records this in the `Run` metadata and
   commits an updated snapshot.
5. Periodic or on-demand status checks query the scheduler for job state, updating
   run metadata accordingly.
6. When the job finishes, the HPC layer collects outputs, updates the `Run.results`
   and attached blobs, and commits another snapshot.

---

## 4. Scheduler Adapter Interface

### 4.1 Goals

The scheduler adapter interface isolates scheduler-specific details (command
names, submission scripts, environment modules) from the rest of the system.

Key requirements:

- Simple, minimal interface;
- Ability to submit jobs, query status, and cancel jobs;
- Extensible to new schedulers without changing callers.

### 4.2 Python Interface

Define a protocol or base class in `chemvcs_py.hpc.adapters`:

```python
from typing import Protocol, Literal, Optional
from pathlib import Path
from dataclasses import dataclass

JobStatus = Literal["pending", "running", "finished", "failed", "unknown"]

@dataclass
class SubmitResult:
    job_id: str
    raw_output: str

class SchedulerAdapter(Protocol):
    name: str

    def submit(self, script_path: Path) -> SubmitResult:
        \"\"\"Submit a job script to the scheduler and return the job ID.\"\"\"

    def status(self, job_id: str) -> JobStatus:
        \"\"\"Query the scheduler for the status of a job.\"\"\"

    def cancel(self, job_id: str) -> None:
        \"\"\"Cancel a job if it is pending or running.\"\"\"
```

Concrete adapters:

- `SlurmAdapter` – uses `sbatch`, `squeue`, `scancel`;
- `PBSSchedulerAdapter` – uses `qsub`, `qstat`, `qdel`;
- Others can be added later.

### 4.3 Implementation Notes

- Adapters call underlying scheduler commands via `subprocess` and parse outputs.
- Error handling must be careful and robust (e.g. timeouts, malformed outputs).
- Adapter configuration (paths to commands, extra flags) can be provided via
  repository or user-level configuration.

---

## 5. Run Lifecycle

The run lifecycle consists of four primary phases:

1. **Create** – define a run and prepare its working directory.
2. **Submit** – submit the run to the scheduler.
3. **Status** – query the scheduler and update the run's state.
4. **Collect** – gather outputs, parse results, and update run metadata.

### 5.1 Run Create

The “create” step associates a domain `Run` object with a working directory
containing input files and a job script.

Responsibilities:

- Generate input files using the domain IO layer (e.g. VASP input deck);
- Create a job script template with scheduler-specific directives;
- Register the run in ChemVCS.

Illustrative function in `chemvcs_py.hpc.run_lifecycle`:

```python
def create_run(
    repo: Repo,
    run: Run,
    structure: Structure,
    base_dir: Path,
    code_handler: RunParser,
    scheduler: SchedulerAdapter,
) -> Run:
    \"\"\"
    Prepare a run directory and register the run in ChemVCS.

    - Creates a run directory under base_dir (e.g. runs/<run_id>/).
    - Uses code_handler to generate input files.
    - Generates a scheduler job script.
    - Stores input files and metadata as blobs/objects in ChemVCS.
    - Returns an updated Run with references to created objects.
    \"\"\"
```

Key outputs:

- A run directory on disk (e.g. `runs/<run-short-id>/`);
- Input files in that directory;
- A job script (e.g. `run.sh`);
- A `Run` object with metadata including:
  - paths or object references to inputs;
  - scheduler type;
  - initial status (e.g. `"planned"`).

### 5.2 Run Submit

Submits the job script to the scheduler using the adapter.

Illustrative function:

```python
def submit_run(
    repo: Repo,
    run: Run,
    run_dir: Path,
    scheduler: SchedulerAdapter,
) -> Run:
    \"\"\"
    Submit a run's job script to the the scheduler.

    - Locates the job script in run_dir.
    - Calls scheduler.submit(script_path).
    - Updates run.status to "running" or "submitted".
    - Records job_id and submission time in run.metadata.
    - Stores updated Run in ChemVCS as a new object/snapshot.
    \"\"\"
```

Key responsibilities:

- Ensure the run directory exists and contains a valid job script;
- Capture the job ID from the scheduler;
- Write updated run metadata back to ChemVCS via the domain layer;
- Optionally create a dedicated snapshot representing "run submitted".

### 5.3 Run Status

Queries the scheduler for the job's status and updates the `Run` accordingly.

Illustrative function:

```python
def update_run_status(
    repo: Repo,
    run: Run,
    scheduler: SchedulerAdapter,
) -> Run:
    \"\"\"
    Query the scheduler for the run's job status and update run.status.

    - Reads job_id from run.metadata.
    - Calls scheduler.status(job_id).
    - Maps scheduler state to Run.status values.
    - Updates run and commits changes to ChemVCS if state changed.
    \"\"\"
```

Possible status transitions:

- `"planned"` → `"submitted"` → `"running"` → `"finished"`;
- `"planned"` / `"submitted"` / `"running"` → `"failed"` (scheduler error or exit code != 0).

The HPC layer should not automatically delete failed runs; failures are important
for provenance and diagnosis.

### 5.4 Run Collect

After a job finishes, the “collect” step gathers outputs and updates results.

Illustrative function:

```python
def collect_run_outputs(
    repo: Repo,
    run: Run,
    run_dir: Path,
    code_handler: RunParser,
) -> Run:
    \"\"\"
    Collect outputs from run_dir, parse results, and update run metadata.

    - Uses code_handler.parse_outputs(run_dir) to build a results dict.
    - Summarises results via code_handler.summarise_results(outputs).
    - Stores output files as blobs/objects in ChemVCS.
    - Updates run.results and status (e.g. "finished" or "failed").
    - Commits an updated Run object as part of a new snapshot.
    \"\"\"
```

Key tasks:

- Store relevant outputs (e.g. OUTCAR, logs) as blobs or folder objects;
- Extract structured results (energies, gaps, convergence info) into `run.results`;
- Update status and commit to ChemVCS.

---

## 6. CLI Integration (Future `chemvcs run` Commands)

While the HPC layer is primarily implemented in Python, some operations may be
surfaced via CLI commands for users who prefer a command-line workflow.

### 6.1 Command Group: `chemvcs run`

Potential subcommands (not necessarily in the MVP of HPC):

- `chemvcs run create` – create a run and prepare a run directory;
- `chemvcs run submit` – submit a prepared run;
- `chemvcs run status` – query status for one or more runs;
- `chemvcs run collect` – collect outputs for finished runs.

Implementation options:

- The Go CLI (`chemvcs`) invokes Python entry points (e.g. via a small wrapper command);
- Or users call the Python functionality directly (`python -m chemvcs_py.run ...`).

The design should avoid hard dependencies from the Go core on Python tools,
but allow optional integration for convenience.

### 6.2 Identifying Runs from CLI

For CLI operations, runs can be identified by:

- Run ID (hash of the run object);
- A human-readable label stored in `run.metadata` (e.g. `name`);
- Directory path (e.g. `runs/systemA-relax/`).

The mapping between these identifiers and the underlying run object is handled
by the Python layer and recorded as metadata in ChemVCS.

---

## 7. Configuration and Environment

### 7.1 Scheduler Configuration

User or repository configuration should specify:

- Which scheduler adapter to use (e.g. `scheduler = "slurm"`);
- Paths to scheduler commands, if not in `PATH`;
- Default partition/queue, account, time limits;
- Module loading or environment setup snippets (if needed).

Configuration sources:

- Repository-level config (`.chemvcs/config` or a dedicated HPC config file);
- User-level config (e.g. `~/.chemvcs/config`);
- Environment variables (e.g. `CHEMVCS_SCHEDULER`, `CHEMVCS_PARTITION`).

### 7.2 Job Script Templates

Job scripts are typically generated from templates that include:

- Scheduler directives (e.g. `#SBATCH --nodes=...`);
- Environment setup commands (modules, virtualenvs);
- A call to the underlying code (e.g. `mpirun vasp_std`);
- Optional hooks to record additional metadata (e.g. git commit of code).

Templates may be stored:

- In the repository (e.g. `hpc_templates/slurm_run.sh.j2`);
- In user-level configuration directories.

The HPC layer uses a templating engine (e.g. Python `string.Template` or Jinja2)
to fill in resource requests and paths.

---

## 8. Error Handling and Robustness

### 8.1 Scheduler Errors

The HPC layer must handle typical scheduler errors, such as:

- Submission failures (invalid job script, unavailable partition);
- Unknown job IDs (expired or not found);
- Transient communication failures (scheduler command not responding).

Strategies:

- Raise domain-specific exceptions (e.g. `SchedulerError`) with detailed messages;
- Log raw scheduler outputs for debugging;
- Store critical scheduler errors in run metadata for provenance.

### 8.2 Partial Failures

Because ChemVCS emphasises provenance, partial failures are data, not noise.
Examples:

- Job submitted but failed immediately:
  - Run status is `"failed"`;
  - Exit code and relevant logs are collected;
  - A snapshot captures this state.

- Output collection partially succeeds:
  - Some files may fail to parse;
  - The system should record which files were collected and which were not;
  - Parsing errors should be captured in metadata.

### 8.3 Safety with Repository State

The HPC layer must respect core repository invariants:

- Changes to run metadata should be committed through regular snapshot creation;
- Blobs and objects must be written before refs are updated;
- Direct manipulation of `.chemvcs` should go through well-tested helpers where possible.

---

## 9. Extensibility

### 9.1 New Schedulers

Supporting a new scheduler involves:

1. Implementing the `SchedulerAdapter` protocol for the scheduler;
2. Providing configuration keys (e.g. `scheduler = "mybatch"`);
3. Registering the adapter in a scheduler registry, e.g.:

   ```python
   SCHEDULERS: Dict[str, SchedulerAdapter] = {}
   ```

4. Testing submission, status, and cancellation flows end-to-end.

### 9.2 Advanced Features

Future features that can be built on the same architecture include:

- Job arrays for parameter sweeps;
- Multi-step workflows with automatic submission of dependent jobs;
- Event-driven status updates (if schedulers support hooks or callbacks).

These do not change the core abstractions but add higher-level logic on top of
the run lifecycle functions.

---

## 10. Summary

The HPC integration layer connects ChemVCS with batch schedulers through a
simple adapter interface and a well-defined run lifecycle:

- `create` prepares run directories and inputs;
- `submit` launches jobs and records scheduler IDs;
- `status` keeps run metadata in sync with scheduler state;
- `collect` captures outputs and derived results into ChemVCS.

By keeping scheduler-specific logic encapsulated in adapters, and by reusing the
domain models (Run, Structure, Workflow) defined in `chemvcs_py`, the design
remains flexible and extensible across different HPC environments and codes.

The Go core continues to provide robust repository management, while the Python
layer orchestrates HPC interactions and domain-specific semantics on top of it.
