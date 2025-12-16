# HPC Workflow Example

This example demonstrates a complete computational chemistry workflow using ChemVCS M6 HPC integration (All Phases Complete).

## Overview

**Goal**: Perform geometry optimization and single-point energy calculation for a water molecule using VASP.

**Workflow**:
1. Create structure from XYZ file
2. Create relaxation Run
3. Submit to HPC scheduler (SLURM/PBS/LSF)
4. Monitor job status with real-time updates
5. Retrieve results (or cancel if needed)
6. Create single-point Run
7. Submit and retrieve

## Files

- `water.xyz` - Water molecule structure
- `vasp_relax.slurm` - Geometry optimization job script
- `vasp_static.slurm` - Single-point calculation job script
- `workflow.py` - Complete Python workflow
- `README.md` - This file

## Prerequisites

- ChemVCS v0.6.0 (M6 Complete)
- Access to HPC cluster (SLURM, PBS, or LSF)
- VASP installed

## Quick Start

### 1. Initialize Repository

```bash
chemvcs init water-calc
cd water-calc
cp /path/to/examples/hpc-workflow/* .
```

### 2. Run Workflow

```python
python workflow.py
```

The script will:
- Create structure and relaxation Run
- Submit to SLURM
- Wait for completion (polling every 60s)
- Retrieve results
- Create and submit single-point Run
- Retrieve final results

### 3. View History

```bash
chemvcs log
```

You should see commits for:
- Initial structure
- Relaxation plan
- Retrieved relaxation results
- Single-point plan
- Retrieved final results

## Detailed Steps

### Step 1: Create Structure

```python
from chemvcs_py import Repo
from chemvcs_py.io import read_xyz

repo = Repo('.')
structure = read_xyz('water.xyz')
obj = structure.to_core_object()
struct_hash = repo.store_object(obj)
repo.commit('Add water molecule structure')
```

### Step 2: Create Relaxation Run

```python
from chemvcs_py.domain import Run

relax_run = Run(
    structure_id=struct_hash,
    code='VASP',
    code_version='6.3.0',
    parameters={
        'IBRION': 2,    # Conjugate gradient
        'NSW': 50,      # Max ionic steps
        'ENCUT': 400,   # Plane wave cutoff
        'EDIFF': 1e-6,  # Electronic convergence
        'ISMEAR': 0,    # Gaussian smearing
        'SIGMA': 0.05
    },
    status='planned'
)

obj = relax_run.to_core_object()
relax_hash = repo.store_object(obj)
repo.commit('Add geometry relaxation plan')
```

### Step 3: Submit to HPC Scheduler

**Option A: Using SLURM**
```python
from chemvcs_py.hpc import SlurmAdapter, JobSubmitter

adapter = SlurmAdapter()
submitter = JobSubmitter(repo, adapter)

relax_job_id = submitter.submit_run(
    relax_run,
    'vasp_relax.slurm',
    job_name='water-relax'
)

print(f"Submitted relaxation: job {relax_job_id}")
```

**Option B: Using PBS/Torque**
```python
from chemvcs_py.hpc import PbsAdapter, JobSubmitter

adapter = PbsAdapter()
submitter = JobSubmitter(repo, adapter)

relax_job_id = submitter.submit_run(
    relax_run,
    'vasp_relax.pbs',
    job_name='water-relax',
    queue='batch'
)
```

**Option C: Using LSF**
```python
from chemvcs_py.hpc import LsfAdapter, JobSubmitter

adapter = LsfAdapter()
submitter = JobSubmitter(repo, adapter)

relax_job_id = submitter.submit_run(
    relax_run,
    'vasp_relax.lsf',
    job_name='water-relax',
    cores=16  # LSF uses cores instead of nodes
)
```

**Option D: Using CLI**
```bash
# Submit using CLI (auto-detects scheduler)
chemvcs submit vasp_relax.slurm --job-name water-relax

# View all jobs
chemvcs jobs

# Check specific job
chemvcs jobs 12345
```

### Step 4: Monitor Status

**Option A: Using Python JobTracker**
```python
from chemvcs_py.hpc import JobTracker
import time

tracker = JobTracker(repo, adapter)

while True:
    status = adapter.get_status(relax_job_id)
    print(f"Job {relax_job_id}: {status.value}")
    
    if status.value in ['completed', 'failed', 'cancelled']:
        break
    
    time.sleep(60)  # Check every minute
```

**Option B: Using JobWatcher (M6 Phase 4)**
```python
from chemvcs_py.hpc import JobWatcher

watcher = JobWatcher(repo, adapter)

# Watch with real-time updates (emoji indicators)
final_status = watcher.watch_job(
    relax_job_id,
    interval=30,     # Check every 30 seconds
    timeout=3600     # Give up after 1 hour
)

print(f"Final status: {final_status}")
```

**Option C: Using CLI watch command**
```bash
# Interactive monitoring with timestamps
chemvcs watch 12345 --interval 30 --timeout 3600

# Output shows real-time updates:
# [2025-12-15 10:30:00] Job 12345: ⏳ PENDING
# [2025-12-15 10:30:30] Job 12345: ⏳ PENDING
# [2025-12-15 10:31:00] Job 12345: ⏳ RUNNING
# [2025-12-15 10:45:00] Job 12345: ✅ COMPLETED
```

### Step 4b: Cancel Job (if needed)

**Using Python**
```python
from chemvcs_py.hpc import JobTracker

tracker = JobTracker(repo, adapter)

# Cancel specific job
success = tracker.cancel_job(relax_job_id)

# Or cancel by run hash
success = tracker.cancel_run(relax_hash)
```

**Using CLI**
```bash
# Cancel by job ID
chemvcs cancel 12345

# Cancel by run hash
chemvcs cancel abc123def456
```

### Step 5: Retrieve Results

```python
from chemvcs_py.hpc import JobRetriever

retriever = JobRetriever(repo, adapter)

relax_run = retriever.retrieve_results(
    job_id=relax_job_id,
    output_files=['OUTCAR', 'CONTCAR', 'vasprun.xml'],
    output_dir='relax_outputs/',
    auto_commit=True
)

print(f"Relaxation energy: {relax_run.get_energy()} eV")
```

### Step 6: Create Single-Point Run

```python
# Read relaxed structure
from chemvcs_py.io import read_poscar

relaxed = read_poscar('relax_outputs/CONTCAR')
obj = relaxed.to_core_object()
relaxed_hash = repo.store_object(obj)

# Create single-point run
static_run = Run(
    structure_id=relaxed_hash,
    code='VASP',
    parameters={
        'IBRION': -1,   # No relaxation
        'NSW': 0,
        'ENCUT': 400,
        'EDIFF': 1e-7,  # Tighter convergence
        'ISMEAR': 0,
        'SIGMA': 0.05
    },
    status='planned'
)

obj = static_run.to_core_object()
static_hash = repo.store_object(obj)
repo.commit('Add single-point calculation plan')
```

### Step 7: Submit and Retrieve

```python
static_job_id = submitter.submit_run(
    static_run,
    'vasp_static.slurm',
    job_name='water-static'
)

# Wait and retrieve (same as steps 4-5)
# ...

print(f"Final energy: {static_run.get_energy()} eV")
```

## Provenance Captured

For each Run, ChemVCS automatically captures:

```python
relax_run.job_id           # '12345'
relax_run.job_system       # 'slurm' / 'pbs' / 'lsf'
relax_run.modules_loaded   # ['vasp/6.3.0', 'intel/2021.4']
relax_run.job_resources    # {'nodes': 2, 'ntasks_per_node': 28, ...}
relax_run.submit_script    # Full script content
relax_run.metadata         # Timestamps for all stages
```

## Viewing Provenance

```python
# Load Run from repository
obj = repo.get_object(relax_hash)
run = Run.from_core_object(obj)

print(f"Job ID: {run.job_id}")
print(f"Scheduler: {run.job_system}")
print(f"Modules: {run.modules_loaded}")
print(f"Submitted: {run.metadata['submitted_at']}")
print(f"Finished: {run.metadata['finished_at']}")
print(f"Resources: {run.job_resources}")
```

## Error Handling

```python
from chemvcs_py.hpc.exceptions import (
    JobSubmissionError,
    InvalidJobStateError,
    JobCancellationError
)

try:
    job_id = submitter.submit_run(run, 'job.slurm')
except JobSubmissionError as e:
    print(f"Submission failed: {e}")
except InvalidJobStateError as e:
    print(f"Invalid state: {e}")

# Cancel with error handling
try:
    tracker.cancel_job(job_id)
except JobCancellationError as e:
    print(f"Cancellation failed: {e}")
```

## Testing Without HPC Cluster

Use mock adapter for local testing:

```python
from unittest.mock import Mock
from chemvcs_py.hpc import JobAdapter, JobStatus

class MockAdapter(JobAdapter):
    def __init__(self):
        self._jobs = {}
    
    @property
    def name(self):
        return "mock"
    
    def submit(self, script_path, **kwargs):
        job_id = f"mock-{len(self._jobs)}"
        self._jobs[job_id] = JobStatus.RUNNING
        return job_id
    
    def get_status(self, job_id):
        return self._jobs.get(job_id, JobStatus.UNKNOWN)
    
    def get_info(self, job_id):
        return Mock(job_id=job_id, status=self.get_status(job_id))
    
    def cancel(self, job_id):
        return True

# Use mock instead of SlurmAdapter
adapter = MockAdapter()
submitter = JobSubmitter(repo, adapter)
```

## Next Steps

- Modify scripts for your HPC system
- Add more complex workflows with dependencies
- Integrate with workflow managers (Snakemake, Nextflow)
- Add automated result parsing

## Troubleshooting

**Job submission fails**: Check scheduler-specific requirements (PBS needs queue, LSF needs cores)

**Module capture empty**: Source module init in job script

**Cannot find completed job**: Job may be purged, check scheduler history (sacct/bjobs -a/qstat -x)

**Watch command hangs**: Set reasonable timeout, check network connectivity

## CLI Quick Reference

```bash
# Submit job
chemvcs submit script.slurm --job-name my-job

# List all jobs
chemvcs jobs

# Watch job (real-time)
chemvcs watch 12345 --interval 30

# Cancel job
chemvcs cancel 12345

# Retrieve results
chemvcs retrieve 12345 --output-dir results/
```

## References

- [HPC User Guide](../../docs/10-hpc-user-guide.md)
- [M6 Design Doc](../../docs/09-hpc-integration-design.md)
- [VASP Manual](https://www.vasp.at/wiki/index.php/The_VASP_Manual)
