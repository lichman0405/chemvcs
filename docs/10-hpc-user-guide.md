# ChemVCS HPC Integration User Guide

**Version**: 0.6.0 (Milestone 6)  
**Date**: 2025-12-15

---

## Overview

ChemVCS M6 integrates version control with HPC job schedulers, enabling you to:

- Submit computational runs directly from your repository
- Track job status across your project history
- Automatically capture computational provenance
- Retrieve results and create atomic commits

This guide covers installation, configuration, and common workflows.

---

## Prerequisites

### Required
- ChemVCS v0.6.0 or later
- Python 3.8+ with chemvcs_py installed
- Access to an HPC cluster (optional for local testing)

### Supported Schedulers
- **SLURM** (fully implemented) - sbatch, squeue, sacct, scancel
- **PBS/Torque** (fully implemented) - qsub, qstat, qdel
- **LSF** (fully implemented) - bsub, bjobs, bkill
- **SGE** (script parsing only, adapter TBD)

---

## Installation

### 1. Install ChemVCS

```bash
# Build Go binary
cd go
go build -o chemvcs ./cmd/chemvcs

# Install Python package
cd ../python
pip install -e .
```

### 2. Verify Installation

```bash
# Check CLI version
chemvcs --version

# Test Python module
python -c "from chemvcs_py.hpc import SlurmAdapter; print('HPC module OK')"
```

---

## Quick Start

This guide shows two ways to use HPC integration: **CLI commands** (simple) and **Python API** (advanced).

### Option 1: Using CLI Commands (Recommended)

**1. Create a Run in the repository** (using Python API or manually)

```python
from chemvcs_py import Repo
from chemvcs_py.domain import Run, RunStatus

repo = Repo('.')
run = Run(
    name="water-relax",
    status=RunStatus.PLANNED,
    structure_hash="abc123...",  # Your structure hash
    command="mpirun vasp_std",
    working_dir="/scratch/user/water"
)
run_hash = repo.add_run(run)
print(f"Run created: {run_hash}")
```

**2. Submit job via CLI**

```bash
chemvcs submit <run-hash> vasp.slurm

# With options:
chemvcs submit <run-hash> vasp.slurm --capture-env=false
```

**3. List tracked jobs**

```bash
# List all jobs
chemvcs jobs

# Filter by status
chemvcs jobs --status=RUNNING
chemvcs jobs --status=COMPLETED

# Verbose output
chemvcs jobs -v
```

**4. Retrieve results**

```bash
# Retrieve all files
chemvcs retrieve <run-hash>

# Retrieve specific patterns
chemvcs retrieve <run-hash> --patterns="*.out,*.log"

# Specify destination
chemvcs retrieve <run-hash> --dest=./results
```

### Option 2: Using Python API (Advanced)

**1. Create a structure**

```python
from chemvcs_py import Repo, Structure
from chemvcs_py.io import read_xyz

repo = Repo('.')
structure = read_xyz('water.xyz')
obj = structure.to_core_object()
struct_hash = repo.store_object(obj)
print(f"Structure: {struct_hash}")
```

**2. Create a Run object**

```python
from chemvcs_py.domain import Run

run = Run(
    structure_id=struct_hash,
    code='VASP',
    code_version='6.3.0',
    parameters={
        'ENCUT': 400,
        'ISMEAR': 0,
        'EDIFF': 1e-6
    },
    status='planned'
)

obj = run.to_core_object()
run_hash = repo.store_object(obj)
repo.commit('Add VASP single-point calculation')
print(f"Run: {run_hash}")
```

**3. Submit to HPC Scheduler**

```python
from chemvcs_py.hpc import SlurmAdapter, PbsAdapter, LsfAdapter, JobSubmitter

# Choose adapter based on your cluster:
# For SLURM:
adapter = SlurmAdapter()
# For PBS/Torque:
# adapter = PbsAdapter()
# For LSF:
# adapter = LsfAdapter()

submitter = JobSubmitter(repo, adapter)

# Submit (captures environment automatically)
job_id = submitter.submit_run(run, 'vasp.slurm')
print(f"Submitted job {job_id}")
```

**4. Check job status**

```python
from chemvcs_py.hpc import JobTracker

tracker = JobTracker(repo, adapter)
jobs = tracker.list_jobs()

for job in jobs:
    print(f"Job {job.job_id}: {job.status.value} ({job.code})")
```

**5. Retrieve results**

```python
from chemvcs_py.hpc import JobRetriever

retriever = JobRetriever(repo, adapter)

# Wait for job to complete, then:
run = retriever.retrieve_results(
    job_id=job_id,
    output_files=['OUTCAR', 'OSZICAR', 'vasprun.xml'],
    auto_commit=True
)

print(f"Energy: {run.get_energy()} eV")
```

---

## Job Scheduler Scripts

### SLURM Script Example

```bash
#!/bin/bash
#SBATCH --job-name=vasp_water
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=28
#SBATCH --time=24:00:00
#SBATCH --partition=batch
#SBATCH --output=vasp.out

# Load modules
module load vasp/6.3.0
module load intel/2021.4

# Run calculation
mpirun -np 56 vasp_std
```

ChemVCS automatically:
- Captures loaded modules (`vasp/6.3.0`, `intel/2021.4`)
- Parses resource requests (2 nodes, 28 tasks/node, 24h walltime)
- Saves complete script snapshot

### PBS/Torque Script Example

```bash
#!/bin/bash
#PBS -N vasp_water
#PBS -l nodes=2:ppn=16
#PBS -l walltime=12:00:00
#PBS -q batch
#PBS -o vasp.out

cd $PBS_O_WORKDIR
module load vasp
mpirun vasp_std
```

### LSF Script Example

```bash
#!/bin/bash
#BSUB -J vasp_water
#BSUB -n 32
#BSUB -W 12:00
#BSUB -q normal
#BSUB -o vasp.out

module load vasp
mpirun -np 32 vasp_std
```

---

## Command Reference

### CLI Commands

#### chemvcs submit

Submit an HPC job for a computational run.

**Syntax:**
```bash
chemvcs submit <run-hash> <script-path> [--capture-env=<bool>]
```

**Arguments:**
- `<run-hash>`: Hash of the Run object to submit (must exist in repository)
- `<script-path>`: Path to the SLURM/PBS submission script

**Options:**
- `--capture-env=<bool>`: Capture environment provenance (default: true)

**Examples:**
```bash
# Submit with automatic environment capture
chemvcs submit a1b2c3d4 vasp_relax.slurm

# Submit without environment capture
chemvcs submit a1b2c3d4 vasp_relax.slurm --capture-env=false
```

**Output:**
```
Submitting HPC job for run a1b2c3d4...
Job submitted successfully!
  Job ID: 12345
  Job System: slurm
  Run: a1b2c3d4
```

#### chemvcs jobs

List all tracked HPC jobs in the repository.

**Syntax:**
```bash
chemvcs jobs [--status=<status>] [-v]
```

**Options:**
- `--status=<status>`: Filter by job status (PENDING, RUNNING, COMPLETED, FAILED, CANCELLED)
- `-v`: Verbose output (show submission time, queue, etc.)

**Examples:**
```bash
# List all jobs
chemvcs jobs

# List only running jobs
chemvcs jobs --status=RUNNING

# List completed jobs with details
chemvcs jobs --status=COMPLETED -v
```

**Output:**
```
Found 3 job(s):

Job ID: 12345
  Run:        a1b2c3d4
  Status:     RUNNING
  System:     slurm
  Queue:      batch

Job ID: 12346
  Run:        e5f6g7h8
  Status:     COMPLETED
  System:     slurm
```

#### chemvcs retrieve

Retrieve results from a completed HPC job.

**Syntax:**
```bash
chemvcs retrieve <run-hash> [--patterns=<patterns>] [--dest=<path>]
```

**Arguments:**
- `<run-hash>`: Hash of the Run object whose results to retrieve

**Options:**
- `--patterns=<patterns>`: Comma-separated glob patterns for files to retrieve (e.g., `*.out,*.log`)
- `--dest=<path>`: Destination directory (default: current directory)

**Examples:**
```bash
# Retrieve all files to current directory
chemvcs retrieve a1b2c3d4

# Retrieve specific file types
chemvcs retrieve a1b2c3d4 --patterns="OUTCAR,OSZICAR,*.xml"

# Retrieve to specific directory
chemvcs retrieve a1b2c3d4 --patterns="*.out" --dest=./results
```

**Output:**
```
Retrieving results for run a1b2c3d4...
Retrieved 3 file(s):
  ./OUTCAR
  ./OSZICAR
  ./vasprun.xml
```

---

### Python API

#### JobSubmitter

```python
from chemvcs_py.hpc import JobSubmitter

submitter = JobSubmitter(repo, adapter)

# Submit run
job_id = submitter.submit_run(
    run,
    script_path='job.slurm',
    capture_env=True,  # Capture provenance (default: True)
    partition='batch',  # SLURM partition
    job_name='my_calc'  # Override job name
)
```

#### JobTracker

```python
from chemvcs_py.hpc import JobTracker

tracker = JobTracker(repo, adapter)

# List all jobs
jobs = tracker.list_jobs()

# Filter by status
running = tracker.list_jobs(status_filter=['running'])

# Limit results
recent = tracker.list_jobs(limit=10)

# Find specific job
job = tracker.get_job_by_id('12345')
if job:
    print(f"Status: {job.status.value}")
```

#### JobRetriever

```python
from chemvcs_py.hpc import JobRetriever

retriever = JobRetriever(repo, adapter)

# Retrieve completed job
run = retriever.retrieve_results(
    job_id='12345',
    output_files=['OUTCAR'],  # Files to copy
    output_dir='results/',    # Destination
    auto_commit=True          # Create commit
)
```

---

## Provenance Tracking

ChemVCS automatically captures:

### 1. Environment Modules

```python
run.modules_loaded
# ['vasp/6.3.0', 'intel/2021.4', 'mkl/2021.4']
```

### 2. Environment Variables

```python
run.environment_vars
# {'OMP_NUM_THREADS': '4', 'MKL_NUM_THREADS': '4', ...}
```

### 3. Resource Requirements

```python
run.job_resources
# {
#   'nodes': 2,
#   'ntasks_per_node': 28,
#   'walltime': '24:00:00',
#   'memory': '64G',
#   'partition': 'batch'
# }
```

### 4. Job Script Snapshot

```python
run.submit_script
# Full script content preserved
```

### 5. Timestamps

```python
run.metadata
# {
#   'submitted_at': '2025-12-15T10:30:00Z',
#   'queued_at': '2025-12-15T10:30:05Z',
#   'started_at': '2025-12-15T10:35:00Z',
#   'finished_at': '2025-12-15T14:20:00Z',
#   'retrieved_at': '2025-12-15T14:25:00Z'
# }
```

---

## Advanced Usage

### Custom Adapters

Create adapter for other schedulers:

```python
from chemvcs_py.hpc import JobAdapter, JobStatus, JobInfo

class CustomAdapter(JobAdapter):
    @property
    def name(self) -> str:
        return "custom"
    
    def submit(self, script_path: str, **kwargs) -> str:
        # Implement submission
        pass
    
    def get_status(self, job_id: str) -> JobStatus:
        # Implement status query
        pass
    
    # Implement other methods...
```

### Selective Environment Capture

```python
from chemvcs_py.hpc.provenance import EnvironmentCapture

# Capture specific variables only
env_vars = EnvironmentCapture.capture_env_vars([
    'VASP_PP_PATH',
    'VASPRUN_COMMAND'
])
```

### Parse Script Without Submission

```python
from chemvcs_py.hpc.provenance import EnvironmentCapture

resources = EnvironmentCapture.parse_slurm_script('job.slurm')
print(f"Will request {resources['nodes']} nodes")
```

---

## Troubleshooting

### Issue: "sbatch command not found"

**Problem**: SLURM not available

**Solution**: 
- Check if on compute node: `which sbatch`
- Use mock adapter for testing (see Testing section)

### Issue: "Job 12345 not found"

**Problem**: Job purged from scheduler

**Solution**:
- Check Run object status: `run.status`
- Jobs may be purged after completion
- Use `sacct` for historical data

### Issue: Module capture returns empty list

**Problem**: `module` command not available

**Solution**:
- Ensure environment modules are installed
- Source module init: `source /etc/profile.d/modules.sh`
- Or manually set `run.modules_loaded`

---

## Testing Without SLURM

Use mock adapter for development/testing:

```python
from unittest.mock import Mock
from chemvcs_py.hpc import JobAdapter, JobStatus

class MockAdapter(JobAdapter):
    @property
    def name(self):
        return "mock"
    
    def submit(self, script_path, **kwargs):
        return "mock-12345"
    
    def get_status(self, job_id):
        return JobStatus.COMPLETED
    
    def get_info(self, job_id):
        return Mock(
            job_id=job_id,
            status=JobStatus.COMPLETED,
            exit_code=0
        )
    
    def cancel(self, job_id):
        return True

# Use mock adapter
adapter = MockAdapter()
submitter = JobSubmitter(repo, adapter)
```

---

## Best Practices

### 1. Always Commit Before Submitting

```python
# Bad
run = Run(...)
submitter.submit_run(run, 'job.slurm')  # Run not in repo

# Good
run = Run(...)
obj = run.to_core_object()
run_hash = repo.store_object(obj)
repo.commit('Add calculation')
submitter.submit_run(run, 'job.slurm')
```

### 2. Use Descriptive Job Names

```python
job_id = submitter.submit_run(
    run,
    'vasp.slurm',
    job_name=f'water-opt-{struct_hash[:8]}'
)
```

### 3. Store Job Scripts in Repository

```
project/
├── structures/
│   └── water.xyz
├── scripts/
│   ├── vasp_relax.slurm
│   └── vasp_static.slurm
└── .chemvcs/
```

### 4. Use Status Filters for Monitoring

```python
# Check for failures
failed = tracker.list_jobs(status_filter=['failed'])
for job in failed:
    print(f"Job {job.job_id} failed!")

# Monitor running jobs
running = tracker.list_jobs(status_filter=['running', 'pending'])
print(f"{len(running)} jobs active")
```

### 5. Automate Retrieval

```bash
# Cron job to auto-retrieve completed jobs
*/10 * * * * cd /path/to/repo && python retrieve_completed.py
```

```python
# retrieve_completed.py
from chemvcs_py import Repo
from chemvcs_py.hpc import JobTracker, JobRetriever, SlurmAdapter

repo = Repo('.')
adapter = SlurmAdapter()
tracker = JobTracker(repo, adapter)
retriever = JobRetriever(repo, adapter)

completed = tracker.list_jobs(status_filter=['completed'])
for job in completed:
    if 'retrieved_at' not in job.run.metadata:
        retriever.retrieve_results(job.job_id, auto_commit=True)
        print(f"Retrieved job {job.job_id}")
```

---

## FAQ

**Q: Can I use ChemVCS without HPC?**

A: Yes! The HPC module is optional. Use local adapters or skip job submission entirely.

**Q: Does it work with Windows?**

A: The Python HPC module works cross-platform. SLURM commands require Linux/Unix cluster.

**Q: What about job arrays?**

A: Not yet supported. Submit individual jobs for now. Arrays planned for future release.

**Q: Can I modify a Run after submission?**

A: No, Run objects are immutable after storage. Create a new Run for different parameters.

**Q: How do I clean up old job records?**

A: Use `chemvcs` garbage collection (planned feature) or manually filter old objects.

---

## Next Steps

- See [examples/hpc-workflow/](../../examples/hpc-workflow/) for complete workflow
- Read [Design Document](../09-hpc-integration-design.md) for architecture details
- Check [FEATURE_STATUS.md](../../FEATURE_STATUS.md) for latest features

---

## Support

- **Issues**: https://github.com/lichman0405/chemvcs/issues
- **Discussions**: https://github.com/lichman0405/chemvcs/discussions
- **Documentation**: https://github.com/lichman0405/chemvcs/tree/main/docs

---

**Last Updated**: 2025-12-15  
**ChemVCS Version**: 0.6.0 (M6)
