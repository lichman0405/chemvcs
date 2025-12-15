# HPC Integration Design (M6)

**Status**: In Progress  
**Date**: 2025-12-15  
**Version**: 0.6.0

---

## 1. Overview

### 1.1 Motivation

The core value proposition of ChemVCS is to integrate version control with computational workflows. While M1-M5 provide solid VCS foundation and domain objects, they don't yet solve the critical problem: **tracking and managing computational jobs on HPC systems**.

This milestone (M6) implements HPC integration, specifically:
- Submitting computational runs to job schedulers (SLURM)
- Tracking job status and lifecycle
- Capturing computational provenance
- Automatically retrieving results and updating version history

### 1.2 Design Goals

1. **Scheduler Agnostic**: Abstract interface supporting multiple job systems (SLURM, PBS, SGE)
2. **Provenance First**: Capture complete computational environment and job metadata
3. **Testable**: Mock-based testing without requiring actual HPC infrastructure
4. **Incremental**: Works locally (no scheduler) and on HPC clusters
5. **Secure**: Proper credential and permission management

---

## 2. Architecture

### 2.1 Component Overview

```
┌─────────────────────────────────────────────────────────────┐
│                     ChemVCS CLI (Go)                        │
│  chemvcs submit | chemvcs jobs | chemvcs retrieve          │
└─────────────────────┬───────────────────────────────────────┘
                      │
                      ▼
┌─────────────────────────────────────────────────────────────┐
│              Python HPC Layer (chemvcs_py.hpc)              │
│                                                             │
│  ┌──────────────┐  ┌────────────────┐  ┌────────────────┐ │
│  │ JobSubmitter │  │ JobTracker     │  │ ProvCapture    │ │
│  └──────────────┘  └────────────────┘  └────────────────┘ │
│           │               │                    │           │
│           └───────────────┴────────────────────┘           │
│                           │                                │
│                           ▼                                │
│              ┌────────────────────────┐                    │
│              │   JobAdapter (ABC)     │                    │
│              └────────────────────────┘                    │
│                           │                                │
│           ┌───────────────┼───────────────┐                │
│           │               │               │                │
│           ▼               ▼               ▼                │
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐    │
│  │SlurmAdapter  │  │ PbsAdapter   │  │LocalAdapter  │    │
│  └──────────────┘  └──────────────┘  └──────────────┘    │
└─────────────────────────────────────────────────────────────┘
                      │
                      ▼
┌─────────────────────────────────────────────────────────────┐
│                 Job Scheduler (SLURM/PBS)                   │
└─────────────────────────────────────────────────────────────┘
```

### 2.2 Data Flow

#### Submit Flow
```
1. User: chemvcs submit <run-id> job.slurm
2. CLI validates run-id exists and is in 'planned' state
3. Python JobSubmitter:
   - Captures environment (modules, resources)
   - Saves job script snapshot
   - Calls adapter.submit()
4. SlurmAdapter executes sbatch, captures job ID
5. Update Run object: status='submitted', job_id, provenance
6. Store updated Run object to repository
7. Return job ID to user
```

#### Query Flow
```
1. User: chemvcs jobs
2. CLI lists all Run objects with job_id != null
3. For each run, query adapter.get_status(job_id)
4. Display: job_id, status, submitted_at, run_id
```

#### Retrieve Flow
```
1. User: chemvcs retrieve <job-id>
2. Find Run object by job_id
3. Check adapter.get_status() == COMPLETED
4. Copy output files to working directory
5. Update Run: status='finished', result data
6. Create commit: "Retrieved results for job <job-id>"
```

---

## 3. Extended Data Model

### 3.1 Run Object Extensions

```python
@dataclass
class Run:
    # Existing fields...
    structure_id: str
    code: str
    parameters: Dict[str, Any]
    status: str  # planned/submitted/running/finished/failed
    
    # New HPC fields
    job_id: Optional[str] = None
    job_system: Optional[str] = None  # "slurm", "pbs", "local"
    queue_name: Optional[str] = None
    submit_script: Optional[str] = None  # Full script content
    
    # Provenance
    modules_loaded: List[str] = field(default_factory=list)
    environment_vars: Dict[str, str] = field(default_factory=dict)
    job_resources: Optional[Dict[str, Any]] = None  # {nodes, cores, walltime}
    
    # Timestamps
    queued_at: Optional[str] = None
    started_at: Optional[str] = None
    
    # New methods
    def mark_queued(self, job_id: str, job_system: str, queue: str = None):
        """Mark as queued in job system"""
        
    def mark_retrieved(self, output_data: Dict[str, Any]):
        """Mark as retrieved with results"""
```

### 3.2 CoreObject Representation

The extended Run serializes to:

```json
{
  "version": "1.0",
  "type": "run",
  "meta": {
    "structure_id": "abc123...",
    "code": "vasp",
    "parameters": {...},
    "status": "submitted",
    "job_id": "12345",
    "job_system": "slurm",
    "queue_name": "batch",
    "modules_loaded": ["vasp/6.3.0", "intel/2021.4"],
    "job_resources": {
      "nodes": 2,
      "ntasks_per_node": 28,
      "walltime": "24:00:00"
    },
    "submit_script": "#!/bin/bash\n#SBATCH --nodes=2\n...",
    "submitted_at": "2025-12-15T10:30:00Z",
    "queued_at": "2025-12-15T10:30:05Z"
  },
  "refs": [...]
}
```

---

## 4. JobAdapter Interface

### 4.1 Abstract Base Class

```python
from abc import ABC, abstractmethod
from enum import Enum
from typing import Dict, Optional, List

class JobStatus(Enum):
    """Standard job status across all schedulers"""
    PENDING = "pending"      # Queued but not running
    RUNNING = "running"      # Currently executing
    COMPLETED = "completed"  # Finished successfully
    FAILED = "failed"        # Terminated with error
    CANCELLED = "cancelled"  # User/admin cancelled
    UNKNOWN = "unknown"      # Cannot determine status

class JobInfo:
    """Job information returned by adapters"""
    job_id: str
    status: JobStatus
    queue: Optional[str]
    nodes: Optional[int]
    start_time: Optional[str]
    end_time: Optional[str]
    exit_code: Optional[int]

class JobAdapter(ABC):
    """Abstract interface for job schedulers"""
    
    @property
    @abstractmethod
    def name(self) -> str:
        """Adapter name: 'slurm', 'pbs', 'local'"""
        pass
    
    @abstractmethod
    def submit(self, script_path: str, **kwargs) -> str:
        """
        Submit a job script.
        
        Args:
            script_path: Path to job script
            **kwargs: Scheduler-specific options
            
        Returns:
            Job ID (string)
            
        Raises:
            JobSubmissionError: If submission fails
        """
        pass
    
    @abstractmethod
    def get_status(self, job_id: str) -> JobStatus:
        """
        Query job status.
        
        Args:
            job_id: Job identifier
            
        Returns:
            JobStatus enum
            
        Raises:
            JobNotFoundError: If job doesn't exist
        """
        pass
    
    @abstractmethod
    def get_info(self, job_id: str) -> JobInfo:
        """
        Get detailed job information.
        
        Args:
            job_id: Job identifier
            
        Returns:
            JobInfo object
        """
        pass
    
    @abstractmethod
    def cancel(self, job_id: str) -> bool:
        """
        Cancel a running/pending job.
        
        Args:
            job_id: Job identifier
            
        Returns:
            True if cancelled successfully
        """
        pass
    
    def validate(self) -> bool:
        """
        Check if scheduler is available.
        
        Returns:
            True if commands are available
        """
        return True
```

### 4.2 SLURM Adapter

```python
class SlurmAdapter(JobAdapter):
    """SLURM workload manager adapter"""
    
    @property
    def name(self) -> str:
        return "slurm"
    
    def submit(self, script_path: str, **kwargs) -> str:
        """Submit via sbatch"""
        cmd = ['sbatch']
        
        # Add optional arguments
        if 'job_name' in kwargs:
            cmd.extend(['--job-name', kwargs['job_name']])
        if 'output' in kwargs:
            cmd.extend(['--output', kwargs['output']])
            
        cmd.append(script_path)
        
        result = subprocess.run(
            cmd, capture_output=True, text=True, check=True
        )
        
        # Parse: "Submitted batch job 12345"
        match = re.search(r'Submitted batch job (\d+)', result.stdout)
        if not match:
            raise JobSubmissionError(f"Failed to parse job ID: {result.stdout}")
        
        return match.group(1)
    
    def get_status(self, job_id: str) -> JobStatus:
        """Query via squeue"""
        result = subprocess.run(
            ['squeue', '-j', job_id, '-h', '-o', '%T'],
            capture_output=True, text=True
        )
        
        if result.returncode != 0:
            # Job not in queue, check sacct for completed jobs
            return self._check_completed(job_id)
        
        status_str = result.stdout.strip()
        return self._parse_slurm_status(status_str)
    
    def _parse_slurm_status(self, status: str) -> JobStatus:
        """Map SLURM status to JobStatus"""
        mapping = {
            'PENDING': JobStatus.PENDING,
            'RUNNING': JobStatus.RUNNING,
            'COMPLETED': JobStatus.COMPLETED,
            'FAILED': JobStatus.FAILED,
            'CANCELLED': JobStatus.CANCELLED,
            'TIMEOUT': JobStatus.FAILED,
        }
        return mapping.get(status, JobStatus.UNKNOWN)
```

---

## 5. CLI Commands

### 5.1 chemvcs submit

```
Usage: chemvcs submit [OPTIONS] <run-id> <script>

Submit a computational run to job scheduler.

Arguments:
  <run-id>   Hash or reference of the Run object
  <script>   Path to job submission script

Options:
  --adapter <name>    Job adapter to use (default: auto-detect)
  --queue <name>      Target queue name
  --job-name <name>   Job name (default: run-id[:8])
  --output <file>     Job output file

Examples:
  chemvcs submit abc123 vasp.slurm
  chemvcs submit @latest job.sh --queue=batch
```

**Implementation**:
- Validates run-id exists and status is 'planned'
- Calls Python: `JobSubmitter.submit_run(run, script_path)`
- Prints: "Submitted job <job-id> for run <run-id>"

### 5.2 chemvcs jobs

```
Usage: chemvcs jobs [OPTIONS]

List tracked computational jobs.

Options:
  --status <status>   Filter by status (pending/running/completed/failed)
  --limit <n>         Show only last n jobs (default: 20)
  --all              Show all jobs including completed

Output Format:
  JOB_ID    RUN_ID      STATUS     SUBMITTED           CODE
  12345     abc123...   RUNNING    2025-12-15 10:30    vasp
  12344     def456...   COMPLETED  2025-12-15 09:15    orca
```

**Implementation**:
- Queries all Run objects with `job_id != null`
- For each, calls `adapter.get_status(job_id)`
- Displays in table format

### 5.3 chemvcs retrieve

```
Usage: chemvcs retrieve [OPTIONS] <job-id>

Retrieve completed job results and create commit.

Arguments:
  <job-id>   Job identifier

Options:
  --output-dir <dir>  Copy outputs to directory (default: working dir)
  --no-commit        Don't create automatic commit

Examples:
  chemvcs retrieve 12345
  chemvcs retrieve 12345 --output-dir=results/
```

**Implementation**:
- Finds Run object by job_id
- Checks status is COMPLETED
- Copies output files to working directory
- Updates Run object with results
- Creates commit: "Retrieved results for job 12345"

---

## 6. Environment Provenance

### 6.1 Module Capture

```python
class EnvironmentCapture:
    """Capture computational environment for provenance"""
    
    @staticmethod
    def capture_modules() -> List[str]:
        """
        Capture loaded environment modules.
        
        Returns:
            List of module names with versions
        """
        try:
            result = subprocess.run(
                ['module', 'list', '-t'],
                capture_output=True, text=True,
                stderr=subprocess.STDOUT  # module outputs to stderr
            )
            
            # Parse output: each line is a module
            modules = []
            for line in result.stdout.strip().split('\n'):
                line = line.strip()
                if line and not line.startswith('Currently'):
                    modules.append(line)
            
            return modules
        except FileNotFoundError:
            # module command not available
            return []
```

### 6.2 Resource Parsing

```python
@staticmethod
def parse_slurm_script(script_path: str) -> Dict[str, Any]:
    """
    Extract resource requirements from SLURM script.
    
    Returns:
        Dictionary with nodes, cores, walltime, memory, etc.
    """
    resources = {}
    
    with open(script_path) as f:
        for line in f:
            line = line.strip()
            if not line.startswith('#SBATCH'):
                continue
            
            # Parse #SBATCH --nodes=4
            if '--nodes=' in line:
                resources['nodes'] = int(line.split('=')[1])
            elif '--ntasks-per-node=' in line:
                resources['ntasks_per_node'] = int(line.split('=')[1])
            elif '--time=' in line:
                resources['walltime'] = line.split('=')[1]
            elif '--mem=' in line:
                resources['memory'] = line.split('=')[1]
    
    return resources
```

### 6.3 Script Snapshot

```python
@staticmethod
def save_script_snapshot(script_path: str) -> str:
    """
    Read and return job script content.
    
    This ensures exact script used for submission is preserved.
    """
    with open(script_path) as f:
        return f.read()
```

---

## 7. Error Handling

### 7.1 Exception Hierarchy

```python
class HpcError(Exception):
    """Base exception for HPC operations"""
    pass

class JobSubmissionError(HpcError):
    """Failed to submit job"""
    pass

class JobNotFoundError(HpcError):
    """Job ID not found in scheduler"""
    pass

class AdapterNotAvailableError(HpcError):
    """Job adapter not available (missing commands)"""
    pass

class InvalidJobStateError(HpcError):
    """Operation invalid for current job state"""
    pass
```

### 7.2 Common Error Scenarios

| Error | Cause | Handling |
|-------|-------|----------|
| Run not in 'planned' state | Trying to submit already-submitted run | Raise error, show current state |
| Job script not found | Invalid path | Clear error message |
| sbatch command fails | No permissions, syntax error | Show sbatch stderr output |
| Job ID not found | Job completed and purged from queue | Check sacct for historical data |
| Network timeout | Cluster unreachable | Retry with exponential backoff |

---

## 8. Testing Strategy

### 8.1 Unit Tests (Mock-based)

**No SLURM Required**: Use `unittest.mock` to simulate subprocess calls.

```python
# test_slurm_adapter.py
class TestSlurmAdapter(unittest.TestCase):
    @mock.patch('subprocess.run')
    def test_submit_success(self, mock_run):
        """Test successful job submission"""
        # Mock sbatch output
        mock_run.return_value = mock.Mock(
            returncode=0,
            stdout="Submitted batch job 12345\n",
            stderr=""
        )
        
        adapter = SlurmAdapter()
        job_id = adapter.submit('test.slurm')
        
        assert job_id == "12345"
        mock_run.assert_called_once()
    
    @mock.patch('subprocess.run')
    def test_get_status_running(self, mock_run):
        """Test status query for running job"""
        mock_run.return_value = mock.Mock(
            returncode=0,
            stdout="RUNNING\n"
        )
        
        adapter = SlurmAdapter()
        status = adapter.get_status('12345')
        
        assert status == JobStatus.RUNNING
```

### 8.2 Integration Tests (Optional)

If SLURM is available:
```python
@pytest.mark.skipif(not has_slurm(), reason="SLURM not available")
def test_real_submission():
    """Test with actual SLURM cluster"""
    # Submit a simple test job
    # Wait for completion
    # Verify results
```

### 8.3 End-to-End Tests

```python
def test_submit_retrieve_workflow(tmp_repo):
    """Test complete workflow: submit → query → retrieve"""
    # 1. Create Run object
    run = Run(structure_id="abc", code="test", status="planned")
    
    # 2. Submit (using mock adapter)
    submitter = JobSubmitter(repo, MockAdapter())
    job_id = submitter.submit_run(run, 'test.sh')
    
    # 3. Verify Run updated
    assert run.job_id == job_id
    assert run.status == "submitted"
    
    # 4. Simulate job completion
    # ...
    
    # 5. Retrieve results
    submitter.retrieve_results(job_id, 'outputs/')
    
    # 6. Verify commit created
    assert repo.has_changes() == False
```

---

## 9. Security Considerations

### 9.1 Credential Management

- **SSH Keys**: Use standard SSH agent for cluster access
- **API Tokens**: If scheduler has REST API, store tokens securely
- **Job Scripts**: Validate script paths to prevent injection

### 9.2 Permission Checks

```python
def submit(self, script_path: str, **kwargs) -> str:
    # Verify script is readable
    if not os.access(script_path, os.R_OK):
        raise PermissionError(f"Cannot read script: {script_path}")
    
    # Verify script is within project directory (prevent arbitrary execution)
    script_abs = os.path.abspath(script_path)
    repo_root = self.repo.get_root()
    if not script_abs.startswith(repo_root):
        raise SecurityError("Job script must be within repository")
    
    # Submit job
    # ...
```

### 9.3 Output Sanitization

- Validate job IDs match expected format (prevent command injection)
- Sanitize paths when copying output files
- Limit output file sizes when storing in repository

---

## 10. Future Extensions

### 10.1 Additional Adapters

- **PBS/Torque**: `qsub`, `qstat`, `qdel`
- **SGE**: `qsub`, `qstat`, `qdel` (different format)
- **LSF**: `bsub`, `bjobs`, `bkill`
- **Local**: Fork process directly (for testing without scheduler)

### 10.2 Advanced Features

- **Job Arrays**: Submit multiple runs as job array
- **Job Dependencies**: Run B after Run A completes
- **Resource Optimization**: Suggest optimal node/core counts
- **Cost Tracking**: Estimate compute costs
- **Auto-resubmit**: Restart failed jobs automatically

### 10.3 Web Dashboard

```
┌─────────────────────────────────────┐
│        ChemVCS Web Dashboard        │
├─────────────────────────────────────┤
│  Active Jobs              [Filter]  │
│                                     │
│  Job 12345  ███████░░░  70%         │
│  RUNNING    2h 15m / 4h             │
│                                     │
│  Job 12344  ██████████  100%        │
│  COMPLETED  3h 45m                  │
└─────────────────────────────────────┘
```

---

## 11. Implementation Checklist

### Phase 1: Core Infrastructure (Week 1-2)
- [ ] Extend Run object with HPC fields
- [ ] Add Run HPC tests (8-10 new tests)
- [ ] Design JobAdapter interface
- [ ] Implement SlurmAdapter
- [ ] Add adapter unit tests (mock-based, 15-20 tests)

### Phase 2: CLI Integration (Week 2-3)
- [ ] Go: Add `submit` command skeleton
- [ ] Go: Add `jobs` command
- [ ] Go: Add `retrieve` command
- [ ] Python: Implement JobSubmitter class
- [ ] Python: Implement JobTracker class
- [ ] Integration tests for CLI (10-15 tests)

### Phase 3: Provenance (Week 3-4)
- [ ] Implement EnvironmentCapture
- [ ] Add module capture tests
- [ ] Add script parsing tests
- [ ] Integrate into submit workflow

### Phase 4: Testing & Documentation (Week 4-5)
- [ ] End-to-end workflow tests
- [ ] Error handling tests
- [ ] Write user guide (docs/10-hpc-user-guide.md)
- [ ] Create example repository
- [ ] Update README and FEATURE_STATUS

### Phase 5: Polish (Week 5-6)
- [ ] Code review and refactoring
- [ ] Performance optimization
- [ ] Documentation improvements
- [ ] Release v0.6.0

---

## 12. Success Metrics

M6 is complete when:

1. ✅ `chemvcs submit` successfully submits jobs (mock-tested)
2. ✅ `chemvcs jobs` displays job status with real-time queries
3. ✅ `chemvcs retrieve` fetches results and creates commits
4. ✅ Run objects contain complete provenance metadata
5. ✅ At least 50 new tests passing (40 Python + 10 Go)
6. ✅ Complete user documentation with examples
7. ✅ Can be extended to other schedulers (PBS, SGE) easily

---

## Appendix A: Example Usage

```bash
# Initialize repository
chemvcs init water-dft
cd water-dft

# Create structure
python -c "
from chemvcs_py import Repo, Structure
from chemvcs_py.io import read_xyz

repo = Repo('.')
structure = read_xyz('water.xyz')
obj = structure.to_core_object()
struct_hash = repo.store_object(obj)
print(f'Structure: {struct_hash}')
"

# Create Run object
python -c "
from chemvcs_py import Repo, Run

repo = Repo('.')
run = Run(
    structure_id='<struct-hash>',
    code='vasp',
    parameters={'encut': 400, 'ismear': 0},
    status='planned'
)
obj = run.to_core_object()
run_hash = repo.store_object(obj)
repo.commit('Add VASP calculation plan')
print(f'Run: {run_hash}')
"

# Submit job
chemvcs submit <run-hash> vasp.slurm
# Output: Submitted job 12345 for run abc123...

# Check status
chemvcs jobs
# JOB_ID  RUN_ID    STATUS   SUBMITTED          CODE
# 12345   abc123... RUNNING  2025-12-15 10:30   vasp

# Wait for completion, then retrieve
chemvcs retrieve 12345
# Output: Retrieved results for job 12345
#         Created commit: abc456...

# View results
chemvcs log
```

---

**End of Design Document**
