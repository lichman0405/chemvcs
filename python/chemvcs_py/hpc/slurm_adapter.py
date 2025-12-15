"""SLURM workload manager adapter."""

import re
import subprocess
import shutil
from typing import Dict, Optional

from .adapter import JobAdapter, JobStatus, JobInfo
from .exceptions import JobSubmissionError, JobNotFoundError, AdapterNotAvailableError


class SlurmAdapter(JobAdapter):
    """
    Adapter for SLURM workload manager.
    
    Uses SLURM commands:
    - sbatch: Submit jobs
    - squeue: Query running/pending jobs
    - sacct: Query completed jobs
    - scancel: Cancel jobs
    """
    
    @property
    def name(self) -> str:
        return "slurm"
    
    def validate(self) -> bool:
        """Check if SLURM commands are available."""
        return shutil.which('sbatch') is not None
    
    def submit(self, script_path: str, **kwargs) -> str:
        """
        Submit job via sbatch.
        
        Args:
            script_path: Path to job script
            **kwargs: Optional arguments:
                - job_name: Job name
                - output: Output file path
                - partition: Partition/queue name
                
        Returns:
            Job ID
            
        Raises:
            JobSubmissionError: If submission fails
            AdapterNotAvailableError: If sbatch not found
        """
        if not self.validate():
            raise AdapterNotAvailableError("sbatch command not found")
        
        cmd = ['sbatch']
        
        # Add optional arguments
        if 'job_name' in kwargs:
            cmd.extend(['--job-name', kwargs['job_name']])
        if 'output' in kwargs:
            cmd.extend(['--output', kwargs['output']])
        if 'partition' in kwargs:
            cmd.extend(['--partition', kwargs['partition']])
        
        cmd.append(script_path)
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
        except subprocess.CalledProcessError as e:
            raise JobSubmissionError(
                f"sbatch failed: {e.stderr}"
            )
        except FileNotFoundError:
            raise AdapterNotAvailableError("sbatch command not found")
        
        # Parse: "Submitted batch job 12345"
        match = re.search(r'Submitted batch job (\d+)', result.stdout)
        if not match:
            raise JobSubmissionError(
                f"Failed to parse job ID from: {result.stdout}"
            )
        
        return match.group(1)
    
    def get_status(self, job_id: str) -> JobStatus:
        """
        Query job status via squeue/sacct.
        
        Args:
            job_id: Job identifier
            
        Returns:
            JobStatus enum
            
        Raises:
            JobNotFoundError: If job doesn't exist
        """
        # First try squeue for running/pending jobs
        try:
            result = subprocess.run(
                ['squeue', '-j', job_id, '-h', '-o', '%T'],
                capture_output=True,
                text=True,
                timeout=10
            )
            
            if result.returncode == 0 and result.stdout.strip():
                status_str = result.stdout.strip()
                return self._parse_slurm_status(status_str)
        except (subprocess.TimeoutExpired, FileNotFoundError):
            pass
        
        # If not in queue, check sacct for completed jobs
        return self._check_completed(job_id)
    
    def _check_completed(self, job_id: str) -> JobStatus:
        """Check sacct for completed job status."""
        try:
            result = subprocess.run(
                ['sacct', '-j', job_id, '-n', '-o', 'State'],
                capture_output=True,
                text=True,
                timeout=10
            )
            
            if result.returncode == 0 and result.stdout.strip():
                # Take first status line
                status_str = result.stdout.strip().split('\n')[0].strip()
                return self._parse_slurm_status(status_str)
            else:
                raise JobNotFoundError(f"Job {job_id} not found")
        except (subprocess.TimeoutExpired, FileNotFoundError):
            raise JobNotFoundError(f"Cannot query job {job_id}")
    
    def _parse_slurm_status(self, status: str) -> JobStatus:
        """Map SLURM status codes to JobStatus."""
        # SLURM status codes: https://slurm.schedmd.com/squeue.html
        mapping = {
            'PENDING': JobStatus.PENDING,
            'RUNNING': JobStatus.RUNNING,
            'COMPLETED': JobStatus.COMPLETED,
            'FAILED': JobStatus.FAILED,
            'CANCELLED': JobStatus.CANCELLED,
            'TIMEOUT': JobStatus.FAILED,
            'NODE_FAIL': JobStatus.FAILED,
            'PREEMPTED': JobStatus.FAILED,
            'OUT_OF_MEMORY': JobStatus.FAILED,
        }
        
        # Handle status with modifiers (e.g., "RUNNING+")
        clean_status = status.split('+')[0].split()[0]
        
        return mapping.get(clean_status, JobStatus.UNKNOWN)
    
    def get_info(self, job_id: str) -> JobInfo:
        """
        Get detailed job information via sacct.
        
        Args:
            job_id: Job identifier
            
        Returns:
            JobInfo object
        """
        try:
            result = subprocess.run(
                ['sacct', '-j', job_id, '-n', '-o',
                 'State,Partition,NNodes,Start,End,ExitCode'],
                capture_output=True,
                text=True,
                timeout=10
            )
            
            if result.returncode != 0 or not result.stdout.strip():
                raise JobNotFoundError(f"Job {job_id} not found")
            
            # Parse first line (job step 0)
            fields = result.stdout.strip().split('\n')[0].split()
            
            status_str = fields[0] if len(fields) > 0 else "UNKNOWN"
            queue = fields[1] if len(fields) > 1 else None
            nodes_str = fields[2] if len(fields) > 2 else None
            start_time = fields[3] if len(fields) > 3 else None
            end_time = fields[4] if len(fields) > 4 else None
            exit_code_str = fields[5] if len(fields) > 5 else None
            
            # Parse exit code (format: "0:0")
            exit_code = None
            if exit_code_str:
                try:
                    exit_code = int(exit_code_str.split(':')[0])
                except (ValueError, IndexError):
                    pass
            
            # Parse node count
            nodes = None
            if nodes_str:
                try:
                    nodes = int(nodes_str)
                except ValueError:
                    pass
            
            return JobInfo(
                job_id=job_id,
                status=self._parse_slurm_status(status_str),
                queue=queue,
                nodes=nodes,
                start_time=start_time,
                end_time=end_time,
                exit_code=exit_code
            )
        except (subprocess.TimeoutExpired, FileNotFoundError):
            raise JobNotFoundError(f"Cannot query job {job_id}")
    
    def cancel(self, job_id: str) -> bool:
        """
        Cancel job via scancel.
        
        Args:
            job_id: Job identifier
            
        Returns:
            True if cancelled successfully
        """
        try:
            result = subprocess.run(
                ['scancel', job_id],
                capture_output=True,
                text=True,
                timeout=10
            )
            
            return result.returncode == 0
        except (subprocess.TimeoutExpired, FileNotFoundError):
            return False
