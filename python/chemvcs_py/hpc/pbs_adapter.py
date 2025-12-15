"""PBS/Torque workload manager adapter."""

import re
import subprocess
import shutil
from typing import Dict, Optional

from .adapter import JobAdapter, JobStatus, JobInfo
from .exceptions import JobSubmissionError, JobNotFoundError, AdapterNotAvailableError


class PbsAdapter(JobAdapter):
    """
    Adapter for PBS/Torque workload manager.
    
    Uses PBS commands:
    - qsub: Submit jobs
    - qstat: Query job status
    - qdel: Cancel jobs
    """
    
    @property
    def name(self) -> str:
        return "pbs"
    
    def validate(self) -> bool:
        """Check if PBS commands are available."""
        return shutil.which('qsub') is not None
    
    def submit(self, script_path: str, **kwargs) -> str:
        """
        Submit job via qsub.
        
        Args:
            script_path: Path to job script
            **kwargs: Optional arguments:
                - job_name: Job name (-N)
                - output: Output file path (-o)
                - queue: Queue name (-q)
                
        Returns:
            Job ID
            
        Raises:
            JobSubmissionError: If submission fails
            AdapterNotAvailableError: If qsub not found
        """
        if not self.validate():
            raise AdapterNotAvailableError("qsub command not found")
        
        cmd = ['qsub']
        
        # Add optional arguments
        if 'job_name' in kwargs:
            cmd.extend(['-N', kwargs['job_name']])
        if 'output' in kwargs:
            cmd.extend(['-o', kwargs['output']])
        if 'queue' in kwargs:
            cmd.extend(['-q', kwargs['queue']])
        
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
                f"qsub failed: {e.stderr}"
            )
        except FileNotFoundError:
            raise AdapterNotAvailableError("qsub command not found")
        
        # Parse job ID from output
        # Format can be: "12345.hostname" or just "12345"
        job_id = result.stdout.strip()
        if not job_id:
            raise JobSubmissionError(
                f"Failed to parse job ID from: {result.stdout}"
            )
        
        return job_id
    
    def get_status(self, job_id: str) -> JobStatus:
        """
        Query job status via qstat.
        
        Args:
            job_id: Job identifier
            
        Returns:
            JobStatus enum
            
        Raises:
            JobNotFoundError: If job doesn't exist
        """
        try:
            result = subprocess.run(
                ['qstat', '-f', job_id],
                capture_output=True,
                text=True,
                timeout=10
            )
            
            if result.returncode != 0:
                # Check if job finished (not in queue)
                if 'Unknown Job Id' in result.stderr or 'job has finished' in result.stderr:
                    # Job completed, try to get final status
                    return self._check_completed(job_id)
                raise JobNotFoundError(f"Job {job_id} not found")
            
            # Parse job_state from qstat output
            # Format: "job_state = R" or similar
            match = re.search(r'job_state\s*=\s*(\w+)', result.stdout)
            if not match:
                return JobStatus.UNKNOWN
            
            status_code = match.group(1)
            return self._parse_pbs_status(status_code)
            
        except (subprocess.TimeoutExpired, FileNotFoundError):
            raise JobNotFoundError(f"Cannot query job {job_id}")
    
    def _check_completed(self, job_id: str) -> JobStatus:
        """
        Try to determine completed job status.
        
        PBS doesn't always keep completed job info, so we assume
        COMPLETED if the job is not in the queue.
        """
        # In a real implementation, you might check job logs or use
        # a job accounting system if available
        return JobStatus.COMPLETED
    
    def _parse_pbs_status(self, status_code: str) -> JobStatus:
        """Map PBS status codes to JobStatus."""
        # PBS status codes:
        # Q = Queued
        # R = Running
        # E = Exiting (completing)
        # C = Completed
        # H = Held
        # W = Waiting
        # S = Suspended
        mapping = {
            'Q': JobStatus.PENDING,
            'W': JobStatus.PENDING,
            'H': JobStatus.PENDING,
            'R': JobStatus.RUNNING,
            'E': JobStatus.RUNNING,  # Still running, cleaning up
            'C': JobStatus.COMPLETED,
            'S': JobStatus.CANCELLED,  # Suspended treated as cancelled
        }
        
        return mapping.get(status_code, JobStatus.UNKNOWN)
    
    def get_info(self, job_id: str) -> JobInfo:
        """
        Get detailed job information via qstat.
        
        Args:
            job_id: Job identifier
            
        Returns:
            JobInfo object
        """
        try:
            result = subprocess.run(
                ['qstat', '-f', job_id],
                capture_output=True,
                text=True,
                timeout=10
            )
            
            if result.returncode != 0:
                raise JobNotFoundError(f"Job {job_id} not found")
            
            output = result.stdout
            
            # Parse various fields from qstat -f output
            status_match = re.search(r'job_state\s*=\s*(\w+)', output)
            queue_match = re.search(r'queue\s*=\s*(\S+)', output)
            nodes_match = re.search(r'Resource_List\.nodes\s*=\s*(\d+)', output)
            start_match = re.search(r'start_time\s*=\s*(.+)', output)
            exit_match = re.search(r'exit_status\s*=\s*(\d+)', output)
            
            status_code = status_match.group(1) if status_match else "U"
            queue = queue_match.group(1) if queue_match else None
            nodes = int(nodes_match.group(1)) if nodes_match else None
            start_time = start_match.group(1).strip() if start_match else None
            exit_code = int(exit_match.group(1)) if exit_match else None
            
            return JobInfo(
                job_id=job_id,
                status=self._parse_pbs_status(status_code),
                queue=queue,
                nodes=nodes,
                start_time=start_time,
                end_time=None,  # PBS qstat doesn't always show end_time
                exit_code=exit_code
            )
        except (subprocess.TimeoutExpired, FileNotFoundError):
            raise JobNotFoundError(f"Cannot query job {job_id}")
    
    def cancel(self, job_id: str) -> bool:
        """
        Cancel job via qdel.
        
        Args:
            job_id: Job identifier
            
        Returns:
            True if cancelled successfully
        """
        try:
            result = subprocess.run(
                ['qdel', job_id],
                capture_output=True,
                text=True,
                timeout=10
            )
            
            return result.returncode == 0
        except (subprocess.TimeoutExpired, FileNotFoundError):
            return False
