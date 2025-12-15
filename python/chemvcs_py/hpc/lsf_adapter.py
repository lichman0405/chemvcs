"""LSF (IBM Spectrum LSF) workload manager adapter."""

import re
import subprocess
import shutil
from typing import Dict, Optional

from .adapter import JobAdapter, JobStatus, JobInfo
from .exceptions import JobSubmissionError, JobNotFoundError, AdapterNotAvailableError


class LsfAdapter(JobAdapter):
    """
    Adapter for LSF (IBM Spectrum LSF) workload manager.
    
    Uses LSF commands:
    - bsub: Submit jobs
    - bjobs: Query job status
    - bkill: Cancel jobs
    """
    
    @property
    def name(self) -> str:
        return "lsf"
    
    def validate(self) -> bool:
        """Check if LSF commands are available."""
        return shutil.which('bsub') is not None
    
    def submit(self, script_path: str, **kwargs) -> str:
        """
        Submit job via bsub.
        
        Args:
            script_path: Path to job script
            **kwargs: Optional arguments:
                - job_name: Job name (-J)
                - output: Output file path (-o)
                - queue: Queue name (-q)
                
        Returns:
            Job ID
            
        Raises:
            JobSubmissionError: If submission fails
            AdapterNotAvailableError: If bsub not found
        """
        if not self.validate():
            raise AdapterNotAvailableError("bsub command not found")
        
        cmd = ['bsub']
        
        # Add optional arguments
        if 'job_name' in kwargs:
            cmd.extend(['-J', kwargs['job_name']])
        if 'output' in kwargs:
            cmd.extend(['-o', kwargs['output']])
        if 'queue' in kwargs:
            cmd.extend(['-q', kwargs['queue']])
        
        # LSF requires input from stdin or file
        cmd.append('<')
        cmd.append(script_path)
        
        try:
            # Use shell=True for input redirection
            result = subprocess.run(
                ' '.join(cmd),
                capture_output=True,
                text=True,
                check=True,
                shell=True
            )
        except subprocess.CalledProcessError as e:
            raise JobSubmissionError(
                f"bsub failed: {e.stderr}"
            )
        except FileNotFoundError:
            raise AdapterNotAvailableError("bsub command not found")
        
        # Parse: "Job <12345> is submitted to queue <normal>."
        match = re.search(r'Job <(\d+)> is submitted', result.stdout)
        if not match:
            raise JobSubmissionError(
                f"Failed to parse job ID from: {result.stdout}"
            )
        
        return match.group(1)
    
    def get_status(self, job_id: str) -> JobStatus:
        """
        Query job status via bjobs.
        
        Args:
            job_id: Job identifier
            
        Returns:
            JobStatus enum
            
        Raises:
            JobNotFoundError: If job doesn't exist
        """
        try:
            result = subprocess.run(
                ['bjobs', '-noheader', '-o', 'stat', job_id],
                capture_output=True,
                text=True,
                timeout=10
            )
            
            if result.returncode != 0:
                # Job not in queue, check if it finished
                if 'not found' in result.stderr.lower():
                    # Assume completed if not in queue
                    return JobStatus.COMPLETED
                raise JobNotFoundError(f"Job {job_id} not found")
            
            status_str = result.stdout.strip()
            if not status_str:
                raise JobNotFoundError(f"Job {job_id} not found")
            
            return self._parse_lsf_status(status_str)
            
        except (subprocess.TimeoutExpired, FileNotFoundError):
            raise JobNotFoundError(f"Cannot query job {job_id}")
    
    def _parse_lsf_status(self, status: str) -> JobStatus:
        """Map LSF status codes to JobStatus."""
        # LSF status codes:
        # PEND = Pending
        # RUN = Running
        # DONE = Completed successfully
        # EXIT = Completed with errors
        # PSUSP = Suspended by user
        # USUSP = Suspended by system
        # SSUSP = Suspended by system
        # WAIT = Waiting for dependency
        # ZOMBI = Zombie job (completed but not reaped)
        mapping = {
            'PEND': JobStatus.PENDING,
            'WAIT': JobStatus.PENDING,
            'RUN': JobStatus.RUNNING,
            'DONE': JobStatus.COMPLETED,
            'EXIT': JobStatus.FAILED,
            'PSUSP': JobStatus.CANCELLED,
            'USUSP': JobStatus.CANCELLED,
            'SSUSP': JobStatus.CANCELLED,
            'ZOMBI': JobStatus.COMPLETED,
        }
        
        return mapping.get(status, JobStatus.UNKNOWN)
    
    def get_info(self, job_id: str) -> JobInfo:
        """
        Get detailed job information via bjobs.
        
        Args:
            job_id: Job identifier
            
        Returns:
            JobInfo object
        """
        try:
            result = subprocess.run(
                ['bjobs', '-l', job_id],
                capture_output=True,
                text=True,
                timeout=10
            )
            
            if result.returncode != 0:
                raise JobNotFoundError(f"Job {job_id} not found")
            
            output = result.stdout
            
            # Parse various fields from bjobs -l output
            status_match = re.search(r'Status\s+<(\w+)>', output)
            queue_match = re.search(r'Queue\s+<(\S+)>', output)
            # LSF shows requested slots/cores, not necessarily nodes
            slots_match = re.search(r'Requested\s+(\d+)\s+Task', output)
            # Match "Mon Dec 15 10:31:00: Started on <exec-host>"
            start_match = re.search(r'(\w{3}\s+\w{3}\s+\d+\s+\d+:\d+:\d+):\s+Started\s+on', output)
            
            status_str = status_match.group(1) if status_match else "UNKN"
            queue = queue_match.group(1) if queue_match else None
            slots = int(slots_match.group(1)) if slots_match else None
            start_time = start_match.group(1) if start_match else None
            
            return JobInfo(
                job_id=job_id,
                status=self._parse_lsf_status(status_str),
                queue=queue,
                nodes=None,  # LSF uses slots/cores, not nodes
                cores=slots,  # Store slots in cores field
                start_time=start_time,
                end_time=None,  # Would need to parse from job completion info
                exit_code=None  # Would need separate query
            )
        except (subprocess.TimeoutExpired, FileNotFoundError):
            raise JobNotFoundError(f"Cannot query job {job_id}")
    
    def cancel(self, job_id: str) -> bool:
        """
        Cancel job via bkill.
        
        Args:
            job_id: Job identifier
            
        Returns:
            True if cancelled successfully
        """
        try:
            result = subprocess.run(
                ['bkill', job_id],
                capture_output=True,
                text=True,
                timeout=10
            )
            
            return result.returncode == 0
        except (subprocess.TimeoutExpired, FileNotFoundError):
            return False
