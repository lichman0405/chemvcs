"""Abstract interface for job schedulers."""

from abc import ABC, abstractmethod
from dataclasses import dataclass
from enum import Enum
from typing import Dict, Optional


class JobStatus(Enum):
    """Standard job status across all schedulers."""
    PENDING = "pending"      # Queued but not running
    RUNNING = "running"      # Currently executing
    COMPLETED = "completed"  # Finished successfully
    FAILED = "failed"        # Terminated with error
    CANCELLED = "cancelled"  # User/admin cancelled
    UNKNOWN = "unknown"      # Cannot determine status


@dataclass
class JobInfo:
    """Job information returned by adapters."""
    job_id: str
    status: JobStatus
    queue: Optional[str] = None
    nodes: Optional[int] = None
    start_time: Optional[str] = None
    end_time: Optional[str] = None
    exit_code: Optional[int] = None


class JobAdapter(ABC):
    """Abstract interface for job schedulers."""
    
    @property
    @abstractmethod
    def name(self) -> str:
        """Adapter name: 'slurm', 'pbs', 'local'."""
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
            
        Raises:
            JobNotFoundError: If job doesn't exist
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
            
        Raises:
            JobNotFoundError: If job doesn't exist
        """
        pass
    
    def validate(self) -> bool:
        """
        Check if scheduler is available.
        
        Returns:
            True if commands are available
        """
        return True
