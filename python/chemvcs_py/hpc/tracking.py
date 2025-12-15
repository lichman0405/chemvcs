"""Job tracking and status queries."""

from typing import List, Optional, Dict
from dataclasses import dataclass

from ..core.repo import Repo
from ..domain.run import Run
from .adapter import JobAdapter, JobStatus
from .exceptions import JobNotFoundError, JobCancellationError


@dataclass
class TrackedJob:
    """Information about a tracked job."""
    job_id: str
    run_id: str
    run: Run
    status: JobStatus
    code: str
    structure_id: str


class JobTracker:
    """Track and query status of submitted jobs."""
    
    def __init__(self, repo: Repo, adapter: JobAdapter):
        """
        Initialize JobTracker.
        
        Args:
            repo: ChemVCS repository
            adapter: Job scheduler adapter
        """
        self.repo = repo
        self.adapter = adapter
    
    def list_jobs(
        self,
        status_filter: Optional[List[str]] = None,
        limit: Optional[int] = None
    ) -> List[TrackedJob]:
        """
        List all tracked jobs.
        
        Args:
            status_filter: Filter by status (e.g., ['running', 'pending'])
            limit: Maximum number of jobs to return
            
        Returns:
            List of TrackedJob objects
        """
        tracked_jobs = []
        
        # Get all Run objects with job_id
        runs = self._get_submitted_runs()
        
        for run in runs:
            if not run.job_id:
                continue
            
            # Query current status from scheduler
            try:
                status = self.adapter.get_status(run.job_id)
            except JobNotFoundError:
                # Job no longer in scheduler, use Run's stored status
                status = self._run_status_to_job_status(run.status)
            
            # Apply status filter
            if status_filter and status.value not in status_filter:
                continue
            
            tracked_job = TrackedJob(
                job_id=run.job_id,
                run_id=run.id or "unknown",
                run=run,
                status=status,
                code=run.code,
                structure_id=run.structure_id
            )
            
            tracked_jobs.append(tracked_job)
            
            # Apply limit
            if limit and len(tracked_jobs) >= limit:
                break
        
        return tracked_jobs
    
    def get_job_by_id(self, job_id: str) -> Optional[TrackedJob]:
        """
        Get tracked job by job ID.
        
        Args:
            job_id: Job identifier
            
        Returns:
            TrackedJob or None if not found
        """
        runs = self._get_submitted_runs()
        
        for run in runs:
            if run.job_id == job_id:
                try:
                    status = self.adapter.get_status(job_id)
                except JobNotFoundError:
                    status = self._run_status_to_job_status(run.status)
                
                return TrackedJob(
                    job_id=job_id,
                    run_id=run.id or "unknown",
                    run=run,
                    status=status,
                    code=run.code,
                    structure_id=run.structure_id
                )
        
        return None
    
    def find_run_by_job_id(self, job_id: str) -> Optional[Run]:
        """
        Find Run object by job ID.
        
        Args:
            job_id: Job identifier
            
        Returns:
            Run object or None if not found
        """
        runs = self._get_submitted_runs()
        
        for run in runs:
            if run.job_id == job_id:
                return run
        
        return None
    
    def _get_submitted_runs(self) -> List[Run]:
        """Get all Run objects that have been submitted."""
        runs = []
        
        # List all objects of type 'run'
        objects = self.repo.list_objects(type_filter='run')
        
        for obj_info in objects:
            try:
                obj = self.repo.get_object(obj_info['hash'])
                run = Run.from_core_object(obj)
                
                # Only include runs with job_id
                if run.job_id:
                    runs.append(run)
            except Exception:
                # Skip objects that can't be parsed
                continue
        
        return runs
    
    def cancel_job(self, job_id: str) -> bool:
        """
        Cancel a running or pending job.
        
        Args:
            job_id: Job identifier
            
        Returns:
            True if cancellation successful
            
        Raises:
            JobNotFoundError: If job or run not found
            JobCancellationError: If cancellation fails
        """
        # Find the associated Run
        run = self.find_run_by_job_id(job_id)
        if not run:
            raise JobNotFoundError(f"No run found for job {job_id}")
        
        # Cancel via adapter
        success = self.adapter.cancel(job_id)
        if not success:
            raise JobCancellationError(f"Failed to cancel job {job_id}")
        
        # Update Run status
        run.status = "cancelled"
        run.metadata['cancelled_at'] = self._get_timestamp()
        
        # Store updated Run
        obj = run.to_core_object()
        self.repo.store_object(obj)
        
        return True
    
    def cancel_run(self, run_hash: str) -> bool:
        """
        Cancel a job by Run hash.
        
        Args:
            run_hash: Hash of the Run object
            
        Returns:
            True if cancellation successful
            
        Raises:
            JobNotFoundError: If run not found or no job_id
            JobCancellationError: If cancellation fails
        """
        # Get Run object
        try:
            obj = self.repo.get_object(run_hash)
            run = Run.from_core_object(obj)
        except Exception as e:
            raise JobNotFoundError(f"Run {run_hash} not found: {e}")
        
        if not run.job_id:
            raise JobNotFoundError(f"Run {run_hash} has no associated job")
        
        return self.cancel_job(run.job_id)
    
    @staticmethod
    def _get_timestamp() -> str:
        """Get current timestamp in ISO format."""
        from datetime import datetime
        return datetime.utcnow().isoformat() + 'Z'
    
    @staticmethod
    def _run_status_to_job_status(run_status: str) -> JobStatus:
        """Convert Run status to JobStatus."""
        mapping = {
            'submitted': JobStatus.PENDING,
            'running': JobStatus.RUNNING,
            'finished': JobStatus.COMPLETED,
            'failed': JobStatus.FAILED,
            'cancelled': JobStatus.CANCELLED,
        }
        return mapping.get(run_status, JobStatus.UNKNOWN)
