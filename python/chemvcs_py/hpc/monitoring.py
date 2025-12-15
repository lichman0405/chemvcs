"""Interactive job monitoring functionality."""

import time
from typing import Optional
from datetime import datetime

from ..core.repo import Repo
from .adapter import JobAdapter, JobStatus
from .tracking import JobTracker, TrackedJob
from .exceptions import JobNotFoundError


class JobWatcher:
    """
    Interactive job monitoring with real-time status updates.
    """
    
    def __init__(self, repo: Repo, adapter: JobAdapter):
        """
        Initialize JobWatcher.
        
        Args:
            repo: ChemVCS repository
            adapter: Job scheduler adapter
        """
        self.repo = repo
        self.adapter = adapter
        self.tracker = JobTracker(repo, adapter)
    
    def watch_job(
        self,
        identifier: str,
        interval: int = 30,
        timeout: Optional[int] = None
    ) -> TrackedJob:
        """
        Watch a job until completion or timeout.
        
        Args:
            identifier: Job ID or run hash
            interval: Polling interval in seconds (default: 30)
            timeout: Maximum watch time in seconds (None for indefinite)
            
        Returns:
            Final TrackedJob state
            
        Raises:
            JobNotFoundError: If job/run not found
        """
        start_time = time.time()
        last_status = None
        
        print(f"Watching job: {identifier}")
        print(f"Polling interval: {interval}s")
        if timeout:
            print(f"Timeout: {timeout}s")
        print("-" * 60)
        
        while True:
            # Get current job info
            try:
                job = self._get_job(identifier)
            except JobNotFoundError as e:
                print(f"\nError: {e}")
                raise
            
            # Print status if changed
            if job.status != last_status:
                timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                status_symbol = self._get_status_symbol(job.status)
                print(f"[{timestamp}] {status_symbol} Status: {job.status.value}")
                
                # Print additional info
                if job.run.job_system:
                    print(f"  System: {job.run.job_system}")
                if job.run.queue_name:
                    print(f"  Queue: {job.run.queue_name}")
                
                last_status = job.status
            
            # Check if job completed
            if job.status in [JobStatus.COMPLETED, JobStatus.FAILED, JobStatus.CANCELLED]:
                print("-" * 60)
                print(f"Job {job.status.value.lower()}")
                return job
            
            # Check timeout
            if timeout and (time.time() - start_time) > timeout:
                print("-" * 60)
                print(f"Watch timeout after {timeout}s (job still {job.status.value})")
                return job
            
            # Wait before next poll
            time.sleep(interval)
    
    def watch_multiple(
        self,
        identifiers: list,
        interval: int = 30,
        stop_on_first: bool = False
    ):
        """
        Watch multiple jobs simultaneously.
        
        Args:
            identifiers: List of job IDs or run hashes
            interval: Polling interval in seconds
            stop_on_first: Stop when first job completes
        """
        print(f"Watching {len(identifiers)} jobs")
        print(f"Polling interval: {interval}s")
        print("-" * 60)
        
        active_jobs = {id: None for id in identifiers}
        
        while active_jobs:
            timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            
            for identifier in list(active_jobs.keys()):
                try:
                    job = self._get_job(identifier)
                    
                    # Check if status changed
                    if active_jobs[identifier] != job.status:
                        status_symbol = self._get_status_symbol(job.status)
                        print(f"[{timestamp}] {identifier[:8]}... {status_symbol} {job.status.value}")
                        active_jobs[identifier] = job.status
                    
                    # Remove completed jobs
                    if job.status in [JobStatus.COMPLETED, JobStatus.FAILED, JobStatus.CANCELLED]:
                        print(f"  → Job {identifier[:8]} {job.status.value.lower()}")
                        del active_jobs[identifier]
                        
                        if stop_on_first:
                            print("-" * 60)
                            print(f"Stopping (first job completed)")
                            return
                
                except JobNotFoundError:
                    print(f"[{timestamp}] {identifier[:8]}... ❌ NOT FOUND")
                    del active_jobs[identifier]
            
            if not active_jobs:
                break
            
            time.sleep(interval)
        
        print("-" * 60)
        print("All jobs completed")
    
    def _get_job(self, identifier: str) -> TrackedJob:
        """Get TrackedJob by identifier (job ID or run hash)."""
        # Try as job ID first
        job = self.tracker.get_job_by_id(identifier)
        if job:
            return job
        
        # Try as run hash
        run = self.tracker.find_run_by_job_id(identifier)
        if run and run.job_id:
            return self.tracker.get_job_by_id(run.job_id)
        
        # Try finding run by hash
        try:
            obj = self.repo.get_object(identifier)
            from ..domain.run import Run
            run = Run.from_core_object(obj)
            if run.job_id:
                return self.tracker.get_job_by_id(run.job_id)
        except:
            pass
        
        raise JobNotFoundError(f"Job or run not found: {identifier}")
    
    @staticmethod
    def _get_status_symbol(status: JobStatus) -> str:
        """Get emoji symbol for job status."""
        symbols = {
            JobStatus.PENDING: "⏳",
            JobStatus.RUNNING: "▶️",
            JobStatus.COMPLETED: "✅",
            JobStatus.FAILED: "❌",
            JobStatus.CANCELLED: "🚫",
            JobStatus.UNKNOWN: "❓",
        }
        return symbols.get(status, "•")
