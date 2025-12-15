"""Job submission logic."""

import os
from typing import Optional, Dict, Any

from ..core.repo import Repo
from ..domain.run import Run
from .adapter import JobAdapter
from .provenance import EnvironmentCapture
from .exceptions import InvalidJobStateError


class JobSubmitter:
    """Handles submitting Run objects to job schedulers."""
    
    def __init__(self, repo: Repo, adapter: JobAdapter):
        """
        Initialize JobSubmitter.
        
        Args:
            repo: ChemVCS repository
            adapter: Job scheduler adapter
        """
        self.repo = repo
        self.adapter = adapter
    
    def submit_run(
        self,
        run: Run,
        script_path: str,
        capture_env: bool = True,
        **submit_kwargs
    ) -> str:
        """
        Submit a Run to job scheduler.
        
        Args:
            run: Run object to submit
            script_path: Path to job submission script
            capture_env: Whether to capture environment provenance
            **submit_kwargs: Additional arguments for adapter.submit()
            
        Returns:
            Job ID from scheduler
            
        Raises:
            InvalidJobStateError: If run is not in 'planned' state
            JobSubmissionError: If submission fails
        """
        # Validate run state
        if run.status != 'planned':
            raise InvalidJobStateError(
                f"Run must be in 'planned' state to submit (current: {run.status})"
            )
        
        # Validate script exists
        if not os.path.exists(script_path):
            raise FileNotFoundError(f"Job script not found: {script_path}")
        
        # Capture environment provenance
        if capture_env:
            self._capture_provenance(run, script_path)
        
        # Submit job
        job_id = self.adapter.submit(script_path, **submit_kwargs)
        
        # Update Run object
        run.mark_submitted()
        run.mark_queued(
            job_id=job_id,
            job_system=self.adapter.name,
            queue=submit_kwargs.get('partition') or submit_kwargs.get('queue')
        )
        
        # Store updated Run to repository
        obj = run.to_core_object()
        hash_id = self.repo.store_object(obj)
        run.id = hash_id
        
        return job_id
    
    def _capture_provenance(self, run: Run, script_path: str) -> None:
        """
        Capture environment provenance and attach to Run.
        
        Args:
            run: Run object to update
            script_path: Path to job script
        """
        # Capture loaded modules
        run.modules_loaded = EnvironmentCapture.capture_modules()
        
        # Capture relevant environment variables
        run.environment_vars = EnvironmentCapture.capture_env_vars()
        
        # Save script snapshot
        run.submit_script = EnvironmentCapture.save_script_snapshot(script_path)
        
        # Parse resource requirements from script
        script_type = EnvironmentCapture.detect_script_type(script_path)
        if script_type == 'slurm':
            run.job_resources = EnvironmentCapture.parse_slurm_script(script_path)
        elif script_type == 'pbs':
            run.job_resources = EnvironmentCapture.parse_pbs_script(script_path)
