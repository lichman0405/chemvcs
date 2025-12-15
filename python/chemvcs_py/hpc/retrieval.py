"""Job result retrieval."""

import os
import shutil
from typing import Optional, Dict, Any, List

from ..core.repo import Repo
from ..domain.run import Run
from .adapter import JobAdapter, JobStatus
from .exceptions import InvalidJobStateError, JobNotFoundError


class JobRetriever:
    """Retrieve completed job results and update repository."""
    
    def __init__(self, repo: Repo, adapter: JobAdapter):
        """
        Initialize JobRetriever.
        
        Args:
            repo: ChemVCS repository
            adapter: Job scheduler adapter
        """
        self.repo = repo
        self.adapter = adapter
    
    def retrieve_results(
        self,
        job_id: str,
        output_files: Optional[List[str]] = None,
        output_dir: Optional[str] = None,
        auto_commit: bool = True
    ) -> Run:
        """
        Retrieve results for a completed job.
        
        Args:
            job_id: Job identifier
            output_files: Optional list of output files to copy
            output_dir: Directory to copy outputs to (default: working dir)
            auto_commit: Whether to create automatic commit
            
        Returns:
            Updated Run object
            
        Raises:
            JobNotFoundError: If job or Run not found
            InvalidJobStateError: If job is not completed
        """
        # Find Run object by job ID
        run = self._find_run_by_job_id(job_id)
        if not run:
            raise JobNotFoundError(f"No Run found for job ID: {job_id}")
        
        # Check job status
        try:
            status = self.adapter.get_status(job_id)
        except JobNotFoundError:
            # Job purged from scheduler, assume completed if Run says so
            if run.status == "finished":
                status = JobStatus.COMPLETED
            else:
                raise
        
        if status != JobStatus.COMPLETED:
            raise InvalidJobStateError(
                f"Job {job_id} is not completed (status: {status.value})"
            )
        
        # Copy output files if specified
        if output_files:
            self._copy_output_files(output_files, output_dir)
        
        # Update Run object
        if run.status != "finished":
            run.mark_finished()
        run.mark_retrieved()
        
        # Store updated Run
        obj = run.to_core_object()
        hash_id = self.repo.store_object(obj)
        run.id = hash_id
        
        # Create automatic commit
        if auto_commit:
            self.repo.commit(f"Retrieved results for job {job_id}")
        
        return run
    
    def parse_output_data(
        self,
        run: Run,
        output_files: Dict[str, str]
    ) -> Dict[str, Any]:
        """
        Parse output files and extract results.
        
        This is a placeholder for code-specific parsing logic.
        Subclasses should implement parsing for specific codes.
        
        Args:
            run: Run object
            output_files: Mapping of file types to file paths
            
        Returns:
            Dictionary of parsed results
        """
        # Placeholder - would need code-specific parsers
        results = {}
        
        # Example: Parse VASP OUTCAR
        if run.code.lower() == 'vasp' and 'outcar' in output_files:
            results.update(self._parse_vasp_outcar(output_files['outcar']))
        
        return results
    
    def _find_run_by_job_id(self, job_id: str) -> Optional[Run]:
        """Find Run object by job ID."""
        objects = self.repo.list_objects(type_filter='run')
        
        for obj_info in objects:
            try:
                obj = self.repo.get_object(obj_info['hash'])
                run = Run.from_core_object(obj)
                
                if run.job_id == job_id:
                    return run
            except Exception:
                continue
        
        return None
    
    def _copy_output_files(
        self,
        output_files: List[str],
        output_dir: Optional[str] = None
    ) -> None:
        """
        Copy output files to destination directory.
        
        Args:
            output_files: List of file paths to copy
            output_dir: Destination directory (default: current directory)
        """
        if output_dir is None:
            output_dir = os.getcwd()
        
        os.makedirs(output_dir, exist_ok=True)
        
        for file_path in output_files:
            if os.path.exists(file_path):
                dest_path = os.path.join(output_dir, os.path.basename(file_path))
                shutil.copy2(file_path, dest_path)
    
    def _parse_vasp_outcar(self, outcar_path: str) -> Dict[str, Any]:
        """
        Parse VASP OUTCAR file (placeholder).
        
        Args:
            outcar_path: Path to OUTCAR file
            
        Returns:
            Dictionary with parsed data
        """
        results = {}
        
        try:
            with open(outcar_path, 'r') as f:
                content = f.read()
            
            # Simple parsing - extract final energy
            for line in content.split('\n'):
                if 'energy  without entropy' in line:
                    parts = line.split()
                    if len(parts) >= 4:
                        results['total_energy'] = float(parts[3])
        except (FileNotFoundError, ValueError):
            pass
        
        return results
