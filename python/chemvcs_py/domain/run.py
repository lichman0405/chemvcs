"""Run domain object for representing computational calculations."""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional
from datetime import datetime

from ..core.objects import CoreObject, Reference
from ..util.errors import ValidationError


@dataclass
class Run:
    """
    Represents a single computational run (e.g., DFT calculation).
    
    Attributes:
        structure_id: Hash of the Structure object used
        code: Computational chemistry code name (e.g., "VASP", "Gaussian", "ORCA")
        code_version: Optional version of the code
        parameters: Input parameters as dict
        status: Run status ("planned", "submitted", "running", "finished", "failed")
        results: Output results as dict (energies, forces, etc.)
        resources: Computational resources requested/used
        metadata: Additional metadata
        id: Optional hash if loaded from repository
        
        # HPC Integration (M6)
        job_id: Job ID from scheduler (SLURM, PBS, etc.)
        job_system: Job scheduler type ("slurm", "pbs", "sge", "local")
        queue_name: Queue/partition name where job was submitted
        submit_script: Full content of the job submission script
        modules_loaded: List of environment modules loaded
        environment_vars: Environment variables captured at submission
        job_resources: Detailed resource requirements from script
    """
    structure_id: str
    code: str
    parameters: Dict[str, Any] = field(default_factory=dict)
    code_version: Optional[str] = None
    status: str = "planned"
    results: Dict[str, Any] = field(default_factory=dict)
    resources: Dict[str, Any] = field(default_factory=dict)
    metadata: Dict[str, Any] = field(default_factory=dict)
    id: Optional[str] = None
    
    # HPC fields (M6)
    job_id: Optional[str] = None
    job_system: Optional[str] = None
    queue_name: Optional[str] = None
    submit_script: Optional[str] = None
    modules_loaded: List[str] = field(default_factory=list)
    environment_vars: Dict[str, str] = field(default_factory=dict)
    job_resources: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        """Validate run data."""
        # Validate status
        valid_statuses = {"planned", "submitted", "running", "finished", "failed", "cancelled"}
        if self.status not in valid_statuses:
            raise ValidationError(
                f"Invalid status '{self.status}'. Must be one of: {valid_statuses}"
            )
        
        # Validate structure_id is not empty
        if not self.structure_id:
            raise ValidationError("structure_id cannot be empty")
        
        # Validate code is not empty
        if not self.code:
            raise ValidationError("code cannot be empty")

    @property
    def is_complete(self) -> bool:
        """Check if run has completed (successfully or failed)."""
        return self.status in {"finished", "failed", "cancelled"}

    @property
    def is_successful(self) -> bool:
        """Check if run completed successfully."""
        return self.status == "finished"

    def get_energy(self, key: str = "total_energy") -> Optional[float]:
        """
        Get energy value from results.
        
        Args:
            key: Energy key in results dict (default: "total_energy")
            
        Returns:
            Energy value or None if not found
        """
        return self.results.get(key)

    def set_result(self, key: str, value: Any) -> None:
        """Set a result value."""
        self.results[key] = value

    def mark_submitted(self, job_id: Optional[str] = None) -> None:
        """
        Mark run as submitted.
        
        Args:
            job_id: Optional job ID from scheduler
        """
        self.status = "submitted"
        if job_id:
            self.metadata["job_id"] = job_id
            self.metadata["submitted_at"] = datetime.now().isoformat()

    def mark_running(self) -> None:
        """Mark run as running."""
        self.status = "running"
        self.metadata["started_at"] = datetime.now().isoformat()

    def mark_finished(self) -> None:
        """Mark run as successfully finished."""
        self.status = "finished"
        self.metadata["finished_at"] = datetime.now().isoformat()

    def mark_failed(self, error: Optional[str] = None) -> None:
        """
        Mark run as failed.
        
        Args:
            error: Optional error message
        """
        self.status = "failed"
        self.metadata["failed_at"] = datetime.now().isoformat()
        if error:
            self.metadata["error"] = error
    
    def mark_queued(self, job_id: str, job_system: str, queue: Optional[str] = None) -> None:
        """
        Mark run as queued in job scheduler.
        
        Args:
            job_id: Job ID from scheduler
            job_system: Scheduler type (e.g., "slurm", "pbs")
            queue: Optional queue/partition name
        """
        self.job_id = job_id
        self.job_system = job_system
        if queue:
            self.queue_name = queue
        self.metadata["queued_at"] = datetime.now().isoformat()
    
    def mark_retrieved(self, output_data: Optional[Dict[str, Any]] = None) -> None:
        """
        Mark run as retrieved (outputs fetched from compute system).
        
        Args:
            output_data: Optional output data to store in results
        """
        self.metadata["retrieved_at"] = datetime.now().isoformat()
        if output_data:
            self.results.update(output_data)

    def to_core_object(self) -> CoreObject:
        """
        Convert Run to CoreObject for storage.
        
        Returns:
            CoreObject with type="run"
        """
        meta = {
            "structure_id": self.structure_id,
            "code": self.code,
            "status": self.status,
            "parameters": self.parameters,
            "results": self.results,
            "resources": self.resources,
        }
        
        if self.code_version:
            meta["code_version"] = self.code_version
        
        # Add HPC fields if present (M6)
        if self.job_id:
            meta["job_id"] = self.job_id
        if self.job_system:
            meta["job_system"] = self.job_system
        if self.queue_name:
            meta["queue_name"] = self.queue_name
        if self.submit_script:
            meta["submit_script"] = self.submit_script
        if self.modules_loaded:
            meta["modules_loaded"] = self.modules_loaded
        if self.environment_vars:
            meta["environment_vars"] = self.environment_vars
        if self.job_resources:
            meta["job_resources"] = self.job_resources
        
        # Merge additional metadata
        meta.update(self.metadata)
        
        # Create reference to structure
        refs = [Reference(kind="object", id=self.structure_id)]
        
        obj = CoreObject(
            version=1,
            type="run",
            meta=meta,
            refs=refs
        )
        
        return obj

    @classmethod
    def from_core_object(cls, obj: CoreObject) -> "Run":
        """
        Create Run from CoreObject.
        
        Args:
            obj: CoreObject with type="run"
            
        Returns:
            Run instance
            
        Raises:
            ValidationError: If object type is not "run"
        """
        if obj.type != "run":
            raise ValidationError(
                f"Expected object type 'run', got '{obj.type}'"
            )
        
        # Extract required fields
        structure_id = obj.get_meta("structure_id")
        code = obj.get_meta("code")
        status = obj.get_meta("status", "planned")
        
        # Extract optional fields
        code_version = obj.get_meta("code_version")
        parameters = obj.get_meta("parameters", {})
        results = obj.get_meta("results", {})
        resources = obj.get_meta("resources", {})
        
        # Extract HPC fields (M6)
        job_id = obj.get_meta("job_id")
        job_system = obj.get_meta("job_system")
        queue_name = obj.get_meta("queue_name")
        submit_script = obj.get_meta("submit_script")
        modules_loaded = obj.get_meta("modules_loaded", [])
        environment_vars = obj.get_meta("environment_vars", {})
        job_resources = obj.get_meta("job_resources", {})
        
        # Extract additional metadata (exclude known fields)
        known_keys = {
            "structure_id", "code", "code_version", "status",
            "parameters", "results", "resources",
            # HPC fields
            "job_id", "job_system", "queue_name", "submit_script",
            "modules_loaded", "environment_vars", "job_resources"
        }
        metadata = {k: v for k, v in obj.meta.items() if k not in known_keys}
        
        return cls(
            structure_id=structure_id,
            code=code,
            code_version=code_version,
            parameters=parameters,
            status=status,
            results=results,
            resources=resources,
            metadata=metadata,
            id=obj.hash,
            # HPC fields
            job_id=job_id,
            job_system=job_system,
            queue_name=queue_name,
            submit_script=submit_script,
            modules_loaded=modules_loaded,
            environment_vars=environment_vars,
            job_resources=job_resources
        )

    def __repr__(self) -> str:
        status_emoji = {
            "planned": "📋",
            "submitted": "📤",
            "running": "▶️",
            "finished": "✅",
            "failed": "❌",
            "cancelled": "🚫"
        }
        emoji = status_emoji.get(self.status, "❓")
        id_str = f", id={self.id[:8]}" if self.id else ""
        return (
            f"Run({emoji} {self.status}, "
            f"code={self.code}, "
            f"structure={self.structure_id[:8]}{id_str})"
        )
