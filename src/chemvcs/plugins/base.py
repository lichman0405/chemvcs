"""Base classes for ChemVCS plugins."""

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Optional


class Plugin(ABC):
    """Base class for all ChemVCS plugins.
    
    Plugins can extend ChemVCS with additional functionality like:
    - Input file validators
    - Custom file parsers
    - Pre/post-commit hooks
    - Custom diff engines
    """
    
    @property
    @abstractmethod
    def name(self) -> str:
        """Plugin name (unique identifier)."""
        pass
    
    @property
    @abstractmethod
    def version(self) -> str:
        """Plugin version."""
        pass
    
    @property
    def description(self) -> str:
        """Plugin description (optional)."""
        return ""
    
    def activate(self) -> None:
        """Called when plugin is loaded.
        
        Override to perform initialization tasks.
        """
        pass
    
    def deactivate(self) -> None:
        """Called when plugin is unloaded.
        
        Override to perform cleanup tasks.
        """
        pass


class ValidatorPlugin(Plugin):
    """Base class for input file validators.
    
    Validators check consistency and correctness of input files.
    They run during `chemvcs add` or `chemvcs commit` operations.
    
    Example:
        class POSCARPOTCARValidator(ValidatorPlugin):
            name = "poscar-potcar"
            version = "1.0.0"
            
            def validate(self, workspace_root: Path, files: list[str]) -> ValidationResult:
                # Check POSCAR-POTCAR element order match
                ...
    """
    
    @property
    def priority(self) -> int:
        """Validation priority (lower runs first).
        
        Default is 50. Use 0-25 for high priority, 75-100 for low priority.
        """
        return 50
    
    @property
    def enabled_by_default(self) -> bool:
        """Whether this validator is enabled by default.
        
        Default is True. Set to False for optional validators.
        """
        return True
    
    @abstractmethod
    def validate(
        self,
        workspace_root: Path,
        files: list[str],
        **kwargs: Any,
    ) -> "ValidationResult":
        """Perform validation on specified files.
        
        Args:
            workspace_root: Root directory of the workspace
            files: List of files to validate (relative paths)
            **kwargs: Additional context (e.g., commit message, author)
        
        Returns:
            ValidationResult object
            
        Raises:
            ValidationError: If validation fails critically
        """
        pass
    
    def can_validate(self, files: list[str]) -> bool:
        """Check if this validator can handle the given files.
        
        Override to implement custom file selection logic.
        Default returns True if any required files are present.
        
        Args:
            files: List of file paths
            
        Returns:
            True if validator should run, False otherwise
        """
        return True


class ValidationResult:
    """Result of a validation operation."""
    
    def __init__(
        self,
        passed: bool,
        message: str = "",
        warnings: Optional[list[str]] = None,
        errors: Optional[list[str]] = None,
        details: Optional[dict[str, Any]] = None,
    ):
        """Initialize validation result.
        
        Args:
            passed: Whether validation passed
            message: Main result message
            warnings: List of warning messages
            errors: List of error messages
            details: Additional details (for debugging)
        """
        self.passed = passed
        self.message = message
        self.warnings = warnings or []
        self.errors = errors or []
        self.details = details or {}
    
    @property
    def has_warnings(self) -> bool:
        """Check if result has warnings."""
        return len(self.warnings) > 0
    
    @property
    def has_errors(self) -> bool:
        """Check if result has errors."""
        return len(self.errors) > 0
    
    def __repr__(self) -> str:
        status = "PASSED" if self.passed else "FAILED"
        return f"ValidationResult(status={status}, warnings={len(self.warnings)}, errors={len(self.errors)})"


class ValidationError(Exception):
    """Raised when validation fails critically."""
    pass
