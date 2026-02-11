"""INCAR-POSCAR consistency validator.

Validates that INCAR parameters are consistent with POSCAR structure.
For example, MAGMOM and LDAUU should match the number of atoms/elements.
"""

from pathlib import Path
from typing import Any

from chemvcs.plugins.base import ValidationResult, ValidatorPlugin


class INCARPOSCARValidator(ValidatorPlugin):
    """Validates INCAR-POSCAR consistency.
    
    Checks:
    - MAGMOM length matches number of atoms
    - LDAUU/LDAUJ length matches number of element types
    - NIONS matches actual atom count (if specified)
    """
    
    @property
    def name(self) -> str:
        return "incar-poscar"
    
    @property
    def version(self) -> str:
        return "0.1.0"
    
    @property
    def description(self) -> str:
        return "Validates INCAR-POSCAR parameter consistency"
    
    @property
    def priority(self) -> int:
        return 20
    
    @property
    def enabled_by_default(self) -> bool:
        return True
    
    def can_validate(self, files: list[str]) -> bool:
        """Check if both INCAR and POSCAR are present."""
        file_names = {Path(f).name for f in files}
        return "INCAR" in file_names and "POSCAR" in file_names
    
    def validate(
        self,
        workspace_root: Path,
        files: list[str],
        **kwargs: Any,
    ) -> ValidationResult:
        """Validate INCAR-POSCAR consistency.
        
        TODO: Implement full validation logic.
        """
        # Placeholder implementation
        return ValidationResult(
            passed=True,
            message="INCAR-POSCAR validation not yet implemented",
            warnings=["This validator is under development"],
        )
