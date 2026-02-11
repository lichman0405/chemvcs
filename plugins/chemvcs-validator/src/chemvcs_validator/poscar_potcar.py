"""POSCAR-POTCAR element order validator.

Ensures that element order in POSCAR matches POTCAR.
This is a critical validation to prevent silently wrong calculations.
"""

from pathlib import Path
from typing import Any

from chemvcs.plugins.base import ValidationError, ValidationResult, ValidatorPlugin

try:
    from pymatgen.io.vasp.inputs import Poscar
except ImportError:
    raise ImportError("pymatgen is required for POSCAR-POTCAR validation")


class POSCARPOTCARValidator(ValidatorPlugin):
    """Validates POSCAR-POTCAR element order consistency.
    
    VASP requires that element order in POSCAR (line 6) exactly matches
    the order in which POTCAR pseudopotential blocks are concatenated.
    
    Example:
        POSCAR: Li Co O
        POTCAR: Li_sv block + Co block + O block  ✓ Correct
        POTCAR: Li_sv block + O block + Co block  ✗ Wrong!
    """
    
    @property
    def name(self) -> str:
        return "poscar-potcar"
    
    @property
    def version(self) -> str:
        return "0.1.0"
    
    @property
    def description(self) -> str:
        return "Validates POSCAR-POTCAR element order consistency"
    
    @property
    def priority(self) -> int:
        return 10  # High priority (run early)
    
    @property
    def enabled_by_default(self) -> bool:
        return True
    
    def can_validate(self, files: list[str]) -> bool:
        """Check if both POSCAR and POTCAR are present."""
        file_names = {Path(f).name for f in files}
        return "POSCAR" in file_names and "POTCAR" in file_names
    
    def validate(
        self,
        workspace_root: Path,
        files: list[str],
        **kwargs: Any,
    ) -> ValidationResult:
        """Validate POSCAR-POTCAR element match.
        
        Args:
            workspace_root: Root directory of the workspace
            files: List of files to validate
            **kwargs: Additional context
        
        Returns:
            ValidationResult
        """
        poscar_path = workspace_root / "POSCAR"
        potcar_path = workspace_root / "POTCAR"
        
        # Check file existence
        if not poscar_path.exists():
            return ValidationResult(
                passed=True,
                message="POSCAR not found, skipping validation",
            )
        
        if not potcar_path.exists():
            return ValidationResult(
                passed=True,
                message="POTCAR not found, skipping validation",
                warnings=["POTCAR missing - recommend adding for completeness"],
            )
        
        try:
            # Extract POSCAR elements
            poscar_elements = self._extract_poscar_elements(poscar_path)
            
            # Extract POTCAR elements
            potcar_elements = self._extract_potcar_elements(potcar_path)
            
            # Compare
            if poscar_elements == potcar_elements:
                return ValidationResult(
                    passed=True,
                    message=f"POSCAR-POTCAR match: {poscar_elements}",
                    details={
                        "poscar_elements": poscar_elements,
                        "potcar_elements": potcar_elements,
                    },
                )
            else:
                # Mismatch!
                error_msg = self._format_mismatch_error(poscar_elements, potcar_elements)
                
                return ValidationResult(
                    passed=False,
                    message="POSCAR-POTCAR element order mismatch",
                    errors=[error_msg],
                    details={
                        "poscar_elements": poscar_elements,
                        "potcar_elements": potcar_elements,
                    },
                )
        
        except Exception as e:
            # Validation failed due to error
            return ValidationResult(
                passed=False,
                message=f"Validation error: {e}",
                errors=[str(e)],
            )
    
    def _extract_poscar_elements(self, poscar_path: Path) -> list[str]:
        """Extract element list from POSCAR.
        
        Args:
            poscar_path: Path to POSCAR file
        
        Returns:
            List of element symbols in order
        
        Raises:
            ValueError: If POSCAR format is invalid
        """
        try:
            poscar = Poscar.from_file(str(poscar_path))
            elements = poscar.site_symbols
            
            if not elements:
                raise ValueError(
                    "VASP 4.x POSCAR format detected (no element symbols). "
                    "Please add element symbols line (VASP 5+ format)."
                )
            
            return elements
        
        except Exception as e:
            raise ValueError(f"Failed to parse POSCAR: {e}") from e
    
    def _extract_potcar_elements(self, potcar_path: Path) -> list[str]:
        """Extract element list from POTCAR.
        
        Args:
            potcar_path: Path to POTCAR file
        
        Returns:
            List of element symbols in order
        
        Raises:
            ValueError: If POTCAR format is invalid
        """
        try:
            with open(potcar_path, "r", encoding="utf-8", errors="ignore") as f:
                content = f.read()
            
            # Find all TITEL lines
            # Format: "   TITEL  = PAW_PBE Li_sv 17Jan2003"
            elements = []
            for line in content.split("\n"):
                if "TITEL" in line and "PAW" in line:
                    parts = line.split()
                    if len(parts) >= 4:
                        # Extract element (may have suffix like Li_sv, Fe_pv)
                        element_raw = parts[3]
                        # Remove suffix
                        element = element_raw.split("_")[0]
                        elements.append(element)
            
            if not elements:
                raise ValueError("No TITEL lines found in POTCAR")
            
            return elements
        
        except Exception as e:
            raise ValueError(f"Failed to parse POTCAR: {e}") from e
    
    def _format_mismatch_error(
        self,
        poscar_elements: list[str],
        potcar_elements: list[str],
    ) -> str:
        """Format a user-friendly error message for element mismatch."""
        msg = f"\nPOSCAR elements (line 6): {poscar_elements}\n"
        msg += f"POTCAR elements (TITEL):  {potcar_elements}\n"
        
        # Detect specific issues
        if set(poscar_elements) != set(potcar_elements):
            missing = set(poscar_elements) - set(potcar_elements)
            extra = set(potcar_elements) - set(poscar_elements)
            
            if missing:
                msg += f"\n❌ Missing in POTCAR: {list(missing)}"
            if extra:
                msg += f"\n❌ Extra in POTCAR: {list(extra)}"
        else:
            # Same elements, wrong order
            msg += "\n⚠️  Same elements but WRONG ORDER!"
            
            # Highlight differences
            for i, (p, q) in enumerate(zip(poscar_elements, potcar_elements)):
                if p != q:
                    msg += f"\n   Position {i+1}: POSCAR={p}, POTCAR={q} ← MISMATCH"
        
        msg += "\n\nHow to fix:"
        msg += "\n  1. Regenerate POTCAR with correct element order:"
        potcar_cmd = " ".join([f"{el}/POTCAR" for el in poscar_elements])
        msg += f"\n     cat {potcar_cmd} > POTCAR"
        msg += "\n  2. Or reorder elements in POSCAR line 6"
        
        return msg
