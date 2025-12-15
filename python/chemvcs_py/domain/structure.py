"""Structure domain object for representing molecular/crystal structures."""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import numpy as np

from ..core.objects import CoreObject, Reference
from ..util.errors import ValidationError


@dataclass
class Structure:
    """
    Represents an atomic structure (molecule or crystal).
    
    Attributes:
        formula: Chemical formula (e.g., "H2O", "Si8")
        positions: N x 3 array of atomic positions (Angstrom)
        species: List of element symbols (length N)
        lattice: Optional 3x3 lattice matrix for periodic systems (Angstrom)
        metadata: Additional metadata dict
        id: Optional hash if loaded from repository
    """
    formula: str
    positions: np.ndarray
    species: List[str]
    lattice: Optional[np.ndarray] = None
    metadata: Dict[str, Any] = field(default_factory=dict)
    id: Optional[str] = None

    def __post_init__(self):
        """Validate structure data."""
        # Convert positions to numpy array if not already
        if not isinstance(self.positions, np.ndarray):
            self.positions = np.array(self.positions)
        
        # Validate positions shape
        if self.positions.ndim != 2 or self.positions.shape[1] != 3:
            raise ValidationError(
                f"Positions must be N x 3 array, got shape {self.positions.shape}"
            )
        
        # Validate species length
        if len(self.species) != len(self.positions):
            raise ValidationError(
                f"Species length ({len(self.species)}) must match "
                f"positions length ({len(self.positions)})"
            )
        
        # Convert lattice to numpy array if present
        if self.lattice is not None:
            if not isinstance(self.lattice, np.ndarray):
                self.lattice = np.array(self.lattice)
            
            if self.lattice.shape != (3, 3):
                raise ValidationError(
                    f"Lattice must be 3x3 matrix, got shape {self.lattice.shape}"
                )

    @property
    def num_atoms(self) -> int:
        """Get number of atoms."""
        return len(self.species)

    @property
    def is_periodic(self) -> bool:
        """Check if structure is periodic."""
        return self.lattice is not None

    def to_core_object(self) -> CoreObject:
        """
        Convert Structure to CoreObject for storage.
        
        Returns:
            CoreObject with type="structure"
        """
        meta = {
            "formula": self.formula,
            "positions": self.positions.tolist(),
            "species": self.species,
        }
        
        if self.lattice is not None:
            meta["lattice"] = self.lattice.tolist()
        
        # Merge additional metadata
        meta.update(self.metadata)
        
        obj = CoreObject(
            version=1,
            type="structure",
            meta=meta,
            refs=[]
        )
        
        return obj

    @classmethod
    def from_core_object(cls, obj: CoreObject) -> "Structure":
        """
        Create Structure from CoreObject.
        
        Args:
            obj: CoreObject with type="structure"
            
        Returns:
            Structure instance
            
        Raises:
            ValidationError: If object type is not "structure"
        """
        if obj.type != "structure":
            raise ValidationError(
                f"Expected object type 'structure', got '{obj.type}'"
            )
        
        # Extract required fields
        formula = obj.get_meta("formula")
        positions = np.array(obj.get_meta("positions"))
        species = obj.get_meta("species")
        
        # Extract optional fields
        lattice_data = obj.get_meta("lattice")
        lattice = np.array(lattice_data) if lattice_data is not None else None
        
        # Extract additional metadata (exclude known fields)
        known_keys = {"formula", "positions", "species", "lattice"}
        metadata = {k: v for k, v in obj.meta.items() if k not in known_keys}
        
        return cls(
            formula=formula,
            positions=positions,
            species=species,
            lattice=lattice,
            metadata=metadata,
            id=obj.hash
        )

    def get_center_of_mass(self) -> np.ndarray:
        """
        Calculate center of mass (assuming all atoms have equal mass).
        
        Returns:
            3D center of mass coordinates
        """
        return np.mean(self.positions, axis=0)

    def translate(self, vector: np.ndarray) -> "Structure":
        """
        Create a new structure translated by vector.
        
        Args:
            vector: 3D translation vector
            
        Returns:
            New Structure instance
        """
        new_positions = self.positions + np.array(vector)
        return Structure(
            formula=self.formula,
            positions=new_positions,
            species=self.species.copy(),
            lattice=self.lattice.copy() if self.lattice is not None else None,
            metadata=self.metadata.copy()
        )

    def __repr__(self) -> str:
        periodic = "periodic" if self.is_periodic else "non-periodic"
        id_str = f", id={self.id[:8]}" if self.id else ""
        return (
            f"Structure(formula='{self.formula}', "
            f"num_atoms={self.num_atoms}, "
            f"{periodic}{id_str})"
        )
