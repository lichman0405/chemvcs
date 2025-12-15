"""XYZ file format parser for molecular structures."""

from pathlib import Path
from typing import Union, Optional

import numpy as np

from ..domain.structure import Structure
from ..util.errors import ParseError


def read_xyz(filepath: Union[str, Path]) -> Structure:
    """
    Read structure from XYZ file.
    
    XYZ format:
        Line 1: Number of atoms
        Line 2: Comment line (optional title)
        Lines 3+: Element X Y Z (in Angstrom)
    
    Args:
        filepath: Path to XYZ file
        
    Returns:
        Structure instance
        
    Raises:
        ParseError: If file format is invalid
    
    Example:
        3
        Water molecule
        O  0.000  0.000  0.000
        H  0.758  0.587  0.000
        H -0.758  0.587  0.000
    """
    filepath = Path(filepath)
    
    if not filepath.exists():
        raise ParseError(f"File not found: {filepath}")
    
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
        
        if len(lines) < 3:
            raise ParseError("XYZ file must have at least 3 lines")
        
        # Parse number of atoms
        try:
            num_atoms = int(lines[0].strip())
        except ValueError as e:
            raise ParseError(f"Invalid atom count on line 1: {lines[0]}") from e
        
        # Get comment line (formula)
        comment = lines[1].strip()
        
        # Parse atom lines
        if len(lines) < 2 + num_atoms:
            raise ParseError(
                f"Expected {num_atoms} atom lines, "
                f"but file has only {len(lines) - 2} lines"
            )
        
        species = []
        positions = []
        
        for i, line in enumerate(lines[2:2+num_atoms], start=3):
            parts = line.strip().split()
            
            if len(parts) < 4:
                raise ParseError(
                    f"Line {i}: Expected 'Element X Y Z', got '{line.strip()}'"
                )
            
            element = parts[0]
            try:
                x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
            except ValueError as e:
                raise ParseError(
                    f"Line {i}: Invalid coordinates '{parts[1:4]}'"
                ) from e
            
            species.append(element)
            positions.append([x, y, z])
        
        positions_array = np.array(positions)
        
        # Generate formula from species
        formula = _generate_formula(species)
        
        # Use comment as title if it's not empty
        metadata = {}
        if comment:
            metadata["title"] = comment
        metadata["source_file"] = str(filepath.name)
        
        return Structure(
            formula=formula,
            positions=positions_array,
            species=species,
            lattice=None,  # XYZ format doesn't include lattice
            metadata=metadata
        )
    
    except ParseError:
        raise
    except Exception as e:
        raise ParseError(f"Failed to parse XYZ file: {e}") from e


def write_xyz(structure: Structure, filepath: Union[str, Path], 
              comment: Optional[str] = None) -> None:
    """
    Write structure to XYZ file.
    
    Args:
        structure: Structure to write
        filepath: Output file path
        comment: Optional comment line (defaults to structure title or formula)
        
    Raises:
        ValueError: If structure has no atoms
    """
    filepath = Path(filepath)
    
    if structure.num_atoms == 0:
        raise ValueError("Cannot write empty structure to XYZ")
    
    with open(filepath, 'w') as f:
        # Line 1: Number of atoms
        f.write(f"{structure.num_atoms}\n")
        
        # Line 2: Comment (use provided comment, title from metadata, or formula)
        if comment is not None:
            title = comment
        else:
            title = structure.metadata.get("title", structure.formula)
        f.write(f"{title}\n")
        
        # Lines 3+: Atom coordinates
        for i, (element, pos) in enumerate(zip(structure.species, structure.positions)):
            f.write(f"{element:2s}  {pos[0]:12.6f}  {pos[1]:12.6f}  {pos[2]:12.6f}\n")


def _generate_formula(species: list) -> str:
    """
    Generate chemical formula from species list.
    
    Example: ['H', 'H', 'O'] -> 'H2O'
    """
    from collections import Counter
    
    counts = Counter(species)
    
    # Sort by element (roughly by electronegativity, but simple alphabetical for now)
    # In real chemistry, you'd want C, H first, then alphabetical
    formula_parts = []
    for element in sorted(counts.keys()):
        count = counts[element]
        if count == 1:
            formula_parts.append(element)
        else:
            formula_parts.append(f"{element}{count}")
    
    return "".join(formula_parts)
