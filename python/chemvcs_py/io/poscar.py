"""POSCAR file format parser for VASP calculations."""

from pathlib import Path
from typing import Union, List, Optional

import numpy as np

from ..domain.structure import Structure
from ..util.errors import ParseError


def read_poscar(filepath: Union[str, Path]) -> Structure:
    """
    Read structure from VASP POSCAR/CONTCAR file.
    
    POSCAR format (VASP 5.x):
        Line 1: Comment line (system name)
        Line 2: Universal scaling factor
        Lines 3-5: Lattice vectors (3x3 matrix)
        Line 6: Element symbols
        Line 7: Number of atoms per element
        Line 8: Coordinate type (Direct/Cartesian)
        Lines 9+: Atomic positions
    
    Args:
        filepath: Path to POSCAR file
        
    Returns:
        Structure instance
        
    Raises:
        ParseError: If file format is invalid
    
    Example POSCAR:
        Cubic BN
        3.57
        0.0 0.5 0.5
        0.5 0.0 0.5
        0.5 0.5 0.0
        B N
        1 1
        Direct
        0.00 0.00 0.00
        0.25 0.25 0.25
    """
    filepath = Path(filepath)
    
    if not filepath.exists():
        raise ParseError(f"File not found: {filepath}")
    
    try:
        with open(filepath, 'r') as f:
            lines = [line.strip() for line in f.readlines()]
        
        if len(lines) < 9:
            raise ParseError("POSCAR file must have at least 9 lines")
        
        # Line 1: Comment/system name
        system_name = lines[0]
        
        # Line 2: Scaling factor
        try:
            scaling_factor = float(lines[1])
        except ValueError as e:
            raise ParseError(f"Invalid scaling factor on line 2: {lines[1]}") from e
        
        # Lines 3-5: Lattice vectors
        lattice_lines = lines[2:5]
        lattice = []
        for i, line in enumerate(lattice_lines, start=3):
            parts = line.split()
            if len(parts) < 3:
                raise ParseError(f"Line {i}: Expected 3 lattice components, got {len(parts)}")
            try:
                vector = [float(x) for x in parts[:3]]
                lattice.append(vector)
            except ValueError as e:
                raise ParseError(f"Line {i}: Invalid lattice vector '{line}'") from e
        
        lattice_array = np.array(lattice) * scaling_factor
        
        # Line 6: Element symbols (VASP 5.x format)
        elements_line = lines[5]
        elements = elements_line.split()
        
        # Line 7: Number of atoms per element
        counts_line = lines[6]
        try:
            counts = [int(x) for x in counts_line.split()]
        except ValueError as e:
            raise ParseError(f"Line 7: Invalid atom counts '{counts_line}'") from e
        
        if len(elements) != len(counts):
            raise ParseError(
                f"Number of elements ({len(elements)}) doesn't match "
                f"number of counts ({len(counts)})"
            )
        
        # Build species list
        species = []
        for element, count in zip(elements, counts):
            species.extend([element] * count)
        
        total_atoms = sum(counts)
        
        # Line 8: Coordinate mode
        coord_mode_line = lines[7].lower()
        if coord_mode_line.startswith('s'):
            # Selective dynamics - skip to next line
            coord_mode_line = lines[8].lower()
            coord_start = 9
        else:
            coord_start = 8
        
        is_direct = coord_mode_line.startswith('d')
        
        # Read atomic positions
        if len(lines) < coord_start + total_atoms:
            raise ParseError(
                f"Expected {total_atoms} position lines, "
                f"but file has only {len(lines) - coord_start} lines"
            )
        
        positions = []
        for i in range(total_atoms):
            line_num = coord_start + i
            parts = lines[line_num].split()
            
            if len(parts) < 3:
                raise ParseError(
                    f"Line {line_num + 1}: Expected 3 coordinates, got {len(parts)}"
                )
            
            try:
                pos = [float(x) for x in parts[:3]]
            except ValueError as e:
                raise ParseError(
                    f"Line {line_num + 1}: Invalid coordinates '{parts[:3]}'"
                ) from e
            
            positions.append(pos)
        
        positions_array = np.array(positions)
        
        # Convert direct (fractional) to Cartesian if needed
        if is_direct:
            # positions_cart = positions_frac @ lattice
            positions_array = positions_array @ lattice_array
        else:
            # Already Cartesian, apply scaling
            positions_array = positions_array * scaling_factor
        
        # Generate formula
        formula = _generate_formula(elements, counts)
        
        metadata = {
            "title": system_name,
            "source_file": str(filepath.name),
            "coordinate_mode": "direct" if is_direct else "cartesian"
        }
        
        return Structure(
            formula=formula,
            positions=positions_array,
            species=species,
            lattice=lattice_array,
            metadata=metadata
        )
    
    except ParseError:
        raise
    except Exception as e:
        raise ParseError(f"Failed to parse POSCAR file: {e}") from e


def write_poscar(structure: Structure, filepath: Union[str, Path], 
                 comment: Optional[str] = None, direct: bool = True) -> None:
    """
    Write structure to VASP POSCAR format.
    
    Args:
        structure: Structure to write
        filepath: Output file path
        comment: Optional comment line (defaults to formula)
        direct: If True, write fractional coordinates; if False, Cartesian
        
    Raises:
        ValueError: If structure is not periodic or has no atoms
    """
    filepath = Path(filepath)
    
    if structure.num_atoms == 0:
        raise ValueError("Cannot write empty structure to POSCAR")
    
    if not structure.is_periodic:
        raise ValueError("POSCAR format requires periodic structure with lattice")
    
    with open(filepath, 'w') as f:
        # Line 1: Comment
        comment_line = comment or structure.metadata.get("title", structure.formula)
        f.write(f"{comment_line}\n")
        
        # Line 2: Scaling factor (use 1.0, include scale in lattice)
        f.write("1.0\n")
        
        # Lines 3-5: Lattice vectors
        for i in range(3):
            f.write(f"  {structure.lattice[i, 0]:16.10f}  "
                   f"{structure.lattice[i, 1]:16.10f}  "
                   f"{structure.lattice[i, 2]:16.10f}\n")
        
        # Lines 6-7: Element symbols and counts
        # Group consecutive same elements
        from collections import OrderedDict
        element_counts = OrderedDict()
        for elem in structure.species:
            if elem not in element_counts:
                element_counts[elem] = 0
            element_counts[elem] += 1
        
        # Line 6: Elements
        f.write("  " + "  ".join(element_counts.keys()) + "\n")
        
        # Line 7: Counts
        f.write("  " + "  ".join(str(c) for c in element_counts.values()) + "\n")
        
        # Line 8: Coordinate mode
        if direct:
            f.write("Direct\n")
            # Convert Cartesian to fractional
            # positions_frac = positions_cart @ inv(lattice)
            lattice_inv = np.linalg.inv(structure.lattice)
            positions = structure.positions @ lattice_inv
        else:
            f.write("Cartesian\n")
            positions = structure.positions
        
        # Lines 9+: Atomic positions
        for pos in positions:
            f.write(f"  {pos[0]:16.10f}  {pos[1]:16.10f}  {pos[2]:16.10f}\n")


def _generate_formula(elements: List[str], counts: List[int]) -> str:
    """Generate chemical formula from elements and counts."""
    formula_parts = []
    for element, count in zip(elements, counts):
        if count == 1:
            formula_parts.append(element)
        else:
            formula_parts.append(f"{element}{count}")
    return "".join(formula_parts)
