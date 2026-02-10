"""VASP file parsers and semantic diff algorithms.

This module provides parsers for INCAR, POSCAR, KPOINTS, OUTCAR and other
VASP files, along with semantic diff capabilities.
"""

from chemvcs.parsers.base_parser import BaseParser, DiffEntry, ParserError
from chemvcs.parsers.incar_parser import IncarParser

__all__ = [
    "BaseParser",
    "DiffEntry",
    "ParserError",
    "IncarParser",
]
