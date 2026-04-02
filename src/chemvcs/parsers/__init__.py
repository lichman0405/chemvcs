"""VASP, LAMMPS, and ORCA file parsers with semantic diff algorithms.

This module provides parsers for INCAR, POSCAR, KPOINTS, OUTCAR (VASP),
LAMMPS input scripts, data files and log files, and ORCA input/output files,
along with semantic diff capabilities.
"""

from chemvcs.parsers.base_parser import BaseParser, DiffEntry, ParserError
from chemvcs.parsers.diff_engine import DiffEngine
from chemvcs.parsers.incar_parser import IncarParser
from chemvcs.parsers.kpoints_parser import KpointsParser
from chemvcs.parsers.lammps_data_parser import LammpsDataParser
from chemvcs.parsers.lammps_input_parser import LammpsInputParser
from chemvcs.parsers.lammps_log_parser import LammpsLogParser
from chemvcs.parsers.orca_input_parser import OrcaInputParser
from chemvcs.parsers.orca_output_parser import OrcaOutputParser
from chemvcs.parsers.outcar_parser import OutcarParser

__all__ = [
    "BaseParser",
    "DiffEntry",
    "DiffEngine",
    "ParserError",
    "IncarParser",
    "KpointsParser",
    "LammpsDataParser",
    "LammpsInputParser",
    "LammpsLogParser",
    "OrcaInputParser",
    "OrcaOutputParser",
    "OutcarParser",
]
