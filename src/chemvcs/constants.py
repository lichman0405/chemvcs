"""Constants used throughout ChemVCS."""

# Version
VERSION = "0.1.0"

# Directory names
CHEMVCS_DIR = ".chemvcs"
OBJECTS_DIR = "objects"
COMMITS_DIR = "commits"
STAGING_DIR = "staging"

# File names
METADATA_DB = "metadata.db"
HEAD_FILE = "HEAD"
STAGING_MANIFEST = "manifest.json"
IGNORE_FILE = ".chemvcsignore"
REMOTES_FILE = "remotes.toml"

# Default ignore patterns
DEFAULT_IGNORE_PATTERNS = [
    "WAVECAR",
    "CHGCAR",
    "*.tmp",
    "CHG",
    "vasprun.xml.bz2",
]

# File size limits (bytes)
MAX_INLINE_SIZE = 100 * 1024 * 1024  # 100 MB
GZIP_THRESHOLD = 200 * 1024 * 1024  # 200 MB

# Hash algorithm
HASH_ALGORITHM = "sha256"
HASH_LENGTH = 64  # SHA-256 produces 64 hex characters

# Exit codes (matching CLI_SPEC.md)
EXIT_SUCCESS = 0
EXIT_USER_ERROR = 1
EXIT_SYSTEM_ERROR = 2
EXIT_DATA_ERROR = 3
EXIT_INTERRUPTED = 130

# VASP file types
VASP_INPUT_FILES = {
    "INCAR": "input_parameters",
    "POSCAR": "structure",
    "KPOINTS": "k_points",
    "POTCAR": "pseudopotential",
}

VASP_OUTPUT_FILES = {
    "OUTCAR": "output_log",
    "OSZICAR": "convergence_log",
    "CONTCAR": "final_structure",
    "vasprun.xml": "full_output",
    "EIGENVAL": "eigenvalues",
    "DOSCAR": "density_of_states",
}

# LAMMPS file types
# Input scripts:  recognised by prefix/suffix (in.*, *.lammps, lammps.in)
# Data files:     recognised by prefix/suffix (data.*, *.lmp, *.data)
# Log files:      recognised by name          (log.lammps, log.*, *.log)
LAMMPS_INPUT_FILES = {
    "LAMMPS_INPUT": "md_input_script",
    "LAMMPS_DATA": "md_structure_topology",
}

LAMMPS_OUTPUT_FILES = {
    "LAMMPS_LOG": "md_thermo_log",
    "LAMMPS_DUMP": "md_trajectory",  # tracked only, not semantically parsed
    "LAMMPS_RESTART": "md_restart",  # binary, tracked only
}

# ORCA file types
# Input scripts:  recognised by suffix (*.inp)
# Output files:   recognised by suffix (*.out, not OUTCAR)
ORCA_INPUT_FILES = {
    "ORCA_INPUT": "orca_input_script",
}

ORCA_OUTPUT_FILES = {
    "ORCA_OUTPUT": "orca_output_log",
}

# LAMMPS files that are typically large and should be ignored by default
LAMMPS_DEFAULT_IGNORE = [
    "dump.*",
    "*.dump",
    "restart.*",
    "*.restart",
    "*.bin",
]

# Database schema version
DB_SCHEMA_VERSION = 1
