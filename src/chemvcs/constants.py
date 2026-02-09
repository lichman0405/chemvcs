"""Constants used throughout ChemVCS."""

from pathlib import Path

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
GZIP_THRESHOLD = 200 * 1024 * 1024   # 200 MB

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

# Database schema version
DB_SCHEMA_VERSION = 1
