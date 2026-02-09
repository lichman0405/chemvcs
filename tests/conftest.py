"""Pytest configuration and shared fixtures."""

import pytest
from pathlib import Path
import tempfile


@pytest.fixture
def temp_dir(tmp_path: Path) -> Path:
    """Create a temporary directory for testing."""
    return tmp_path


@pytest.fixture
def temp_repo(tmp_path: Path) -> Path:
    """Create a temporary test repository with sample VASP files."""
    repo = tmp_path / "test_repo"
    repo.mkdir()
    
    # Create sample INCAR
    (repo / "INCAR").write_text(
        "ENCUT = 520\n"
        "ISMEAR = 0\n"
        "SIGMA = 0.05\n"
        "LDAUU = 3.5 0 0\n"
    )
    
    # Create sample POSCAR (simplified)
    (repo / "POSCAR").write_text(
        "Li4 Co4 O8\n"
        "1.0\n"
        "2.830000   0.000000   0.000000\n"
        "-1.415000   2.450794   0.000000\n"
        "0.000000   0.000000  14.050000\n"
        "Li Co O\n"
        "4 4 8\n"
        "direct\n"
        "0.000000  0.000000  0.000000 Li\n"
        "0.500000  0.500000  0.500000 Li\n"
        "0.250000  0.250000  0.250000 Li\n"
        "0.750000  0.750000  0.750000 Li\n"
        "0.000000  0.500000  0.500000 Co\n"
        "0.500000  0.000000  0.000000 Co\n"
        "0.250000  0.750000  0.750000 Co\n"
        "0.750000  0.250000  0.250000 Co\n"
        "0.125000  0.125000  0.125000 O\n"
        "0.375000  0.375000  0.375000 O\n"
        "0.625000  0.625000  0.625000 O\n"
        "0.875000  0.875000  0.875000 O\n"
        "0.125000  0.625000  0.625000 O\n"
        "0.375000  0.875000  0.875000 O\n"
        "0.625000  0.125000  0.125000 O\n"
        "0.875000  0.375000  0.375000 O\n"
    )
    
    # Create sample KPOINTS
    (repo / "KPOINTS").write_text(
        "Automatic\n"
        "0\n"
        "Gamma\n"
        "4 4 4\n"
    )
    
    return repo
