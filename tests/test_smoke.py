"""Basic smoke tests to verify project setup."""

import pytest
from chemvcs import __version__


def test_version() -> None:
    """Test that version is correctly defined."""
    assert __version__ == "0.1.0"


def test_import_storage() -> None:
    """Test that storage module can be imported."""
    from chemvcs import storage  # noqa: F401


def test_import_parsers() -> None:
    """Test that parsers module can be imported."""
    from chemvcs import parsers  # noqa: F401


def test_import_cli() -> None:
    """Test that cli module can be imported."""
    from chemvcs.cli import main  # noqa: F401


def test_temp_repo_fixture(temp_repo) -> None:
    """Test that temp_repo fixture creates VASP files."""
    assert (temp_repo / "INCAR").exists()
    assert (temp_repo / "POSCAR").exists()
    assert (temp_repo / "KPOINTS").exists()
    
    # Check INCAR content
    incar_content = (temp_repo / "INCAR").read_text()
    assert "ENCUT" in incar_content
    assert "520" in incar_content
