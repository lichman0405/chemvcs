"""Tests for semantic diff engine."""

import pytest

from chemvcs.parsers.diff_engine import DiffEngine


class TestDiffEngine:
    """Test semantic diff engine."""
    
    def test_get_file_type_incar(self) -> None:
        """Test file type detection for INCAR."""
        engine = DiffEngine()
        
        assert engine.get_file_type("INCAR") == "INCAR"
        assert engine.get_file_type("incar") == "INCAR"
        assert engine.get_file_type("path/to/INCAR") == "INCAR"
    
    def test_get_file_type_kpoints(self) -> None:
        """Test file type detection for KPOINTS."""
        engine = DiffEngine()
        
        assert engine.get_file_type("KPOINTS") == "KPOINTS"
        assert engine.get_file_type("kpoints") == "KPOINTS"
    
    def test_get_file_type_unknown(self) -> None:
        """Test file type detection for unknown files."""
        engine = DiffEngine()
        
        assert engine.get_file_type("unknown.txt") is None
        assert engine.get_file_type("README.md") is None
    
    def test_can_parse_supported(self) -> None:
        """Test can_parse for supported files."""
        engine = DiffEngine()
        
        assert engine.can_parse("INCAR") is True
        assert engine.can_parse("KPOINTS") is True
    
    def test_can_parse_outcar(self) -> None:
        """Test can_parse for OUTCAR."""
        engine = DiffEngine()

        assert engine.can_parse("OUTCAR") is True
        assert engine.get_file_type("OUTCAR") == "OUTCAR"
        assert engine.get_file_type("path/to/OUTCAR") == "OUTCAR"

    def test_can_parse_unsupported(self) -> None:
        """Test can_parse for unsupported files."""
        engine = DiffEngine()

        assert engine.can_parse("unknown.txt") is False
        assert engine.can_parse("POSCAR") is False  # Not yet implemented
    
    def test_diff_files_incar(self) -> None:
        """Test semantic diff for INCAR files."""
        engine = DiffEngine()
        
        old = "ENCUT = 520\nPREC = Accurate"
        new = "ENCUT = 600\nPREC = Accurate"
        
        diff = engine.diff_files(old, new, "INCAR")
        
        assert diff is not None
        assert len(diff) == 1
        assert diff[0].path == "ENCUT"
        assert diff[0].old_value == 520
        assert diff[0].new_value == 600
    
    def test_diff_files_kpoints(self) -> None:
        """Test semantic diff for KPOINTS files."""
        engine = DiffEngine()
        
        old = "Test\n0\nGamma\n4 4 4"
        new = "Test\n0\nGamma\n8 8 8"
        
        diff = engine.diff_files(old, new, "KPOINTS")
        
        assert diff is not None
        assert len(diff) == 1
        assert diff[0].path == "grid"
    
    def test_diff_files_unsupported(self) -> None:
        """Test diff for unsupported file types."""
        engine = DiffEngine()
        
        diff = engine.diff_files("old", "new", "unknown.txt")
        
        assert diff is None
    
    def test_diff_files_no_changes(self) -> None:
        """Test diff with no changes."""
        engine = DiffEngine()
        
        content = "ENCUT = 520"
        diff = engine.diff_files(content, content, "INCAR")
        
        assert diff is not None
        assert len(diff) == 0
    
    def test_format_diff(self) -> None:
        """Test diff formatting."""
        engine = DiffEngine()
        
        old = "ENCUT = 520"
        new = "ENCUT = 600"
        
        diff = engine.diff_files(old, new, "INCAR")
        formatted = engine.format_diff(diff)
        
        assert "ENCUT" in formatted
        assert "520" in formatted
        assert "600" in formatted
    
    def test_summarize_diff_no_changes(self) -> None:
        """Test diff summary with no changes."""
        engine = DiffEngine()
        
        summary = engine.summarize_diff([])
        
        assert summary["total_changes"] == 0
        assert summary["by_type"] == {}
        assert summary["by_significance"] == {}
    
    def test_summarize_diff_with_changes(self) -> None:
        """Test diff summary with changes."""
        engine = DiffEngine()
        
        old = "ENCUT = 520"
        new = "ENCUT = 600\nPREC = Accurate"
        
        diff = engine.diff_files(old, new, "INCAR")
        summary = engine.summarize_diff(diff)
        
        assert summary["total_changes"] == 2
        assert summary["by_type"]["modified"] == 1
        assert summary["by_type"]["added"] == 1
        assert summary["has_critical_changes"] is True  # ENCUT is critical
    
    def test_validate_file_valid_incar(self) -> None:
        """Test validation of valid INCAR."""
        engine = DiffEngine()
        
        content = "ENCUT = 520\nISMEAR = 0\nSIGMA = 0.05"
        is_valid, errors = engine.validate_file(content, "INCAR")
        
        assert is_valid is True
        assert len(errors) == 0
    
    def test_validate_file_invalid_incar(self) -> None:
        """Test validation of invalid INCAR."""
        engine = DiffEngine()
        
        content = "ENCUT = -520"
        is_valid, errors = engine.validate_file(content, "INCAR")
        
        assert is_valid is False
        assert len(errors) > 0
    
    def test_validate_file_unsupported(self) -> None:
        """Test validation of unsupported file."""
        engine = DiffEngine()
        
        is_valid, errors = engine.validate_file("content", "unknown.txt")
        
        assert is_valid is True  # Unknown types pass validation
        assert len(errors) == 0
    
    def test_diff_files_parse_error_fallback(self) -> None:
        """Test that parse errors return None (fallback to text diff)."""
        engine = DiffEngine()
        
        # Malformed KPOINTS
        malformed = "Too short"
        valid = "Test\n0\nGamma\n4 4 4"
        
        diff = engine.diff_files(malformed, valid, "KPOINTS")
        
        # Should return None to indicate fallback needed
        assert diff is None
