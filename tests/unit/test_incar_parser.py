"""Tests for INCAR parser."""

import pytest

from chemvcs.parsers.incar_parser import IncarParser
from chemvcs.parsers.base_parser import ParserError


class TestIncarParser:
    """Test INCAR file parser."""
    
    def test_parse_simple_incar(self) -> None:
        """Test parsing simple INCAR file."""
        content = """
ENCUT = 520
PREC = Accurate
ISMEAR = 0
SIGMA = 0.05
"""
        parser = IncarParser()
        data = parser.parse(content)
        
        assert data["ENCUT"] == 520
        assert data["PREC"] == "Accurate"
        assert data["ISMEAR"] == 0
        assert data["SIGMA"] == 0.05
    
    def test_parse_with_comments(self) -> None:
        """Test parsing with comments."""
        content = """
# Energy cutoff
ENCUT = 520  # in eV
PREC = Accurate  ! precision
! This is a comment line
ISMEAR = 0
"""
        parser = IncarParser()
        data = parser.parse(content)
        
        assert data["ENCUT"] == 520
        assert data["PREC"] == "Accurate"
        assert data["ISMEAR"] == 0
    
    def test_parse_semicolon_separated(self) -> None:
        """Test parsing semicolon-separated assignments."""
        content = "ISMEAR = 0; SIGMA = 0.05; NSW = 100"
        parser = IncarParser()
        data = parser.parse(content)
        
        assert data["ISMEAR"] == 0
        assert data["SIGMA"] == 0.05
        assert data["NSW"] == 100
    
    def test_parse_boolean_values(self) -> None:
        """Test parsing boolean-like values."""
        content = """
LDAU = .TRUE.
LWAVE = .FALSE.
LCHARG = T
LORBIT = F
"""
        parser = IncarParser()
        data = parser.parse(content)
        
        assert data["LDAU"] is True
        assert data["LWAVE"] is False
        assert data["LCHARG"] is True
        assert data["LORBIT"] is False
    
    def test_parse_multi_value(self) -> None:
        """Test parsing multi-value parameters."""
        content = """
MAGMOM = 5.0 -5.0 5.0 -5.0
LDAUL = 2 -1 -1
LDAUU = 4.0 0.0 0.0
"""
        parser = IncarParser()
        data = parser.parse(content)
        
        assert data["MAGMOM"] == [5.0, -5.0, 5.0, -5.0]
        assert data["LDAUL"] == [2, -1, -1]
        assert data["LDAUU"] == [4.0, 0.0, 0.0]
    
    def test_parse_case_insensitive_keys(self) -> None:
        """Test that keys are normalized to uppercase."""
        content = """
encut = 520
Prec = Accurate
IsMear = 0
"""
        parser = IncarParser()
        data = parser.parse(content)
        
        assert "ENCUT" in data
        assert "PREC" in data
        assert "ISMEAR" in data
    
    def test_parse_empty_lines(self) -> None:
        """Test parsing with empty lines."""
        content = """

ENCUT = 520

PREC = Accurate

"""
        parser = IncarParser()
        data = parser.parse(content)
        
        assert len(data) == 2
        assert data["ENCUT"] == 520
    
    def test_diff_no_changes(self) -> None:
        """Test diff with no changes."""
        content = "ENCUT = 520\nPREC = Accurate"
        parser = IncarParser()
        
        data1 = parser.parse(content)
        data2 = parser.parse(content)
        
        diff = parser.diff(data1, data2)
        assert len(diff) == 0
    
    def test_diff_added_parameter(self) -> None:
        """Test diff with added parameter."""
        parser = IncarParser()
        
        old = parser.parse("ENCUT = 520")
        new = parser.parse("ENCUT = 520\nPREC = Accurate")
        
        diff = parser.diff(old, new)
        
        assert len(diff) == 1
        assert diff[0].path == "PREC"
        assert diff[0].change_type == "added"
        assert diff[0].new_value == "Accurate"
    
    def test_diff_deleted_parameter(self) -> None:
        """Test diff with deleted parameter."""
        parser = IncarParser()
        
        old = parser.parse("ENCUT = 520\nPREC = Accurate")
        new = parser.parse("ENCUT = 520")
        
        diff = parser.diff(old, new)
        
        assert len(diff) == 1
        assert diff[0].path == "PREC"
        assert diff[0].change_type == "deleted"
        assert diff[0].old_value == "Accurate"
    
    def test_diff_modified_parameter(self) -> None:
        """Test diff with modified parameter."""
        parser = IncarParser()
        
        old = parser.parse("ENCUT = 520")
        new = parser.parse("ENCUT = 600")
        
        diff = parser.diff(old, new)
        
        assert len(diff) == 1
        assert diff[0].path == "ENCUT"
        assert diff[0].change_type == "modified"
        assert diff[0].old_value == 520
        assert diff[0].new_value == 600
    
    def test_diff_significance_critical(self) -> None:
        """Test that critical parameters are marked correctly."""
        parser = IncarParser()
        
        old = parser.parse("ENCUT = 520")
        new = parser.parse("ENCUT = 600")
        
        diff = parser.diff(old, new)
        
        assert diff[0].significance == "critical"
    
    def test_diff_significance_major(self) -> None:
        """Test that major parameters are marked correctly."""
        parser = IncarParser()
        
        old = parser.parse("LWAVE = .TRUE.")
        new = parser.parse("LWAVE = .FALSE.")
        
        diff = parser.diff(old, new)
        
        assert diff[0].significance == "major"
    
    def test_diff_significance_minor(self) -> None:
        """Test that minor parameters are marked correctly."""
        parser = IncarParser()
        
        old = parser.parse("SYSTEM = test1")
        new = parser.parse("SYSTEM = test2")
        
        diff = parser.diff(old, new)
        
        assert diff[0].significance == "minor"
    
    def test_format_diff_default(self) -> None:
        """Test default diff formatting."""
        parser = IncarParser()
        
        old = parser.parse("ENCUT = 520")
        new = parser.parse("ENCUT = 600\nPREC = Accurate")
        
        diff = parser.diff(old, new)
        formatted = parser.format_diff(diff)
        
        assert "ENCUT" in formatted
        assert "520" in formatted
        assert "600" in formatted
        assert "PREC" in formatted
    
    def test_format_diff_compact(self) -> None:
        """Test compact diff formatting."""
        parser = IncarParser()
        
        old = parser.parse("ENCUT = 520")
        new = parser.parse("ENCUT = 600")
        
        diff = parser.diff(old, new)
        formatted = parser.format_diff(diff, style="compact")
        
        assert "~" in formatted
        assert "ENCUT" in formatted
    
    def test_validate_valid_incar(self) -> None:
        """Test validation of valid INCAR."""
        parser = IncarParser()
        data = parser.parse("ENCUT = 520\nISMEAR = 0\nSIGMA = 0.05")
        
        is_valid, errors = parser.validate(data)
        
        assert is_valid
        assert len(errors) == 0
    
    def test_validate_negative_encut(self) -> None:
        """Test validation catches negative ENCUT."""
        parser = IncarParser()
        data = {"ENCUT": -520}
        
        is_valid, errors = parser.validate(data)
        
        assert not is_valid
        assert any("ENCUT" in err for err in errors)
    
    def test_validate_invalid_ismear(self) -> None:
        """Test validation catches invalid ISMEAR."""
        parser = IncarParser()
        data = {"ISMEAR": 10}
        
        is_valid, errors = parser.validate(data)
        
        assert not is_valid
        assert any("ISMEAR" in err for err in errors)
    
    def test_validate_negative_sigma(self) -> None:
        """Test validation catches negative SIGMA."""
        parser = IncarParser()
        data = {"SIGMA": -0.05}
        
        is_valid, errors = parser.validate(data)
        
        assert not is_valid
        assert any("SIGMA" in err for err in errors)
    
    def test_validate_large_sigma_warning(self) -> None:
        """Test validation warns about large SIGMA with ISMEAR=0."""
        parser = IncarParser()
        data = {"ISMEAR": 0, "SIGMA": 0.5}
        
        is_valid, errors = parser.validate(data)
        
        assert not is_valid
        assert any("SIGMA" in err for err in errors)
    
    def test_parse_scientific_notation(self) -> None:
        """Test parsing scientific notation."""
        content = "EDIFF = 1e-6\nEDIFFG = -1.5e-2"
        parser = IncarParser()
        data = parser.parse(content)
        
        assert data["EDIFF"] == 1e-6
        assert data["EDIFFG"] == -1.5e-2
    
    def test_parse_string_values(self) -> None:
        """Test parsing string values."""
        content = """
SYSTEM = Fe bulk calculation
PREC = Accurate
ALGO = Fast
"""
        parser = IncarParser()
        data = parser.parse(content)
        
        assert data["SYSTEM"] == "Fe bulk calculation"
        assert data["PREC"] == "Accurate"
        assert data["ALGO"] == "Fast"
