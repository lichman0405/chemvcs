"""Tests for KPOINTS parser."""

import pytest

from chemvcs.parsers.kpoints_parser import KpointsParser
from chemvcs.parsers.base_parser import ParserError


class TestKpointsParser:
    """Test KPOINTS file parser."""
    
    def test_parse_automatic_gamma(self) -> None:
        """Test parsing automatic Gamma-centered mesh."""
        content = """Automatic mesh
0
Gamma
4 4 4
"""
        parser = KpointsParser()
        data = parser.parse(content)
        
        assert data["type"] == "automatic"
        assert data["gamma_centered"] is True
        assert data["grid"] == [4, 4, 4]
        assert data["shift"] == [0.0, 0.0, 0.0]
    
    def test_parse_automatic_monkhorst(self) -> None:
        """Test parsing automatic Monkhorst-Pack mesh."""
        content = """Automatic mesh
0
Monkhorst-Pack
8 8 8
"""
        parser = KpointsParser()
        data = parser.parse(content)
        
        assert data["type"] == "automatic"
        assert data["gamma_centered"] is False
        assert data["grid"] == [8, 8, 8]
    
    def test_parse_automatic_with_shift(self) -> None:
        """Test parsing automatic mesh with shift."""
        content = """Automatic mesh
0
Gamma
4 4 4
0.5 0.5 0.5
"""
        parser = KpointsParser()
        data = parser.parse(content)
        
        assert data["shift"] == [0.5, 0.5, 0.5]
    
    def test_parse_explicit_kpoints(self) -> None:
        """Test parsing explicit k-points."""
        content = """Explicit k-points
4
Reciprocal
0.0 0.0 0.0  1.0
0.5 0.0 0.0  1.0
0.5 0.5 0.0  1.0
0.5 0.5 0.5  1.0
"""
        parser = KpointsParser()
        data = parser.parse(content)
        
        assert data["type"] == "explicit"
        assert data["num_kpoints"] == 4
        assert data["cartesian"] is False
        assert len(data["kpoints"]) == 4
        assert data["kpoints"][0]["coords"] == [0.0, 0.0, 0.0]
        assert data["kpoints"][0]["weight"] == 1.0
    
    def test_parse_explicit_cartesian(self) -> None:
        """Test parsing explicit Cartesian k-points."""
        content = """Explicit k-points
2
Cartesian
0.0 0.0 0.0  0.5
0.1 0.1 0.1  0.5
"""
        parser = KpointsParser()
        data = parser.parse(content)
        
        assert data["cartesian"] is True
        assert data["kpoints"][1]["coords"] == [0.1, 0.1, 0.1]
        assert data["kpoints"][1]["weight"] == 0.5
    
    def test_parse_line_mode(self) -> None:
        """Test parsing line-mode for band structure."""
        content = """k-points for band structure
40
Line-mode
Reciprocal
0.0 0.0 0.0  ! Gamma
0.5 0.0 0.0  ! X

0.5 0.0 0.0  ! X
0.5 0.5 0.0  ! M
"""
        parser = KpointsParser()
        data = parser.parse(content)
        
        assert data["type"] == "line"
        assert data["num_points"] == 40
        assert data["reciprocal"] is True
        assert len(data["segments"]) == 2
        assert data["segments"][0]["start"]["label"] == "Gamma"
        assert data["segments"][0]["end"]["label"] == "X"
    
    def test_parse_comment_preserved(self) -> None:
        """Test that comment is preserved."""
        content = """My calculation
0
Gamma
4 4 4
"""
        parser = KpointsParser()
        data = parser.parse(content)
        
        assert data["comment"] == "My calculation"
    
    def test_parse_too_short_fails(self) -> None:
        """Test that file with too few lines fails."""
        content = """Too short
0
"""
        parser = KpointsParser()
        
        with pytest.raises(ParserError):
            parser.parse(content)
    
    def test_parse_invalid_grid_fails(self) -> None:
        """Test that invalid grid specification fails."""
        content = """Test
0
Gamma
not a number
"""
        parser = KpointsParser()
        
        with pytest.raises(ParserError):
            parser.parse(content)
    
    def test_diff_no_changes(self) -> None:
        """Test diff with no changes."""
        content = """Test
0
Gamma
4 4 4
"""
        parser = KpointsParser()
        
        data1 = parser.parse(content)
        data2 = parser.parse(content)
        
        diff = parser.diff(data1, data2)
        assert len(diff) == 0
    
    def test_diff_type_change(self) -> None:
        """Test diff when k-point type changes."""
        parser = KpointsParser()
        
        old = parser.parse("A\n0\nGamma\n4 4 4")
        new = parser.parse("A\n4\nReciprocal\n0 0 0 1\n0.5 0 0 1\n0.5 0.5 0 1\n0.5 0.5 0.5 1")
        
        diff = parser.diff(old, new)
        
        assert len(diff) == 1
        assert diff[0].path == "type"
        assert diff[0].change_type == "modified"
        assert diff[0].significance == "critical"
    
    def test_diff_grid_change(self) -> None:
        """Test diff when grid changes."""
        parser = KpointsParser()
        
        old = parser.parse("A\n0\nGamma\n4 4 4")
        new = parser.parse("A\n0\nGamma\n8 8 8")
        
        diff = parser.diff(old, new)
        
        assert len(diff) == 1
        assert diff[0].path == "grid"
        assert diff[0].old_value == [4, 4, 4]
        assert diff[0].new_value == [8, 8, 8]
        assert diff[0].significance == "critical"
    
    def test_diff_gamma_to_monkhorst(self) -> None:
        """Test diff when switching Gamma to Monkhorst-Pack."""
        parser = KpointsParser()
        
        old = parser.parse("A\n0\nGamma\n4 4 4")
        new = parser.parse("A\n0\nMonkhorst-Pack\n4 4 4")
        
        diff = parser.diff(old, new)
        
        assert len(diff) == 1
        assert diff[0].path == "gamma_centered"
        assert diff[0].significance == "major"
    
    def test_diff_shift_change(self) -> None:
        """Test diff when shift changes."""
        parser = KpointsParser()
        
        old = parser.parse("A\n0\nGamma\n4 4 4\n0.0 0.0 0.0")
        new = parser.parse("A\n0\nGamma\n4 4 4\n0.5 0.5 0.5")
        
        diff = parser.diff(old, new)
        
        # Should have shift change
        shift_diff = [d for d in diff if d.path == "shift"]
        assert len(shift_diff) == 1
        assert shift_diff[0].significance == "major"
    
    def test_diff_explicit_num_kpoints(self) -> None:
        """Test diff when number of explicit k-points changes."""
        parser = KpointsParser()
        
        old = parser.parse("A\n2\nReciprocal\n0 0 0 1\n0.5 0 0 1")
        new = parser.parse("A\n3\nReciprocal\n0 0 0 1\n0.5 0 0 1\n0.5 0.5 0 1")
        
        diff = parser.diff(old, new)
        
        num_diff = [d for d in diff if d.path == "num_kpoints"]
        assert len(num_diff) == 1
        assert num_diff[0].significance == "critical"
    
    def test_parse_case_insensitive_style(self) -> None:
        """Test that generation style is case-insensitive."""
        content1 = "A\n0\ngamma\n4 4 4"
        content2 = "A\n0\nGAMMA\n4 4 4"
        content3 = "A\n0\nGamma\n4 4 4"
        
        parser = KpointsParser()
        
        data1 = parser.parse(content1)
        data2 = parser.parse(content2)
        data3 = parser.parse(content3)
        
        assert data1["gamma_centered"] is True
        assert data2["gamma_centered"] is True
        assert data3["gamma_centered"] is True
    
    def test_parse_abbreviated_style(self) -> None:
        """Test parsing abbreviated style (G for Gamma, M for Monkhorst)."""
        parser = KpointsParser()
        
        gamma = parser.parse("A\n0\nG\n4 4 4")
        monkhorst = parser.parse("A\n0\nM\n4 4 4")
        
        assert gamma["gamma_centered"] is True
        assert monkhorst["gamma_centered"] is False
