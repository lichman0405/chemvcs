"""Tests for XYZ file parser."""

import pytest
import numpy as np
from pathlib import Path
import tempfile

from chemvcs_py.io.xyz import read_xyz, write_xyz
from chemvcs_py.domain.structure import Structure
from chemvcs_py.util.errors import ParseError


class TestXYZParser:
    """Test cases for XYZ file I/O."""
    
    @pytest.fixture
    def water_xyz_content(self):
        """Sample water molecule XYZ file."""
        return """3
Water molecule
O      0.000000    0.000000    0.000000
H      0.000000    0.757000    0.586000
H      0.000000   -0.757000    0.586000
"""
    
    @pytest.fixture
    def water_xyz_file(self, water_xyz_content):
        """Create temporary XYZ file."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.xyz', delete=False) as f:
            f.write(water_xyz_content)
            temp_path = f.name
        yield temp_path
        Path(temp_path).unlink()
    
    def test_read_xyz_water(self, water_xyz_file):
        """Test reading water molecule from XYZ file."""
        struct = read_xyz(water_xyz_file)
        
        assert struct.formula == "H2O"
        assert struct.num_atoms == 3
        assert struct.species == ["O", "H", "H"]
        
        # Check positions
        expected_positions = np.array([
            [0.0, 0.0, 0.0],
            [0.0, 0.757, 0.586],
            [0.0, -0.757, 0.586]
        ])
        np.testing.assert_array_almost_equal(struct.positions, expected_positions)
        
        # Check metadata
        assert struct.metadata["title"] == "Water molecule"
        assert "source_file" in struct.metadata
    
    def test_read_xyz_file_not_found(self):
        """Test error handling for missing file."""
        with pytest.raises(ParseError, match="File not found"):
            read_xyz("nonexistent.xyz")
    
    def test_read_xyz_invalid_atom_count(self):
        """Test error handling for invalid atom count."""
        content = """abc
Test
H 0 0 0
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.xyz', delete=False) as f:
            f.write(content)
            temp_path = f.name
        
        try:
            with pytest.raises(ParseError, match="Invalid atom count"):
                read_xyz(temp_path)
        finally:
            Path(temp_path).unlink()
    
    def test_read_xyz_too_few_lines(self):
        """Test error handling for truncated file."""
        content = """3
Test molecule
H 0 0 0
"""  # Only 1 atom, but header says 3
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.xyz', delete=False) as f:
            f.write(content)
            temp_path = f.name
        
        try:
            with pytest.raises(ParseError, match="Expected 3 atom lines"):
                read_xyz(temp_path)
        finally:
            Path(temp_path).unlink()
    
    def test_read_xyz_invalid_coordinates(self):
        """Test error handling for non-numeric coordinates."""
        content = """1
Test
H abc def ghi
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.xyz', delete=False) as f:
            f.write(content)
            temp_path = f.name
        
        try:
            with pytest.raises(ParseError, match="Invalid coordinates"):
                read_xyz(temp_path)
        finally:
            Path(temp_path).unlink()
    
    def test_write_xyz(self):
        """Test writing structure to XYZ file."""
        struct = Structure(
            formula="H2O",
            positions=np.array([
                [0.0, 0.0, 0.0],
                [0.0, 0.757, 0.586],
                [0.0, -0.757, 0.586]
            ]),
            species=["O", "H", "H"]
        )
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.xyz', delete=False) as f:
            temp_path = f.name
        
        try:
            write_xyz(struct, temp_path)
            
            # Read it back
            reconstructed = read_xyz(temp_path)
            
            assert reconstructed.formula == struct.formula
            assert reconstructed.num_atoms == struct.num_atoms
            assert reconstructed.species == struct.species
            np.testing.assert_array_almost_equal(reconstructed.positions, struct.positions)
        finally:
            Path(temp_path).unlink()
    
    def test_write_xyz_with_comment(self):
        """Test writing XYZ with custom comment."""
        struct = Structure(
            formula="H2",
            positions=np.array([[0, 0, 0], [0, 0, 1]]),
            species=["H", "H"]
        )
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.xyz', delete=False) as f:
            temp_path = f.name
        
        try:
            write_xyz(struct, temp_path, comment="Custom comment")
            
            # Check comment line
            with open(temp_path, 'r') as f:
                lines = f.readlines()
            assert lines[1].strip() == "Custom comment"
        finally:
            Path(temp_path).unlink()
    
    def test_round_trip(self):
        """Test that write->read is lossless."""
        original = Structure(
            formula="CH4",
            positions=np.array([
                [0.0, 0.0, 0.0],
                [0.629, 0.629, 0.629],
                [-0.629, -0.629, 0.629],
                [-0.629, 0.629, -0.629],
                [0.629, -0.629, -0.629]
            ]),
            species=["C", "H", "H", "H", "H"],
            metadata={"method": "test"}
        )
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.xyz', delete=False) as f:
            temp_path = f.name
        
        try:
            write_xyz(original, temp_path, comment="Methane")
            reconstructed = read_xyz(temp_path)
            
            assert reconstructed.formula == original.formula
            assert reconstructed.num_atoms == original.num_atoms
            assert reconstructed.species == original.species
            np.testing.assert_array_almost_equal(reconstructed.positions, original.positions, decimal=5)
        finally:
            Path(temp_path).unlink()
