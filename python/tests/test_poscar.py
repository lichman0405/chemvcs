"""Tests for POSCAR file parser."""

import pytest
import numpy as np
from pathlib import Path
import tempfile

from chemvcs_py.io.poscar import read_poscar, write_poscar
from chemvcs_py.domain.structure import Structure
from chemvcs_py.util.errors import ParseError


class TestPOSCARParser:
    """Test cases for POSCAR file I/O."""
    
    @pytest.fixture
    def nacl_poscar_content(self):
        """Sample NaCl POSCAR in direct coordinates."""
        return """NaCl rocksalt structure
5.64
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
Na Cl
1 1
Direct
0.0 0.0 0.0
0.5 0.5 0.5
"""
    
    @pytest.fixture
    def nacl_poscar_file(self, nacl_poscar_content):
        """Create temporary POSCAR file."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='', delete=False, 
                                        prefix='POSCAR_') as f:
            f.write(nacl_poscar_content)
            temp_path = f.name
        yield temp_path
        Path(temp_path).unlink()
    
    def test_read_poscar_direct(self, nacl_poscar_file):
        """Test reading POSCAR with direct coordinates."""
        struct = read_poscar(nacl_poscar_file)
        
        assert struct.formula == "NaCl"
        assert struct.num_atoms == 2
        assert struct.species == ["Na", "Cl"]
        assert struct.is_periodic
        
        # Check lattice (5.64 * identity)
        expected_lattice = np.eye(3) * 5.64
        np.testing.assert_array_almost_equal(struct.lattice, expected_lattice)
        
        # Check positions (converted to Cartesian)
        expected_positions = np.array([
            [0.0, 0.0, 0.0],
            [2.82, 2.82, 2.82]  # 0.5 * 5.64
        ])
        np.testing.assert_array_almost_equal(struct.positions, expected_positions)
        
        # Check metadata
        assert struct.metadata["title"] == "NaCl rocksalt structure"
        assert struct.metadata["coordinate_mode"] == "direct"
    
    @pytest.fixture
    def cartesian_poscar_content(self):
        """Sample POSCAR in Cartesian coordinates."""
        return """H2 molecule
1.0
10.0 0.0 0.0
0.0 10.0 0.0
0.0 0.0 10.0
H
2
Cartesian
0.0 0.0 0.0
0.0 0.0 0.74
"""
    
    def test_read_poscar_cartesian(self, cartesian_poscar_content):
        """Test reading POSCAR with Cartesian coordinates."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='', delete=False) as f:
            f.write(cartesian_poscar_content)
            temp_path = f.name
        
        try:
            struct = read_poscar(temp_path)
            
            assert struct.formula == "H2"
            assert struct.num_atoms == 2
            assert struct.species == ["H", "H"]
            
            # Positions should be as-is (Cartesian)
            expected = np.array([[0, 0, 0], [0, 0, 0.74]])
            np.testing.assert_array_almost_equal(struct.positions, expected)
            
            assert struct.metadata["coordinate_mode"] == "cartesian"
        finally:
            Path(temp_path).unlink()
    
    def test_read_poscar_file_not_found(self):
        """Test error handling for missing file."""
        with pytest.raises(ParseError, match="File not found"):
            read_poscar("nonexistent_POSCAR")
    
    def test_read_poscar_too_few_lines(self):
        """Test error handling for truncated file."""
        content = """Title
1.0
1 0 0
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='', delete=False) as f:
            f.write(content)
            temp_path = f.name
        
        try:
            with pytest.raises(ParseError, match="at least 9 lines"):
                read_poscar(temp_path)
        finally:
            Path(temp_path).unlink()
    
    def test_read_poscar_invalid_scaling(self):
        """Test error handling for invalid scaling factor."""
        content = """Title
abc
1 0 0
0 1 0
0 0 1
H
1
Direct
0 0 0
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='', delete=False) as f:
            f.write(content)
            temp_path = f.name
        
        try:
            with pytest.raises(ParseError, match="Invalid scaling factor"):
                read_poscar(temp_path)
        finally:
            Path(temp_path).unlink()
    
    def test_read_poscar_element_count_mismatch(self):
        """Test error handling for element/count mismatch."""
        content = """Title
1.0
1 0 0
0 1 0
0 0 1
Na Cl
1
Direct
0 0 0
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='', delete=False) as f:
            f.write(content)
            temp_path = f.name
        
        try:
            with pytest.raises(ParseError, match="doesn't match"):
                read_poscar(temp_path)
        finally:
            Path(temp_path).unlink()
    
    def test_write_poscar_direct(self):
        """Test writing POSCAR with direct coordinates."""
        struct = Structure(
            formula="NaCl",
            positions=np.array([
                [0.0, 0.0, 0.0],
                [2.82, 2.82, 2.82]
            ]),
            species=["Na", "Cl"],
            lattice=np.array([
                [5.64, 0, 0],
                [0, 5.64, 0],
                [0, 0, 5.64]
            ])
        )
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='', delete=False) as f:
            temp_path = f.name
        
        try:
            write_poscar(struct, temp_path, direct=True)
            
            # Read it back
            reconstructed = read_poscar(temp_path)
            
            assert reconstructed.formula == struct.formula
            assert reconstructed.num_atoms == struct.num_atoms
            assert reconstructed.species == struct.species
            np.testing.assert_array_almost_equal(reconstructed.lattice, struct.lattice)
            np.testing.assert_array_almost_equal(reconstructed.positions, struct.positions, decimal=5)
        finally:
            Path(temp_path).unlink()
    
    def test_write_poscar_cartesian(self):
        """Test writing POSCAR with Cartesian coordinates."""
        struct = Structure(
            formula="H2",
            positions=np.array([[0, 0, 0], [0, 0, 0.74]]),
            species=["H", "H"],
            lattice=np.eye(3) * 10
        )
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='', delete=False) as f:
            temp_path = f.name
        
        try:
            write_poscar(struct, temp_path, direct=False, comment="H2 test")
            
            with open(temp_path, 'r') as f:
                lines = f.readlines()
            
            assert lines[0].strip() == "H2 test"
            assert "Cartesian" in lines[7]
        finally:
            Path(temp_path).unlink()
    
    def test_write_poscar_requires_lattice(self):
        """Test that writing requires periodic structure."""
        struct = Structure(
            formula="H2",
            positions=np.array([[0, 0, 0], [0, 0, 1]]),
            species=["H", "H"]
        )
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='', delete=False) as f:
            temp_path = f.name
        
        try:
            with pytest.raises(ValueError, match="periodic structure"):
                write_poscar(struct, temp_path)
        finally:
            Path(temp_path).unlink()
    
    def test_round_trip_direct(self):
        """Test that write->read with direct coords is lossless."""
        original = Structure(
            formula="Si2",
            positions=np.array([
                [0.0, 0.0, 0.0],
                [1.357, 1.357, 1.357]
            ]),
            species=["Si", "Si"],
            lattice=np.array([
                [5.43, 0, 0],
                [0, 5.43, 0],
                [0, 0, 5.43]
            ]),
            metadata={"source": "test"}
        )
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='', delete=False) as f:
            temp_path = f.name
        
        try:
            write_poscar(original, temp_path, direct=True)
            reconstructed = read_poscar(temp_path)
            
            assert reconstructed.formula == original.formula
            assert reconstructed.species == original.species
            np.testing.assert_array_almost_equal(reconstructed.lattice, original.lattice, decimal=8)
            np.testing.assert_array_almost_equal(reconstructed.positions, original.positions, decimal=8)
        finally:
            Path(temp_path).unlink()
