"""Tests for Structure domain object."""

import pytest
import numpy as np

from chemvcs_py.domain.structure import Structure
from chemvcs_py.core.objects import CoreObject


class TestStructure:
    """Test cases for Structure class."""
    
    def test_create_molecular_structure(self):
        """Test creating a molecular structure (no lattice)."""
        positions = np.array([
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0]
        ])
        species = ["H", "H"]
        
        struct = Structure(
            formula="H2",
            positions=positions,
            species=species
        )
        
        assert struct.formula == "H2"
        assert struct.num_atoms == 2
        assert not struct.is_periodic
        np.testing.assert_array_equal(struct.positions, positions)
        assert struct.species == species
        assert struct.lattice is None
    
    def test_create_periodic_structure(self):
        """Test creating a periodic structure (with lattice)."""
        positions = np.array([
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.5]
        ])
        species = ["Na", "Cl"]
        lattice = np.array([
            [5.64, 0.0, 0.0],
            [0.0, 5.64, 0.0],
            [0.0, 0.0, 5.64]
        ])
        
        struct = Structure(
            formula="NaCl",
            positions=positions,
            species=species,
            lattice=lattice
        )
        
        assert struct.formula == "NaCl"
        assert struct.num_atoms == 2
        assert struct.is_periodic
        np.testing.assert_array_equal(struct.lattice, lattice)
    
    def test_validation_position_shape(self):
        """Test that positions must be Nx3 array."""
        from chemvcs_py.util.errors import ValidationError
        with pytest.raises(ValidationError, match="must be N x 3"):
            Structure(
                formula="H2O",
                positions=np.array([[0, 0], [1, 1]]),  # Wrong shape
                species=["H", "H"]
            )
    
    def test_validation_species_length(self):
        """Test that species list must match number of atoms."""
        from chemvcs_py.util.errors import ValidationError
        positions = np.array([[0, 0, 0], [1, 1, 1]])
        
        with pytest.raises(ValidationError, match="Species length"):
            Structure(
                formula="H2O",
                positions=positions,
                species=["H"]  # Too few
            )
    
    def test_validation_lattice_shape(self):
        """Test that lattice must be 3x3 if provided."""
        from chemvcs_py.util.errors import ValidationError
        positions = np.array([[0, 0, 0]])
        species = ["H"]
        
        with pytest.raises(ValidationError, match="must be 3x3"):
            Structure(
                formula="H",
                positions=positions,
                species=species,
                lattice=np.array([[1, 0], [0, 1]])  # Wrong shape
            )
    
    def test_translate(self):
        """Test structure translation."""
        struct = Structure(
            formula="H2",
            positions=np.array([[0, 0, 0], [0, 0, 1]]),
            species=["H", "H"]
        )
        
        displacement = np.array([1.0, 2.0, 3.0])
        struct.translate(displacement)
        
        expected = np.array([[1, 2, 3], [1, 2, 4]], dtype=float)
        np.testing.assert_array_almost_equal(struct.positions, expected)
    
    def test_center_of_mass(self):
        """Test center of mass calculation."""
        # Water molecule
        struct = Structure(
            formula="H2O",
            positions=np.array([
                [0.0, 0.0, 0.0],    # O
                [0.0, 0.757, 0.586],  # H
                [0.0, -0.757, 0.586]  # H
            ]),
            species=["O", "H", "H"]
        )
        
        com = struct.get_center_of_mass()
        
        # For water, COM should be close to oxygen
        assert com[0] == pytest.approx(0.0)
        assert com[1] == pytest.approx(0.0, abs=0.1)
        assert com[2] > 0  # Shifted toward hydrogens
    
    def test_to_core_object(self):
        """Test conversion to CoreObject."""
        struct = Structure(
            formula="H2",
            positions=np.array([[0, 0, 0], [0, 0, 1]]),
            species=["H", "H"],
            metadata={"source": "test"}
        )
        
        core_obj = struct.to_core_object()
        
        assert isinstance(core_obj, CoreObject)
        assert core_obj.type == "structure"
        assert core_obj.meta["formula"] == "H2"
        assert core_obj.meta["num_atoms"] == 2
        assert core_obj.meta["is_periodic"] is False
        assert "positions" in core_obj.meta
        assert "species" in core_obj.meta
    
    def test_from_core_object(self):
        """Test reconstruction from CoreObject."""
        original = Structure(
            formula="NaCl",
            positions=np.array([[0, 0, 0], [0.5, 0.5, 0.5]]),
            species=["Na", "Cl"],
            lattice=np.array([[5, 0, 0], [0, 5, 0], [0, 0, 5]]),
            metadata={"tag": "test"}
        )
        
        core_obj = original.to_core_object()
        reconstructed = Structure.from_core_object(core_obj)
        
        assert reconstructed.formula == original.formula
        assert reconstructed.num_atoms == original.num_atoms
        np.testing.assert_array_almost_equal(reconstructed.positions, original.positions)
        assert reconstructed.species == original.species
        np.testing.assert_array_almost_equal(reconstructed.lattice, original.lattice)
        assert reconstructed.metadata["tag"] == "test"
    
    def test_round_trip_conversion(self):
        """Test that to/from CoreObject is lossless."""
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
            metadata={"method": "generated"}
        )
        
        core_obj = original.to_core_object()
        reconstructed = Structure.from_core_object(core_obj)
        
        # Check all attributes match
        assert reconstructed.formula == original.formula
        assert reconstructed.num_atoms == original.num_atoms
        assert reconstructed.species == original.species
        np.testing.assert_array_almost_equal(reconstructed.positions, original.positions)
        assert reconstructed.metadata == original.metadata
