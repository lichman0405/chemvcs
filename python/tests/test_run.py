"""Tests for Run domain object."""

import pytest
from datetime import datetime

from chemvcs_py.domain.run import Run
from chemvcs_py.core.objects import CoreObject
from chemvcs_py.util.errors import ValidationError


class TestRun:
    """Test cases for Run class."""
    
    def test_create_basic_run(self):
        """Test creating a basic run."""
        run = Run(
            structure_id="struct_001",
            code="ORCA",
            status="planned"
        )
        
        assert run.structure_id == "struct_001"
        assert run.code == "ORCA"
        assert run.status == "planned"
        assert run.parameters == {}
        assert run.results == {}
        assert run.resources == {}
    
    def test_create_run_with_parameters(self):
        """Test creating run with calculation parameters."""
        params = {
            "method": "DFT",
            "functional": "B3LYP",
            "basis": "def2-TZVP"
        }
        
        run = Run(
            structure_id="struct_001",
            code="ORCA",
            code_version="5.0.3",
            parameters=params,
            status="planned"
        )
        
        assert run.code_version == "5.0.3"
        assert run.parameters == params
    
    def test_validation_empty_structure_id(self):
        """Test that structure_id cannot be empty."""
        with pytest.raises(ValidationError, match="structure_id cannot be empty"):
            Run(
                structure_id="",
                code="ORCA",
                status="planned"
            )
    
    def test_validation_empty_code(self):
        """Test that code cannot be empty."""
        with pytest.raises(ValidationError, match="code cannot be empty"):
            Run(
                structure_id="struct_001",
                code="",
                status="planned"
            )
    
    def test_validation_invalid_status(self):
        """Test that invalid status is rejected."""
        with pytest.raises(ValidationError, match="Invalid status"):
            Run(
                structure_id="struct_001",
                code="ORCA",
                status="invalid_status"
            )
    
    def test_status_lifecycle_planned_to_submitted(self):
        """Test status transition from planned to submitted."""
        run = Run(
            structure_id="struct_001",
            code="ORCA",
            status="planned"
        )
        
        run.mark_submitted(job_id="12345")
        
        assert run.status == "submitted"
        assert run.metadata["job_id"] == "12345"
        assert "submitted_at" in run.metadata
    
    def test_status_lifecycle_submitted_to_running(self):
        """Test status transition from submitted to running."""
        run = Run(
            structure_id="struct_001",
            code="ORCA",
            status="submitted"
        )
        
        run.mark_running()
        
        assert run.status == "running"
        assert "started_at" in run.metadata
    
    def test_status_lifecycle_running_to_finished(self):
        """Test status transition from running to finished."""
        run = Run(
            structure_id="struct_001",
            code="ORCA",
            status="running"
        )
        
        run.mark_finished()
        
        assert run.status == "finished"
        assert "finished_at" in run.metadata
        assert run.is_complete
        assert run.is_successful
    
    def test_status_lifecycle_running_to_failed(self):
        """Test status transition from running to failed."""
        run = Run(
            structure_id="struct_001",
            code="ORCA",
            status="running"
        )
        
        run.mark_failed(error="Convergence failure")
        
        assert run.status == "failed"
        assert "failed_at" in run.metadata
        assert run.metadata["error"] == "Convergence failure"
        assert run.is_complete
        assert not run.is_successful
    
    def test_full_lifecycle(self):
        """Test complete run lifecycle."""
        run = Run(
            structure_id="struct_001",
            code="ORCA",
            parameters={"method": "DFT"}
        )
        
        # Initial state
        assert run.status == "planned"
        assert not run.is_complete
        
        # Submit
        run.mark_submitted(job_id="job_12345")
        assert run.status == "submitted"
        
        # Start running
        run.mark_running()
        assert run.status == "running"
        
        # Finish successfully
        run.mark_finished()
        assert run.status == "finished"
        assert run.is_complete
        assert run.is_successful
    
    def test_set_and_get_results(self):
        """Test setting and retrieving results."""
        run = Run(
            structure_id="struct_001",
            code="ORCA",
            status="finished"
        )
        
        run.set_result("energy", -76.4321)
        run.set_result("dipole", [0.1, 0.2, 0.3])
        run.set_result("iterations", 25)
        
        assert run.results["energy"] == -76.4321
        assert run.results["dipole"] == [0.1, 0.2, 0.3]
        assert run.results["iterations"] == 25
    
    def test_get_energy_default_key(self):
        """Test getting energy with default key."""
        run = Run(
            structure_id="struct_001",
            code="ORCA",
            status="finished"
        )
        
        run.set_result("total_energy", -152.7893)
        
        assert run.get_energy() == -152.7893
    
    def test_get_energy_custom_key(self):
        """Test getting energy with custom key."""
        run = Run(
            structure_id="struct_001",
            code="ORCA",
            status="finished"
        )
        
        run.set_result("scf_energy", -150.1234)
        
        assert run.get_energy("scf_energy") == -150.1234
    
    def test_get_energy_not_found(self):
        """Test getting non-existent energy returns None."""
        run = Run(
            structure_id="struct_001",
            code="ORCA",
            status="finished"
        )
        
        assert run.get_energy() is None
        assert run.get_energy("nonexistent") is None
    
    def test_to_core_object(self):
        """Test conversion to CoreObject."""
        run = Run(
            structure_id="struct_001",
            code="ORCA",
            code_version="5.0.3",
            parameters={"method": "DFT"},
            status="finished",
            results={"energy": -76.4321},
            metadata={"description": "Test run"}
        )
        
        core_obj = run.to_core_object()
        
        assert isinstance(core_obj, CoreObject)
        assert core_obj.type == "run"
        assert core_obj.meta["structure_id"] == "struct_001"
        assert core_obj.meta["code"] == "ORCA"
        assert core_obj.meta["code_version"] == "5.0.3"
        assert core_obj.meta["status"] == "finished"
        assert core_obj.meta["parameters"] == {"method": "DFT"}
        assert core_obj.meta["results"] == {"energy": -76.4321}
        
        # Check reference to structure
        assert len(core_obj.refs) == 1
        assert core_obj.refs[0].kind == "object"
        assert core_obj.refs[0].id == "struct_001"
    
    def test_from_core_object(self):
        """Test reconstruction from CoreObject."""
        original = Run(
            structure_id="struct_002",
            code="Gaussian",
            code_version="16",
            parameters={"basis": "6-31G*"},
            status="running",
            results={"iterations": 10},
            resources={"nodes": 4},
            metadata={"queue": "normal"}
        )
        
        core_obj = original.to_core_object()
        reconstructed = Run.from_core_object(core_obj)
        
        assert reconstructed.structure_id == original.structure_id
        assert reconstructed.code == original.code
        assert reconstructed.code_version == original.code_version
        assert reconstructed.parameters == original.parameters
        assert reconstructed.status == original.status
        assert reconstructed.results == original.results
        assert reconstructed.resources == original.resources
        assert reconstructed.metadata["queue"] == "normal"
    
    def test_round_trip_conversion(self):
        """Test that to/from CoreObject is lossless."""
        original = Run(
            structure_id="struct_003",
            code="VASP",
            code_version="6.3.0",
            parameters={
                "ENCUT": 500,
                "ISMEAR": 0,
                "SIGMA": 0.05
            },
            status="finished",
            results={
                "energy": -10.5678,
                "forces": [[0.1, 0.2, 0.3]]
            },
            resources={
                "cores": 16,
                "walltime": "2:00:00"
            },
            metadata={"project": "test_project"}
        )
        
        core_obj = original.to_core_object()
        reconstructed = Run.from_core_object(core_obj)
        
        # Check all attributes match
        assert reconstructed.structure_id == original.structure_id
        assert reconstructed.code == original.code
        assert reconstructed.code_version == original.code_version
        assert reconstructed.parameters == original.parameters
        assert reconstructed.status == original.status
        assert reconstructed.results == original.results
        assert reconstructed.resources == original.resources
        assert reconstructed.metadata == original.metadata
    
    def test_repr(self):
        """Test string representation."""
        run = Run(
            structure_id="struct_001",
            code="ORCA",
            status="running"
        )
        
        repr_str = repr(run)
        assert "Run" in repr_str
        assert "running" in repr_str
        assert "ORCA" in repr_str
        assert "struct_" in repr_str
