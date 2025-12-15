"""Tests for environment provenance capture."""

import pytest
import os
import tempfile
from unittest import mock

from chemvcs_py.hpc.provenance import EnvironmentCapture


class TestEnvironmentCapture:
    """Test cases for EnvironmentCapture."""
    
    @mock.patch('subprocess.run')
    def test_capture_modules_success(self, mock_run):
        """Test capturing loaded modules."""
        mock_run.return_value = mock.Mock(
            returncode=0,
            stdout="vasp/6.3.0\nintel/2021.4\nmkl/2021.4\n"
        )
        
        modules = EnvironmentCapture.capture_modules()
        
        assert modules == ['vasp/6.3.0', 'intel/2021.4', 'mkl/2021.4']
    
    @mock.patch('subprocess.run')
    def test_capture_modules_with_header(self, mock_run):
        """Test capturing modules with header lines."""
        mock_run.return_value = mock.Mock(
            returncode=0,
            stdout="Currently Loaded Modulefiles:\nvasp/6.3.0\nintel/2021.4\n"
        )
        
        modules = EnvironmentCapture.capture_modules()
        
        # Should skip "Currently" line
        assert modules == ['vasp/6.3.0', 'intel/2021.4']
    
    @mock.patch('subprocess.run')
    def test_capture_modules_command_not_found(self, mock_run):
        """Test when module command doesn't exist."""
        mock_run.side_effect = FileNotFoundError()
        
        modules = EnvironmentCapture.capture_modules()
        
        assert modules == []
    
    @mock.patch('subprocess.run')
    def test_capture_modules_timeout(self, mock_run):
        """Test when module command times out."""
        import subprocess
        mock_run.side_effect = subprocess.TimeoutExpired(
            cmd=['module', 'list'],
            timeout=5
        )
        
        modules = EnvironmentCapture.capture_modules()
        
        assert modules == []
    
    def test_capture_env_vars_specific(self):
        """Test capturing specific environment variables."""
        os.environ['TEST_VAR1'] = 'value1'
        os.environ['TEST_VAR2'] = 'value2'
        
        env_vars = EnvironmentCapture.capture_env_vars(['TEST_VAR1', 'TEST_VAR2', 'NONEXISTENT'])
        
        assert env_vars == {
            'TEST_VAR1': 'value1',
            'TEST_VAR2': 'value2'
        }
        
        # Cleanup
        del os.environ['TEST_VAR1']
        del os.environ['TEST_VAR2']
    
    def test_capture_env_vars_default(self):
        """Test capturing default environment variables."""
        os.environ['OMP_NUM_THREADS'] = '4'
        
        env_vars = EnvironmentCapture.capture_env_vars()
        
        assert 'OMP_NUM_THREADS' in env_vars
        assert env_vars['OMP_NUM_THREADS'] == '4'
        
        # Cleanup
        del os.environ['OMP_NUM_THREADS']
    
    def test_parse_slurm_script(self):
        """Test parsing SLURM script directives."""
        script_content = """#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=28
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --partition=batch
#SBATCH --job-name=my_job

module load vasp/6.3.0
mpirun vasp_std
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.slurm', delete=False) as f:
            f.write(script_content)
            script_path = f.name
        
        try:
            resources = EnvironmentCapture.parse_slurm_script(script_path)
            
            assert resources['nodes'] == 4
            assert resources['ntasks_per_node'] == 28
            assert resources['walltime'] == '24:00:00'
            assert resources['memory'] == '64G'
            assert resources['partition'] == 'batch'
            assert resources['job_name'] == 'my_job'
        finally:
            os.unlink(script_path)
    
    def test_parse_slurm_script_minimal(self):
        """Test parsing minimal SLURM script."""
        script_content = """#!/bin/bash
#SBATCH --nodes=2
#SBATCH --time=1:00:00

echo "Running job"
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.slurm', delete=False) as f:
            f.write(script_content)
            script_path = f.name
        
        try:
            resources = EnvironmentCapture.parse_slurm_script(script_path)
            
            assert resources['nodes'] == 2
            assert resources['walltime'] == '1:00:00'
            assert 'ntasks_per_node' not in resources
        finally:
            os.unlink(script_path)
    
    def test_parse_slurm_script_file_not_found(self):
        """Test parsing non-existent script."""
        resources = EnvironmentCapture.parse_slurm_script('/nonexistent/script.slurm')
        
        assert resources == {}
    
    def test_parse_pbs_script(self):
        """Test parsing PBS script directives."""
        script_content = """#!/bin/bash
#PBS -l nodes=2:ppn=16
#PBS -l walltime=12:00:00
#PBS -l mem=32gb
#PBS -q batch
#PBS -N my_pbs_job

cd $PBS_O_WORKDIR
mpirun my_program
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pbs', delete=False) as f:
            f.write(script_content)
            script_path = f.name
        
        try:
            resources = EnvironmentCapture.parse_pbs_script(script_path)
            
            assert resources['nodes'] == 2
            assert resources['ppn'] == 16
            assert resources['walltime'] == '12:00:00'
            assert resources['memory'] == '32gb'
            assert resources['queue'] == 'batch'
            assert resources['job_name'] == 'my_pbs_job'
        finally:
            os.unlink(script_path)
    
    def test_save_script_snapshot(self):
        """Test saving script snapshot."""
        script_content = """#!/bin/bash
#SBATCH --nodes=2
echo "Test"
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.slurm', delete=False) as f:
            f.write(script_content)
            script_path = f.name
        
        try:
            snapshot = EnvironmentCapture.save_script_snapshot(script_path)
            
            assert snapshot == script_content
        finally:
            os.unlink(script_path)
    
    def test_save_script_snapshot_not_found(self):
        """Test snapshot of non-existent script."""
        snapshot = EnvironmentCapture.save_script_snapshot('/nonexistent/script.slurm')
        
        assert snapshot == ""
    
    def test_detect_script_type_slurm(self):
        """Test detecting SLURM script."""
        script_content = """#!/bin/bash
#SBATCH --nodes=2
echo "Test"
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.slurm', delete=False) as f:
            f.write(script_content)
            script_path = f.name
        
        try:
            script_type = EnvironmentCapture.detect_script_type(script_path)
            assert script_type == 'slurm'
        finally:
            os.unlink(script_path)
    
    def test_detect_script_type_pbs(self):
        """Test detecting PBS script."""
        script_content = """#!/bin/bash
#PBS -l nodes=2
echo "Test"
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pbs', delete=False) as f:
            f.write(script_content)
            script_path = f.name
        
        try:
            script_type = EnvironmentCapture.detect_script_type(script_path)
            assert script_type == 'pbs'
        finally:
            os.unlink(script_path)
    
    def test_detect_script_type_sge(self):
        """Test detecting SGE script."""
        script_content = """#!/bin/bash
#$ -pe mpi 16
echo "Test"
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sge', delete=False) as f:
            f.write(script_content)
            script_path = f.name
        
        try:
            script_type = EnvironmentCapture.detect_script_type(script_path)
            assert script_type == 'sge'
        finally:
            os.unlink(script_path)
    
    def test_detect_script_type_unknown(self):
        """Test detecting unknown script type."""
        script_content = """#!/bin/bash
echo "Plain bash script"
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sh', delete=False) as f:
            f.write(script_content)
            script_path = f.name
        
        try:
            script_type = EnvironmentCapture.detect_script_type(script_path)
            assert script_type is None
        finally:
            os.unlink(script_path)
