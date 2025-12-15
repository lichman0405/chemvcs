"""Tests for SLURM adapter."""

import pytest
from unittest import mock
import subprocess

from chemvcs_py.hpc.slurm_adapter import SlurmAdapter
from chemvcs_py.hpc.adapter import JobStatus
from chemvcs_py.hpc.exceptions import (
    JobSubmissionError,
    JobNotFoundError,
    AdapterNotAvailableError
)


class TestSlurmAdapter:
    """Test cases for SlurmAdapter."""
    
    def test_adapter_name(self):
        """Test adapter name property."""
        adapter = SlurmAdapter()
        assert adapter.name == "slurm"
    
    @mock.patch('shutil.which')
    def test_validate_available(self, mock_which):
        """Test validate when SLURM is available."""
        mock_which.return_value = '/usr/bin/sbatch'
        
        adapter = SlurmAdapter()
        assert adapter.validate() is True
    
    @mock.patch('shutil.which')
    def test_validate_not_available(self, mock_which):
        """Test validate when SLURM is not available."""
        mock_which.return_value = None
        
        adapter = SlurmAdapter()
        assert adapter.validate() is False
    
    @mock.patch('shutil.which')
    @mock.patch('subprocess.run')
    def test_submit_success(self, mock_run, mock_which):
        """Test successful job submission."""
        mock_which.return_value = '/usr/bin/sbatch'
        mock_run.return_value = mock.Mock(
            returncode=0,
            stdout="Submitted batch job 12345\n",
            stderr=""
        )
        
        adapter = SlurmAdapter()
        job_id = adapter.submit('test.slurm')
        
        assert job_id == "12345"
        mock_run.assert_called_once()
        args = mock_run.call_args[0][0]
        assert args[0] == 'sbatch'
        assert args[-1] == 'test.slurm'
    
    @mock.patch('shutil.which')
    @mock.patch('subprocess.run')
    def test_submit_with_options(self, mock_run, mock_which):
        """Test submission with optional arguments."""
        mock_which.return_value = '/usr/bin/sbatch'
        mock_run.return_value = mock.Mock(
            returncode=0,
            stdout="Submitted batch job 12345\n",
            stderr=""
        )
        
        adapter = SlurmAdapter()
        job_id = adapter.submit(
            'test.slurm',
            job_name='my_job',
            output='output.log',
            partition='batch'
        )
        
        assert job_id == "12345"
        args = mock_run.call_args[0][0]
        assert '--job-name' in args
        assert 'my_job' in args
        assert '--output' in args
        assert 'output.log' in args
        assert '--partition' in args
        assert 'batch' in args
    
    @mock.patch('shutil.which')
    @mock.patch('subprocess.run')
    def test_submit_failure(self, mock_run, mock_which):
        """Test failed job submission."""
        mock_which.return_value = '/usr/bin/sbatch'
        mock_run.side_effect = subprocess.CalledProcessError(
            returncode=1,
            cmd=['sbatch', 'test.slurm'],
            stderr="Error: Invalid script"
        )
        
        adapter = SlurmAdapter()
        with pytest.raises(JobSubmissionError, match="sbatch failed"):
            adapter.submit('test.slurm')
    
    @mock.patch('shutil.which')
    def test_submit_not_available(self, mock_which):
        """Test submission when SLURM not available."""
        mock_which.return_value = None
        
        adapter = SlurmAdapter()
        with pytest.raises(AdapterNotAvailableError, match="sbatch command not found"):
            adapter.submit('test.slurm')
    
    @mock.patch('shutil.which')
    @mock.patch('subprocess.run')
    def test_submit_invalid_output(self, mock_run, mock_which):
        """Test submission with unparseable output."""
        mock_which.return_value = '/usr/bin/sbatch'
        mock_run.return_value = mock.Mock(
            returncode=0,
            stdout="Some unexpected output\n",
            stderr=""
        )
        
        adapter = SlurmAdapter()
        with pytest.raises(JobSubmissionError, match="Failed to parse job ID"):
            adapter.submit('test.slurm')
    
    @mock.patch('subprocess.run')
    def test_get_status_pending(self, mock_run):
        """Test status query for pending job."""
        mock_run.return_value = mock.Mock(
            returncode=0,
            stdout="PENDING\n",
            stderr=""
        )
        
        adapter = SlurmAdapter()
        status = adapter.get_status('12345')
        
        assert status == JobStatus.PENDING
    
    @mock.patch('subprocess.run')
    def test_get_status_running(self, mock_run):
        """Test status query for running job."""
        mock_run.return_value = mock.Mock(
            returncode=0,
            stdout="RUNNING\n",
            stderr=""
        )
        
        adapter = SlurmAdapter()
        status = adapter.get_status('12345')
        
        assert status == JobStatus.RUNNING
    
    @mock.patch('subprocess.run')
    def test_get_status_completed(self, mock_run):
        """Test status query for completed job (from sacct)."""
        # First call: squeue returns nothing (job not in queue)
        # Second call: sacct returns COMPLETED
        mock_run.side_effect = [
            mock.Mock(returncode=1, stdout="", stderr=""),  # squeue
            mock.Mock(returncode=0, stdout="COMPLETED\n", stderr="")  # sacct
        ]
        
        adapter = SlurmAdapter()
        status = adapter.get_status('12345')
        
        assert status == JobStatus.COMPLETED
        assert mock_run.call_count == 2
    
    @mock.patch('subprocess.run')
    def test_get_status_failed(self, mock_run):
        """Test status query for failed job."""
        mock_run.side_effect = [
            mock.Mock(returncode=1, stdout="", stderr=""),
            mock.Mock(returncode=0, stdout="FAILED\n", stderr="")
        ]
        
        adapter = SlurmAdapter()
        status = adapter.get_status('12345')
        
        assert status == JobStatus.FAILED
    
    @mock.patch('subprocess.run')
    def test_get_status_timeout(self, mock_run):
        """Test status query with timeout."""
        mock_run.side_effect = [
            mock.Mock(returncode=1, stdout="", stderr=""),
            mock.Mock(returncode=0, stdout="TIMEOUT\n", stderr="")
        ]
        
        adapter = SlurmAdapter()
        status = adapter.get_status('12345')
        
        assert status == JobStatus.FAILED  # TIMEOUT maps to FAILED
    
    @mock.patch('subprocess.run')
    def test_get_status_cancelled(self, mock_run):
        """Test status query for cancelled job."""
        mock_run.side_effect = [
            mock.Mock(returncode=1, stdout="", stderr=""),
            mock.Mock(returncode=0, stdout="CANCELLED\n", stderr="")
        ]
        
        adapter = SlurmAdapter()
        status = adapter.get_status('12345')
        
        assert status == JobStatus.CANCELLED
    
    @mock.patch('subprocess.run')
    def test_get_status_not_found(self, mock_run):
        """Test status query for non-existent job."""
        mock_run.side_effect = [
            mock.Mock(returncode=1, stdout="", stderr=""),
            mock.Mock(returncode=1, stdout="", stderr="")
        ]
        
        adapter = SlurmAdapter()
        with pytest.raises(JobNotFoundError, match="Job 12345 not found"):
            adapter.get_status('12345')
    
    @mock.patch('subprocess.run')
    def test_get_info_success(self, mock_run):
        """Test getting detailed job information."""
        # sacct output format: State Partition NNodes Start End ExitCode
        mock_run.return_value = mock.Mock(
            returncode=0,
            stdout="COMPLETED batch 2 2025-12-15T10:30:00 2025-12-15T12:30:00 0:0\n",
            stderr=""
        )
        
        adapter = SlurmAdapter()
        info = adapter.get_info('12345')
        
        assert info.job_id == "12345"
        assert info.status == JobStatus.COMPLETED
        assert info.queue == "batch"
        assert info.nodes == 2
        assert info.start_time == "2025-12-15T10:30:00"
        assert info.end_time == "2025-12-15T12:30:00"
        assert info.exit_code == 0
    
    @mock.patch('subprocess.run')
    def test_get_info_not_found(self, mock_run):
        """Test get_info for non-existent job."""
        mock_run.return_value = mock.Mock(
            returncode=1,
            stdout="",
            stderr=""
        )
        
        adapter = SlurmAdapter()
        with pytest.raises(JobNotFoundError, match="Job 12345 not found"):
            adapter.get_info('12345')
    
    @mock.patch('subprocess.run')
    def test_cancel_success(self, mock_run):
        """Test successful job cancellation."""
        mock_run.return_value = mock.Mock(
            returncode=0,
            stdout="",
            stderr=""
        )
        
        adapter = SlurmAdapter()
        result = adapter.cancel('12345')
        
        assert result is True
        mock_run.assert_called_once()
        args = mock_run.call_args[0][0]
        assert args == ['scancel', '12345']
    
    @mock.patch('subprocess.run')
    def test_cancel_failure(self, mock_run):
        """Test failed job cancellation."""
        mock_run.return_value = mock.Mock(
            returncode=1,
            stdout="",
            stderr="Error: Job not found"
        )
        
        adapter = SlurmAdapter()
        result = adapter.cancel('12345')
        
        assert result is False
    
    @mock.patch('subprocess.run')
    def test_cancel_timeout(self, mock_run):
        """Test cancellation with timeout."""
        mock_run.side_effect = subprocess.TimeoutExpired(
            cmd=['scancel', '12345'],
            timeout=10
        )
        
        adapter = SlurmAdapter()
        result = adapter.cancel('12345')
        
        assert result is False
    
    def test_parse_status_with_modifiers(self):
        """Test parsing status with modifiers like RUNNING+."""
        adapter = SlurmAdapter()
        
        assert adapter._parse_slurm_status("RUNNING+") == JobStatus.RUNNING
        assert adapter._parse_slurm_status("PENDING extra") == JobStatus.PENDING
        assert adapter._parse_slurm_status("COMPLETED+") == JobStatus.COMPLETED
