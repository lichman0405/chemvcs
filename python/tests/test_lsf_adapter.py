"""Tests for LSF adapter."""

import pytest
from unittest.mock import patch, MagicMock
import subprocess

from chemvcs_py.hpc import LsfAdapter, JobStatus, JobNotFoundError, JobSubmissionError, AdapterNotAvailableError


@pytest.fixture
def lsf_adapter():
    """Create LsfAdapter instance."""
    return LsfAdapter()


class TestLsfSubmit:
    """Test LSF job submission."""
    
    def test_submit_success(self, lsf_adapter):
        """Test successful job submission."""
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = "Job <12345> is submitted to queue <normal>.\n"
        mock_result.stderr = ""
        
        with patch('subprocess.run', return_value=mock_result):
            with patch('shutil.which', return_value='/usr/bin/bsub'):
                job_id = lsf_adapter.submit('/path/to/script.sh')
        
        assert job_id == "12345"
    
    def test_submit_with_options(self, lsf_adapter):
        """Test submission with optional parameters."""
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = "Job <12345> is submitted to queue <batch>.\n"
        
        with patch('subprocess.run', return_value=mock_result) as mock_run:
            with patch('shutil.which', return_value='/usr/bin/bsub'):
                job_id = lsf_adapter.submit(
                    '/path/to/script.sh',
                    job_name='test_job',
                    queue='batch',
                    output='output.log'
                )
        
        # Verify command includes options (as a shell command string)
        call_args = mock_run.call_args[0][0]
        assert '-J test_job' in call_args
        assert '-q batch' in call_args
        assert '-o output.log' in call_args
    
    def test_submit_failure(self, lsf_adapter):
        """Test submission failure."""
        with patch('subprocess.run', side_effect=subprocess.CalledProcessError(1, 'bsub', stderr='Permission denied')):
            with patch('shutil.which', return_value='/usr/bin/bsub'):
                with pytest.raises(JobSubmissionError, match='bsub failed'):
                    lsf_adapter.submit('/path/to/script.sh')
    
    def test_submit_parse_failure(self, lsf_adapter):
        """Test when job ID cannot be parsed."""
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = "Unexpected output format\n"
        
        with patch('subprocess.run', return_value=mock_result):
            with patch('shutil.which', return_value='/usr/bin/bsub'):
                with pytest.raises(JobSubmissionError, match='Failed to parse job ID'):
                    lsf_adapter.submit('/path/to/script.sh')
    
    def test_submit_bsub_not_found(self, lsf_adapter):
        """Test when bsub command not available."""
        with patch('shutil.which', return_value=None):
            with pytest.raises(AdapterNotAvailableError, match='bsub command not found'):
                lsf_adapter.submit('/path/to/script.sh')


class TestLsfStatus:
    """Test LSF job status queries."""
    
    def test_status_running(self, lsf_adapter):
        """Test querying running job."""
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = "RUN\n"
        
        with patch('subprocess.run', return_value=mock_result):
            status = lsf_adapter.get_status('12345')
        
        assert status == JobStatus.RUNNING
    
    def test_status_pending(self, lsf_adapter):
        """Test querying pending job."""
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = "PEND\n"
        
        with patch('subprocess.run', return_value=mock_result):
            status = lsf_adapter.get_status('12345')
        
        assert status == JobStatus.PENDING
    
    def test_status_done(self, lsf_adapter):
        """Test querying completed job."""
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = "DONE\n"
        
        with patch('subprocess.run', return_value=mock_result):
            status = lsf_adapter.get_status('12345')
        
        assert status == JobStatus.COMPLETED
    
    def test_status_exit(self, lsf_adapter):
        """Test querying failed job."""
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = "EXIT\n"
        
        with patch('subprocess.run', return_value=mock_result):
            status = lsf_adapter.get_status('12345')
        
        assert status == JobStatus.FAILED
    
    def test_status_suspended(self, lsf_adapter):
        """Test querying suspended job."""
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = "USUSP\n"
        
        with patch('subprocess.run', return_value=mock_result):
            status = lsf_adapter.get_status('12345')
        
        assert status == JobStatus.CANCELLED
    
    def test_status_not_in_queue(self, lsf_adapter):
        """Test querying job not in queue (completed)."""
        mock_result = MagicMock()
        mock_result.returncode = 1
        mock_result.stderr = "Job <12345> is not found"
        
        with patch('subprocess.run', return_value=mock_result):
            status = lsf_adapter.get_status('12345')
        
        # Jobs not in queue are assumed completed
        assert status == JobStatus.COMPLETED
    
    def test_status_command_failed(self, lsf_adapter):
        """Test when bjobs command fails."""
        mock_result = MagicMock()
        mock_result.returncode = 1
        mock_result.stderr = "Connection error"
        
        with patch('subprocess.run', return_value=mock_result):
            with pytest.raises(JobNotFoundError):
                lsf_adapter.get_status('12345')


class TestLsfInfo:
    """Test LSF job info queries."""
    
    def test_get_info_success(self, lsf_adapter):
        """Test getting detailed job information."""
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = """Job <12345>, User <user>, Project <default>, Status <RUN>
                     Queue <normal>, Command <script.sh>
 Mon Dec 15 10:30:00: Submitted from host <submit-host>
 Mon Dec 15 10:31:00: Started on <exec-host>
                     Requested 4 Task(s) on 1 Processor(s)
"""
        
        with patch('subprocess.run', return_value=mock_result):
            info = lsf_adapter.get_info('12345')
        
        assert info.job_id == '12345'
        assert info.status == JobStatus.RUNNING
        assert info.queue == 'normal'
        assert info.cores == 4
        assert 'Mon Dec 15 10:31:00' in info.start_time
    
    def test_get_info_minimal(self, lsf_adapter):
        """Test getting info with minimal output."""
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = "Job <12345>, Status <PEND>"
        
        with patch('subprocess.run', return_value=mock_result):
            info = lsf_adapter.get_info('12345')
        
        assert info.job_id == '12345'
        assert info.status == JobStatus.PENDING
        assert info.queue is None
    
    def test_get_info_not_found(self, lsf_adapter):
        """Test info for non-existent job."""
        mock_result = MagicMock()
        mock_result.returncode = 1
        
        with patch('subprocess.run', return_value=mock_result):
            with pytest.raises(JobNotFoundError):
                lsf_adapter.get_info('99999')


class TestLsfCancel:
    """Test LSF job cancellation."""
    
    def test_cancel_success(self, lsf_adapter):
        """Test successful job cancellation."""
        mock_result = MagicMock()
        mock_result.returncode = 0
        
        with patch('subprocess.run', return_value=mock_result):
            result = lsf_adapter.cancel('12345')
        
        assert result is True
    
    def test_cancel_failure(self, lsf_adapter):
        """Test failed job cancellation."""
        mock_result = MagicMock()
        mock_result.returncode = 1
        
        with patch('subprocess.run', return_value=mock_result):
            result = lsf_adapter.cancel('12345')
        
        assert result is False
    
    def test_cancel_timeout(self, lsf_adapter):
        """Test cancellation timeout."""
        with patch('subprocess.run', side_effect=subprocess.TimeoutExpired('bkill', 10)):
            result = lsf_adapter.cancel('12345')
        
        assert result is False


class TestLsfValidate:
    """Test LSF adapter validation."""
    
    def test_validate_available(self, lsf_adapter):
        """Test validation when bsub is available."""
        with patch('shutil.which', return_value='/usr/bin/bsub'):
            assert lsf_adapter.validate() is True
    
    def test_validate_not_available(self, lsf_adapter):
        """Test validation when bsub is not available."""
        with patch('shutil.which', return_value=None):
            assert lsf_adapter.validate() is False
    
    def test_name(self, lsf_adapter):
        """Test adapter name."""
        assert lsf_adapter.name == "lsf"
