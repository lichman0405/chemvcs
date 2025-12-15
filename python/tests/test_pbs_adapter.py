"""Tests for PBS adapter."""

import pytest
from unittest.mock import patch, MagicMock
import subprocess

from chemvcs_py.hpc import PbsAdapter, JobStatus, JobNotFoundError, JobSubmissionError, AdapterNotAvailableError


@pytest.fixture
def pbs_adapter():
    """Create PbsAdapter instance."""
    return PbsAdapter()


class TestPbsSubmit:
    """Test PBS job submission."""
    
    def test_submit_success(self, pbs_adapter):
        """Test successful job submission."""
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = "12345.hostname.domain\n"
        mock_result.stderr = ""
        
        with patch('subprocess.run', return_value=mock_result):
            with patch('shutil.which', return_value='/usr/bin/qsub'):
                job_id = pbs_adapter.submit('/path/to/script.sh')
        
        assert job_id == "12345.hostname.domain"
    
    def test_submit_with_options(self, pbs_adapter):
        """Test submission with optional parameters."""
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = "12345.host\n"
        
        with patch('subprocess.run', return_value=mock_result) as mock_run:
            with patch('shutil.which', return_value='/usr/bin/qsub'):
                job_id = pbs_adapter.submit(
                    '/path/to/script.sh',
                    job_name='test_job',
                    queue='batch',
                    output='output.log'
                )
        
        # Verify command includes options
        call_args = mock_run.call_args[0][0]
        assert '-N' in call_args
        assert 'test_job' in call_args
        assert '-q' in call_args
        assert 'batch' in call_args
        assert '-o' in call_args
        assert 'output.log' in call_args
    
    def test_submit_failure(self, pbs_adapter):
        """Test submission failure."""
        with patch('subprocess.run', side_effect=subprocess.CalledProcessError(1, 'qsub', stderr='Permission denied')):
            with patch('shutil.which', return_value='/usr/bin/qsub'):
                with pytest.raises(JobSubmissionError, match='qsub failed'):
                    pbs_adapter.submit('/path/to/script.sh')
    
    def test_submit_qsub_not_found(self, pbs_adapter):
        """Test when qsub command not available."""
        with patch('shutil.which', return_value=None):
            with pytest.raises(AdapterNotAvailableError, match='qsub command not found'):
                pbs_adapter.submit('/path/to/script.sh')


class TestPbsStatus:
    """Test PBS job status queries."""
    
    def test_status_running(self, pbs_adapter):
        """Test querying running job."""
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = "Job Id: 12345.host\n    job_state = R\n"
        
        with patch('subprocess.run', return_value=mock_result):
            status = pbs_adapter.get_status('12345.host')
        
        assert status == JobStatus.RUNNING
    
    def test_status_pending(self, pbs_adapter):
        """Test querying pending job."""
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = "Job Id: 12345.host\n    job_state = Q\n"
        
        with patch('subprocess.run', return_value=mock_result):
            status = pbs_adapter.get_status('12345.host')
        
        assert status == JobStatus.PENDING
    
    def test_status_completed(self, pbs_adapter):
        """Test querying completed job (not in queue)."""
        mock_result = MagicMock()
        mock_result.returncode = 1
        mock_result.stderr = "qstat: Unknown Job Id 12345.host"
        
        with patch('subprocess.run', return_value=mock_result):
            status = pbs_adapter.get_status('12345.host')
        
        # Completed jobs return COMPLETED
        assert status == JobStatus.COMPLETED
    
    def test_status_held(self, pbs_adapter):
        """Test querying held job."""
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = "Job Id: 12345.host\n    job_state = H\n"
        
        with patch('subprocess.run', return_value=mock_result):
            status = pbs_adapter.get_status('12345.host')
        
        assert status == JobStatus.PENDING
    
    def test_status_not_found(self, pbs_adapter):
        """Test querying non-existent job."""
        mock_result = MagicMock()
        mock_result.returncode = 1
        mock_result.stderr = "qstat: 12345 not found"
        
        with patch('subprocess.run', return_value=mock_result):
            with pytest.raises(JobNotFoundError):
                pbs_adapter.get_status('12345')


class TestPbsInfo:
    """Test PBS job info queries."""
    
    def test_get_info_success(self, pbs_adapter):
        """Test getting detailed job information."""
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = """Job Id: 12345.host
    job_state = R
    queue = batch
    Resource_List.nodes = 4
    start_time = Mon Dec 15 10:30:00 2025
    exit_status = 0
"""
        
        with patch('subprocess.run', return_value=mock_result):
            info = pbs_adapter.get_info('12345.host')
        
        assert info.job_id == '12345.host'
        assert info.status == JobStatus.RUNNING
        assert info.queue == 'batch'
        assert info.nodes == 4
        assert info.start_time == 'Mon Dec 15 10:30:00 2025'
        assert info.exit_code == 0
    
    def test_get_info_minimal(self, pbs_adapter):
        """Test getting info with minimal output."""
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = "Job Id: 12345.host\n    job_state = Q\n"
        
        with patch('subprocess.run', return_value=mock_result):
            info = pbs_adapter.get_info('12345.host')
        
        assert info.job_id == '12345.host'
        assert info.status == JobStatus.PENDING
        assert info.queue is None
        assert info.nodes is None
    
    def test_get_info_not_found(self, pbs_adapter):
        """Test info for non-existent job."""
        mock_result = MagicMock()
        mock_result.returncode = 1
        
        with patch('subprocess.run', return_value=mock_result):
            with pytest.raises(JobNotFoundError):
                pbs_adapter.get_info('99999')


class TestPbsCancel:
    """Test PBS job cancellation."""
    
    def test_cancel_success(self, pbs_adapter):
        """Test successful job cancellation."""
        mock_result = MagicMock()
        mock_result.returncode = 0
        
        with patch('subprocess.run', return_value=mock_result):
            result = pbs_adapter.cancel('12345.host')
        
        assert result is True
    
    def test_cancel_failure(self, pbs_adapter):
        """Test failed job cancellation."""
        mock_result = MagicMock()
        mock_result.returncode = 1
        
        with patch('subprocess.run', return_value=mock_result):
            result = pbs_adapter.cancel('12345.host')
        
        assert result is False
    
    def test_cancel_timeout(self, pbs_adapter):
        """Test cancellation timeout."""
        with patch('subprocess.run', side_effect=subprocess.TimeoutExpired('qdel', 10)):
            result = pbs_adapter.cancel('12345.host')
        
        assert result is False


class TestPbsValidate:
    """Test PBS adapter validation."""
    
    def test_validate_available(self, pbs_adapter):
        """Test validation when qsub is available."""
        with patch('shutil.which', return_value='/usr/bin/qsub'):
            assert pbs_adapter.validate() is True
    
    def test_validate_not_available(self, pbs_adapter):
        """Test validation when qsub is not available."""
        with patch('shutil.which', return_value=None):
            assert pbs_adapter.validate() is False
    
    def test_name(self, pbs_adapter):
        """Test adapter name."""
        assert pbs_adapter.name == "pbs"
