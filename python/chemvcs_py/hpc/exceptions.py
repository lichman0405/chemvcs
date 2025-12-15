"""HPC-related exceptions."""


class HpcError(Exception):
    """Base exception for HPC operations."""
    pass


class JobSubmissionError(HpcError):
    """Failed to submit job to scheduler."""
    pass


class JobNotFoundError(HpcError):
    """Job ID not found in scheduler."""
    pass


class AdapterNotAvailableError(HpcError):
    """Job adapter not available (missing commands or environment)."""
    pass


class InvalidJobStateError(HpcError):
    """Operation invalid for current job state."""
    pass


class JobCancellationError(HpcError):
    """Failed to cancel job."""
    pass
