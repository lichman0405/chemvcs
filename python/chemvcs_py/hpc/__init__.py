"""HPC integration module for ChemVCS.

This module provides integration with HPC job schedulers (SLURM, PBS, etc.)
for submitting computational runs and tracking their execution.
"""

from .adapter import JobAdapter, JobStatus, JobInfo
from .slurm_adapter import SlurmAdapter
from .exceptions import (
    HpcError,
    JobSubmissionError,
    JobNotFoundError,
    AdapterNotAvailableError,
    InvalidJobStateError
)
from .provenance import EnvironmentCapture
from .submission import JobSubmitter
from .tracking import JobTracker, TrackedJob
from .retrieval import JobRetriever

__all__ = [
    "JobAdapter",
    "JobStatus",
    "JobInfo",
    "SlurmAdapter",
    "HpcError",
    "JobSubmissionError",
    "JobNotFoundError",
    "AdapterNotAvailableError",
    "InvalidJobStateError",
    "EnvironmentCapture",
    "JobSubmitter",
    "JobTracker",
    "TrackedJob",
    "JobRetriever",
]
