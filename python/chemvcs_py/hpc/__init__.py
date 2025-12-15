"""HPC integration module for ChemVCS.

This module provides integration with HPC job schedulers (SLURM, PBS, etc.)
for submitting computational runs and tracking their execution.
"""

from .adapter import JobAdapter, JobStatus, JobInfo
from .slurm_adapter import SlurmAdapter
from .pbs_adapter import PbsAdapter
from .lsf_adapter import LsfAdapter
from .exceptions import (
    HpcError,
    JobSubmissionError,
    JobNotFoundError,
    AdapterNotAvailableError,
    InvalidJobStateError,
    JobCancellationError,
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
    "PbsAdapter",
    "LsfAdapter",
    "HpcError",
    "JobSubmissionError",
    "JobNotFoundError",
    "AdapterNotAvailableError",
    "InvalidJobStateError",
    "JobCancellationError",
    "EnvironmentCapture",
    "JobSubmitter",
    "JobTracker",
    "TrackedJob",
    "JobRetriever",
]
