"""Core engine layer for ChemVCS.

This module provides the core business logic for version control operations,
including staging, commit building, and history management.
"""

from chemvcs.core.staging import StagingManager, StagingError

__all__ = [
    "StagingManager",
    "StagingError",
]
