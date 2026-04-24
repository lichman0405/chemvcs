"""Remote collaboration module for ChemVCS.

Provides push/pull operations over SSH + rsync.
"""

from __future__ import annotations

from chemvcs.remote.config import (
    RemoteAlreadyExistsError,
    RemoteNotFoundError,
    add_remote,
    get_remote_url,
    list_remotes,
    remove_remote,
)
from chemvcs.remote.manager import (
    RemoteConnectionError,
    RemoteManager,
    RemoteNotInitializedError,
)

__all__: list[str] = [
    "RemoteManager",
    "RemoteConnectionError",
    "RemoteNotInitializedError",
    "RemoteNotFoundError",
    "RemoteAlreadyExistsError",
    "list_remotes",
    "add_remote",
    "remove_remote",
    "get_remote_url",
]
