"""Remote transport layer — SSH and rsync wrappers.

All network operations are performed via subprocess calls to system
``ssh`` and ``rsync`` binaries.  No third-party SSH library is required.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

from chemvcs.constants import CHEMVCS_DIR, COMMITS_DIR, OBJECTS_DIR


class RemoteConnectionError(Exception):
    """Raised when SSH connection to the remote host fails."""


class RemoteNotInitializedError(Exception):
    """Raised when the remote repository exists but is not a ChemVCS repo."""


class RemoteManager:
    """Manage a remote ChemVCS repository over SSH + rsync.

    Parameters
    ----------
    url:
        Remote URL in the form ``user@host:/path/to/repo``.
    """

    def __init__(self, url: str) -> None:
        if ":" not in url:
            raise ValueError(f"Invalid remote URL format: {url!r} (expected user@host:path)")
        self.user_host, self.remote_path = url.rsplit(":", 1)
        if not self.user_host or not self.remote_path:
            raise ValueError(f"Invalid remote URL format: {url!r} (expected user@host:path)")
        if "@" not in self.user_host:
            # user defaults to the current local user if not specified
            self.user_host = f"{os.environ.get('USER', 'unknown')}@{self.user_host}"

    # ------------------------------------------------------------------ #
    # SSH helpers                                                         #
    # ------------------------------------------------------------------ #

    def check_connection(self) -> bool:
        """Test that SSH to the remote host works.

        Returns:
            ``True`` when ``ssh user@host 'echo ok'`` succeeds.
        """
        try:
            result = subprocess.run(
                ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=10",
                 self.user_host, "echo", "ok"],
                capture_output=True,
                text=True,
                timeout=15,
            )
            return result.returncode == 0 and "ok" in result.stdout
        except (subprocess.TimeoutExpired, OSError):
            return False

    def get_remote_head(self) -> str | None:
        """Read the remote ``.chemvcs/HEAD`` file.

        Returns:
            Commit hash string or ``None`` if the remote has no commits yet.
        """
        head_path = f"{self.remote_path}/{CHEMVCS_DIR}/HEAD"
        try:
            result = subprocess.run(
                ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=10",
                 self.user_host, "cat", head_path],
                capture_output=True,
                text=True,
                timeout=15,
            )
            if result.returncode == 0 and result.stdout.strip():
                return result.stdout.strip()
            return None
        except (subprocess.TimeoutExpired, OSError):
            return None

    def init_remote_if_needed(self) -> None:
        """Create the ChemVCS directory structure on the remote if it is missing."""
        base = f"{self.remote_path}/{CHEMVCS_DIR}"
        subprocess.run(
            ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=10",
             self.user_host, "mkdir", "-p", f"{base}/{OBJECTS_DIR}", f"{base}/{COMMITS_DIR}"],
            capture_output=True,
            timeout=15,
            check=False,
        )

    def update_remote_head(self, commit_hash: str) -> None:
        """Write *commit_hash* into the remote ``.chemvcs/HEAD`` file."""
        head_path = f"{self.remote_path}/{CHEMVCS_DIR}/HEAD"
        subprocess.run(
            ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=10",
             self.user_host, "sh", "-c",
             f"printf '%s' {commit_hash} > {head_path}"],
            capture_output=True,
            timeout=15,
            check=True,
        )

    # ------------------------------------------------------------------ #
    # rsync helpers                                                       #
    # ------------------------------------------------------------------ #

    def push_objects(self, local_chemvcs_dir: Path) -> None:
        """Rsync local ``objects/`` → remote, skipping blobs that already exist."""
        src = f"{local_chemvcs_dir}/{OBJECTS_DIR}/"
        dst = f"{self.user_host}:{self.remote_path}/{CHEMVCS_DIR}/{OBJECTS_DIR}/"
        self._rsync(src, dst)

    def push_commits(self, local_chemvcs_dir: Path) -> None:
        """Rsync local ``commits/`` → remote, skipping files that already exist."""
        src = f"{local_chemvcs_dir}/{COMMITS_DIR}/"
        dst = f"{self.user_host}:{self.remote_path}/{CHEMVCS_DIR}/{COMMITS_DIR}/"
        self._rsync(src, dst)

    def fetch_objects(self, local_chemvcs_dir: Path) -> None:
        """Rsync remote ``objects/`` → local, skipping blobs that already exist."""
        src = f"{self.user_host}:{self.remote_path}/{CHEMVCS_DIR}/{OBJECTS_DIR}/"
        dst = f"{local_chemvcs_dir}/{OBJECTS_DIR}/"
        self._rsync(src, dst)

    def fetch_commits(self, local_chemvcs_dir: Path) -> None:
        """Rsync remote ``commits/`` → local, skipping files that already exist."""
        src = f"{self.user_host}:{self.remote_path}/{CHEMVCS_DIR}/{COMMITS_DIR}/"
        dst = f"{local_chemvcs_dir}/{COMMITS_DIR}/"
        self._rsync(src, dst)

    @staticmethod
    def _rsync(src: str, dst: str) -> None:
        """Run rsync with ``--ignore-existing`` to avoid re-transferring blobs."""
        subprocess.run(
            ["rsync", "-avz", "--ignore-existing", src, dst],
            check=True,
            timeout=300,
        )
