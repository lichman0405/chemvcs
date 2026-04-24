"""Unit tests for remote/manager.py.

Uses ``unittest.mock.patch("subprocess.run")`` to verify SSH/rsync
command construction without making real network calls.
"""

from __future__ import annotations

import subprocess
from pathlib import Path
from unittest.mock import patch

import pytest

from chemvcs.remote.manager import RemoteManager


class TestRemoteManagerInit:
    def test_parses_user_host_path(self) -> None:
        mgr = RemoteManager("alice@server.example.com:/data/repo")
        assert mgr.user_host == "alice@server.example.com"
        assert mgr.remote_path == "/data/repo"

    def test_defaults_to_current_user(self) -> None:
        mgr = RemoteManager("server.example.com:/data/repo")
        assert mgr.user_host.endswith("@server.example.com")
        assert mgr.remote_path == "/data/repo"

    def test_rejects_invalid_url(self) -> None:
        with pytest.raises(ValueError, match="Invalid remote URL"):
            RemoteManager("no-colon-here")


class TestCheckConnection:
    def test_success(self) -> None:
        mgr = RemoteManager("user@host:/path")
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = subprocess.CompletedProcess(
                args=[], returncode=0, stdout="ok\n", stderr=""
            )
            assert mgr.check_connection() is True
            mock_run.assert_called_once()
            args = mock_run.call_args[0][0]
            assert "ssh" in args
            assert "user@host" in args

    def test_failure(self) -> None:
        mgr = RemoteManager("user@host:/path")
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = subprocess.CompletedProcess(
                args=[], returncode=255, stdout="", stderr="Permission denied"
            )
            assert mgr.check_connection() is False

    def test_timeout(self) -> None:
        mgr = RemoteManager("user@host:/path")
        with patch("subprocess.run", side_effect=subprocess.TimeoutExpired("ssh", 15)):
            assert mgr.check_connection() is False


class TestGetRemoteHead:
    def test_returns_hash(self) -> None:
        mgr = RemoteManager("user@host:/path")
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = subprocess.CompletedProcess(
                args=[],
                returncode=0,
                stdout="abc123def456\n",
                stderr="",
            )
            assert mgr.get_remote_head() == "abc123def456"

    def test_returns_none_when_empty(self) -> None:
        mgr = RemoteManager("user@host:/path")
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = subprocess.CompletedProcess(
                args=[], returncode=0, stdout="\n", stderr=""
            )
            assert mgr.get_remote_head() is None

    def test_returns_none_when_missing(self) -> None:
        mgr = RemoteManager("user@host:/path")
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = subprocess.CompletedProcess(
                args=[], returncode=1, stdout="", stderr="cat: No such file"
            )
            assert mgr.get_remote_head() is None


class TestInitRemoteIfNeeded:
    def test_creates_directories(self) -> None:
        mgr = RemoteManager("user@host:/path")
        with patch("subprocess.run") as mock_run:
            mgr.init_remote_if_needed()
            mock_run.assert_called_once()
            args = mock_run.call_args[0][0]
            assert "mkdir" in args
            assert "-p" in args


class TestUpdateRemoteHead:
    def test_writes_head_file(self) -> None:
        mgr = RemoteManager("user@host:/path")
        with patch("subprocess.run") as mock_run:
            mgr.update_remote_head("abc123")
            mock_run.assert_called_once()
            args = mock_run.call_args[0][0]
            assert "abc123" in str(args)


class TestRsyncOperations:
    def test_push_objects_calls_rsync(self) -> None:
        mgr = RemoteManager("user@host:/path")
        local = Path("/tmp/dot_chemvcs")
        with patch("subprocess.run") as mock_run:
            mgr.push_objects(local)
            mock_run.assert_called_once()
            call_args = mock_run.call_args[0][0]
            assert call_args[0] == "rsync"
            assert "--ignore-existing" in call_args
            # Source: local objects/
            assert str(local / "objects") in call_args[3]
            # Dest: remote objects/
            assert "objects" in call_args[4]

    def test_push_commits_calls_rsync(self) -> None:
        mgr = RemoteManager("user@host:/path")
        local = Path("/tmp/dot_chemvcs")
        with patch("subprocess.run") as mock_run:
            mgr.push_commits(local)
            call_args = mock_run.call_args[0][0]
            assert "rsync" in str(call_args)
            assert "--ignore-existing" in call_args

    def test_fetch_objects_calls_rsync(self) -> None:
        mgr = RemoteManager("user@host:/path")
        local = Path("/tmp/dot_chemvcs")
        with patch("subprocess.run") as mock_run:
            mgr.fetch_objects(local)
            call_args = mock_run.call_args[0][0]
            assert "rsync" in str(call_args)
            assert "--ignore-existing" in call_args

    def test_fetch_commits_calls_rsync(self) -> None:
        mgr = RemoteManager("user@host:/path")
        local = Path("/tmp/dot_chemvcs")
        with patch("subprocess.run") as mock_run:
            mgr.fetch_commits(local)
            call_args = mock_run.call_args[0][0]
            assert "rsync" in str(call_args)
            assert "--ignore-existing" in call_args
