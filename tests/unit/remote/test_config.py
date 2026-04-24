"""Unit tests for remote/config.py."""

from __future__ import annotations

from pathlib import Path

import pytest

from chemvcs.constants import CHEMVCS_DIR
from chemvcs.remote.config import (
    RemoteAlreadyExistsError,
    RemoteNotFoundError,
    _remotes_path,
    add_remote,
    get_remote_url,
    list_remotes,
    remove_remote,
)


class TestListRemotes:
    def test_empty_when_no_file(self, tmp_path: Path) -> None:
        chemvcs = tmp_path / CHEMVCS_DIR
        chemvcs.mkdir()
        assert list_remotes(chemvcs) == {}

    def test_returns_mapping(self, tmp_path: Path) -> None:
        chemvcs = tmp_path / CHEMVCS_DIR
        chemvcs.mkdir()
        add_remote(chemvcs, "origin", "user@host:/path/repo")
        result = list_remotes(chemvcs)
        assert result == {"origin": "user@host:/path/repo"}

    def test_multiple_remotes(self, tmp_path: Path) -> None:
        chemvcs = tmp_path / CHEMVCS_DIR
        chemvcs.mkdir()
        add_remote(chemvcs, "origin", "user@host:/path/repo")
        add_remote(chemvcs, "backup", "user@host2:/backup/repo")
        result = list_remotes(chemvcs)
        assert result == {
            "origin": "user@host:/path/repo",
            "backup": "user@host2:/backup/repo",
        }


class TestAddRemote:
    def test_add_new_remote(self, tmp_path: Path) -> None:
        chemvcs = tmp_path / CHEMVCS_DIR
        chemvcs.mkdir()
        add_remote(chemvcs, "origin", "user@host:/path/repo")
        assert (chemvcs / "remotes.toml").exists()
        assert get_remote_url(chemvcs, "origin") == "user@host:/path/repo"

    def test_duplicate_raises(self, tmp_path: Path) -> None:
        chemvcs = tmp_path / CHEMVCS_DIR
        chemvcs.mkdir()
        add_remote(chemvcs, "origin", "user@host:/path/repo")
        with pytest.raises(RemoteAlreadyExistsError, match="origin"):
            add_remote(chemvcs, "origin", "other@host:/other")


class TestRemoveRemote:
    def test_remove_existing(self, tmp_path: Path) -> None:
        chemvcs = tmp_path / CHEMVCS_DIR
        chemvcs.mkdir()
        add_remote(chemvcs, "origin", "user@host:/path/repo")
        add_remote(chemvcs, "backup", "user@host2:/backup")
        remove_remote(chemvcs, "origin")
        assert list_remotes(chemvcs) == {"backup": "user@host2:/backup"}

    def test_remove_nonexistent_raises(self, tmp_path: Path) -> None:
        chemvcs = tmp_path / CHEMVCS_DIR
        chemvcs.mkdir()
        with pytest.raises(RemoteNotFoundError, match="ghost"):
            remove_remote(chemvcs, "ghost")

    def test_remove_last_remote(self, tmp_path: Path) -> None:
        chemvcs = tmp_path / CHEMVCS_DIR
        chemvcs.mkdir()
        add_remote(chemvcs, "origin", "user@host:/path/repo")
        remove_remote(chemvcs, "origin")
        assert list_remotes(chemvcs) == {}


class TestGetRemoteURL:
    def test_returns_url(self, tmp_path: Path) -> None:
        chemvcs = tmp_path / CHEMVCS_DIR
        chemvcs.mkdir()
        add_remote(chemvcs, "origin", "user@host:/path/repo")
        assert get_remote_url(chemvcs, "origin") == "user@host:/path/repo"

    def test_nonexistent_raises(self, tmp_path: Path) -> None:
        chemvcs = tmp_path / CHEMVCS_DIR
        chemvcs.mkdir()
        with pytest.raises(RemoteNotFoundError, match="ghost"):
            get_remote_url(chemvcs, "ghost")


class TestRemotesPath:
    def test_returns_correct_path(self, tmp_path: Path) -> None:
        chemvcs = tmp_path / CHEMVCS_DIR
        chemvcs.mkdir()
        assert _remotes_path(chemvcs) == chemvcs / "remotes.toml"
