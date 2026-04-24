"""Integration tests for push/pull workflow.

Simulates two local directories: one "remote" and one "local client",
using SSH/rsync mocks that replace network operations with local file copies.
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path
from typing import Any
from unittest.mock import patch

from typer.testing import CliRunner

from chemvcs.cli.main import app
from chemvcs.constants import COMMITS_DIR, OBJECTS_DIR
from chemvcs.remote.config import add_remote
from chemvcs.storage import ObjectStore


def _make_commit(
    repo_dir: Path,
    message: str,
    content: dict[str, bytes],
    parent: str | None = None,
) -> str:
    """Create a minimal commit file and return its hash."""
    import hashlib

    chemvcs = repo_dir / ".chemvcs"
    chemvcs.mkdir(parents=True, exist_ok=True)
    (chemvcs / OBJECTS_DIR).mkdir(parents=True, exist_ok=True)
    (chemvcs / COMMITS_DIR).mkdir(parents=True, exist_ok=True)

    store = ObjectStore(chemvcs)
    files: list[dict[str, Any]] = []
    for fname, fcontent in content.items():
        blob_hash = store.write_blob(fcontent)
        files.append(
            {
                "path": fname,
                "blob_hash": blob_hash,
                "file_type": "UNKNOWN",
                "size_bytes": len(fcontent),
                "is_reference": False,
            }
        )

    commit_obj: dict[str, Any] = {
        "parent": parent,
        "timestamp": "2026-04-24T00:00:00+00:00",
        "author": "test@localhost",
        "message": message,
        "files": files,
    }
    canonical = json.dumps(commit_obj, sort_keys=True, separators=(",", ":"))
    commit_hash = hashlib.sha256(canonical.encode()).hexdigest()
    commit_obj["hash"] = commit_hash

    commit_path = chemvcs / COMMITS_DIR / commit_hash
    commit_path.write_text(json.dumps(commit_obj, indent=2), encoding="utf-8")

    (chemvcs / "HEAD").write_text(commit_hash, encoding="utf-8")

    return commit_hash


def _build_remote_repo(remote_root: Path) -> tuple[str, str]:
    """Set up a remote repo with two commits. Returns (root_hash, tip_hash)."""
    c1 = _make_commit(remote_root, "Initial commit", {"README.md": b"# My Calc"})
    c2 = _make_commit(
        remote_root,
        "Add INCAR",
        {"README.md": b"# My Calc", "INCAR": b"ENCUT=520"},
        parent=c1,
    )
    return c1, c2


# -------------------- mock helpers --------------------


def _make_ssh_mock(remote_root: Path):
    """Build a callable that mimics SSH for the given remote_root."""

    def _fake_ssh(args, **kwargs):
        import subprocess as sp

        # Find the remote command after the ssh flags and user@host
        # args layout: ssh [-o flag]... user@host cmd [cmd_args...]
        cmd_idx = None
        for i, a in enumerate(args):
            if a in ("echo", "cat", "mkdir", "sh"):
                cmd_idx = i
                break

        if cmd_idx is None:
            return sp.CompletedProcess(args, 0, "", "")

        cmd = args[cmd_idx]

        if cmd == "echo":
            return sp.CompletedProcess(args, 0, "ok\n", "")

        if cmd == "cat":
            head_path = args[cmd_idx + 1]
            filepath = Path(head_path)
            if filepath.exists():
                return sp.CompletedProcess(
                    args, 0, filepath.read_text(encoding="utf-8"), ""
                )
            return sp.CompletedProcess(args, 1, "", "")

        if cmd == "mkdir":
            for p in args[cmd_idx + 1 :]:
                Path(p).mkdir(parents=True, exist_ok=True)
            return sp.CompletedProcess(args, 0, "", "")

        if cmd == "sh":
            actual_cmd = args[cmd_idx + 2]  # skip "-c"
            if ">" in actual_cmd:
                _, outpath = actual_cmd.split(">", 1)
                outpath = outpath.strip()
                content = actual_cmd.split("printf '%s' ")[1].split(" >")[0]
                Path(outpath).parent.mkdir(parents=True, exist_ok=True)
                Path(outpath).write_text(content, encoding="utf-8")
            return sp.CompletedProcess(args, 0, "", "")

        return sp.CompletedProcess(args, 0, "", "")

    return _fake_ssh


def _make_rsync_mock(remote_root: Path):
    """Build a callable that mimics rsync using local copies."""

    def _fake_rsync(args, **kwargs):
        import subprocess as sp

        # rsync -avz --ignore-existing src dst
        src = args[3].rstrip("/")
        dst = args[4].rstrip("/")

        # Resolve remote paths: "user@host:/path" → "/path"
        if ":" in src:
            src = src.split(":", 1)[1]
        if ":" in dst:
            dst = dst.split(":", 1)[1]

        src_path = Path(src)
        dst_path = Path(dst)

        if not src_path.exists():
            return sp.CompletedProcess(args, 0, "", "")

        dst_path.mkdir(parents=True, exist_ok=True)
        for item in src_path.iterdir():
            dst_item = dst_path / item.name
            if not dst_item.exists():
                if item.is_dir():
                    shutil.copytree(item, dst_item, dirs_exist_ok=True)
                else:
                    shutil.copy2(item, dst_item)
        return sp.CompletedProcess(args, 0, "", "")

    return _fake_rsync


def _dispatch_mock(remote_root: Path):
    """Return a single `side_effect` that routes SSH/rsync calls."""
    ssh = _make_ssh_mock(remote_root)
    rsync = _make_rsync_mock(remote_root)

    def _dispatch(args, **kwargs):
        import subprocess as sp

        if args[0] == "ssh":
            return ssh(args, **kwargs)
        if args[0] == "rsync":
            return rsync(args, **kwargs)
        return sp.CompletedProcess(args, 0, "", "")

    return _dispatch


# -------------------- tests --------------------


class TestPushPullWorkflow:
    def test_full_push_pull_cycle(self, tmp_path: Path) -> None:
        """Simulate push → pull with two local directories."""
        remote_root = tmp_path / "remote_repo"
        remote_root.mkdir()
        _build_remote_repo(remote_root)

        local_root = tmp_path / "local_client"
        local_chemvcs = local_root / ".chemvcs"
        local_chemvcs.mkdir(parents=True)
        (local_chemvcs / OBJECTS_DIR).mkdir(parents=True, exist_ok=True)
        (local_chemvcs / COMMITS_DIR).mkdir(parents=True, exist_ok=True)
        add_remote(local_chemvcs, "origin", f"user@host:{remote_root}")

        dispatch = _dispatch_mock(remote_root)

        import chemvcs.remote.manager as mgr_module

        with patch.object(mgr_module.subprocess, "run", side_effect=dispatch):
            runner = CliRunner()
            import os

            old_cwd = os.getcwd()
            try:
                os.chdir(local_root)
                result = runner.invoke(app, ["pull", "origin"])
                assert result.exit_code == 0, f"pull failed: {result.output}"
                assert "Pulled" in result.output

                # HEAD should now match remote
                head = (
                    (local_root / ".chemvcs/HEAD")
                    .read_text(encoding="utf-8")
                    .strip()
                )
                remote_head = (
                    (remote_root / ".chemvcs/HEAD")
                    .read_text(encoding="utf-8")
                    .strip()
                )
                assert head == remote_head

                # Already up to date
                result2 = runner.invoke(app, ["pull", "origin"])
                assert "Already up to date" in result2.output

            finally:
                os.chdir(old_cwd)

    def test_push_fast_forward_accepted(self, tmp_path: Path) -> None:
        """Push a new commit on top of remote HEAD (fast-forward)."""
        remote_root = tmp_path / "remote_repo"
        remote_root.mkdir()
        remote_chemvcs = remote_root / ".chemvcs"
        remote_chemvcs.mkdir(parents=True)
        (remote_chemvcs / OBJECTS_DIR).mkdir(parents=True, exist_ok=True)
        (remote_chemvcs / COMMITS_DIR).mkdir(parents=True, exist_ok=True)
        # No HEAD file → simulates empty remote (first push)

        local_root = tmp_path / "local_client"
        local_chemvcs = local_root / ".chemvcs"
        local_chemvcs.mkdir(parents=True)
        (local_chemvcs / OBJECTS_DIR).mkdir(parents=True, exist_ok=True)
        (local_chemvcs / COMMITS_DIR).mkdir(parents=True, exist_ok=True)
        add_remote(local_chemvcs, "origin", f"user@host:{remote_root}")

        c1 = _make_commit(
            local_root, "Initial commit", {"README.md": b"# Local"}
        )
        (local_chemvcs / "HEAD").write_text(c1, encoding="utf-8")

        dispatch = _dispatch_mock(remote_root)

        import chemvcs.remote.manager as mgr_module

        with patch.object(mgr_module.subprocess, "run", side_effect=dispatch):
            runner = CliRunner()
            import os

            old_cwd = os.getcwd()
            try:
                os.chdir(local_root)
                result = runner.invoke(app, ["push", "origin"])
                assert result.exit_code == 0, f"push failed: {result.output}"
            finally:
                os.chdir(old_cwd)

    def test_remote_empty_on_pull(self, tmp_path: Path) -> None:
        """Pull from a remote with no commits."""
        remote_root = tmp_path / "remote_repo"
        remote_root.mkdir()
        remote_chemvcs = remote_root / ".chemvcs"
        remote_chemvcs.mkdir(parents=True)

        local_root = tmp_path / "local_client"
        local_chemvcs = local_root / ".chemvcs"
        local_chemvcs.mkdir(parents=True)
        (local_chemvcs / OBJECTS_DIR).mkdir(parents=True, exist_ok=True)
        (local_chemvcs / COMMITS_DIR).mkdir(parents=True, exist_ok=True)
        add_remote(local_chemvcs, "origin", f"user@host:{remote_root}")

        dispatch = _dispatch_mock(remote_root)

        import chemvcs.remote.manager as mgr_module

        with patch.object(mgr_module.subprocess, "run", side_effect=dispatch):
            runner = CliRunner()
            import os

            old_cwd = os.getcwd()
            try:
                os.chdir(local_root)
                result = runner.invoke(app, ["pull", "origin"])
                assert result.exit_code != 0
                assert "empty" in result.output.lower()
            finally:
                os.chdir(old_cwd)

    def test_remote_commands(self, tmp_path: Path) -> None:
        """Test remote add/list/remove commands."""
        local_root = tmp_path / "local_client"
        local_chemvcs = local_root / ".chemvcs"
        local_chemvcs.mkdir(parents=True)

        runner = CliRunner()
        import os

        old_cwd = os.getcwd()
        try:
            os.chdir(local_root)

            # add
            result = runner.invoke(
                app, ["remote", "add", "origin", "user@host:/path/repo"]
            )
            assert result.exit_code == 0
            assert "origin" in result.output

            # list
            result = runner.invoke(app, ["remote", "list"])
            assert result.exit_code == 0
            assert "origin" in result.output
            assert "user@host:/path/repo" in result.output

            # duplicate
            result = runner.invoke(
                app, ["remote", "add", "origin", "other@host:/other"]
            )
            assert result.exit_code != 0

            # remove
            result = runner.invoke(app, ["remote", "remove", "origin"])
            assert result.exit_code == 0

            # list empty
            result = runner.invoke(app, ["remote", "list"])
            assert "No remotes" in result.output

        finally:
            os.chdir(old_cwd)
