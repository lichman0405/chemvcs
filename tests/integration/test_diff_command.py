"""Integration tests for 'chemvcs diff' command."""

import os
from pathlib import Path

import pytest
from typer.testing import CliRunner

from chemvcs.cli.main import app

runner = CliRunner()


def _init_and_commit(tmp_path: Path, files: dict, message: str) -> str:
    """Helper: init repo, write files, add, commit. Returns commit hash."""
    runner.invoke(app, ["init", "--quiet"])
    for name, content in files.items():
        (tmp_path / name).write_text(content, encoding="utf-8")
        runner.invoke(app, ["add", name])
    result = runner.invoke(app, ["commit", "-m", message])
    assert result.exit_code == 0, f"commit failed:\n{result.stdout}"
    for line in result.stdout.splitlines():
        if "Committed" in line:
            return line.strip().split()[-1]
    raise AssertionError("Could not extract commit hash from output")


class TestDiffCommand:
    """Tests for the 'chemvcs diff' command."""

    def test_diff_no_commits(self, tmp_path: Path) -> None:
        """diff should fail gracefully when there are no commits yet."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        try:
            runner.invoke(app, ["init", "--quiet"])
            result = runner.invoke(app, ["diff"])
            assert result.exit_code == 1
            assert "No commits yet" in result.stdout
        finally:
            os.chdir(original_cwd)

    def test_diff_working_tree_no_changes(self, tmp_path: Path) -> None:
        """diff against working tree shows 'No changes' when files are unmodified."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        try:
            _init_and_commit(tmp_path, {"INCAR": "ENCUT = 520\nISMEAR = 0\n"}, "init")
            result = runner.invoke(app, ["diff"])
            assert result.exit_code == 0
            assert "No changes" in result.stdout
        finally:
            os.chdir(original_cwd)

    def test_diff_working_tree_detects_semantic_change(self, tmp_path: Path) -> None:
        """diff against working tree shows semantic changes for INCAR edits."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        try:
            _init_and_commit(
                tmp_path,
                {"INCAR": "ENCUT = 520\nISMEAR = 0\nSIGMA = 0.05\n"},
                "init",
            )
            # Modify INCAR on disk (change ENCUT — critical change)
            (tmp_path / "INCAR").write_text(
                "ENCUT = 600\nISMEAR = 0\nSIGMA = 0.05\n", encoding="utf-8"
            )
            result = runner.invoke(app, ["diff"])
            assert result.exit_code == 0
            assert "INCAR" in result.stdout
            assert "modified" in result.stdout
        finally:
            os.chdir(original_cwd)

    def test_diff_cosmetic_change_not_reported(self, tmp_path: Path) -> None:
        """Cosmetic change (comment/whitespace) should NOT appear as modified."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        try:
            _init_and_commit(
                tmp_path,
                {"INCAR": "ENCUT = 520\nISMEAR = 0\n"},
                "init",
            )
            # Only add a comment — semantically identical
            (tmp_path / "INCAR").write_text(
                "# This is a comment\nENCUT = 520\nISMEAR = 0\n", encoding="utf-8"
            )
            result = runner.invoke(app, ["diff"])
            assert result.exit_code == 0
            assert "No changes" in result.stdout
        finally:
            os.chdir(original_cwd)

    def test_diff_two_commits(self, tmp_path: Path) -> None:
        """diff between two explicit commit hashes shows changes."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        try:
            hash1 = _init_and_commit(
                tmp_path, {"INCAR": "ENCUT = 520\n"}, "first"
            )
            (tmp_path / "INCAR").write_text("ENCUT = 600\n", encoding="utf-8")
            runner.invoke(app, ["add", "INCAR"])
            result_c2 = runner.invoke(app, ["commit", "-m", "second"])
            assert result_c2.exit_code == 0
            hash2 = None
            for line in result_c2.stdout.splitlines():
                if "Committed" in line:
                    hash2 = line.strip().split()[-1]
                    break
            assert hash2 is not None

            result = runner.invoke(app, ["diff", hash2, hash1])
            assert result.exit_code == 0
            assert "INCAR" in result.stdout
        finally:
            os.chdir(original_cwd)

    def test_diff_specific_file(self, tmp_path: Path) -> None:
        """diff --file filters output to the named file only."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        try:
            _init_and_commit(
                tmp_path,
                {
                    "INCAR": "ENCUT = 520\n",
                    "KPOINTS": "Automatic\n0\nGamma\n4 4 4\n",
                },
                "init",
            )
            (tmp_path / "INCAR").write_text("ENCUT = 600\n", encoding="utf-8")
            (tmp_path / "KPOINTS").write_text("Automatic\n0\nGamma\n6 6 6\n", encoding="utf-8")

            result = runner.invoke(app, ["diff", "--file", "INCAR"])
            assert result.exit_code == 0
            assert "INCAR" in result.stdout
            assert "KPOINTS" not in result.stdout
        finally:
            os.chdir(original_cwd)

    def test_diff_unsupported_file_shows_contents_differ(self, tmp_path: Path) -> None:
        """Unsupported file types fall back to 'contents differ' message."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        try:
            _init_and_commit(
                tmp_path, {"notes.txt": "hello world\n"}, "init"
            )
            (tmp_path / "notes.txt").write_text("goodbye world\n", encoding="utf-8")
            result = runner.invoke(app, ["diff"])
            assert result.exit_code == 0
            assert "notes.txt" in result.stdout
            # The file is reported as modified (not parseable)
            assert "modified" in result.stdout
        finally:
            os.chdir(original_cwd)

    def test_diff_not_initialized(self, tmp_path: Path) -> None:
        """diff in non-initialized directory exits with error."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        try:
            result = runner.invoke(app, ["diff"])
            assert result.exit_code == 1
            assert "Not a ChemVCS repository" in result.stdout
        finally:
            os.chdir(original_cwd)

    def test_diff_summary_flag(self, tmp_path: Path) -> None:
        """diff --summary shows counts only, not full diff details."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        try:
            _init_and_commit(
                tmp_path,
                {"INCAR": "ENCUT = 520\nISMEAR = 0\n"},
                "init",
            )
            (tmp_path / "INCAR").write_text("ENCUT = 600\nISMEAR = 1\n", encoding="utf-8")
            result = runner.invoke(app, ["diff", "--summary"])
            assert result.exit_code == 0
            assert "INCAR" in result.stdout
        finally:
            os.chdir(original_cwd)
