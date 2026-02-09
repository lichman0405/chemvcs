"""Integration tests for chemvcs status command."""

import os
from pathlib import Path

import pytest
from typer.testing import CliRunner

from chemvcs.cli.main import app
from chemvcs.constants import CHEMVCS_DIR

runner = CliRunner()


class TestStatusCommand:
    """Test chemvcs status command."""

    def test_status_empty_repo(self, tmp_path: Path) -> None:
        """Test status on newly initialized repo."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        
        try:
            runner.invoke(app, ["init", "--quiet"])
            
            result = runner.invoke(app, ["status"])
            
            assert result.exit_code == 0
            assert "no commits yet" in result.stdout.lower()
            assert "no files staged" in result.stdout.lower() or "nothing to commit" in result.stdout.lower()
            
        finally:
            os.chdir(original_cwd)

    def test_status_with_staged_files(self, tmp_path: Path) -> None:
        """Test status showing staged files."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        
        try:
            runner.invoke(app, ["init", "--quiet"])
            
            # Create and stage files
            (tmp_path / "INCAR").write_text("ENCUT = 520\n")
            (tmp_path / "POSCAR").write_text("Fe\n1.0\n")
            
            runner.invoke(app, ["add", "INCAR", "POSCAR"])
            
            result = runner.invoke(app, ["status"])
            
            assert result.exit_code == 0
            assert "changes to be committed" in result.stdout.lower()
            assert "INCAR" in result.stdout
            assert "POSCAR" in result.stdout
            
        finally:
            os.chdir(original_cwd)

    def test_status_after_commit(self, tmp_path: Path) -> None:
        """Test status after making a commit."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        
        try:
            runner.invoke(app, ["init", "--quiet"])
            
            # Create, stage, and commit
            (tmp_path / "test.txt").write_text("test")
            runner.invoke(app, ["add", "test.txt"])
            runner.invoke(app, ["commit", "-m", "First commit"])
            
            result = runner.invoke(app, ["status"])
            
            assert result.exit_code == 0
            assert "nothing to commit" in result.stdout.lower() or "working tree clean" in result.stdout.lower()
            
            # Should show HEAD commit
            assert "HEAD:" in result.stdout
            
        finally:
            os.chdir(original_cwd)

    def test_status_short_format(self, tmp_path: Path) -> None:
        """Test status --short flag."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        
        try:
            runner.invoke(app, ["init", "--quiet"])
            
            # Stage a file
            (tmp_path / "file.txt").write_text("content")
            runner.invoke(app, ["add", "file.txt"])
            
            result = runner.invoke(app, ["status", "--short"])
            
            assert result.exit_code == 0
            # Short format should just show "A  file.txt"
            assert "A  file.txt" in result.stdout
            # Should not have verbose descriptions
            assert "changes to be committed" not in result.stdout.lower()
            
        finally:
            os.chdir(original_cwd)

    def test_status_shows_file_details(self, tmp_path: Path) -> None:
        """Test that status shows file type and size."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        
        try:
            runner.invoke(app, ["init", "--quiet"])
            
            # Create an INCAR file
            (tmp_path / "INCAR").write_text("ENCUT = 520\n")
            runner.invoke(app, ["add", "INCAR"])
            
            result = runner.invoke(app, ["status"])
            
            assert result.exit_code == 0
            assert "INCAR" in result.stdout
            # Should show file type
            assert "INCAR" in result.stdout  # file type INCAR
            # Should show size
            assert "B" in result.stdout or "KB" in result.stdout
            
        finally:
            os.chdir(original_cwd)

    def test_status_not_initialized(self, tmp_path: Path) -> None:
        """Test status in non-initialized directory."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        
        try:
            result = runner.invoke(app, ["status"])
            
            assert result.exit_code == 1
            assert "not a chemvcs repository" in result.stdout.lower()
            
        finally:
            os.chdir(original_cwd)

    def test_status_multiple_files_sorted(self, tmp_path: Path) -> None:
        """Test that status shows files in sorted order."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        
        try:
            runner.invoke(app, ["init", "--quiet"])
            
            # Create files in non-alphabetical order
            (tmp_path / "zebra.txt").write_text("z")
            (tmp_path / "apple.txt").write_text("a")
            (tmp_path / "mango.txt").write_text("m")
            
            runner.invoke(app, ["add", "zebra.txt", "apple.txt", "mango.txt"])
            
            result = runner.invoke(app, ["status"])
            
            assert result.exit_code == 0
            
            # Check that files appear in sorted order
            lines = result.stdout.split("\n")
            file_lines = [l for l in lines if any(f in l for f in ["apple", "mango", "zebra"])]
            
            # Extract file names
            found_files = []
            for line in file_lines:
                if "apple" in line:
                    found_files.append("apple")
                elif "mango" in line:
                    found_files.append("mango")
                elif "zebra" in line:
                    found_files.append("zebra")
            
            assert found_files == ["apple", "mango", "zebra"]
            
        finally:
            os.chdir(original_cwd)
