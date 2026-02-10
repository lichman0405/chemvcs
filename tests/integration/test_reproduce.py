"""Integration tests for chemvcs reproduce command."""

import os
from pathlib import Path

import pytest
from typer.testing import CliRunner

from chemvcs.cli.main import app

runner = CliRunner()


class TestReproduceCommand:
    """Test chemvcs reproduce command."""

    def test_reproduce_basic(self, tmp_path: Path) -> None:
        """Test basic reproduce functionality."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        
        try:
            runner.invoke(app, ["init", "--quiet"])
            
            # Create and commit a file
            (tmp_path / "test.txt").write_text("test content")
            runner.invoke(app, ["add", "test.txt"])
            commit_result = runner.invoke(app, ["commit", "-m", "Test commit"])
            
            # Extract commit hash from output (format: "> Committed <hash>")
            for line in commit_result.stdout.split("\n"):
                if "> Committed" in line:
                    commit_hash = line.split()[-1]
                    break
            
            # Reproduce the commit
            result = runner.invoke(app, ["reproduce", commit_hash])
            
            assert result.exit_code == 0
            assert "Reproducing commit" in result.stdout
            assert f"reproduce_{commit_hash[:7]}" in result.stdout
            
            # Check that file was reproduced
            reproduced_dir = tmp_path / f"reproduce_{commit_hash[:7]}"
            assert reproduced_dir.exists()
            assert (reproduced_dir / "test.txt").exists()
            assert (reproduced_dir / "test.txt").read_text() == "test content"
            
        finally:
            os.chdir(original_cwd)

    def test_reproduce_multiple_files(self, tmp_path: Path) -> None:
        """Test reproducing commit with multiple files."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        
        try:
            runner.invoke(app, ["init", "--quiet"])
            
            # Create multiple files
            (tmp_path / "file1.txt").write_text("content1")
            (tmp_path / "file2.txt").write_text("content2")
            (tmp_path / "file3.txt").write_text("content3")
            
            runner.invoke(app, ["add", "file1.txt", "file2.txt", "file3.txt"])
            commit_result = runner.invoke(app, ["commit", "-m", "Multi-file commit"])
            
            # Extract commit hash
            for line in commit_result.stdout.split("\n"):
                if "> Committed" in line:
                    commit_hash = line.split()[-1]
                    break
            
            # Reproduce
            result = runner.invoke(app, ["reproduce", commit_hash])
            
            assert result.exit_code == 0
            
            # Check all files reproduced
            reproduced_dir = tmp_path / f"reproduce_{commit_hash[:7]}"
            assert (reproduced_dir / "file1.txt").read_text() == "content1"
            assert (reproduced_dir / "file2.txt").read_text() == "content2"
            assert (reproduced_dir / "file3.txt").read_text() == "content3"
            
        finally:
            os.chdir(original_cwd)

    def test_reproduce_custom_output_dir(self, tmp_path: Path) -> None:
        """Test reproduce with custom output directory."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        
        try:
            runner.invoke(app, ["init", "--quiet"])
            
            (tmp_path / "file.txt").write_text("content")
            runner.invoke(app, ["add", "file.txt"])
            commit_result = runner.invoke(app, ["commit", "-m", "Test"])
            
            # Extract commit hash
            for line in commit_result.stdout.split("\n"):
                if "> Committed" in line:
                    commit_hash = line.split()[-1]
                    break
            
            # Reproduce to custom directory
            custom_dir = "custom_output"
            result = runner.invoke(app, ["reproduce", commit_hash, "--output-dir", custom_dir])
            
            assert result.exit_code == 0
            assert custom_dir in result.stdout
            
            # Check file in custom directory
            assert (tmp_path / custom_dir / "file.txt").exists()
            assert (tmp_path / custom_dir / "file.txt").read_text() == "content"
            
        finally:
            os.chdir(original_cwd)

    def test_reproduce_short_hash(self, tmp_path: Path) -> None:
        """Test reproducing with short hash prefix."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        
        try:
            runner.invoke(app, ["init", "--quiet"])
            
            (tmp_path / "file.txt").write_text("content")
            runner.invoke(app, ["add", "file.txt"])
            commit_result = runner.invoke(app, ["commit", "-m", "Test"])
            
            # Extract commit hash
            for line in commit_result.stdout.split("\n"):
                if "> Committed" in line:
                    full_hash = line.split()[-1]
                    short_hash = full_hash[:7]
                    break
            
            # Reproduce with short hash
            result = runner.invoke(app, ["reproduce", short_hash])
            
            assert result.exit_code == 0
            
            # Check file reproduced
            reproduced_dir = tmp_path / f"reproduce_{short_hash}"
            assert (reproduced_dir / "file.txt").exists()
            
        finally:
            os.chdir(original_cwd)

    def test_reproduce_nonexistent_commit(self, tmp_path: Path) -> None:
        """Test reproducing non-existent commit fails."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        
        try:
            runner.invoke(app, ["init", "--quiet"])
            
            # Try to reproduce fake commit
            result = runner.invoke(app, ["reproduce", "nonexistent"])
            
            assert result.exit_code == 1
            assert "not found" in result.stdout.lower()
            
        finally:
            os.chdir(original_cwd)

    def test_reproduce_not_initialized(self, tmp_path: Path) -> None:
        """Test reproduce in non-initialized repo."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        
        try:
            result = runner.invoke(app, ["reproduce", "somehash"])
            
            assert result.exit_code == 1
            assert "not a chemvcs repository" in result.stdout.lower()
            
        finally:
            os.chdir(original_cwd)

    def test_reproduce_preserves_directory_structure(self, tmp_path: Path) -> None:
        """Test that reproduce preserves directory structure."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        
        try:
            runner.invoke(app, ["init", "--quiet"])
            
            # Create files in subdirectories
            (tmp_path / "subdir").mkdir()
            (tmp_path / "subdir" / "file.txt").write_text("nested content")
            (tmp_path / "root.txt").write_text("root content")
            
            runner.invoke(app, ["add", "subdir/file.txt", "root.txt"])
            commit_result = runner.invoke(app, ["commit", "-m", "Nested files"])
            
            # Extract commit hash
            for line in commit_result.stdout.split("\n"):
                if "> Committed" in line:
                    commit_hash = line.split()[-1]
                    break
            
            # Reproduce
            result = runner.invoke(app, ["reproduce", commit_hash])
            
            assert result.exit_code == 0
            
            # Check directory structure preserved
            reproduced_dir = tmp_path / f"reproduce_{commit_hash[:7]}"
            assert (reproduced_dir / "subdir" / "file.txt").exists()
            assert (reproduced_dir / "subdir" / "file.txt").read_text() == "nested content"
            assert (reproduced_dir / "root.txt").read_text() == "root content"
            
        finally:
            os.chdir(original_cwd)

    def test_reproduce_shows_file_count(self, tmp_path: Path) -> None:
        """Test that reproduce shows file count in output."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        
        try:
            runner.invoke(app, ["init", "--quiet"])
            
            # Create 3 files
            for i in range(3):
                (tmp_path / f"file{i}.txt").write_text(f"content{i}")
            
            runner.invoke(app, ["add", "file0.txt", "file1.txt", "file2.txt"])
            commit_result = runner.invoke(app, ["commit", "-m", "Three files"])
            
            # Extract commit hash
            for line in commit_result.stdout.split("\n"):
                if "> Committed" in line:
                    commit_hash = line.split()[-1]
                    break
            
            # Reproduce
            result = runner.invoke(app, ["reproduce", commit_hash])
            
            assert result.exit_code == 0
            assert "3 file(s)" in result.stdout
            assert "Reproduced successfully" in result.stdout
            
        finally:
            os.chdir(original_cwd)
