"""Integration tests for chemvcs log command."""

import os
from pathlib import Path

import pytest
from typer.testing import CliRunner

from chemvcs.cli.main import app

runner = CliRunner()


class TestLogCommand:
    """Test chemvcs log command."""

    def test_log_empty_repo(self, tmp_path: Path) -> None:
        """Test log on newly initialized repo with no commits."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        
        try:
            runner.invoke(app, ["init", "--quiet"])
            
            result = runner.invoke(app, ["log"])
            
            assert result.exit_code == 0
            assert "no commits yet" in result.stdout.lower()
            
        finally:
            os.chdir(original_cwd)

    def test_log_single_commit(self, tmp_path: Path) -> None:
        """Test log with a single commit."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        
        try:
            runner.invoke(app, ["init", "--quiet"])
            
            # Create and commit
            (tmp_path / "file.txt").write_text("content")
            runner.invoke(app, ["add", "file.txt"])
            runner.invoke(app, ["commit", "-m", "First commit"])
            
            result = runner.invoke(app, ["log"])
            
            assert result.exit_code == 0
            assert "commit" in result.stdout.lower()
            assert "First commit" in result.stdout
            assert "Author:" in result.stdout
            assert "Date:" in result.stdout
            assert "(root commit)" in result.stdout
            
        finally:
            os.chdir(original_cwd)

    def test_log_multiple_commits(self, tmp_path: Path) -> None:
        """Test log with multiple commits."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        
        try:
            runner.invoke(app, ["init", "--quiet"])
            
            # First commit
            (tmp_path / "file1.txt").write_text("content1")
            runner.invoke(app, ["add", "file1.txt"])
            runner.invoke(app, ["commit", "-m", "First commit"])
            
            # Second commit
            (tmp_path / "file2.txt").write_text("content2")
            runner.invoke(app, ["add", "file2.txt"])
            runner.invoke(app, ["commit", "-m", "Second commit"])
            
            # Third commit
            (tmp_path / "file3.txt").write_text("content3")
            runner.invoke(app, ["add", "file3.txt"])
            runner.invoke(app, ["commit", "-m", "Third commit"])
            
            result = runner.invoke(app, ["log"])
            
            assert result.exit_code == 0
            # Should show all commits in reverse chronological order
            assert "First commit" in result.stdout
            assert "Second commit" in result.stdout
            assert "Third commit" in result.stdout
            
            # Third commit should appear before First commit
            assert result.stdout.index("Third commit") < result.stdout.index("First commit")
            
        finally:
            os.chdir(original_cwd)

    def test_log_max_count(self, tmp_path: Path) -> None:
        """Test log --max-count option."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        
        try:
            runner.invoke(app, ["init", "--quiet"])
            
            # Create 5 commits
            for i in range(1, 6):
                (tmp_path / f"file{i}.txt").write_text(f"content{i}")
                runner.invoke(app, ["add", f"file{i}.txt"])
                runner.invoke(app, ["commit", "-m", f"Commit {i}"])
            
            # Get only 2 commits
            result = runner.invoke(app, ["log", "--max-count", "2"])
            
            assert result.exit_code == 0
            # Should show commits 5 and 4 (most recent)
            assert "Commit 5" in result.stdout
            assert "Commit 4" in result.stdout
            # Should not show older commits
            assert "Commit 3" not in result.stdout
            assert "Commit 1" not in result.stdout
            
        finally:
            os.chdir(original_cwd)

    def test_log_oneline_format(self, tmp_path: Path) -> None:
        """Test log --oneline option."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        
        try:
            runner.invoke(app, ["init", "--quiet"])
            
            # Create commits
            (tmp_path / "file.txt").write_text("content")
            runner.invoke(app, ["add", "file.txt"])
            runner.invoke(app, ["commit", "-m", "Test commit"])
            
            result = runner.invoke(app, ["log", "--oneline"])
            
            assert result.exit_code == 0
            # Oneline should be compact
            assert "Test commit" in result.stdout
            # Should NOT have detailed formatting
            assert "Author:" not in result.stdout
            assert "Date:" not in result.stdout
            
        finally:
            os.chdir(original_cwd)

    def test_log_n_shorthand(self, tmp_path: Path) -> None:
        """Test log -n (shorthand for --max-count)."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        
        try:
            runner.invoke(app, ["init", "--quiet"])
            
            # Create 3 commits
            for i in range(1, 4):
                (tmp_path / f"file{i}.txt").write_text(f"content{i}")
                runner.invoke(app, ["add", f"file{i}.txt"])
                runner.invoke(app, ["commit", "-m", f"Commit {i}"])
            
            result = runner.invoke(app, ["log", "-n", "1"])
            
            assert result.exit_code == 0
            assert "Commit 3" in result.stdout
            assert "Commit 2" not in result.stdout
            
        finally:
            os.chdir(original_cwd)

    def test_log_not_initialized(self, tmp_path: Path) -> None:
        """Test log in non-initialized directory."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        
        try:
            result = runner.invoke(app, ["log"])
            
            assert result.exit_code == 1
            assert "not a chemvcs repository" in result.stdout.lower()
            
        finally:
            os.chdir(original_cwd)

    def test_log_multiline_message(self, tmp_path: Path) -> None:
        """Test log with multiline commit message."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        
        try:
            runner.invoke(app, ["init", "--quiet"])
            
            # Commit with multiline message
            (tmp_path / "file.txt").write_text("content")
            runner.invoke(app, ["add", "file.txt"])
            multiline_msg = "First line\n\nDetailed description\nwith multiple lines"
            runner.invoke(app, ["commit", "-m", multiline_msg])
            
            # Default format should show full message
            result = runner.invoke(app, ["log"])
            assert result.exit_code == 0
            assert "First line" in result.stdout
            assert "Detailed description" in result.stdout
            
            # Oneline should show only first line
            result_oneline = runner.invoke(app, ["log", "--oneline"])
            assert result_oneline.exit_code == 0
            assert "First line" in result_oneline.stdout
            # Second paragraph should not appear
            assert "Detailed description" not in result_oneline.stdout
            
        finally:
            os.chdir(original_cwd)
