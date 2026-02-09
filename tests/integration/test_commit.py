"""Integration tests for chemvcs commit command."""

import os
from pathlib import Path

import pytest
from typer.testing import CliRunner

from chemvcs.cli.main import app
from chemvcs.constants import CHEMVCS_DIR

runner = CliRunner()


class TestCommitCommand:
    """Test chemvcs commit command."""

    def test_commit_basic(self, tmp_path: Path) -> None:
        """Test basic commit workflow."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        
        try:
            # Initialize repo
            runner.invoke(app, ["init", "--quiet"])
            
            # Create and add a file
            test_file = tmp_path / "INCAR"
            test_file.write_text("ENCUT = 520\n")
            
            result_add = runner.invoke(app, ["add", "INCAR"])
            assert result_add.exit_code == 0
            
            # Commit
            result = runner.invoke(app, ["commit", "-m", "Initial commit"])
            
            assert result.exit_code == 0
            assert "Committed" in result.stdout
            assert "Initial commit" in result.stdout
            
            # Verify HEAD was updated
            head_file = tmp_path / CHEMVCS_DIR / "HEAD"
            assert head_file.exists()
            commit_hash = head_file.read_text().strip()
            assert len(commit_hash) == 64
            
            # Verify staging area was cleared
            index_file = tmp_path / CHEMVCS_DIR / "index"
            import json
            index_data = json.loads(index_file.read_text())
            assert len(index_data["entries"]) == 0
            
        finally:
            os.chdir(original_cwd)

    def test_commit_requires_message(self, tmp_path: Path) -> None:
        """Test that commit requires a message."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        
        try:
            runner.invoke(app, ["init", "--quiet"])
            
            result = runner.invoke(app, ["commit"])
            
            assert result.exit_code == 1
            assert "message is required" in result.stdout.lower()
            
        finally:
            os.chdir(original_cwd)

    def test_commit_empty_staging_fails(self, tmp_path: Path) -> None:
        """Test that commit fails on empty staging area."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        
        try:
            runner.invoke(app, ["init", "--quiet"])
            
            result = runner.invoke(app, ["commit", "-m", "Empty"])
            
            assert result.exit_code == 1
            assert "nothing to commit" in result.stdout.lower()
            
        finally:
            os.chdir(original_cwd)

    def test_commit_allow_empty(self, tmp_path: Path) -> None:
        """Test --allow-empty flag."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        
        try:
            runner.invoke(app, ["init", "--quiet"])
            
            result = runner.invoke(app, ["commit", "-m", "Empty milestone", "--allow-empty"])
            
            assert result.exit_code == 0
            assert "Committed" in result.stdout
            
        finally:
            os.chdir(original_cwd)

    def test_commit_multiple_files(self, tmp_path: Path) -> None:
        """Test committing multiple files."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        
        try:
            runner.invoke(app, ["init", "--quiet"])
            
            # Create multiple files
            (tmp_path / "INCAR").write_text("ENCUT = 520\n")
            (tmp_path / "POSCAR").write_text("Fe\n1.0\n1 0 0\n0 1 0\n0 0 1\n")
            (tmp_path / "KPOINTS").write_text("Auto\n0\nG\n4 4 4\n")
            
            runner.invoke(app, ["add", "INCAR", "POSCAR", "KPOINTS"])
            
            result = runner.invoke(app, ["commit", "-m", "Add VASP inputs"])
            
            assert result.exit_code == 0
            assert "Committing 3 file(s)" in result.stdout
            
        finally:
            os.chdir(original_cwd)

    def test_commit_custom_author(self, tmp_path: Path) -> None:
        """Test --author flag."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        
        try:
            runner.invoke(app, ["init", "--quiet"])
            
            test_file = tmp_path / "test.txt"
            test_file.write_text("test")
            
            runner.invoke(app, ["add", "test.txt"])
            
            result = runner.invoke(
                app,
                ["commit", "-m", "Test", "--author", "testuser@testhost"]
            )
            
            assert result.exit_code == 0
            assert "testuser@testhost" in result.stdout
            
        finally:
            os.chdir(original_cwd)

    def test_commit_chain(self, tmp_path: Path) -> None:
        """Test creating a chain of commits."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        
        try:
            runner.invoke(app, ["init", "--quiet"])
            
            # First commit
            (tmp_path / "file1.txt").write_text("v1")
            runner.invoke(app, ["add", "file1.txt"])
            result1 = runner.invoke(app, ["commit", "-m", "First"])
            assert result1.exit_code == 0
            assert "(root commit)" in result1.stdout
            
            # Second commit
            (tmp_path / "file2.txt").write_text("v2")
            runner.invoke(app, ["add", "file2.txt"])
            result2 = runner.invoke(app, ["commit", "-m", "Second"])
            assert result2.exit_code == 0
            assert "(root commit)" not in result2.stdout
            
            # Verify HEAD points to second commit
            head_file = tmp_path / CHEMVCS_DIR / "HEAD"
            second_hash = head_file.read_text().strip()
            
            # Verify commit parent chain
            from chemvcs.storage import CommitBuilder, MetadataDB, ObjectStore
            
            db = MetadataDB(tmp_path / CHEMVCS_DIR)
            db.open()
            object_store = ObjectStore(tmp_path / CHEMVCS_DIR)
            builder = CommitBuilder(tmp_path / CHEMVCS_DIR, object_store, db)
            
            second_commit = builder.read_commit(second_hash)
            assert second_commit["parent"] is not None
            
            first_commit = builder.read_commit(second_commit["parent"])
            assert first_commit["parent"] is None
            
            db.close()
            
        finally:
            os.chdir(original_cwd)

    def test_commit_not_initialized(self, tmp_path: Path) -> None:
        """Test commit in non-initialized directory."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        
        try:
            result = runner.invoke(app, ["commit", "-m", "Test"])
            
            assert result.exit_code == 1
            assert "not a chemvcs repository" in result.stdout.lower()
            
        finally:
            os.chdir(original_cwd)
