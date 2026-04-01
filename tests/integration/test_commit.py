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


class TestCommitAllFlag:
    """Tests for 'chemvcs commit --all' auto-staging behaviour."""

    def test_commit_all_auto_stages_tracked_files(self, tmp_path: Path) -> None:
        """--all re-stages files from the previous commit that exist on disk."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        try:
            runner.invoke(app, ["init", "--quiet"])

            # First commit — stage and commit INCAR manually
            (tmp_path / "INCAR").write_text("ENCUT = 520\n", encoding="utf-8")
            runner.invoke(app, ["add", "INCAR"])
            runner.invoke(app, ["commit", "-m", "Initial"])

            # Modify file on disk WITHOUT calling 'add'
            (tmp_path / "INCAR").write_text("ENCUT = 600\n", encoding="utf-8")

            # commit --all should auto-stage the modified tracked file
            result = runner.invoke(app, ["commit", "--all", "-m", "Updated via --all"])
            assert result.exit_code == 0
            assert "Committed" in result.stdout

            # Verify the new commit contains the updated content
            from chemvcs.constants import CHEMVCS_DIR
            from chemvcs.storage import CommitBuilder, MetadataDB, ObjectStore

            chemvcs_dir = tmp_path / CHEMVCS_DIR
            db = MetadataDB(chemvcs_dir)
            db.open()
            obj_store = ObjectStore(chemvcs_dir)
            builder = CommitBuilder(chemvcs_dir, obj_store, db)

            head_hash = (chemvcs_dir / "HEAD").read_text().strip()
            commit = builder.read_commit(head_hash)
            incar_entry = next(f for f in commit["files"] if f["path"] == "INCAR")
            content = obj_store.read_blob(incar_entry["blob_hash"]).decode()
            assert "600" in content
            db.close()

        finally:
            os.chdir(original_cwd)

    def test_commit_all_on_first_commit_gives_hint(self, tmp_path: Path) -> None:
        """--all on first commit (no HEAD) prints a hint and does not crash."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        try:
            runner.invoke(app, ["init", "--quiet"])
            # Stage a file manually so the commit is not empty
            (tmp_path / "INCAR").write_text("ENCUT = 520\n", encoding="utf-8")
            runner.invoke(app, ["add", "INCAR"])
            result = runner.invoke(app, ["commit", "--all", "-m", "First"])
            assert result.exit_code == 0
            # Should mention that there are no tracked files to auto-stage
            assert "no tracked files" in result.stdout.lower() or "no previous commit" in result.stdout.lower()
        finally:
            os.chdir(original_cwd)

    def test_commit_all_skips_deleted_tracked_file(self, tmp_path: Path) -> None:
        """--all skips files that were tracked but no longer exist on disk."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        try:
            runner.invoke(app, ["init", "--quiet"])

            (tmp_path / "INCAR").write_text("ENCUT = 520\n", encoding="utf-8")
            (tmp_path / "KPOINTS").write_text("Automatic\n0\nGamma\n4 4 4\n", encoding="utf-8")
            runner.invoke(app, ["add", "INCAR", "KPOINTS"])
            runner.invoke(app, ["commit", "-m", "Initial"])

            # Delete INCAR from disk
            (tmp_path / "INCAR").unlink()
            # Modify KPOINTS
            (tmp_path / "KPOINTS").write_text("Automatic\n0\nGamma\n6 6 6\n", encoding="utf-8")

            result = runner.invoke(app, ["commit", "--all", "-m", "After deletion"])
            assert result.exit_code == 0
            # Should mention skipping the missing file
            assert "skipped" in result.stdout.lower()
        finally:
            os.chdir(original_cwd)
