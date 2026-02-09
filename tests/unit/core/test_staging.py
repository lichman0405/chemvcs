"""Unit tests for StagingManager."""

import json
from pathlib import Path

import pytest

from chemvcs.core.staging import StagingError, StagingManager
from chemvcs.storage import ObjectStore


@pytest.fixture
def workspace(tmp_path: Path) -> Path:
    """Create a workspace with initialized .chemvcs."""
    workspace_root = tmp_path / "workspace"
    workspace_root.mkdir()
    
    chemvcs_dir = workspace_root / ".chemvcs"
    chemvcs_dir.mkdir()
    (chemvcs_dir / "objects").mkdir()
    
    return workspace_root


@pytest.fixture
def object_store(workspace: Path) -> ObjectStore:
    """Create ObjectStore instance."""
    return ObjectStore(workspace / ".chemvcs")


@pytest.fixture
def staging(workspace: Path, object_store: ObjectStore) -> StagingManager:
    """Create StagingManager instance."""
    return StagingManager(workspace, object_store)


class TestStagingManagerInit:
    """Test StagingManager initialization."""

    def test_init_valid_workspace(self, workspace: Path, object_store: ObjectStore) -> None:
        """Test initialization with valid workspace."""
        staging = StagingManager(workspace, object_store)
        
        assert staging.workspace_root == workspace
        assert staging.index_path == workspace / ".chemvcs" / "index"

    def test_init_no_chemvcs(self, tmp_path: Path, object_store: ObjectStore) -> None:
        """Test initialization fails without .chemvcs directory."""
        with pytest.raises(StagingError, match="Not a ChemVCS repository"):
            StagingManager(tmp_path, object_store)


class TestAdd:
    """Test adding files to staging area."""

    def test_add_single_file(self, staging: StagingManager, workspace: Path) -> None:
        """Test adding a single file."""
        test_file = workspace / "test.txt"
        test_file.write_text("Hello, world!")
        
        stats = staging.add([test_file])
        
        assert "test.txt" in stats["added"]
        assert len(stats["errors"]) == 0
        
        # Verify file is in index
        staged = staging.get_staged_files()
        assert "test.txt" in staged
        assert staged["test.txt"]["size_bytes"] == 13

    def test_add_multiple_files(self, staging: StagingManager, workspace: Path) -> None:
        """Test adding multiple files at once."""
        file1 = workspace / "file1.txt"
        file2 = workspace / "file2.txt"
        file1.write_text("content1")
        file2.write_text("content2")
        
        stats = staging.add([file1, file2])
        
        assert "file1.txt" in stats["added"]
        assert "file2.txt" in stats["added"]
        
        staged = staging.get_staged_files()
        assert len(staged) == 2

    def test_add_update_existing(self, staging: StagingManager, workspace: Path) -> None:
        """Test updating an already-staged file."""
        test_file = workspace / "test.txt"
        test_file.write_text("original")
        
        # First add
        stats1 = staging.add([test_file])
        assert "test.txt" in stats1["added"]
        
        # Modify and re-add
        test_file.write_text("modified")
        stats2 = staging.add([test_file])
        
        assert "test.txt" in stats2["updated"]
        assert "test.txt" not in stats2["added"]

    def test_add_directory(self, staging: StagingManager, workspace: Path) -> None:
        """Test adding a directory recursively."""
        subdir = workspace / "subdir"
        subdir.mkdir()
        (subdir / "file1.txt").write_text("content1")
        (subdir / "file2.txt").write_text("content2")
        
        stats = staging.add([subdir])
        
        assert "subdir/file1.txt" in stats["added"]
        assert "subdir/file2.txt" in stats["added"]

    def test_add_nonexistent_file(self, staging: StagingManager, workspace: Path) -> None:
        """Test adding a file that doesn't exist."""
        nonexistent = workspace / "does_not_exist.txt"
        
        stats = staging.add([nonexistent])
        
        assert len(stats["added"]) == 0
        assert any("not found" in err for err in stats["errors"])

    def test_add_ignores_chemvcs(self, staging: StagingManager, workspace: Path) -> None:
        """Test that .chemvcs directory is automatically ignored."""
        chemvcs_file = workspace / ".chemvcs" / "test.txt"
        chemvcs_file.write_text("should be ignored")
        
        stats = staging.add([chemvcs_file])
        
        # Should be silently skipped
        assert len(stats["added"]) == 0
        assert len(stats["errors"]) == 0

    def test_add_with_ignore_patterns(
        self, staging: StagingManager, workspace: Path
    ) -> None:
        """Test that .chemvcsignore patterns are respected."""
        # Create .chemvcsignore
        ignore_file = workspace / ".chemvcsignore"
        ignore_file.write_text("*.tmp\n__pycache__/\n")
        
        # Create files
        normal = workspace / "normal.txt"
        ignored = workspace / "test.tmp"
        normal.write_text("normal")
        ignored.write_text("ignored")
        
        stats = staging.add([normal, ignored])
        
        assert "normal.txt" in stats["added"]
        assert "test.tmp" in stats["ignored"]

    def test_add_force_overrides_ignore(
        self, staging: StagingManager, workspace: Path
    ) -> None:
        """Test that --force overrides ignore patterns."""
        ignore_file = workspace / ".chemvcsignore"
        ignore_file.write_text("*.tmp\n")
        
        ignored = workspace / "test.tmp"
        ignored.write_text("content")
        
        stats = staging.add([ignored], force=True)
        
        assert "test.tmp" in stats["added"]
        assert len(stats["ignored"]) == 0

    def test_add_detects_vasp_files(
        self, staging: StagingManager, workspace: Path
    ) -> None:
        """Test that VASP file types are detected correctly."""
        incar = workspace / "INCAR"
        poscar = workspace / "POSCAR"
        incar.write_text("ENCUT = 520")
        poscar.write_text("Fe\n1.0\n...")
        
        staging.add([incar, poscar])
        
        staged = staging.get_staged_files()
        assert staged["INCAR"]["file_type"] == "INCAR"
        assert staged["POSCAR"]["file_type"] == "POSCAR"


class TestRemove:
    """Test removing files from staging area."""

    def test_remove_staged_file(self, staging: StagingManager, workspace: Path) -> None:
        """Test removing a file from staging."""
        test_file = workspace / "test.txt"
        test_file.write_text("content")
        
        # Add then remove
        staging.add([test_file])
        stats = staging.remove([test_file])
        
        assert "test.txt" in stats["removed"]
        assert staging.is_empty()

    def test_remove_not_staged(self, staging: StagingManager, workspace: Path) -> None:
        """Test removing a file that wasn't staged."""
        test_file = workspace / "test.txt"
        test_file.write_text("content")
        
        stats = staging.remove([test_file])
        
        assert "test.txt" in stats["not_staged"]


class TestGetStagedFiles:
    """Test retrieving staged files."""

    def test_get_staged_empty(self, staging: StagingManager) -> None:
        """Test getting staged files when empty."""
        staged = staging.get_staged_files()
        assert staged == {}

    def test_get_staged_with_files(
        self, staging: StagingManager, workspace: Path
    ) -> None:
        """Test getting staged files."""
        file1 = workspace / "file1.txt"
        file2 = workspace / "file2.txt"
        file1.write_text("content1")
        file2.write_text("content2")
        
        staging.add([file1, file2])
        
        staged = staging.get_staged_files()
        assert len(staged) == 2
        assert "file1.txt" in staged
        assert "file2.txt" in staged


class TestClear:
    """Test clearing staging area."""

    def test_clear(self, staging: StagingManager, workspace: Path) -> None:
        """Test clearing all staged files."""
        file1 = workspace / "file1.txt"
        file1.write_text("content")
        
        staging.add([file1])
        assert not staging.is_empty()
        
        staging.clear()
        assert staging.is_empty()


class TestIsEmpty:
    """Test checking if staging area is empty."""

    def test_is_empty_initially(self, staging: StagingManager) -> None:
        """Test that staging area is empty initially."""
        assert staging.is_empty()

    def test_is_empty_after_add(self, staging: StagingManager, workspace: Path) -> None:
        """Test that staging area is not empty after adding."""
        file1 = workspace / "file1.txt"
        file1.write_text("content")
        
        staging.add([file1])
        assert not staging.is_empty()


class TestIndexPersistence:
    """Test index file persistence."""

    def test_index_persists_across_instances(
        self, workspace: Path, object_store: ObjectStore
    ) -> None:
        """Test that index is saved and loaded correctly."""
        file1 = workspace / "file1.txt"
        file1.write_text("content")
        
        # Create first instance and add file
        staging1 = StagingManager(workspace, object_store)
        staging1.add([file1])
        
        # Create second instance
        staging2 = StagingManager(workspace, object_store)
        
        # Should still have the file
        staged = staging2.get_staged_files()
        assert "file1.txt" in staged

    def test_index_atomic_write(self, staging: StagingManager, workspace: Path) -> None:
        """Test that index updates are atomic."""
        file1 = workspace / "file1.txt"
        file1.write_text("content")
        
        staging.add([file1])
        
        # Index file should exist
        assert staging.index_path.exists()
        
        # Should be valid JSON
        with open(staging.index_path, "r") as f:
            index = json.load(f)
        
        assert index["version"] == 1
        assert "file1.txt" in index["entries"]


class TestPathResolution:
    """Test path resolution and validation."""

    def test_relative_path_resolution(
        self, staging: StagingManager, workspace: Path
    ) -> None:
        """Test that relative paths are resolved correctly."""
        import os
        
        subdir = workspace / "subdir"
        subdir.mkdir()
        file1 = subdir / "file1.txt"
        file1.write_text("content")
        
        # Change to workspace directory
        original_cwd = Path.cwd()
        os.chdir(workspace)
        
        try:
            # Add using relative path
            stats = staging.add([Path("subdir/file1.txt")])
            
            assert "subdir/file1.txt" in stats["added"]
        finally:
            os.chdir(original_cwd)

    def test_outside_workspace_rejected(
        self, staging: StagingManager, tmp_path: Path
    ) -> None:
        """Test that paths outside workspace are rejected."""
        outside_file = tmp_path / "outside.txt"
        outside_file.write_text("content")
        
        stats = staging.add([outside_file])
        
        assert any("outside workspace" in err for err in stats["errors"])
