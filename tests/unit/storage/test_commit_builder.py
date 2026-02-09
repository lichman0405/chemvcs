"""Unit tests for CommitBuilder."""

import json
from pathlib import Path

import pytest

from chemvcs.storage.commit_builder import CommitBuilder, CommitBuilderError
from chemvcs.storage.metadata_db import MetadataDB
from chemvcs.storage.object_store import ObjectStore


@pytest.fixture
def chemvcs_dir(tmp_path: Path) -> Path:
    """Create a temporary .chemvcs directory."""
    chemvcs = tmp_path / ".chemvcs"
    chemvcs.mkdir()
    return chemvcs


@pytest.fixture
def object_store(chemvcs_dir: Path) -> ObjectStore:
    """Create ObjectStore instance."""
    return ObjectStore(chemvcs_dir)


@pytest.fixture
def metadata_db(chemvcs_dir: Path) -> MetadataDB:
    """Create and initialize MetadataDB instance."""
    db = MetadataDB(chemvcs_dir)
    db.open()
    db.init_schema()
    yield db
    db.close()


@pytest.fixture
def commit_builder(
    chemvcs_dir: Path, object_store: ObjectStore, metadata_db: MetadataDB
) -> CommitBuilder:
    """Create CommitBuilder instance."""
    return CommitBuilder(chemvcs_dir, object_store, metadata_db)


class TestCommitBuilderInit:
    """Test CommitBuilder initialization."""

    def test_init_creates_commits_dir(
        self, chemvcs_dir: Path, object_store: ObjectStore, metadata_db: MetadataDB
    ) -> None:
        """Test that CommitBuilder creates commits directory."""
        builder = CommitBuilder(chemvcs_dir, object_store, metadata_db)
        
        assert builder.commits_dir.exists()
        assert builder.commits_dir.is_dir()
        assert builder.commits_dir == chemvcs_dir / "commits"

    def test_init_with_existing_commits_dir(
        self, chemvcs_dir: Path, object_store: ObjectStore, metadata_db: MetadataDB
    ) -> None:
        """Test init with existing commits directory."""
        commits_dir = chemvcs_dir / "commits"
        commits_dir.mkdir()
        
        builder = CommitBuilder(chemvcs_dir, object_store, metadata_db)
        
        assert builder.commits_dir.exists()


class TestCreateCommit:
    """Test commit creation."""

    def test_create_commit_basic(self, commit_builder: CommitBuilder) -> None:
        """Test basic commit creation with minimal data."""
        files = [
            {
                "path": "INCAR",
                "content": b"ENCUT = 520\nISMEAR = 0\n",
                "file_type": "INCAR",
            }
        ]
        
        commit_hash = commit_builder.create_commit(
            files=files,
            message="Initial commit",
            author="test@host",
        )
        
        assert isinstance(commit_hash, str)
        assert len(commit_hash) == 64
        assert all(c in "0123456789abcdef" for c in commit_hash)

    def test_create_commit_with_parent(self, commit_builder: CommitBuilder) -> None:
        """Test creating commit with parent."""
        # Create parent commit
        parent_hash = commit_builder.create_commit(
            files=[{"path": "file1.txt", "content": b"parent", "file_type": "text"}],
            message="Parent commit",
            author="test@host",
        )
        
        # Create child commit
        child_hash = commit_builder.create_commit(
            files=[{"path": "file2.txt", "content": b"child", "file_type": "text"}],
            message="Child commit",
            author="test@host",
            parent_hash=parent_hash,
        )
        
        assert child_hash != parent_hash
        
        # Verify parent relationship in commit object
        child_obj = commit_builder.read_commit(child_hash)
        assert child_obj["parent"] == parent_hash

    def test_create_commit_with_metadata(self, commit_builder: CommitBuilder) -> None:
        """Test creating commit with semantic and output metadata."""
        files = [
            {
                "path": "POSCAR",
                "content": b"Fe\n1.0\n1 0 0\n0 1 0\n0 0 1\nFe\n1\nDirect\n0 0 0\n",
                "file_type": "POSCAR",
            }
        ]
        
        semantic_data = {
            "formula": "Fe1",
            "space_group": "Pm-3m",
            "lattice_parameters": {"a": 1.0, "b": 1.0, "c": 1.0},
        }
        
        output_data = {
            "final_energy": -8.532,
            "converged": True,
        }
        
        environment = {
            "VASP_VERSION": "6.3.0",
            "HOSTNAME": "node01",
        }
        
        commit_hash = commit_builder.create_commit(
            files=files,
            message="Commit with metadata",
            author="test@host",
            semantic_data=semantic_data,
            output_data=output_data,
            environment=environment,
        )
        
        # Verify metadata is stored
        commit_obj = commit_builder.read_commit(commit_hash)
        assert commit_obj["semantic_summary"] == semantic_data
        assert commit_obj["output_summary"] == output_data
        assert commit_obj["environment"] == environment

    def test_create_commit_multiple_files(self, commit_builder: CommitBuilder) -> None:
        """Test creating commit with multiple files."""
        files = [
            {"path": "INCAR", "content": b"ENCUT = 520\n", "file_type": "INCAR"},
            {"path": "POSCAR", "content": b"Fe\n...", "file_type": "POSCAR"},
            {"path": "KPOINTS", "content": b"Auto\n0\nG\n4 4 4\n", "file_type": "KPOINTS"},
        ]
        
        commit_hash = commit_builder.create_commit(
            files=files,
            message="Multi-file commit",
            author="test@host",
        )
        
        commit_obj = commit_builder.read_commit(commit_hash)
        assert len(commit_obj["files"]) == 3
        
        paths = [f["path"] for f in commit_obj["files"]]
        assert "INCAR" in paths
        assert "POSCAR" in paths
        assert "KPOINTS" in paths

    def test_create_commit_reference_file(self, commit_builder: CommitBuilder) -> None:
        """Test creating commit with reference file (POTCAR)."""
        files = [
            {
                "path": "POTCAR",
                "content": b"PAW_PBE Fe 06Sep2000\n...",
                "file_type": "POTCAR",
                "is_reference": True,
            }
        ]
        
        commit_hash = commit_builder.create_commit(
            files=files,
            message="Commit with POTCAR reference",
            author="test@host",
        )
        
        commit_obj = commit_builder.read_commit(commit_hash)
        assert commit_obj["files"][0]["is_reference"] is True

    def test_create_commit_blob_deduplication(
        self, commit_builder: CommitBuilder
    ) -> None:
        """Test that identical file content reuses same blob."""
        content = b"Same content"
        
        # Create two commits with same file content
        hash1 = commit_builder.create_commit(
            files=[{"path": "file1.txt", "content": content, "file_type": "text"}],
            message="First",
            author="test@host",
        )
        
        hash2 = commit_builder.create_commit(
            files=[{"path": "file2.txt", "content": content, "file_type": "text"}],
            message="Second",
            author="test@host",
        )
        
        # Verify they reference the same blob
        obj1 = commit_builder.read_commit(hash1)
        obj2 = commit_builder.read_commit(hash2)
        
        blob_hash1 = obj1["files"][0]["blob_hash"]
        blob_hash2 = obj2["files"][0]["blob_hash"]
        
        assert blob_hash1 == blob_hash2


class TestReadCommit:
    """Test commit reading."""

    def test_read_commit_exists(self, commit_builder: CommitBuilder) -> None:
        """Test reading an existing commit."""
        files = [{"path": "test.txt", "content": b"test", "file_type": "text"}]
        
        commit_hash = commit_builder.create_commit(
            files=files,
            message="Test commit",
            author="test@host",
        )
        
        commit_obj = commit_builder.read_commit(commit_hash)
        
        assert commit_obj["hash"] == commit_hash
        assert commit_obj["message"] == "Test commit"
        assert commit_obj["author"] == "test@host"
        assert len(commit_obj["files"]) == 1

    def test_read_commit_short_hash(self, commit_builder: CommitBuilder) -> None:
        """Test reading commit with short hash."""
        commit_hash = commit_builder.create_commit(
            files=[{"path": "file.txt", "content": b"data", "file_type": "text"}],
            message="Test",
            author="test@host",
        )
        
        short_hash = commit_hash[:7]
        commit_obj = commit_builder.read_commit(short_hash)
        
        assert commit_obj["hash"] == commit_hash

    def test_read_commit_not_found(self, commit_builder: CommitBuilder) -> None:
        """Test reading non-existent commit raises error."""
        with pytest.raises(CommitBuilderError, match="not found"):
            commit_builder.read_commit("nonexistent123")

    def test_read_commit_hash_mismatch(
        self, commit_builder: CommitBuilder, chemvcs_dir: Path
    ) -> None:
        """Test reading commit with corrupted hash raises error."""
        # Create valid commit
        commit_hash = commit_builder.create_commit(
            files=[{"path": "file.txt", "content": b"test", "file_type": "text"}],
            message="Test",
            author="test@host",
        )
        
        # Corrupt the commit file by changing hash field
        commit_path = chemvcs_dir / "commits" / commit_hash
        with open(commit_path, "r", encoding="utf-8") as f:
            commit_obj = json.load(f)
        
        commit_obj["hash"] = "0" * 64  # Wrong hash
        
        with open(commit_path, "w", encoding="utf-8") as f:
            json.dump(commit_obj, f)
        
        with pytest.raises(CommitBuilderError, match="hash mismatch"):
            commit_builder.read_commit(commit_hash)


class TestCommitExists:
    """Test commit existence checking."""

    def test_commit_exists_true(self, commit_builder: CommitBuilder) -> None:
        """Test checking existing commit."""
        commit_hash = commit_builder.create_commit(
            files=[{"path": "file.txt", "content": b"test", "file_type": "text"}],
            message="Test",
            author="test@host",
        )
        
        assert commit_builder.commit_exists(commit_hash) is True

    def test_commit_exists_false(self, commit_builder: CommitBuilder) -> None:
        """Test checking non-existent commit."""
        assert commit_builder.commit_exists("a" * 64) is False


class TestCommitHash:
    """Test commit hash computation."""

    def test_commit_hash_deterministic(self, commit_builder: CommitBuilder) -> None:
        """Test that same commit data produces same hash."""
        files = [{"path": "test.txt", "content": b"test content", "file_type": "text"}]
        
        hash1 = commit_builder.create_commit(
            files=files,
            message="Test message",
            author="test@host",
        )
        
        # Create another commit with same files/message
        # Note: timestamp will differ, so hashes will differ
        # Instead, test _compute_commit_hash directly
        
        commit_obj = {
            "parent": None,
            "timestamp": "2026-02-09T10:00:00Z",
            "author": "test@host",
            "message": "Test",
            "files": [{"path": "f.txt", "blob_hash": "abc123"}],
        }
        
        hash1 = commit_builder._compute_commit_hash(commit_obj)
        hash2 = commit_builder._compute_commit_hash(commit_obj)
        
        assert hash1 == hash2

    def test_commit_hash_different_content(self, commit_builder: CommitBuilder) -> None:
        """Test that different commit data produces different hash."""
        obj1 = {
            "parent": None,
            "timestamp": "2026-02-09T10:00:00Z",
            "author": "test@host",
            "message": "Message 1",
            "files": [],
        }
        
        obj2 = {
            "parent": None,
            "timestamp": "2026-02-09T10:00:00Z",
            "author": "test@host",
            "message": "Message 2",
            "files": [],
        }
        
        hash1 = commit_builder._compute_commit_hash(obj1)
        hash2 = commit_builder._compute_commit_hash(obj2)
        
        assert hash1 != hash2

    def test_commit_hash_ignores_hash_field(
        self, commit_builder: CommitBuilder
    ) -> None:
        """Test that hash field is excluded from hash computation."""
        obj_without_hash = {
            "parent": None,
            "timestamp": "2026-02-09T10:00:00Z",
            "author": "test@host",
            "message": "Test",
            "files": [],
        }
        
        obj_with_hash = {
            **obj_without_hash,
            "hash": "should_be_ignored",
        }
        
        hash1 = commit_builder._compute_commit_hash(obj_without_hash)
        hash2 = commit_builder._compute_commit_hash(obj_with_hash)
        
        assert hash1 == hash2


class TestCommitPersistence:
    """Test commit persistence and file format."""

    def test_commit_file_created(
        self, commit_builder: CommitBuilder, chemvcs_dir: Path
    ) -> None:
        """Test that commit file is created on disk."""
        commit_hash = commit_builder.create_commit(
            files=[{"path": "file.txt", "content": b"data", "file_type": "text"}],
            message="Test",
            author="test@host",
        )
        
        commit_path = chemvcs_dir / "commits" / commit_hash
        assert commit_path.exists()
        assert commit_path.is_file()

    def test_commit_file_valid_json(
        self, commit_builder: CommitBuilder, chemvcs_dir: Path
    ) -> None:
        """Test that commit file is valid JSON."""
        commit_hash = commit_builder.create_commit(
            files=[{"path": "file.txt", "content": b"data", "file_type": "text"}],
            message="Test",
            author="test@host",
        )
        
        commit_path = chemvcs_dir / "commits" / commit_hash
        
        with open(commit_path, "r", encoding="utf-8") as f:
            commit_obj = json.load(f)
        
        assert isinstance(commit_obj, dict)
        assert commit_obj["hash"] == commit_hash


class TestDatabaseIntegration:
    """Test integration with MetadataDB."""

    def test_commit_indexed_in_db(
        self, commit_builder: CommitBuilder, metadata_db: MetadataDB
    ) -> None:
        """Test that commit is indexed in database."""
        commit_hash = commit_builder.create_commit(
            files=[{"path": "file.txt", "content": b"data", "file_type": "text"}],
            message="Test commit",
            author="test@host",
        )
        
        # Query database
        db_commit = metadata_db.get_commit_by_hash(commit_hash)
        
        assert db_commit is not None
        assert db_commit["commit_hash"] == commit_hash
        assert db_commit["message"] == "Test commit"
        assert db_commit["author"] == "test@host"

    def test_files_indexed_in_db(
        self, commit_builder: CommitBuilder, metadata_db: MetadataDB
    ) -> None:
        """Test that files are indexed in database."""
        files = [
            {"path": "INCAR", "content": b"ENCUT = 520\n", "file_type": "INCAR"},
            {"path": "POSCAR", "content": b"Fe\n...", "file_type": "POSCAR"},
        ]
        
        commit_hash = commit_builder.create_commit(
            files=files,
            message="Multi-file",
            author="test@host",
        )
        
        # Get commit ID from database
        db_commit = metadata_db.get_commit_by_hash(commit_hash)
        commit_id = db_commit["id"]
        
        # Query files
        db_files = metadata_db.get_files_for_commit(commit_id)
        
        assert len(db_files) == 2
        paths = [f["path"] for f in db_files]
        assert "INCAR" in paths
        assert "POSCAR" in paths


class TestObjectStoreIntegration:
    """Test integration with ObjectStore."""

    def test_file_content_stored_in_blobs(
        self, commit_builder: CommitBuilder, object_store: ObjectStore
    ) -> None:
        """Test that file content is stored as blobs."""
        content = b"Test file content"
        
        commit_hash = commit_builder.create_commit(
            files=[{"path": "test.txt", "content": content, "file_type": "text"}],
            message="Test",
            author="test@host",
        )
        
        # Get blob hash from commit
        commit_obj = commit_builder.read_commit(commit_hash)
        blob_hash = commit_obj["files"][0]["blob_hash"]
        
        # Verify blob exists and contains correct content
        assert object_store.blob_exists(blob_hash)
        retrieved_content = object_store.read_blob(blob_hash)
        assert retrieved_content == content
