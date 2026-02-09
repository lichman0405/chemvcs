"""Unit tests for MetadataDB."""

import sqlite3
from pathlib import Path

import pytest

from chemvcs.storage.metadata_db import DatabaseError, MetadataDB
from chemvcs.constants import DB_SCHEMA_VERSION


@pytest.fixture
def chemvcs_dir(tmp_path: Path) -> Path:
    """Create a temporary .chemvcs directory."""
    chemvcs = tmp_path / ".chemvcs"
    chemvcs.mkdir()
    return chemvcs


@pytest.fixture
def db(chemvcs_dir: Path) -> MetadataDB:
    """Create and initialize a MetadataDB instance."""
    metadata_db = MetadataDB(chemvcs_dir)
    metadata_db.open()
    metadata_db.init_schema()
    return metadata_db


class TestMetadataDBInit:
    """Test database initialization."""

    def test_init_creates_db_path(self, chemvcs_dir: Path) -> None:
        """Test that MetadataDB sets correct path."""
        db = MetadataDB(chemvcs_dir)
        assert db.db_path == chemvcs_dir / "metadata.db"
        assert db.conn is None

    def test_open_creates_connection(self, chemvcs_dir: Path) -> None:
        """Test opening database connection."""
        db = MetadataDB(chemvcs_dir)
        db.open()
        
        assert db.conn is not None
        assert isinstance(db.conn, sqlite3.Connection)
        
        db.close()

    def test_open_idempotent(self, chemvcs_dir: Path) -> None:
        """Test that calling open multiple times is safe."""
        db = MetadataDB(chemvcs_dir)
        db.open()
        conn1 = db.conn
        
        db.open()  # Second call
        conn2 = db.conn
        
        assert conn1 is conn2
        db.close()

    def test_context_manager(self, chemvcs_dir: Path) -> None:
        """Test using database as context manager."""
        db = MetadataDB(chemvcs_dir)
        
        with db:
            assert db.conn is not None
        
        assert db.conn is None  # Closed after context

    def test_init_schema_creates_tables(self, chemvcs_dir: Path) -> None:
        """Test schema initialization creates all tables."""
        db = MetadataDB(chemvcs_dir)
        db.open()
        db.init_schema()
        
        cursor = db.conn.cursor()  # type: ignore
        cursor.execute(
            "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name"
        )
        tables = [row[0] for row in cursor.fetchall()]
        
        expected_tables = [
            "commits",
            "environment",
            "files",
            "metadata",
            "output_summary",
            "semantic_summary",
        ]
        
        for table in expected_tables:
            assert table in tables
        
        db.close()

    def test_init_schema_sets_version(self, db: MetadataDB) -> None:
        """Test that schema version is set correctly."""
        version = db.get_schema_version()
        assert version == DB_SCHEMA_VERSION

    def test_init_schema_idempotent(self, db: MetadataDB) -> None:
        """Test that calling init_schema multiple times is safe."""
        # First init already done by fixture
        db.init_schema()  # Second call should not error
        
        version = db.get_schema_version()
        assert version == DB_SCHEMA_VERSION


class TestInsertCommit:
    """Test commit insertion."""

    def test_insert_commit_basic(self, db: MetadataDB) -> None:
        """Test basic commit insertion."""
        commit_id = db.insert_commit(
            commit_hash="a" * 64,
            parent_hash=None,
            timestamp="2026-02-09T10:00:00Z",
            author="test@host",
            message="Initial commit",
        )
        
        assert isinstance(commit_id, int)
        assert commit_id > 0

    def test_insert_commit_with_parent(self, db: MetadataDB) -> None:
        """Test inserting commit with parent."""
        parent_id = db.insert_commit(
            "a" * 64, None, "2026-02-09T10:00:00Z", "user@host", "First"
        )
        
        child_id = db.insert_commit(
            "b" * 64, "a" * 64, "2026-02-09T10:05:00Z", "user@host", "Second"
        )
        
        assert child_id > parent_id

    def test_insert_commit_duplicate_hash(self, db: MetadataDB) -> None:
        """Test that duplicate commit hash raises error."""
        db.insert_commit(
            "c" * 64, None, "2026-02-09T10:00:00Z", "user@host", "Test"
        )
        
        with pytest.raises(DatabaseError, match="already exists"):
            db.insert_commit(
                "c" * 64, None, "2026-02-09T10:01:00Z", "user@host", "Duplicate"
            )

    def test_insert_commit_db_not_open(self, chemvcs_dir: Path) -> None:
        """Test that operation fails if DB not open."""
        db = MetadataDB(chemvcs_dir)
        # Don't call db.open()
        
        with pytest.raises(DatabaseError, match="not open"):
            db.insert_commit("d" * 64, None, "2026-02-09T10:00:00Z", "user@host", "Test")


class TestInsertFile:
    """Test file insertion."""

    def test_insert_file_basic(self, db: MetadataDB) -> None:
        """Test basic file insertion."""
        commit_id = db.insert_commit(
            "e" * 64, None, "2026-02-09T10:00:00Z", "user@host", "Test"
        )
        
        file_id = db.insert_file(
            commit_id=commit_id,
            path="INCAR",
            blob_hash="f" * 64,
            is_reference=False,
            file_type="INCAR",
            size_bytes=1024,
        )
        
        assert isinstance(file_id, int)
        assert file_id > 0

    def test_insert_file_reference(self, db: MetadataDB) -> None:
        """Test inserting POTCAR as reference."""
        commit_id = db.insert_commit(
            "g" * 64, None, "2026-02-09T10:00:00Z", "user@host", "Test"
        )
        
        file_id = db.insert_file(
            commit_id=commit_id,
            path="POTCAR",
            blob_hash="h" * 64,
            is_reference=True,
            file_type="POTCAR",
        )
        
        assert file_id > 0

    def test_insert_file_foreign_key(self, db: MetadataDB) -> None:
        """Test that foreign key constraint is enforced."""
        with pytest.raises(DatabaseError, match="Failed to insert file"):
            db.insert_file(
                commit_id=99999,  # Non-existent commit
                path="INCAR",
                blob_hash="i" * 64,
            )


class TestGetCommitByHash:
    """Test commit retrieval by hash."""

    def test_get_commit_by_full_hash(self, db: MetadataDB) -> None:
        """Test retrieving commit by full hash."""
        commit_hash = "j" * 64
        db.insert_commit(
            commit_hash, None, "2026-02-09T10:00:00Z", "user@host", "Test message"
        )
        
        commit = db.get_commit_by_hash(commit_hash)
        
        assert commit is not None
        assert commit["commit_hash"] == commit_hash
        assert commit["message"] == "Test message"
        assert commit["author"] == "user@host"

    def test_get_commit_by_short_hash(self, db: MetadataDB) -> None:
        """Test retrieving commit by short hash (7 chars)."""
        commit_hash = "abc1234" + "k" * 57
        db.insert_commit(
            commit_hash, None, "2026-02-09T10:00:00Z", "user@host", "Short hash test"
        )
        
        commit = db.get_commit_by_hash("abc1234")
        
        assert commit is not None
        assert commit["commit_hash"] == commit_hash

    def test_get_commit_not_found(self, db: MetadataDB) -> None:
        """Test that non-existent commit returns None."""
        commit = db.get_commit_by_hash("nonexistent123")
        assert commit is None


class TestGetLatestCommit:
    """Test latest commit retrieval."""

    def test_get_latest_commit_empty_db(self, db: MetadataDB) -> None:
        """Test getting latest commit from empty database."""
        latest = db.get_latest_commit()
        assert latest is None

    def test_get_latest_commit_single(self, db: MetadataDB) -> None:
        """Test getting latest commit with one commit."""
        commit_hash = "m" * 64
        db.insert_commit(
            commit_hash, None, "2026-02-09T10:00:00Z", "user@host", "Only one"
        )
        
        latest = db.get_latest_commit()
        
        assert latest is not None
        assert latest["commit_hash"] == commit_hash

    def test_get_latest_commit_multiple(self, db: MetadataDB) -> None:
        """Test that latest commit is actually the most recent."""
        db.insert_commit(
            "n" * 64, None, "2026-02-09T10:00:00Z", "user@host", "First"
        )
        db.insert_commit(
            "o" * 64, "n" * 64, "2026-02-09T10:05:00Z", "user@host", "Second"
        )
        latest_hash = "p" * 64
        db.insert_commit(
            latest_hash, "o" * 64, "2026-02-09T10:10:00Z", "user@host", "Third"
        )
        
        latest = db.get_latest_commit()
        
        assert latest is not None
        assert latest["commit_hash"] == latest_hash


class TestGetFilesForCommit:
    """Test file listing for commits."""

    def test_get_files_empty(self, db: MetadataDB) -> None:
        """Test getting files for commit with no files."""
        commit_id = db.insert_commit(
            "q" * 64, None, "2026-02-09T10:00:00Z", "user@host", "Empty"
        )
        
        files = db.get_files_for_commit(commit_id)
        assert files == []

    def test_get_files_multiple(self, db: MetadataDB) -> None:
        """Test getting multiple files for commit."""
        commit_id = db.insert_commit(
            "r" * 64, None, "2026-02-09T10:00:00Z", "user@host", "Multi-file"
        )
        
        db.insert_file(commit_id, "INCAR", "s" * 64, file_type="INCAR")
        db.insert_file(commit_id, "POSCAR", "t" * 64, file_type="POSCAR")
        db.insert_file(commit_id, "KPOINTS", "u" * 64, file_type="KPOINTS")
        
        files = db.get_files_for_commit(commit_id)
        
        assert len(files) == 3
        paths = [f["path"] for f in files]
        assert "INCAR" in paths
        assert "POSCAR" in paths
        assert "KPOINTS" in paths


class TestGetCommitHistory:
    """Test commit history retrieval."""

    def test_get_history_empty(self, db: MetadataDB) -> None:
        """Test getting history from empty database."""
        history = db.get_commit_history()
        assert history == []

    def test_get_history_order(self, db: MetadataDB) -> None:
        """Test that history is in reverse chronological order."""
        hash1 = "v" * 64
        hash2 = "w" * 64
        hash3 = "x" * 64
        
        db.insert_commit(hash1, None, "2026-02-09T10:00:00Z", "user@host", "First")
        db.insert_commit(hash2, hash1, "2026-02-09T10:05:00Z", "user@host", "Second")
        db.insert_commit(hash3, hash2, "2026-02-09T10:10:00Z", "user@host", "Third")
        
        history = db.get_commit_history()
        
        assert len(history) == 3
        assert history[0]["commit_hash"] == hash3  # Most recent first
        assert history[1]["commit_hash"] == hash2
        assert history[2]["commit_hash"] == hash1

    def test_get_history_with_limit(self, db: MetadataDB) -> None:
        """Test limiting history results."""
        for i in range(10):
            db.insert_commit(
                chr(ord("a") + i) * 64,
                None,
                f"2026-02-09T10:{i:02d}:00Z",
                "user@host",
                f"Commit {i}",
            )
        
        history = db.get_commit_history(limit=5)
        
        assert len(history) == 5


class TestSchemaVersion:
    """Test schema version management."""

    def test_get_schema_version(self, db: MetadataDB) -> None:
        """Test retrieving schema version."""
        version = db.get_schema_version()
        assert version == DB_SCHEMA_VERSION
        assert isinstance(version, int)


class TestForeignKeys:
    """Test foreign key constraints."""

    def test_cascade_delete_files(self, db: MetadataDB) -> None:
        """Test that deleting commit cascades to files."""
        commit_id = db.insert_commit(
            "y" * 64, None, "2026-02-09T10:00:00Z", "user@host", "Delete test"
        )
        
        db.insert_file(commit_id, "INCAR", "z" * 64)
        db.insert_file(commit_id, "POSCAR", "A" * 64)
        
        # Verify files exist
        files_before = db.get_files_for_commit(commit_id)
        assert len(files_before) == 2
        
        # Delete commit
        cursor = db.conn.cursor()  # type: ignore
        cursor.execute("DELETE FROM commits WHERE id = ?", (commit_id,))
        db.conn.commit()  # type: ignore
        
        # Verify files are gone
        files_after = db.get_files_for_commit(commit_id)
        assert len(files_after) == 0


class TestWALMode:
    """Test WAL mode detection and fallback."""

    def test_wal_detection(self, chemvcs_dir: Path) -> None:
        """Test WAL mode detection."""
        db = MetadataDB(chemvcs_dir)
        db.open()
        
        # Just verify it doesn't crash
        # Actual WAL support depends on filesystem
        assert db._wal_mode_supported is not None
        assert isinstance(db._wal_mode_supported, bool)
        
        db.close()
