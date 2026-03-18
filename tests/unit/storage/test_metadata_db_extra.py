"""Extra coverage tests for metadata_db.py.

Targets uncovered lines:
   84      – WAL mode disabled → uses DELETE journal mode
   89-90   – open() sqlite3.Error path
  116-117  – _detect_wal_support sqlite3.Error path
   129     – init_schema when conn is None
  244-246  – insert_commit duplicate hash (IntegrityError)
  289-291  – insert_file when conn is None (and sqlite3 error)
   316     – insert_file when conn is None
   344     – get_latest_commit when conn is None
  367-368  – get_files_for_commit when conn is None
   377     – get_latest_commit raises when conn is None
  387-388  – get_latest_commit sqlite3.Error
   400     – get_commit_history conn is None
  411-412  – get_commit_history raises
   429     – get_commit_history error handling
  436-446  – get_commit_history with start_hash parameter
  458-459  – get_commit_history sqlite3.Error
   468     – get_schema_version raises if conn is None
   475     – get_schema_version returns 0 when no version row
  478-479  – get_schema_version sqlite3.Error
"""

import sqlite3
from pathlib import Path
from unittest.mock import patch, MagicMock

import pytest

from chemvcs.storage.metadata_db import DatabaseError, MetadataDB


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def chemvcs_dir(tmp_path: Path) -> Path:
    d = tmp_path / ".chemvcs"
    d.mkdir()
    return d


@pytest.fixture
def db(chemvcs_dir: Path) -> MetadataDB:
    m = MetadataDB(chemvcs_dir)
    m.open()
    m.init_schema()
    return m


# ---------------------------------------------------------------------------
# WAL / open edge cases
# ---------------------------------------------------------------------------

class TestWalAndOpen:
    def test_open_delete_journal_mode(self, chemvcs_dir: Path) -> None:
        """Force non-WAL mode to cover the else-branch (journal_mode=DELETE)."""
        db = MetadataDB(chemvcs_dir)
        # Set _wal_mode_supported=False before open so the else-branch is taken
        db._wal_mode_supported = False
        db.open()
        assert db.conn is not None
        db.close()

    def test_open_sqlite_error_raises_database_error(self, chemvcs_dir: Path) -> None:
        """If sqlite3.connect() fails inside open(), DatabaseError is raised."""
        db = MetadataDB(chemvcs_dir)
        with patch(
            "chemvcs.storage.metadata_db.sqlite3.connect",
            side_effect=sqlite3.OperationalError("no such file"),
        ):
            with pytest.raises(DatabaseError, match="Failed to open database"):
                db.open()

    def test_detect_wal_sqlite_error(self, chemvcs_dir: Path) -> None:
        """If WAL detection raises sqlite3.Error, wal is set to False."""
        db = MetadataDB(chemvcs_dir)
        # Replace conn with a MagicMock so execute() raises sqlite3.Error
        mock_conn = MagicMock()
        mock_conn.execute.side_effect = sqlite3.OperationalError("no WAL")
        db.conn = mock_conn

        db._detect_wal_support()

        assert db._wal_mode_supported is False


# ---------------------------------------------------------------------------
# "conn is None" error paths for all public methods
# ---------------------------------------------------------------------------

class TestConnNoneErrors:
    def test_init_schema_no_conn(self, chemvcs_dir: Path) -> None:
        db = MetadataDB(chemvcs_dir)
        with pytest.raises(DatabaseError, match="not open"):
            db.init_schema()

    def test_insert_commit_no_conn(self, chemvcs_dir: Path) -> None:
        db = MetadataDB(chemvcs_dir)
        with pytest.raises(DatabaseError, match="not open"):
            db.insert_commit("a" * 64, None, "2026-01-01T00:00:00Z", "u", "msg")

    def test_insert_file_no_conn(self, chemvcs_dir: Path) -> None:
        db = MetadataDB(chemvcs_dir)
        with pytest.raises(DatabaseError, match="not open"):
            db.insert_file(1, "INCAR", "b" * 64)

    def test_get_commit_by_hash_no_conn(self, chemvcs_dir: Path) -> None:
        db = MetadataDB(chemvcs_dir)
        with pytest.raises(DatabaseError, match="not open"):
            db.get_commit_by_hash("abc1234")

    def test_get_latest_commit_no_conn(self, chemvcs_dir: Path) -> None:
        db = MetadataDB(chemvcs_dir)
        with pytest.raises(DatabaseError, match="not open"):
            db.get_latest_commit()

    def test_get_files_for_commit_no_conn(self, chemvcs_dir: Path) -> None:
        db = MetadataDB(chemvcs_dir)
        with pytest.raises(DatabaseError, match="not open"):
            db.get_files_for_commit(1)

    def test_get_commit_history_no_conn(self, chemvcs_dir: Path) -> None:
        db = MetadataDB(chemvcs_dir)
        with pytest.raises(DatabaseError, match="not open"):
            db.get_commit_history()

    def test_get_schema_version_no_conn(self, chemvcs_dir: Path) -> None:
        db = MetadataDB(chemvcs_dir)
        with pytest.raises(DatabaseError, match="not open"):
            db.get_schema_version()


# ---------------------------------------------------------------------------
# get_commit_history – start_hash pagination
# ---------------------------------------------------------------------------

class TestGetCommitHistoryStartHash:
    def _populate(self, db: MetadataDB) -> list[str]:
        hashes = [chr(ord("a") + i) * 64 for i in range(5)]
        for i, h in enumerate(hashes):
            db.insert_commit(
                h,
                None,
                f"2026-01-01T10:{i:02d}:00Z",
                "user@host",
                f"Commit {i}",
            )
        return hashes

    def test_start_hash_returns_subset(self, db: MetadataDB) -> None:
        """History starting from a given hash includes that commit
        and all earlier ones."""
        hashes = self._populate(db)
        # Start from the 3rd commit (index 2; timestamp 10:02)
        history = db.get_commit_history(start_hash=hashes[2])
        # Should include commits at 10:00, 10:01, and 10:02
        assert len(history) == 3
        returned_hashes = {c["commit_hash"] for c in history}
        assert hashes[2] in returned_hashes
        assert hashes[0] in returned_hashes

    def test_start_hash_not_found_returns_empty(self, db: MetadataDB) -> None:
        """Unknown start_hash returns empty list."""
        self._populate(db)
        result = db.get_commit_history(start_hash="z" * 64)
        assert result == []

    def test_start_hash_with_limit(self, db: MetadataDB) -> None:
        """Combining start_hash and limit."""
        hashes = self._populate(db)
        history = db.get_commit_history(start_hash=hashes[4], limit=2)
        assert len(history) == 2

    def test_limit_no_start_hash(self, db: MetadataDB) -> None:
        """Limit without start_hash applies to full history."""
        self._populate(db)
        history = db.get_commit_history(limit=3)
        assert len(history) == 3


# ---------------------------------------------------------------------------
# get_schema_version – edge cases
# ---------------------------------------------------------------------------

class TestGetSchemaVersionExtra:
    def test_returns_zero_when_no_version_row(self, chemvcs_dir: Path) -> None:
        """Should return 0 if schema_version key is absent."""
        db = MetadataDB(chemvcs_dir)
        db.open()
        # Manually create metadata table without inserting schema_version
        db.conn.execute(  # type: ignore
            "CREATE TABLE IF NOT EXISTS metadata (key TEXT UNIQUE, value TEXT)"
        )
        db.conn.commit()  # type: ignore

        version = db.get_schema_version()
        assert version == 0
        db.close()

    def test_raises_database_error_when_not_open(self, chemvcs_dir: Path) -> None:
        db = MetadataDB(chemvcs_dir)
        with pytest.raises(DatabaseError):
            db.get_schema_version()


# ---------------------------------------------------------------------------
# init_schema – sqlite3.Error path
# ---------------------------------------------------------------------------

class TestInitSchemaSqliteError:
    def test_sqlite_error_raises_database_error(self, chemvcs_dir: Path) -> None:
        db = MetadataDB(chemvcs_dir)
        db.open()
        mock_conn = MagicMock()
        cursor_mock = MagicMock()
        cursor_mock.execute.side_effect = sqlite3.OperationalError("table locked")
        mock_conn.cursor.return_value = cursor_mock
        db.conn = mock_conn
        with pytest.raises(DatabaseError, match="Failed to initialize schema"):
            db.init_schema()


# ---------------------------------------------------------------------------
# insert_commit – generic sqlite3.Error path
# ---------------------------------------------------------------------------

class TestInsertCommitSqliteError:
    def test_generic_sqlite_error_raises_database_error(self, db: MetadataDB) -> None:
        mock_conn = MagicMock()
        cursor_mock = MagicMock()
        cursor_mock.execute.side_effect = sqlite3.OperationalError("disk full")
        mock_conn.cursor.return_value = cursor_mock
        db.conn = mock_conn
        with pytest.raises(DatabaseError, match="Failed to insert commit"):
            db.insert_commit("a" * 64, None, "2026-01-01T00:00:00Z", "u", "msg")


# ---------------------------------------------------------------------------
# get_commit_by_hash – sqlite3.Error path
# ---------------------------------------------------------------------------

class TestGetCommitByHashSqliteError:
    def test_sqlite_error_raises_database_error(self, db: MetadataDB) -> None:
        mock_conn = MagicMock()
        mock_conn.cursor.side_effect = sqlite3.OperationalError("disk error")
        db.conn = mock_conn
        with pytest.raises(DatabaseError, match="Failed to query commit"):
            db.get_commit_by_hash("abc1234")


# ---------------------------------------------------------------------------
# get_files_for_commit – sqlite3.Error path
# ---------------------------------------------------------------------------

class TestGetFilesForCommitSqliteError:
    def test_sqlite_error_raises_database_error(self, db: MetadataDB) -> None:
        mock_conn = MagicMock()
        mock_conn.cursor.side_effect = sqlite3.OperationalError("disk error")
        db.conn = mock_conn
        with pytest.raises(DatabaseError, match="Failed to get files"):
            db.get_files_for_commit(1)


# ---------------------------------------------------------------------------
# insert_commit – duplicate hash (IntegrityError path)
# ---------------------------------------------------------------------------

class TestInsertCommitDuplicateHash:
    def test_duplicate_hash_raises_database_error(self, db: MetadataDB) -> None:
        db.insert_commit("c" * 64, None, "2026-01-01T00:00:00Z", "u", "First")
        with pytest.raises(DatabaseError, match="already exists"):
            db.insert_commit("c" * 64, None, "2026-01-01T00:01:00Z", "u", "Dup")


# ---------------------------------------------------------------------------
# insert_file – sqlite3 error on foreign-key violation
# ---------------------------------------------------------------------------

class TestInsertFileSqliteError:
    def test_insert_file_invalid_commit_id_raises(self, db: MetadataDB) -> None:
        """Foreign key violation → DatabaseError."""
        with pytest.raises(DatabaseError, match="Failed to insert file"):
            db.insert_file(
                commit_id=99999,
                path="INCAR",
                blob_hash="d" * 64,
            )


# ---------------------------------------------------------------------------
# get_latest_commit – sqlite3.Error path
# ---------------------------------------------------------------------------

class TestGetLatestCommitSqliteError:
    def test_sqlite_error_raises_database_error(self, db: MetadataDB) -> None:
        mock_conn = MagicMock()
        mock_conn.cursor.side_effect = sqlite3.OperationalError("simulated")
        db.conn = mock_conn
        with pytest.raises(DatabaseError, match="Failed to get latest commit"):
            db.get_latest_commit()


# ---------------------------------------------------------------------------
# get_commit_history – sqlite3.Error path
# ---------------------------------------------------------------------------

class TestGetCommitHistorySqliteError:
    def test_sqlite_error_raises_database_error(self, db: MetadataDB) -> None:
        mock_conn = MagicMock()
        mock_conn.cursor.side_effect = sqlite3.OperationalError("simulated")
        db.conn = mock_conn
        with pytest.raises(DatabaseError, match="Failed to get history"):
            db.get_commit_history()


# ---------------------------------------------------------------------------
# get_schema_version – sqlite3.Error path
# ---------------------------------------------------------------------------

class TestGetSchemaVersionSqliteError:
    def test_sqlite_error_raises_database_error(self, db: MetadataDB) -> None:
        mock_conn = MagicMock()
        mock_conn.cursor.side_effect = sqlite3.OperationalError("simulated")
        db.conn = mock_conn
        with pytest.raises(DatabaseError, match="Failed to get schema version"):
            db.get_schema_version()
