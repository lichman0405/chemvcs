"""SQLite metadata database for ChemVCS.

This module manages commit metadata, file references, and semantic summaries
in a SQLite database. The database serves as a rebuildable index - the true
source of truth is the commits/ directory.
"""

import json
import sqlite3
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from chemvcs.constants import DB_SCHEMA_VERSION, METADATA_DB


class DatabaseError(Exception):
    """Base exception for database errors."""

    pass


class MetadataDB:
    """SQLite database manager for ChemVCS metadata.

    Manages commit history, file references, semantic summaries, and environment
    information. The database is designed to be rebuildable from commit JSON files.

    Schema Tables:
        - commits: Commit records with hash, parent, timestamp, message
        - files: File list for each commit
        - semantic_summary: Parsed INCAR/POSCAR parameters
        - output_summary: OUTCAR energy and convergence status
        - environment: Runtime environment (VASP version, modules)
        - metadata: Schema version and configuration

    Attributes:
        db_path: Path to the SQLite database file
        conn: Active database connection (if open)

    Example:
        >>> db = MetadataDB(Path(".chemvcs"))
        >>> db.init_schema()
        >>> db.insert_commit("abc123...", None, "2026-02-09T10:00:00Z", "user@host", "Initial")
    """

    def __init__(self, chemvcs_dir: Path) -> None:
        """Initialize database manager.

        Args:
            chemvcs_dir: Path to .chemvcs directory
        """
        self.chemvcs_dir = Path(chemvcs_dir)
        self.db_path = self.chemvcs_dir / METADATA_DB
        self.conn: Optional[sqlite3.Connection] = None
        self._wal_mode_supported: Optional[bool] = None

    def open(self) -> None:
        """Open database connection and configure journal mode.

        Attempts to use WAL mode for better concurrency. Falls back to
        DELETE mode if WAL is not supported (e.g., on NFS).

        Raises:
            DatabaseError: If connection fails
        """
        if self.conn is not None:
            return  # Already open

        try:
            self.conn = sqlite3.connect(
                str(self.db_path),
                timeout=30.0,  # Wait up to 30s for locks
                check_same_thread=False,
            )
            self.conn.row_factory = sqlite3.Row  # Access columns by name

            # Try WAL mode first
            if self._wal_mode_supported is None:
                self._detect_wal_support()

            if self._wal_mode_supported:
                self.conn.execute("PRAGMA journal_mode=WAL")
            else:
                self.conn.execute("PRAGMA journal_mode=DELETE")

            # Enable foreign keys
            self.conn.execute("PRAGMA foreign_keys=ON")

        except sqlite3.Error as e:
            raise DatabaseError(f"Failed to open database: {e}") from e

    def close(self) -> None:
        """Close database connection."""
        if self.conn is not None:
            self.conn.close()
            self.conn = None

    def __enter__(self) -> "MetadataDB":
        """Context manager entry."""
        self.open()
        return self

    def __exit__(self, *args: Any) -> None:
        """Context manager exit."""
        self.close()

    def _detect_wal_support(self) -> None:
        """Detect if WAL journal mode is supported.

        WAL may not work on network filesystems like NFS or some Lustre configs.
        """
        try:
            cursor = self.conn.execute("PRAGMA journal_mode=WAL")  # type: ignore
            result = cursor.fetchone()
            self._wal_mode_supported = result[0].upper() == "WAL"
        except sqlite3.Error:
            self._wal_mode_supported = False

    def init_schema(self) -> None:
        """Initialize database schema.

        Creates all tables and indices. Safe to call on existing database
        (uses IF NOT EXISTS).

        Raises:
            DatabaseError: If schema creation fails
        """
        if self.conn is None:
            raise DatabaseError("Database not open")

        try:
            cursor = self.conn.cursor()

            # Metadata table (schema version, etc.)
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS metadata (
                    key TEXT PRIMARY KEY,
                    value TEXT NOT NULL
                )
            """)

            # Commits table
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS commits (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    commit_hash TEXT UNIQUE NOT NULL,
                    parent_hash TEXT,
                    timestamp TEXT NOT NULL,
                    author TEXT NOT NULL,
                    message TEXT NOT NULL,
                    created_at DATETIME DEFAULT CURRENT_TIMESTAMP
                )
            """)
            cursor.execute("""
                CREATE INDEX IF NOT EXISTS idx_commits_hash 
                ON commits(commit_hash)
            """)
            cursor.execute("""
                CREATE INDEX IF NOT EXISTS idx_commits_parent 
                ON commits(parent_hash)
            """)

            # Files table
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS files (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    commit_id INTEGER NOT NULL,
                    path TEXT NOT NULL,
                    blob_hash TEXT NOT NULL,
                    is_reference BOOLEAN NOT NULL DEFAULT 0,
                    file_type TEXT,
                    size_bytes INTEGER,
                    FOREIGN KEY (commit_id) REFERENCES commits(id) ON DELETE CASCADE
                )
            """)
            cursor.execute("""
                CREATE INDEX IF NOT EXISTS idx_files_commit 
                ON files(commit_id)
            """)
            cursor.execute("""
                CREATE INDEX IF NOT EXISTS idx_files_path 
                ON files(path)
            """)

            # Semantic summary table
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS semantic_summary (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    file_id INTEGER NOT NULL,
                    summary_json TEXT NOT NULL,
                    FOREIGN KEY (file_id) REFERENCES files(id) ON DELETE CASCADE
                )
            """)
            cursor.execute("""
                CREATE INDEX IF NOT EXISTS idx_semantic_file 
                ON semantic_summary(file_id)
            """)

            # Output summary table
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS output_summary (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    commit_id INTEGER NOT NULL,
                    total_energy_eV REAL,
                    is_converged BOOLEAN,
                    ionic_steps INTEGER,
                    max_force REAL,
                    warnings TEXT,
                    summary_json TEXT,
                    FOREIGN KEY (commit_id) REFERENCES commits(id) ON DELETE CASCADE
                )
            """)
            cursor.execute("""
                CREATE INDEX IF NOT EXISTS idx_output_commit 
                ON output_summary(commit_id)
            """)

            # Environment table
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS environment (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    commit_id INTEGER NOT NULL,
                    hostname TEXT,
                    vasp_version TEXT,
                    modules TEXT,
                    python_version TEXT,
                    env_json TEXT,
                    FOREIGN KEY (commit_id) REFERENCES commits(id) ON DELETE CASCADE
                )
            """)
            cursor.execute("""
                CREATE INDEX IF NOT EXISTS idx_env_commit 
                ON environment(commit_id)
            """)

            # Set schema version
            cursor.execute(
                "INSERT OR REPLACE INTO metadata (key, value) VALUES (?, ?)",
                ("schema_version", str(DB_SCHEMA_VERSION)),
            )

            self.conn.commit()

        except sqlite3.Error as e:
            self.conn.rollback()
            raise DatabaseError(f"Failed to initialize schema: {e}") from e

    def insert_commit(
        self,
        commit_hash: str,
        parent_hash: Optional[str],
        timestamp: str,
        author: str,
        message: str,
    ) -> int:
        """Insert a new commit record.

        Args:
            commit_hash: SHA-256 hash of commit (64 hex chars)
            parent_hash: Parent commit hash (None for initial commit)
            timestamp: ISO 8601 timestamp
            author: Author identifier (e.g., "user@hostname")
            message: Commit message

        Returns:
            Database row ID of inserted commit

        Raises:
            DatabaseError: If insert fails (e.g., duplicate hash)
        """
        if self.conn is None:
            raise DatabaseError("Database not open")

        try:
            cursor = self.conn.cursor()
            cursor.execute(
                """
                INSERT INTO commits (commit_hash, parent_hash, timestamp, author, message)
                VALUES (?, ?, ?, ?, ?)
                """,
                (commit_hash, parent_hash, timestamp, author, message),
            )
            self.conn.commit()
            return cursor.lastrowid  # type: ignore

        except sqlite3.IntegrityError as e:
            self.conn.rollback()
            raise DatabaseError(f"Commit hash already exists: {commit_hash}") from e
        except sqlite3.Error as e:
            self.conn.rollback()
            raise DatabaseError(f"Failed to insert commit: {e}") from e

    def insert_file(
        self,
        commit_id: int,
        path: str,
        blob_hash: str,
        is_reference: bool = False,
        file_type: Optional[str] = None,
        size_bytes: Optional[int] = None,
    ) -> int:
        """Insert a file record for a commit.

        Args:
            commit_id: Database ID of parent commit
            path: Relative file path
            blob_hash: SHA-256 hash of file content
            is_reference: True for POTCAR reference-only tracking
            file_type: File type (INCAR, POSCAR, etc.)
            size_bytes: Original file size

        Returns:
            Database row ID of inserted file
        """
        if self.conn is None:
            raise DatabaseError("Database not open")

        try:
            cursor = self.conn.cursor()
            cursor.execute(
                """
                INSERT INTO files (commit_id, path, blob_hash, is_reference, file_type, size_bytes)
                VALUES (?, ?, ?, ?, ?, ?)
                """,
                (commit_id, path, blob_hash, is_reference, file_type, size_bytes),
            )
            self.conn.commit()
            return cursor.lastrowid  # type: ignore

        except sqlite3.Error as e:
            self.conn.rollback()
            raise DatabaseError(f"Failed to insert file: {e}") from e

    def get_commit_by_hash(self, commit_hash: str) -> Optional[Dict[str, Any]]:
        """Retrieve a commit by its hash.

        Args:
            commit_hash: Full or short commit hash (min 7 chars)

        Returns:
            Dictionary with commit data, or None if not found
        """
        if self.conn is None:
            raise DatabaseError("Database not open")

        try:
            cursor = self.conn.cursor()

            # Support short hashes (7+ chars)
            if len(commit_hash) < 64:
                cursor.execute(
                    "SELECT * FROM commits WHERE commit_hash LIKE ? || '%' LIMIT 1",
                    (commit_hash,),
                )
            else:
                cursor.execute(
                    "SELECT * FROM commits WHERE commit_hash = ?",
                    (commit_hash,),
                )

            row = cursor.fetchone()
            if row is None:
                return None

            return dict(row)

        except sqlite3.Error as e:
            raise DatabaseError(f"Failed to query commit: {e}") from e

    def get_latest_commit(self) -> Optional[Dict[str, Any]]:
        """Get the most recent commit.

        Returns:
            Dictionary with latest commit data, or None if empty
        """
        if self.conn is None:
            raise DatabaseError("Database not open")

        try:
            cursor = self.conn.cursor()
            cursor.execute("SELECT * FROM commits ORDER BY timestamp DESC LIMIT 1")
            row = cursor.fetchone()
            if row is None:
                return None
            return dict(row)

        except sqlite3.Error as e:
            raise DatabaseError(f"Failed to get latest commit: {e}") from e

    def get_files_for_commit(self, commit_id: int) -> List[Dict[str, Any]]:
        """Get all files for a commit.

        Args:
            commit_id: Database ID of commit

        Returns:
            List of file dictionaries
        """
        if self.conn is None:
            raise DatabaseError("Database not open")

        try:
            cursor = self.conn.cursor()
            cursor.execute(
                "SELECT * FROM files WHERE commit_id = ? ORDER BY path",
                (commit_id,),
            )
            rows = cursor.fetchall()
            return [dict(row) for row in rows]

        except sqlite3.Error as e:
            raise DatabaseError(f"Failed to get files: {e}") from e

    def get_commit_history(
        self,
        limit: Optional[int] = None,
        start_hash: Optional[str] = None,
    ) -> List[Dict[str, Any]]:
        """Get commit history in reverse chronological order.

        Args:
            limit: Maximum number of commits to return
            start_hash: Start from this commit (for pagination)

        Returns:
            List of commit dictionaries
        """
        if self.conn is None:
            raise DatabaseError("Database not open")

        try:
            cursor = self.conn.cursor()

            if start_hash:
                # Find start commit's timestamp
                cursor.execute(
                    "SELECT timestamp FROM commits WHERE commit_hash = ?",
                    (start_hash,),
                )
                row = cursor.fetchone()
                if row is None:
                    return []
                start_time = row[0]

                query = "SELECT * FROM commits WHERE timestamp <= ? ORDER BY timestamp DESC"
                params: Tuple[Any, ...] = (start_time,)
            else:
                query = "SELECT * FROM commits ORDER BY timestamp DESC"
                params = ()

            if limit:
                query += f" LIMIT {limit}"

            cursor.execute(query, params)
            rows = cursor.fetchall()
            return [dict(row) for row in rows]

        except sqlite3.Error as e:
            raise DatabaseError(f"Failed to get history: {e}") from e

    def get_schema_version(self) -> int:
        """Get the database schema version.

        Returns:
            Schema version number
        """
        if self.conn is None:
            raise DatabaseError("Database not open")

        try:
            cursor = self.conn.cursor()
            cursor.execute("SELECT value FROM metadata WHERE key = 'schema_version'")
            row = cursor.fetchone()
            if row is None:
                return 0
            return int(row[0])

        except sqlite3.Error as e:
            raise DatabaseError(f"Failed to get schema version: {e}") from e
