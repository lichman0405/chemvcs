"""Commit object builder and serializer.

This module handles the creation and serialization of commit objects,
which represent snapshots of the computational workspace at a point in time.
"""

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

from chemvcs.constants import COMMITS_DIR
from chemvcs.storage.metadata_db import MetadataDB
from chemvcs.storage.object_store import ObjectStore


class CommitBuilderError(Exception):
    """Exception raised during commit building."""


class CommitBuilder:
    """Builder for creating and persisting commit objects.
    
    A commit represents a snapshot of tracked files with semantic metadata.
    Commits are stored as JSON files in `.chemvcs/commits/<hash>`.
    
    Attributes:
        chemvcs_dir: Path to .chemvcs directory
        object_store: ObjectStore instance for blob management
        metadata_db: MetadataDB instance for indexing
    """

    def __init__(
        self,
        chemvcs_dir: Path,
        object_store: ObjectStore,
        metadata_db: MetadataDB,
    ):
        """Initialize CommitBuilder.
        
        Args:
            chemvcs_dir: Path to .chemvcs directory
            object_store: ObjectStore for blob storage
            metadata_db: MetadataDB for indexing
        """
        self.chemvcs_dir = Path(chemvcs_dir)
        self.commits_dir = self.chemvcs_dir / COMMITS_DIR
        self.object_store = object_store
        self.metadata_db = metadata_db
        
        # Ensure commits directory exists
        self.commits_dir.mkdir(parents=True, exist_ok=True)

    def create_commit(
        self,
        files: List[Dict[str, Any]],
        message: str,
        author: str,
        parent_hash: Optional[str] = None,
        semantic_data: Optional[Dict[str, Any]] = None,
        output_data: Optional[Dict[str, Any]] = None,
        environment: Optional[Dict[str, str]] = None,
    ) -> str:
        """Create a new commit object.
        
        Args:
            files: List of file dictionaries with keys:
                - path: File path relative to workspace root
                - content: File content as bytes
                - file_type: Type identifier (INCAR, POSCAR, etc.)
                - is_reference: Whether file is stored as reference (POTCAR)
                - size_bytes: Optional file size
            message: Commit message
            author: Author identifier (e.g., "user@hostname")
            parent_hash: Hash of parent commit, or None for first commit
            semantic_data: Semantic summary from VASP parsers
            output_data: Output summary (energies, forces, etc.)
            environment: Environment variables/settings
        
        Returns:
            Commit hash (SHA-256 hex string)
        
        Raises:
            CommitBuilderError: If commit creation fails
        """
        try:
            # Generate timestamp
            timestamp = datetime.now(timezone.utc).isoformat()
            
            # Process files - store blobs and collect metadata
            processed_files = []
            for file_entry in files:
                path = file_entry["path"]
                content = file_entry["content"]
                file_type = file_entry.get("file_type", "unknown")
                is_reference = file_entry.get("is_reference", False)
                size_bytes = file_entry.get("size_bytes", len(content))
                
                # Write blob to object store
                blob_hash = self.object_store.write_blob(content)
                
                processed_files.append({
                    "path": path,
                    "blob_hash": blob_hash,
                    "file_type": file_type,
                    "is_reference": is_reference,
                    "size_bytes": size_bytes,
                })
            
            # Build commit object
            commit_obj = {
                "parent": parent_hash,
                "timestamp": timestamp,
                "author": author,
                "message": message,
                "files": processed_files,
            }
            
            # Add optional metadata
            if semantic_data:
                commit_obj["semantic_summary"] = semantic_data
            if output_data:
                commit_obj["output_summary"] = output_data
            if environment:
                commit_obj["environment"] = environment
            
            # Compute commit hash
            commit_hash = self._compute_commit_hash(commit_obj)
            commit_obj["hash"] = commit_hash
            
            # Persist commit to disk
            self._write_commit_file(commit_hash, commit_obj)
            
            # Index in metadata database
            self._index_commit(commit_obj)
            
            return commit_hash
            
        except Exception as e:
            raise CommitBuilderError(f"Failed to create commit: {e}") from e

    def read_commit(self, commit_hash: str) -> Dict[str, Any]:
        """Read a commit object from disk.
        
        Args:
            commit_hash: Commit hash (full or short)
        
        Returns:
            Commit object dictionary
        
        Raises:
            CommitBuilderError: If commit not found or corrupted
        """
        try:
            # Resolve short hash via database
            if len(commit_hash) < 64:
                db_commit = self.metadata_db.get_commit_by_hash(commit_hash)
                if db_commit is None:
                    raise CommitBuilderError(f"Commit not found: {commit_hash}")
                commit_hash = db_commit["commit_hash"]
            
            commit_path = self.commits_dir / commit_hash
            
            if not commit_path.exists():
                raise CommitBuilderError(f"Commit file not found: {commit_hash}")
            
            with open(commit_path, "r", encoding="utf-8") as f:
                commit_obj = json.load(f)
            
            # Verify hash integrity
            stored_hash = commit_obj.get("hash")
            if stored_hash != commit_hash:
                raise CommitBuilderError(
                    f"Commit hash mismatch: expected {commit_hash}, got {stored_hash}"
                )
            
            return commit_obj
            
        except CommitBuilderError:
            raise
        except Exception as e:
            raise CommitBuilderError(f"Failed to read commit: {e}") from e

    def commit_exists(self, commit_hash: str) -> bool:
        """Check if a commit exists.
        
        Args:
            commit_hash: Commit hash to check
        
        Returns:
            True if commit exists, False otherwise
        """
        commit_path = self.commits_dir / commit_hash
        return commit_path.exists()

    def _compute_commit_hash(self, commit_obj: Dict[str, Any]) -> str:
        """Compute SHA-256 hash of commit object.
        
        Hash is computed over the canonical JSON representation
        (sorted keys, no whitespace) of the commit object WITHOUT
        the 'hash' field itself.
        
        Args:
            commit_obj: Commit object dictionary (without 'hash' field)
        
        Returns:
            SHA-256 hash as hex string
        """
        # Create copy without 'hash' field
        obj_for_hash = {k: v for k, v in commit_obj.items() if k != "hash"}
        
        # Canonical JSON: sorted keys, no whitespace
        canonical_json = json.dumps(
            obj_for_hash,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
        )
        
        # Compute hash
        hasher = hashlib.sha256()
        hasher.update(canonical_json.encode("utf-8"))
        return hasher.hexdigest()

    def _write_commit_file(self, commit_hash: str, commit_obj: Dict[str, Any]) -> None:
        """Write commit object to disk as JSON.
        
        Args:
            commit_hash: Commit hash (filename)
            commit_obj: Commit object to serialize
        
        Raises:
            CommitBuilderError: If write fails
        """
        commit_path = self.commits_dir / commit_hash
        
        # Pretty-print JSON for human readability
        json_str = json.dumps(commit_obj, indent=2, ensure_ascii=False)
        
        try:
            # Atomic write via temp file
            import tempfile
            import os
            
            fd, tmp_path = tempfile.mkstemp(
                dir=self.commits_dir,
                prefix=".tmp_commit_",
                suffix=".json",
            )
            
            try:
                with os.fdopen(fd, "w", encoding="utf-8") as f:
                    f.write(json_str)
                    f.flush()
                    os.fsync(f.fileno())
                
                # Atomic rename
                os.replace(tmp_path, commit_path)
                
            except Exception:
                # Clean up temp file on error
                try:
                    os.unlink(tmp_path)
                except OSError:
                    pass
                raise
                
        except Exception as e:
            raise CommitBuilderError(f"Failed to write commit file: {e}") from e

    def _index_commit(self, commit_obj: Dict[str, Any]) -> None:
        """Index commit in metadata database.
        
        Args:
            commit_obj: Commit object to index
        
        Raises:
            CommitBuilderError: If indexing fails
        """
        try:
            # Insert commit record
            commit_id = self.metadata_db.insert_commit(
                commit_hash=commit_obj["hash"],
                parent_hash=commit_obj.get("parent"),
                timestamp=commit_obj["timestamp"],
                author=commit_obj["author"],
                message=commit_obj["message"],
            )
            
            # Insert file records
            for file_entry in commit_obj["files"]:
                self.metadata_db.insert_file(
                    commit_id=commit_id,
                    path=file_entry["path"],
                    blob_hash=file_entry["blob_hash"],
                    is_reference=file_entry.get("is_reference", False),
                    file_type=file_entry.get("file_type"),
                    size_bytes=file_entry.get("size_bytes"),
                )
            
            # TODO: Index semantic_summary, output_summary, environment
            # when those tables are fully integrated
            
        except Exception as e:
            raise CommitBuilderError(f"Failed to index commit: {e}") from e
