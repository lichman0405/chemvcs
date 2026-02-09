"""Staging area management for ChemVCS.

The staging area (index) tracks which files should be included in the next commit.
Unlike Git, we use a simpler JSON-based index format.
"""

import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

from chemvcs.constants import CHEMVCS_DIR
from chemvcs.storage import ObjectStore


class StagingError(Exception):
    """Exception raised during staging operations."""


class StagingManager:
    """Manager for the staging area (index).
    
    The staging area tracks files to be included in the next commit.
    Index format (JSON):
    {
        "version": 1,
        "entries": {
            "relative/path/to/file": {
                "blob_hash": "sha256...",
                "size_bytes": 1024,
                "mtime": "2026-02-09T10:00:00Z",
                "file_type": "INCAR"  // detected type
            },
            ...
        }
    }
    
    Attributes:
        workspace_root: Root directory of the workspace
        index_path: Path to the index file (.chemvcs/index)
        object_store: ObjectStore instance for blob storage
    """

    def __init__(self, workspace_root: Path, object_store: ObjectStore):
        """Initialize StagingManager.
        
        Args:
            workspace_root: Root directory of workspace
            object_store: ObjectStore for blob management
        """
        self.workspace_root = Path(workspace_root)
        self.chemvcs_dir = self.workspace_root / CHEMVCS_DIR
        self.index_path = self.chemvcs_dir / "index"
        self.object_store = object_store
        
        # Ensure .chemvcs exists
        if not self.chemvcs_dir.exists():
            raise StagingError(
                f"Not a ChemVCS repository (no .chemvcs/ found in {workspace_root})"
            )

    def add(self, paths: List[Path], force: bool = False) -> Dict[str, Any]:
        """Add files to the staging area.
        
        Args:
            paths: List of paths (absolute or relative) to add
            force: Override .chemvcsignore rules
        
        Returns:
            Dictionary with statistics:
            {
                "added": ["file1", "file2"],
                "updated": ["file3"],
                "ignored": ["file4"],
                "errors": ["file5: reason"]
            }
        
        Raises:
            StagingError: If staging fails critically
        """
        # Load current index
        index = self._load_index()
        
        # Load ignore patterns
        ignore_patterns = self._load_ignore_patterns() if not force else []
        
        stats = {
            "added": [],
            "updated": [],
            "ignored": [],
            "errors": [],
        }
        
        for path in paths:
            try:
                # Normalize to absolute path
                abs_path = self._resolve_path(path)
                
                # Check if path exists
                if not abs_path.exists():
                    stats["errors"].append(f"{path}: file not found")
                    continue
                
                # Handle directory recursively
                if abs_path.is_dir():
                    dir_stats = self._add_directory(abs_path, ignore_patterns, force)
                    # Merge stats
                    for key in stats:
                        stats[key].extend(dir_stats[key])
                    continue
                
                # Skip .chemvcs directory itself
                if self._is_chemvcs_path(abs_path):
                    continue
                
                # Get relative path (use POSIX format for cross-platform compatibility)
                rel_path = abs_path.relative_to(self.workspace_root)
                rel_path_str = rel_path.as_posix()
                
                # Check ignore patterns
                if self._should_ignore(rel_path, ignore_patterns):
                    stats["ignored"].append(rel_path_str)
                    continue
                
                # Read file content
                content = abs_path.read_bytes()
                
                # Write blob to object store
                blob_hash = self.object_store.write_blob(content)
                
                # Detect file type
                file_type = self._detect_file_type(abs_path)
                
                # Get file metadata
                mtime = abs_path.stat().st_mtime
                import datetime
                mtime_iso = datetime.datetime.fromtimestamp(
                    mtime, tz=datetime.timezone.utc
                ).isoformat()
                
                # Check if file already staged
                was_staged = rel_path_str in index["entries"]
                
                # Update index
                index["entries"][rel_path_str] = {
                    "blob_hash": blob_hash,
                    "size_bytes": len(content),
                    "mtime": mtime_iso,
                    "file_type": file_type,
                }
                
                # Update stats
                if was_staged:
                    stats["updated"].append(rel_path_str)
                else:
                    stats["added"].append(rel_path_str)
                    
            except Exception as e:
                stats["errors"].append(f"{path}: {e}")
        
        # Save updated index
        self._save_index(index)
        
        return stats

    def remove(self, paths: List[Path]) -> Dict[str, Any]:
        """Remove files from the staging area.
        
        Args:
            paths: List of paths to remove from staging
        
        Returns:
            Dictionary with statistics:
            {
                "removed": ["file1", "file2"],
                "not_staged": ["file3"],
                "errors": ["file4: reason"]
            }
        """
        index = self._load_index()
        
        stats = {
            "removed": [],
            "not_staged": [],
            "errors": [],
        }
        
        for path in paths:
            try:
                abs_path = self._resolve_path(path)
                rel_path = abs_path.relative_to(self.workspace_root)
                rel_path_str = rel_path.as_posix()
                
                if rel_path_str in index["entries"]:
                    del index["entries"][rel_path_str]
                    stats["removed"].append(rel_path_str)
                else:
                    stats["not_staged"].append(rel_path_str)
                    
            except Exception as e:
                stats["errors"].append(f"{path}: {e}")
        
        self._save_index(index)
        return stats

    def get_staged_files(self) -> Dict[str, Dict[str, Any]]:
        """Get all currently staged files.
        
        Returns:
            Dictionary mapping relative paths to file metadata
        """
        index = self._load_index()
        return index["entries"]

    def clear(self) -> None:
        """Clear all staged files."""
        self._save_index({"version": 1, "entries": {}})

    def is_empty(self) -> bool:
        """Check if staging area is empty."""
        index = self._load_index()
        return len(index["entries"]) == 0

    def _load_index(self) -> Dict[str, Any]:
        """Load index from disk."""
        if not self.index_path.exists():
            return {"version": 1, "entries": {}}
        
        try:
            with open(self.index_path, "r", encoding="utf-8") as f:
                index = json.load(f)
            
            # Validate version
            if index.get("version") != 1:
                raise StagingError(f"Unsupported index version: {index.get('version')}")
            
            return index
            
        except json.JSONDecodeError as e:
            raise StagingError(f"Corrupted index file: {e}") from e

    def _save_index(self, index: Dict[str, Any]) -> None:
        """Save index to disk."""
        # Atomic write via temp file
        import tempfile
        import os
        
        fd, tmp_path = tempfile.mkstemp(
            dir=self.chemvcs_dir,
            prefix=".tmp_index_",
            suffix=".json",
        )
        
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as f:
                json.dump(index, f, indent=2, ensure_ascii=False)
                f.flush()
                os.fsync(f.fileno())
            
            # Atomic rename
            os.replace(tmp_path, self.index_path)
            
        except Exception:
            # Clean up temp file on error
            try:
                os.unlink(tmp_path)
            except OSError:
                pass
            raise

    def _resolve_path(self, path: Path) -> Path:
        """Resolve path to absolute path within workspace."""
        if path.is_absolute():
            abs_path = path
        else:
            abs_path = (self.workspace_root / path).resolve()
        
        # Verify path is within workspace
        try:
            abs_path.relative_to(self.workspace_root)
        except ValueError:
            raise StagingError(
                f"Path {path} is outside workspace root {self.workspace_root}"
            )
        
        return abs_path

    def _is_chemvcs_path(self, abs_path: Path) -> bool:
        """Check if path is within .chemvcs directory."""
        try:
            abs_path.relative_to(self.chemvcs_dir)
            return True
        except ValueError:
            return False

    def _add_directory(
        self, dir_path: Path, ignore_patterns: List[str], force: bool
    ) -> Dict[str, Any]:
        """Recursively add all files in directory."""
        stats = {
            "added": [],
            "updated": [],
            "ignored": [],
            "errors": [],
        }
        
        # Get all files recursively
        try:
            for item in dir_path.rglob("*"):
                if item.is_file():
                    # Skip .chemvcs
                    if self._is_chemvcs_path(item):
                        continue
                    
                    # Delegate to add() to reuse logic
                    item_stats = self.add([item], force=force)
                    for key in stats:
                        stats[key].extend(item_stats[key])
                        
        except Exception as e:
            stats["errors"].append(f"{dir_path}: {e}")
        
        return stats

    def _load_ignore_patterns(self) ->List[str]:
        """Load patterns from .chemvcsignore file."""
        ignore_file = self.workspace_root / ".chemvcsignore"
        
        if not ignore_file.exists():
            return []
        
        patterns = []
        try:
            content = ignore_file.read_text(encoding="utf-8")
            for line in content.splitlines():
                line = line.strip()
                # Skip empty lines and comments
                if line and not line.startswith("#"):
                    patterns.append(line)
        except Exception:
            # If can't read ignore file, just continue without it
            pass
        
        return patterns

    def _should_ignore(self, rel_path: Path, patterns: List[str]) -> bool:
        """Check if path matches any ignore pattern."""
        import fnmatch
        
        path_str = rel_path.as_posix()
        
        for pattern in patterns:
            # Support directory patterns ending with /
            if pattern.endswith("/"):
                # Match if path starts with this directory
                if path_str.startswith(pattern.rstrip("/")):
                    return True
            else:
                # fnmatch for glob patterns
                if fnmatch.fnmatch(path_str, pattern):
                    return True
                # Also check individual path components
                if fnmatch.fnmatch(rel_path.name, pattern):
                    return True
        
        return False

    def _detect_file_type(self, path: Path) -> str:
        """Detect file type based on name and content.
        
        For now, simple name-based detection. Will be enhanced
        with VASP parsers in Phase 1d.
        """
        name = path.name.upper()
        
        # VASP input files
        if name == "INCAR":
            return "INCAR"
        elif name == "POSCAR" or name == "CONTCAR":
            return "POSCAR"
        elif name == "KPOINTS":
            return "KPOINTS"
        elif name == "POTCAR":
            return "POTCAR"
        
        # VASP output files
        elif name == "OUTCAR":
            return "OUTCAR"
        elif name == "VASPRUN.XML":
            return "VASPRUN"
        elif name == "OSZICAR":
            return "OSZICAR"
        
        # Other
        elif path.suffix == ".json":
            return "JSON"
        elif path.suffix in [".txt", ".log"]:
            return "TEXT"
        else:
            return "UNKNOWN"
