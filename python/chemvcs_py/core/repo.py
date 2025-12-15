"""Repository discovery and CLI integration."""

import json
import os
import subprocess
from pathlib import Path
from typing import List, Optional

from ..util.errors import RepositoryNotFoundError, ObjectNotFoundError, CLIError
from .objects import CoreObject


class Repo:
    """
    Interface to a ChemVCS repository.
    
    This class provides methods to interact with the Go-based ChemVCS core
    via CLI commands, handling repository discovery, object inspection,
    and command execution.
    """

    def __init__(self, root: Optional[Path] = None):
        """
        Initialize repository interface.
        
        Args:
            root: Optional path to repository. If None, searches upwards
                  from current directory.
        
        Raises:
            RepositoryNotFoundError: If no repository found
        """
        if root is None:
            self._root = self._find_repository(Path.cwd())
        else:
            self._root = Path(root).resolve()
            if not self._is_repository(self._root):
                raise RepositoryNotFoundError(f"No repository at {self._root}")

    @property
    def root(self) -> Path:
        """Get repository root path."""
        return self._root

    @property
    def chemvcs_dir(self) -> Path:
        """Get .chemvcs directory path."""
        return self._root / ".chemvcs"

    def _find_repository(self, start_path: Path) -> Path:
        """
        Search upwards for a .chemvcs directory.
        
        Args:
            start_path: Starting directory for search
            
        Returns:
            Path to repository root
            
        Raises:
            RepositoryNotFoundError: If no repository found
        """
        current = start_path.resolve()
        
        while True:
            if self._is_repository(current):
                return current
            
            parent = current.parent
            if parent == current:  # Reached filesystem root
                raise RepositoryNotFoundError(
                    "Not a chemvcs repository (or any of the parent directories)"
                )
            current = parent

    def _is_repository(self, path: Path) -> bool:
        """Check if path contains a ChemVCS repository."""
        chemvcs_dir = path / ".chemvcs"
        return chemvcs_dir.is_dir()

    def _run_chemvcs(self, args: List[str], check: bool = True) -> subprocess.CompletedProcess:
        """
        Run chemvcs CLI command.
        
        Args:
            args: Command arguments
            check: Whether to raise exception on non-zero exit
            
        Returns:
            CompletedProcess result
            
        Raises:
            CLIError: If command fails and check=True
        """
        # Find chemvcs executable
        # For now, assume it's in PATH or in go/ directory
        chemvcs_paths = [
            "chemvcs",
            "chemvcs.exe",
            str(self._root.parent / "go" / "chemvcs.exe"),
            str(self._root.parent / "go" / "chemvcs"),
        ]
        
        chemvcs_cmd = None
        for path in chemvcs_paths:
            try:
                result = subprocess.run(
                    [path, "version"],
                    capture_output=True,
                    cwd=self._root,
                    timeout=5
                )
                if result.returncode == 0:
                    chemvcs_cmd = path
                    break
            except (FileNotFoundError, subprocess.TimeoutExpired):
                continue
        
        if chemvcs_cmd is None:
            raise CLIError("chemvcs executable not found in PATH or ../go/")
        
        try:
            result = subprocess.run(
                [chemvcs_cmd] + args,
                capture_output=True,
                text=True,
                cwd=self._root,
                check=check,
                timeout=30
            )
            return result
        except subprocess.CalledProcessError as e:
            raise CLIError(f"chemvcs command failed: {e.stderr}") from e
        except subprocess.TimeoutExpired as e:
            raise CLIError(f"chemvcs command timed out: {' '.join(args)}") from e

    def list_objects(self, type_filter: Optional[str] = None) -> List[dict]:
        """
        List all objects in the repository.
        
        Args:
            type_filter: Optional type to filter by (e.g., "file", "folder", "structure")
            
        Returns:
            List of dicts with 'Hash' and 'Type' keys
            
        Raises:
            CLIError: If command fails
        """
        args = ["list-objects", "--format=json"]
        if type_filter:
            args.extend(["--type", type_filter])
        
        result = self._run_chemvcs(args)
        
        # Handle empty result
        if not result.stdout.strip():
            return []
        
        return json.loads(result.stdout)

    def get_object(self, hash: str) -> CoreObject:
        """
        Get an object by its hash.
        
        Args:
            hash: Object hash (full or prefix)
            
        Returns:
            CoreObject instance
            
        Raises:
            ObjectNotFoundError: If object not found
            CLIError: If command fails
        """
        args = ["inspect-object", "--format=json", hash]
        
        try:
            result = self._run_chemvcs(args)
        except CLIError as e:
            if "not found" in str(e).lower():
                raise ObjectNotFoundError(f"Object not found: {hash}") from e
            raise
        
        data = json.loads(result.stdout)
        return CoreObject.from_dict(data, hash=hash)

    def commit(self, message: str) -> str:
        """
        Create a new commit.
        
        Args:
            message: Commit message
            
        Returns:
            Snapshot hash
            
        Raises:
            CLIError: If commit fails
        """
        result = self._run_chemvcs(["commit", "-m", message])
        
        # Parse output like "[main abc123] message"
        output = result.stdout.strip()
        if "[" in output and "]" in output:
            parts = output.split("]")[0].split()
            if len(parts) >= 2:
                return parts[1]  # The hash
        
        raise CLIError(f"Failed to parse commit output: {output}")

    def get_current_branch(self) -> str:
        """
        Get current branch name.
        
        Returns:
            Branch name
        """
        # Read HEAD file
        head_file = self.chemvcs_dir / "HEAD"
        if not head_file.exists():
            raise CLIError("HEAD file not found")
        
        content = head_file.read_text().strip()
        if content.startswith("ref: refs/heads/"):
            return content[16:]  # Remove "ref: refs/heads/"
        return "detached"
