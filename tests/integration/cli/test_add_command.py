"""Integration tests for 'chemvcs add' command."""

import subprocess
from pathlib import Path

import pytest


class TestAddCommand:
    """Test suite for 'chemvcs add' command."""

    def test_add_not_initialized(self, tmp_path):
        """Test 'add' command in uninitialized directory."""
        # Change to temp directory
        result = subprocess.run(
            ["chemvcs", "add", "file.txt"],
            cwd=tmp_path,
            capture_output=True,
            text=True,
        )
        
        assert result.returncode == 1
        assert "Not a ChemVCS repository" in result.stdout

    def test_add_single_file(self, initialized_repo):
        """Test adding a single file."""
        workspace = initialized_repo
        
        # Create a test file
        test_file = workspace / "INCAR"
        test_file.write_text("ENCUT = 400\nNSW = 100\n")
        
        # Run add command
        result = subprocess.run(
            ["chemvcs", "add", "INCAR"],
            cwd=workspace,
            capture_output=True,
            text=True,
        )
        
        assert result.returncode == 0
        assert "INCAR" in result.stdout
        assert "staged for commit" in result.stdout
        
        # Verify staging area
        import json
        index_path = workspace / ".chemvcs" / "index"
        assert index_path.exists()
        
        with open(index_path) as f:
            index = json.load(f)
        
        assert "INCAR" in index["entries"]
        assert index["entries"]["INCAR"]["file_type"] == "INCAR"

    def test_add_multiple_files(self, initialized_repo):
        """Test adding multiple files at once."""
        workspace = initialized_repo
        
        # Create test files
        (workspace / "INCAR").write_text("ENCUT = 400\n")
        (workspace / "POSCAR").write_text("Test Structure\n1.0\n")
        (workspace / "KPOINTS").write_text("Automatic\n0\nGamma\n4 4 4\n")
        
        # Run add command
        result = subprocess.run(
            ["chemvcs", "add", "INCAR", "POSCAR", "KPOINTS"],
            cwd=workspace,
            capture_output=True,
            text=True,
        )
        
        assert result.returncode == 0
        assert "INCAR" in result.stdout
        assert "POSCAR" in result.stdout
        assert "KPOINTS" in result.stdout
        assert "3 file(s) staged" in result.stdout

    def test_add_directory(self, initialized_repo):
        """Test adding a directory recursively."""
        workspace = initialized_repo
        
        # Create test directory with files
        subdir = workspace / "inputs"
        subdir.mkdir()
        (subdir / "INCAR").write_text("ENCUT = 400\n")
        (subdir / "POSCAR").write_text("Structure\n")
        
        # Run add command
        result = subprocess.run(
            ["chemvcs", "add", "inputs"],
            cwd=workspace,
            capture_output=True,
            text=True,
        )
        
        assert result.returncode == 0
        assert "inputs/INCAR" in result.stdout
        assert "inputs/POSCAR" in result.stdout

    def test_add_with_ignore(self, initialized_repo):
        """Test adding files with .chemvcsignore."""
        workspace = initialized_repo
        
        # Create .chemvcsignore
        ignore_file = workspace / ".chemvcsignore"
        ignore_file.write_text("*.tmp\nWAVECAR\n")
        
        # Create test files
        (workspace / "INCAR").write_text("ENCUT = 400\n")
        (workspace / "data.tmp").write_text("Temporary\n")
        (workspace / "WAVECAR").write_text("Large wave function\n")
        
        # Run add command
        result = subprocess.run(
            ["chemvcs", "add", "INCAR", "data.tmp", "WAVECAR"],
            cwd=workspace,
            capture_output=True,
            text=True,
        )
        
        assert result.returncode == 0
        assert "INCAR" in result.stdout
        assert "data.tmp" in result.stdout
        assert "WAVECAR" in result.stdout
        assert "Ignored" in result.stdout

    def test_add_with_force(self, initialized_repo):
        """Test adding ignored files with --force."""
        workspace = initialized_repo
        
        # Create .chemvcsignore
        ignore_file = workspace / ".chemvcsignore"
        ignore_file.write_text("*.tmp\n")
        
        # Create test file
        (workspace / "data.tmp").write_text("Temporary\n")
        
        # Run add command with --force
        result = subprocess.run(
            ["chemvcs", "add", "--force", "data.tmp"],
            cwd=workspace,
            capture_output=True,
            text=True,
        )
        
        assert result.returncode == 0
        assert "data.tmp" in result.stdout
        assert "staged for commit" in result.stdout

    def test_add_update_existing(self, initialized_repo):
        """Test updating an already staged file."""
        workspace = initialized_repo
        
        # Create and add file
        test_file = workspace / "INCAR"
        test_file.write_text("ENCUT = 400\n")
        
        subprocess.run(
            ["chemvcs", "add", "INCAR"],
            cwd=workspace,
            capture_output=True,
        )
        
        # Modify file
        test_file.write_text("ENCUT = 500\n")
        
        # Add again
        result = subprocess.run(
            ["chemvcs", "add", "INCAR"],
            cwd=workspace,
            capture_output=True,
            text=True,
        )
        
        assert result.returncode == 0
        assert "Updated" in result.stdout
        assert "INCAR" in result.stdout

    def test_add_nonexistent_file(self, initialized_repo):
        """Test adding a file that doesn't exist."""
        workspace = initialized_repo
        
        result = subprocess.run(
            ["chemvcs", "add", "nonexistent.txt"],
            cwd=workspace,
            capture_output=True,
            text=True,
        )
        
        assert result.returncode == 1
        assert "Errors" in result.stdout
        assert "nonexistent.txt" in result.stdout

    def test_add_vasp_file_types(self, initialized_repo):
        """Test that VASP file types are correctly detected."""
        workspace = initialized_repo
        
        # Create various VASP files
        (workspace / "INCAR").write_text("ENCUT = 400\n")
        (workspace / "POSCAR").write_text("Structure\n")
        (workspace / "KPOINTS").write_text("Gamma\n")
        (workspace / "OUTCAR").write_text("Output\n")
        
        # Add files
        result = subprocess.run(
            ["chemvcs", "add", "INCAR", "POSCAR", "KPOINTS", "OUTCAR"],
            cwd=workspace,
            capture_output=True,
            text=True,
        )
        
        assert result.returncode == 0
        assert "INCAR" in result.stdout
        assert "POSCAR" in result.stdout
        assert "KPOINTS" in result.stdout
        assert "OUTCAR" in result.stdout

    def test_add_displays_file_sizes(self, initialized_repo):
        """Test that file sizes are displayed correctly."""
        workspace = initialized_repo
        
        # Create file with known size
        test_file = workspace / "INCAR"
        content = "ENCUT = 400\n" * 100  # ~1.2 KB
        test_file.write_text(content)
        
        result = subprocess.run(
            ["chemvcs", "add", "INCAR"],
            cwd=workspace,
            capture_output=True,
            text=True,
        )
        
        assert result.returncode == 0
        # Should show file size in output
        assert "KB" in result.stdout or "B" in result.stdout
