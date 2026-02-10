"""Test semantic diff integration with commit and diff commands."""

import os
import subprocess
from pathlib import Path

import pytest


def run_chemvcs(args, cwd):
    """Helper to run chemvcs command with proper encoding."""
    env = os.environ.copy()
    env["PYTHONIOENCODING"] = "utf-8"
    
    result = subprocess.run(
        ["chemvcs"] + args,
        cwd=cwd,
        capture_output=True,
        text=True,
        encoding="utf-8",
        env=env,
    )
    return result


@pytest.fixture
def initialized_repo(tmp_path):
    """Create a temporary initialized repository."""
    repo_dir = tmp_path / "test_repo"
    repo_dir.mkdir()
    
    # Initialize repo
    result = run_chemvcs(["init"], repo_dir)
    assert result.returncode == 0, f"Init failed: {result.stderr}"
    
    return repo_dir


class TestSemanticDiffIntegration:
    """Test semantic diff integration across commands."""
    
    def test_commit_shows_semantic_diff_for_incar_changes(self, initialized_repo):
        """Test that commit shows semantic diff when INCAR changes."""
        workspace = initialized_repo
        
        # Create initial INCAR
        incar_v1 = workspace / "INCAR"
        incar_v1.write_text(
            "ENCUT = 520\n"
            "PREC = Accurate\n"
            "ISMEAR = 0\n"
            "SIGMA = 0.05\n"
        )
        
        # First commit
        result = run_chemvcs(["add", "INCAR"], workspace)
        assert result.returncode == 0
        
        result = run_chemvcs(["commit", "-m", "Initial INCAR"], workspace)
        assert result.returncode == 0
        
        # Modify INCAR with critical and major changes
        incar_v2 = workspace / "INCAR"
        incar_v2.write_text(
            "ENCUT = 600\n"  # Critical change
            "PREC = Accurate\n"
            "ISMEAR = 1\n"   # Critical change
            "SIGMA = 0.1\n"  # Critical change
            "LWAVE = .FALSE.\n"  # Major change (new)
        )
        
        result = run_chemvcs(["add", "INCAR"], workspace)
        assert result.returncode == 0
        
        result = run_chemvcs(["commit", "-m", "Update INCAR parameters"], workspace)
        assert result.returncode == 0
        
        # Verify semantic diff output in commit
        assert "Semantic Changes:" in result.stdout
        assert "INCAR" in result.stdout
        assert "critical" in result.stdout.lower()
        assert "ENCUT" in result.stdout or "change" in result.stdout.lower()
    
    def test_commit_shows_semantic_diff_for_kpoints_changes(self, initialized_repo):
        """Test that commit shows semantic diff when KPOINTS changes."""
        workspace = initialized_repo
        
        # Create initial KPOINTS
        kpoints_v1 = workspace / "KPOINTS"
        kpoints_v1.write_text(
            "Automatic mesh\n"
            "0\n"
            "Gamma\n"
            "4 4 4\n"
            "0 0 0\n"
        )
        
        # First commit
        result = run_chemvcs(["add", "KPOINTS"], workspace)
        assert result.returncode == 0
        
        result = run_chemvcs(["commit", "-m", "Initial KPOINTS"], workspace)
        assert result.returncode == 0
        
        # Modify KPOINTS - change grid (critical)
        kpoints_v2 = workspace / "KPOINTS"
        kpoints_v2.write_text(
            "Automatic mesh\n"
            "0\n"
            "Gamma\n"
            "8 8 8\n"  # Critical change
            "0 0 0\n"
        )
        
        result = run_chemvcs(["add", "KPOINTS"], workspace)
        assert result.returncode == 0
        
        result = run_chemvcs(["commit", "-m", "Increase k-point density"], workspace)
        assert result.returncode == 0
        
        # Verify semantic diff output
        assert "Semantic Changes:" in result.stdout
        assert "KPOINTS" in result.stdout
        assert "critical" in result.stdout.lower()
    
    def test_diff_command_semantic_format(self, initialized_repo):
        """Test diff command with semantic format."""
        workspace = initialized_repo
        
        # Create and commit INCAR v1
        incar = workspace / "INCAR"
        incar.write_text("ENCUT = 520\nISMEAR = 0\n")
        
        result = run_chemvcs(["add", "INCAR"], workspace)
        assert result.returncode == 0
        
        result = run_chemvcs(["commit", "-m", "Initial"], workspace)
        assert result.returncode == 0
        
        # Create and commit INCAR v2
        incar.write_text("ENCUT = 600\nISMEAR = 1\n")
        
        result = run_chemvcs(["add", "INCAR"], workspace)
        assert result.returncode == 0
        
        result = run_chemvcs(["commit", "-m", "Update"], workspace)
        assert result.returncode == 0
        
        # Run diff command
        result = run_chemvcs(["diff"], workspace)
        
        # Should show comparison between HEAD and parent
        assert result.returncode == 0
        assert "Comparing:" in result.stdout
        assert "INCAR" in result.stdout
        assert "modified" in result.stdout
    
    def test_diff_command_summary_flag(self, initialized_repo):
        """Test diff command with --summary flag."""
        workspace = initialized_repo
        
        # Create two commits with INCAR changes
        incar = workspace / "INCAR"
        incar.write_text("ENCUT = 520\n")
        
        result = run_chemvcs(["add", "INCAR"], workspace)
        assert result.returncode == 0
        result = run_chemvcs(["commit", "-m", "v1"], workspace)
        assert result.returncode == 0
        
        incar.write_text("ENCUT = 600\nPREC = Accurate\n")
        result = run_chemvcs(["add", "INCAR"], workspace)
        assert result.returncode == 0
        result = run_chemvcs(["commit", "-m", "v2"], workspace)
        assert result.returncode == 0
        
        # Run diff with summary
        result = run_chemvcs(["diff", "--summary"], workspace)
        
        assert result.returncode == 0
        assert "Total:" in result.stdout or "change" in result.stdout.lower()
    
    def test_diff_command_json_format(self, initialized_repo):
        """Test diff command with JSON format."""
        workspace = initialized_repo
        
        # Create two commits
        incar = workspace / "INCAR"
        incar.write_text("ENCUT = 520\n")
        
        result = run_chemvcs(["add", "INCAR"], workspace)
        assert result.returncode == 0
        result = run_chemvcs(["commit", "-m", "v1"], workspace)
        assert result.returncode == 0
        
        incar.write_text("ENCUT = 600\n")
        result = run_chemvcs(["add", "INCAR"], workspace)
        assert result.returncode == 0
        result = run_chemvcs(["commit", "-m", "v2"], workspace)
        assert result.returncode == 0
        
        # Run diff with JSON format
        result = run_chemvcs(["diff", "--format", "json"], workspace)
        
        assert result.returncode == 0
        # Should contain JSON-like output
        assert "{" in result.stdout and "}" in result.stdout
    
    def test_diff_command_file_filter(self, initialized_repo):
        """Test diff command with --file filter."""
        workspace = initialized_repo
        
        # Create multiple files
        incar = workspace / "INCAR"
        incar.write_text("ENCUT = 520\n")
        kpoints = workspace / "KPOINTS"
        kpoints.write_text("Automatic\n0\nGamma\n4 4 4\n0 0 0\n")
        
        result = run_chemvcs(["add", "INCAR", "KPOINTS"], workspace)
        assert result.returncode == 0
        result = run_chemvcs(["commit", "-m", "v1"], workspace)
        assert result.returncode == 0
        
        # Modify both
        incar.write_text("ENCUT = 600\n")
        kpoints.write_text("Automatic\n0\nGamma\n8 8 8\n0 0 0\n")
        
        result = run_chemvcs(["add", "INCAR", "KPOINTS"], workspace)
        assert result.returncode == 0
        result = run_chemvcs(["commit", "-m", "v2"], workspace)
        assert result.returncode == 0
        
        # Diff only INCAR
        result = run_chemvcs(["diff", "--file", "INCAR"], workspace)
        
        assert result.returncode == 0
        assert "INCAR" in result.stdout
        # Should not show KPOINTS
        assert "KPOINTS" not in result.stdout or result.stdout.count("KPOINTS") == 0
    
    def test_commit_no_semantic_diff_for_first_commit(self, initialized_repo):
        """Test that first commit doesn't show semantic diff (no parent)."""
        workspace = initialized_repo
        
        # Create INCAR
        incar = workspace / "INCAR"
        incar.write_text("ENCUT = 520\n")
        
        result = run_chemvcs(["add", "INCAR"], workspace)
        assert result.returncode == 0
        
        result = run_chemvcs(["commit", "-m", "Initial commit"], workspace)
        assert result.returncode == 0
        
        # Should not show semantic diff for first commit
        assert "Semantic Changes:" not in result.stdout
    
    def test_commit_semantic_diff_skips_non_vasp_files(self, initialized_repo):
        """Test that semantic diff only processes VASP files."""
        workspace = initialized_repo
        
        # Create INCAR and a regular file
        incar = workspace / "INCAR"
        incar.write_text("ENCUT = 520\n")
        readme = workspace / "README.md"
        readme.write_text("# Project\n")
        
        result = run_chemvcs(["add", "INCAR", "README.md"], workspace)
        assert result.returncode == 0
        result = run_chemvcs(["commit", "-m", "v1"], workspace)
        assert result.returncode == 0
        
        # Modify both
        incar.write_text("ENCUT = 600\n")
        readme.write_text("# Project\nUpdated\n")
        
        result = run_chemvcs(["add", "INCAR", "README.md"], workspace)
        assert result.returncode == 0
        result = run_chemvcs(["commit", "-m", "v2"], workspace)
        assert result.returncode == 0
        
        # Should show semantic diff only for INCAR
        assert "Semantic Changes:" in result.stdout
        assert "INCAR" in result.stdout
        # README.md should not appear in semantic changes section
        # (it may appear in file list, but not in semantic analysis)
