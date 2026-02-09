"""Fixtures for integration tests."""

import subprocess
from pathlib import Path

import pytest


@pytest.fixture
def initialized_repo(tmp_path):
    """Create a temporary directory with initialized ChemVCS repository.
    
    Returns:
        Path: Path to the workspace root
    """
    workspace = tmp_path / "test_workspace"
    workspace.mkdir()
    
    # Initialize repository
    result = subprocess.run(
        ["chemvcs", "init", "--quiet"],
        cwd=workspace,
        capture_output=True,
        text=True,
    )
    
    if result.returncode != 0:
        raise RuntimeError(f"Failed to initialize repo: {result.stdout}\n{result.stderr}")
    
    return workspace
