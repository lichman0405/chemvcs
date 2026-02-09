"""Unit tests for chemvcs init command."""

import os
from pathlib import Path

import pytest
from typer.testing import CliRunner

from chemvcs.cli.main import app
from chemvcs.constants import CHEMVCS_DIR

runner = CliRunner()


class TestInitCommand:
    """Test chemvcs init command."""

    def test_init_creates_directory_structure(self, tmp_path: Path) -> None:
        """Test that init creates required directories."""
        # Change to tmp_path directory
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        
        try:
            result = runner.invoke(app, ["init", "--quiet"])
            
            assert result.exit_code == 0
            
            chemvcs_dir = tmp_path / CHEMVCS_DIR
            assert chemvcs_dir.exists()
            assert (chemvcs_dir / "objects").exists()
            assert (chemvcs_dir / "commits").exists()
            assert (chemvcs_dir / "metadata.db").exists()
        finally:
            os.chdir(original_cwd)

    def test_init_creates_chemvcsignore(self, tmp_path: Path) -> None:
        """Test that init creates default .chemvcsignore file."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        
        try:
            result = runner.invoke(app, ["init", "--quiet"])
            
            assert result.exit_code == 0
            
            ignore_file = tmp_path / ".chemvcsignore"
            assert ignore_file.exists()
            
            content = ignore_file.read_text(encoding="utf-8")
            assert "*.tmp" in content
            assert "__pycache__/" in content
        finally:
            os.chdir(original_cwd)

    def test_init_does_not_overwrite_existing_ignore(self, tmp_path: Path) -> None:
        """Test that init preserves existing .chemvcsignore."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        
        try:
            ignore_file = tmp_path / ".chemvcsignore"
            custom_content = "# My custom ignore rules\n*.custom\n"
            ignore_file.write_text(custom_content, encoding="utf-8")
            
            result = runner.invoke(app, ["init", "--quiet"])
            
            assert result.exit_code == 0
            assert ignore_file.read_text(encoding="utf-8") == custom_content
        finally:
            os.chdir(original_cwd)

    def test_init_fails_when_already_initialized(self, tmp_path: Path) -> None:
        """Test that init fails on already initialized directory."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        
        try:
            # First init
            result1 = runner.invoke(app, ["init", "--quiet"])
            assert result1.exit_code == 0
            
            # Second init without --force
            result2 = runner.invoke(app, ["init"])
            assert result2.exit_code == 1
            assert "already exists" in result2.stdout
        finally:
            os.chdir(original_cwd)

    def test_init_force_reinitializes(self, tmp_path: Path) -> None:
        """Test that --force flag reinitializes repository."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        
        try:
            # First init
            result1 = runner.invoke(app, ["init", "--quiet"])
            assert result1.exit_code == 0
            
            # Create a test file in .chemvcs to verify deletion
            test_file = tmp_path / CHEMVCS_DIR / "test_marker.txt"
            test_file.write_text("marker")
            assert test_file.exists()
            
            # Force reinit
            result2 = runner.invoke(app, ["init", "--force", "--quiet"])
            assert result2.exit_code == 0
            
            # Verify test file was removed
            assert not test_file.exists()
            
            # Verify structure was recreated
            chemvcs_dir = tmp_path / CHEMVCS_DIR
            assert chemvcs_dir.exists()
            assert (chemvcs_dir / "metadata.db").exists()
        finally:
            os.chdir(original_cwd)

    def test_init_quiet_mode_no_output(self, tmp_path: Path) -> None:
        """Test that --quiet suppresses output."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        
        try:
            result = runner.invoke(app, ["init", "--quiet"])
            
            assert result.exit_code == 0
            # In quiet mode, stdout should be minimal or empty
            assert "ChemVCS Initialized" not in result.stdout
        finally:
            os.chdir(original_cwd)

    def test_init_normal_mode_shows_panel(self, tmp_path: Path) -> None:
        """Test that normal mode shows success panel."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        
        try:
            result = runner.invoke(app, ["init"])
            
            assert result.exit_code == 0
            assert "ChemVCS Initialized" in result.stdout
            assert "Repository root:" in result.stdout
            assert "Next steps:" in result.stdout
        finally:
            os.chdir(original_cwd)

    def test_init_database_schema_initialized(self, tmp_path: Path) -> None:
        """Test that database schema is correctly initialized."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)
        
        try:
            result = runner.invoke(app, ["init", "--quiet"])
            assert result.exit_code == 0
            
            # Verify database exists and has correct schema
            from chemvcs.storage import MetadataDB
            
            db = MetadataDB(tmp_path / CHEMVCS_DIR)
            db.open()
            
            # Check schema version
            version = db.get_schema_version()
            assert version > 0
            
            # Check tables exist
            cursor = db.conn.cursor()  # type: ignore
            cursor.execute(
                "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name"
            )
            tables = [row[0] for row in cursor.fetchall()]
            
            expected_tables = ["commits", "environment", "files", "metadata", "output_summary", "semantic_summary"]
            for table in expected_tables:
                assert table in tables
            
            db.close()
        finally:
            os.chdir(original_cwd)
