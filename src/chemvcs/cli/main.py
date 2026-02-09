"""Main CLI entry point for ChemVCS."""

import os
import socket
import sys
from pathlib import Path
from typing import Optional

import typer
from rich.console import Console
from rich.panel import Panel

from chemvcs.constants import CHEMVCS_DIR, COMMITS_DIR, OBJECTS_DIR
from chemvcs.storage import MetadataDB

console = Console()
app = typer.Typer(
    name="chemvcs",
    help="Version control for computational materials science calculations",
    add_completion=False,
)


@app.command()
def version() -> None:
    """Show ChemVCS version."""
    from chemvcs import __version__
    typer.echo(f"ChemVCS version {__version__}")


@app.command()
def init(
    force: bool = typer.Option(
        False,
        "--force",
        help="Overwrite existing .chemvcs/ directory (dangerous!)",
    ),
    quiet: bool = typer.Option(
        False,
        "--quiet",
        "-q",
        help="Suppress output except errors",
    ),
) -> None:
    """Initialize a ChemVCS repository in the current directory."""
    workspace_root = Path.cwd()
    chemvcs_dir = workspace_root / CHEMVCS_DIR
    
    # Check if already initialized
    if chemvcs_dir.exists():
        if not force:
            console.print(
                f"[bold red]Error:[/bold red] ChemVCS repository already exists in {workspace_root}",
                style="red",
            )
            console.print(
                f"  .chemvcs/ directory found at: {chemvcs_dir}",
                style="dim",
            )
            console.print(
                "\nUse [bold]--force[/bold] to reinitialize (will delete existing data!)",
                style="yellow",
            )
            raise typer.Exit(1)
        
        # Force mode: remove existing directory
        if not quiet:
            console.print(
                f"[yellow]⚠️  Removing existing .chemvcs/ directory...[/yellow]"
            )
        import shutil
        shutil.rmtree(chemvcs_dir)
    
    try:
        # Create directory structure
        chemvcs_dir.mkdir(exist_ok=True)
        (chemvcs_dir / OBJECTS_DIR).mkdir(exist_ok=True)
        (chemvcs_dir / COMMITS_DIR).mkdir(exist_ok=True)
        
        # Initialize metadata database
        db = MetadataDB(chemvcs_dir)
        db.open()
        db.init_schema()
        db.close()
        
        # Create default .chemvcsignore
        chemvcsignore_path = workspace_root / ".chemvcsignore"
        if not chemvcsignore_path.exists():
            default_ignore_content = """# ChemVCS Ignore Rules
# 
# Patterns follow .gitignore syntax
# Lines starting with # are comments

# Temporary files
*.tmp
*.swp
*~
.DS_Store
Thumbs.db

# Large output files (consider storing separately)
# WAVECAR
# CHG
# CHGCAR

# Compiled Python
__pycache__/
*.pyc
*.pyo

# Virtual environments
venv/
.venv/
env/

# IDE files
.vscode/
.idea/
*.sublime-*

# OS files
.Trash-*/
"""
            chemvcsignore_path.write_text(default_ignore_content, encoding="utf-8")
        
        if not quiet:
            # Success message with fancy panel
            hostname = socket.gethostname()
            username = os.getenv("USER") or os.getenv("USERNAME") or "unknown"
            author_id = f"{username}@{hostname}"
            
            success_message = f"""[bold green]✓[/bold green] Initialized ChemVCS repository

[dim]Repository root:[/dim] {workspace_root}
[dim]Storage location:[/dim] {chemvcs_dir}
[dim]Default author:[/dim] {author_id}

[bold]Next steps:[/bold]
  1. Add VASP input files: [cyan]chemvcs add INCAR POSCAR KPOINTS POTCAR[/cyan]
  2. Create initial commit: [cyan]chemvcs commit -m "Initial setup"[/cyan]
  3. Run your calculation
  4. Track outputs: [cyan]chemvcs add OUTCAR vasprun.xml[/cyan]
  5. Commit results: [cyan]chemvcs commit -m "Completed SCF calculation"[/cyan]
"""
            console.print(Panel(success_message, border_style="green", title="ChemVCS Initialized"))
        
    except Exception as e:
        console.print(
            f"[bold red]Error:[/bold red] Failed to initialize repository: {e}",
            style="red",
        )
        # Clean up partial initialization
        if chemvcs_dir.exists():
            import shutil
            shutil.rmtree(chemvcs_dir)
        raise typer.Exit(1)


@app.command()
def add(
    paths: list[str] = typer.Argument(..., help="Files or directories to add"),
    force: bool = typer.Option(
        False,
        "--force",
        help="Override .chemvcsignore rules",
    ),
    dry_run: bool = typer.Option(
        False,
        "--dry-run",
        help="Show what would be added without modifying staging area",
    ),
) -> None:
    """Add files to the staging area."""
    typer.echo("🚧 chemvcs add - Not implemented yet")
    typer.echo(f"Would add: {', '.join(paths)}")
    raise typer.Exit(1)


@app.command()
def commit(
    message: Optional[str] = typer.Option(
        None,
        "--message",
        "-m",
        help="Commit message (required)",
    ),
    all: bool = typer.Option(  # noqa: A002
        False,
        "--all",
        "-a",
        help="Automatically stage all tracked files",
    ),
    allow_empty: bool = typer.Option(
        False,
        "--allow-empty",
        help="Allow empty commits (for milestones)",
    ),
    author: Optional[str] = typer.Option(
        None,
        "--author",
        help="Override default author (format: name@host)",
    ),
) -> None:
    """Commit staged files as a snapshot."""
    typer.echo("🚧 chemvcs commit - Not implemented yet")
    if message:
        typer.echo(f"Message: {message}")
    raise typer.Exit(1)


@app.command()
def log(
    max_count: Optional[int] = typer.Option(
        None,
        "--max-count",
        "-n",
        help="Limit number of commits to show",
    ),
    oneline: bool = typer.Option(
        False,
        "--oneline",
        help="Show each commit on a single line",
    ),
    format: str = typer.Option(  # noqa: A002
        "default",
        "--format",
        help="Output format: default, oneline, json",
    ),
) -> None:
    """Show commit history."""
    typer.echo("🚧 chemvcs log - Not implemented yet")
    raise typer.Exit(1)


@app.command()
def diff(
    rev1: Optional[str] = typer.Argument(None, help="Source revision (default: HEAD)"),
    rev2: Optional[str] = typer.Argument(None, help="Target revision (default: working dir)"),
    format: str = typer.Option(  # noqa: A002
        "semantic",
        "--format",
        help="Output format: semantic, unified, json",
    ),
    summary: bool = typer.Option(
        False,
        "--summary",
        help="Show only change summary",
    ),
) -> None:
    """Compare two revisions or working directory."""
    typer.echo("🚧 chemvcs diff - Not implemented yet")
    raise typer.Exit(1)


@app.command()
def reproduce(
    revision: str = typer.Argument(..., help="Commit hash or tag to reproduce"),
    output_dir: Optional[str] = typer.Option(
        None,
        "--output-dir",
        "-o",
        help="Output directory (default: reproduce_<rev>/)",
    ),
    verify_potcar: bool = typer.Option(
        True,
        "--verify-potcar/--no-verify-potcar",
        help="Verify POTCAR file hashes",
    ),
    verify_env: bool = typer.Option(
        True,
        "--verify-env/--no-verify-env",
        help="Compare environment info",
    ),
) -> None:
    """Restore files from a historical commit."""
    typer.echo("🚧 chemvcs reproduce - Not implemented yet")
    typer.echo(f"Would reproduce commit: {revision}")
    raise typer.Exit(1)


@app.command()
def status(
    short: bool = typer.Option(
        False,
        "--short",
        help="Show short format output",
    ),
) -> None:
    """Show working directory and staging area status."""
    typer.echo("🚧 chemvcs status - Not implemented yet")
    raise typer.Exit(1)


def main() -> None:
    """Entry point for the CLI."""
    app()


if __name__ == "__main__":
    main()
