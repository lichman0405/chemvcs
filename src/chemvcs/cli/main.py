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
from chemvcs.core import StagingManager
from chemvcs.storage import CommitBuilder, MetadataDB, ObjectStore

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
    workspace_root = Path.cwd()
    chemvcs_dir = workspace_root / CHEMVCS_DIR
    
    # Check if repository is initialized
    if not chemvcs_dir.exists():
        console.print(
            "[bold red]Error:[/bold red] Not a ChemVCS repository",
            style="red",
        )
        console.print(
            f"  No .chemvcs/ directory found in {workspace_root}",
            style="dim",
        )
        console.print(
            "\nRun [bold]chemvcs init[/bold] to initialize a repository",
            style="yellow",
        )
        raise typer.Exit(1)
    
    try:
        # Initialize storage
        object_store = ObjectStore(chemvcs_dir)
        staging_manager = StagingManager(workspace_root, object_store)
        
        # Convert string paths to Path objects
        path_objects = [Path(p) for p in paths]
        
        # Dry run mode: just show what would be added
        if dry_run:
            console.print("[bold yellow]Dry run mode - no changes will be made[/bold yellow]\n")
        
        # Add files to staging
        console.print("[bold]Staging files...[/bold]\n")
        stats = staging_manager.add(path_objects, force=force)
        
        # Display added files
        if stats["added"]:
            console.print("[bold green]Added:[/bold green]")
            for path_str in stats["added"]:
                # Get file info from staging
                staged_files = staging_manager.get_staged_files()
                file_info = staged_files.get(path_str, {})
                file_type = file_info.get("file_type", "UNKNOWN")
                size_bytes = file_info.get("size_bytes", 0)
                
                # Format size
                if size_bytes < 1024:
                    size_str = f"{size_bytes} B"
                elif size_bytes < 1024 * 1024:
                    size_str = f"{size_bytes / 1024:.1f} KB"
                else:
                    size_str = f"{size_bytes / (1024 * 1024):.2f} MB"
                
                console.print(f"  [green]+[/green] {path_str}  [dim]({size_str}, {file_type})[/dim]")
        
        # Display updated files
        if stats["updated"]:
            console.print("\n[bold yellow]Updated:[/bold yellow]")
            for path_str in stats["updated"]:
                staged_files = staging_manager.get_staged_files()
                file_info = staged_files.get(path_str, {})
                file_type = file_info.get("file_type", "UNKNOWN")
                size_bytes = file_info.get("size_bytes", 0)
                
                if size_bytes < 1024:
                    size_str = f"{size_bytes} B"
                elif size_bytes < 1024 * 1024:
                    size_str = f"{size_bytes / 1024:.1f} KB"
                else:
                    size_str = f"{size_bytes / (1024 * 1024):.2f} MB"
                
                console.print(f"  [yellow]*[/yellow] {path_str}  [dim]({size_str}, {file_type})[/dim]")
        
        # Display ignored files
        if stats["ignored"]:
            console.print("\n[bold dim]Ignored:[/bold dim]")
            for path_str in stats["ignored"]:
                console.print(f"  [dim]-[/dim] {path_str}  [dim](.chemvcsignore)[/dim]")
        
        # Display errors
        if stats["errors"]:
            console.print("\n[bold red]Errors:[/bold red]")
            for error in stats["errors"]:
                console.print(f"  [red]x[/red] {error}")
        
        # Summary
        total_added = len(stats["added"]) + len(stats["updated"])
        if total_added > 0:
            console.print(f"\n[bold green]>[/bold green] {total_added} file(s) staged for commit")
        elif stats["ignored"]:
            console.print("\n[yellow]No files staged. All files were ignored.[/yellow]")
            console.print("  Use [bold]--force[/bold] to override .chemvcsignore rules")
        else:
            console.print("\n[yellow]No files found to add[/yellow]")
        
        if stats["errors"]:
            raise typer.Exit(1)
            
    except Exception as e:
        console.print(
            f"[bold red]Error:[/bold red] {e}",
            style="red",
        )
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
    workspace_root = Path.cwd()
    chemvcs_dir = workspace_root / CHEMVCS_DIR
    
    # Check if repository is initialized
    if not chemvcs_dir.exists():
        console.print(
            "[bold red]Error:[/bold red] Not a ChemVCS repository",
            style="red",
        )
        console.print(
            f"  No .chemvcs/ directory found in {workspace_root}",
            style="dim",
        )
        console.print(
            "\nRun [bold]chemvcs init[/bold] to initialize a repository",
            style="yellow",
        )
        raise typer.Exit(1)
    
    # Require message
    if not message:
        console.print(
            "[bold red]Error:[/bold red] Commit message is required",
            style="red",
        )
        console.print(
            "  Use [bold]-m \"your message\"[/bold] to provide a commit message",
            style="yellow",
        )
        raise typer.Exit(1)
    
    try:
        # Initialize storage
        object_store = ObjectStore(chemvcs_dir)
        staging_manager = StagingManager(workspace_root, object_store)
        metadata_db = MetadataDB(chemvcs_dir)
        metadata_db.open()
        
        # Handle -a/--all flag (auto-stage tracked files)
        if all:
            console.print("[dim]Auto-staging tracked files...[/dim]")
            # TODO: Implement tracking logic in Phase 1c later
            # For now, just note that this is not implemented
            console.print("[yellow]Warning: --all flag not yet implemented[/yellow]")
        
        # Check if staging area is empty
        staged_files = staging_manager.get_staged_files()
        if not staged_files and not allow_empty:
            console.print(
                "[bold yellow]Warning:[/bold yellow] Nothing to commit (staging area is empty)",
                style="yellow",
            )
            console.print(
                "  Use [bold]chemvcs add <files>[/bold] to stage files",
                style="dim",
            )
            console.print(
                "  Or use [bold]--allow-empty[/bold] to create an empty commit",
                style="dim",
            )
            raise typer.Exit(1)
        
        # Determine author
        if not author:
            hostname = socket.gethostname()
            username = os.getenv("USER") or os.getenv("USERNAME") or "unknown"
            author = f"{username}@{hostname}"
        
        # Get parent commit (HEAD)
        head_file = chemvcs_dir / "HEAD"
        parent_id = None
        if head_file.exists():
            parent_content = head_file.read_text(encoding="utf-8").strip()
            if parent_content:
                parent_id = parent_content
        
        # Build commit
        console.print(f"\n[bold]Committing {len(staged_files)} file(s)...[/bold]\n")
        
        commit_builder = CommitBuilder(chemvcs_dir, object_store, metadata_db)
        
        # Prepare file list for commit
        files_for_commit = []
        for rel_path, file_info in staged_files.items():
            blob_hash = file_info["blob_hash"]
            file_type = file_info.get("file_type", "unknown")
            size_bytes = file_info.get("size_bytes", 0)
            
            console.print(
                f"  [{len(files_for_commit) + 1}/{len(staged_files)}] "
                f"{rel_path}  [dim]-> {blob_hash[:8]}[/dim]"
            )
            
            # Read content from blob store
            content = object_store.read_blob(blob_hash)
            
            files_for_commit.append({
                "path": rel_path,
                "content": content,
                "file_type": file_type,
                "size_bytes": size_bytes,
                "is_reference": False,
            })
        
        # Create commit
        commit_hash = commit_builder.create_commit(
            files=files_for_commit,
            message=message,
            author=author,
            parent_hash=parent_id,
        )
        
        # Update HEAD
        head_file.write_text(commit_hash, encoding="utf-8")
        
        # Clear staging area
        staging_manager.clear()
        
        # Success output
        console.print(f"\n[bold green]>[/bold green] Committed [bold cyan]{commit_hash[:7]}[/bold cyan]")
        console.print(f"  [dim]Author:[/dim]  {author}")
        
        # Get commit timestamp from metadata
        commit_info = metadata_db.get_commit_by_hash(commit_hash)
        if commit_info:
            from datetime import datetime as dt
            timestamp = dt.fromisoformat(commit_info["timestamp"])
            console.print(f"  [dim]Date:[/dim]    {timestamp.strftime('%Y-%m-%d %H:%M:%S')}")
        
        console.print(f"  [dim]Parent:[/dim]  {parent_id[:7] if parent_id else '(root commit)'}")
        console.print(f"\n  {message}")
        
        metadata_db.close()
        
    except Exception as e:
        console.print(
            f"[bold red]Error:[/bold red] {e}",
            style="red",
        )
        if 'metadata_db' in locals():
            metadata_db.close()
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
