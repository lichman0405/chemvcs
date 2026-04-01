"""Main CLI entry point for ChemVCS."""

import os
import socket
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import typer
from rich.console import Console
from rich.panel import Panel

from chemvcs.constants import CHEMVCS_DIR, COMMITS_DIR, OBJECTS_DIR
from chemvcs.core import StagingManager
from chemvcs.parsers import DiffEngine
from chemvcs.plugins.manager import PluginManager
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

# ---- VASP large output files (consider storing separately) ----
# WAVECAR
# CHG
# CHGCAR

# ---- LAMMPS large output files --------------------------------
# Trajectory dump files can be many GB – ignore by default.
# Use `chemvcs add --force dump.myrun` to track a specific file.
dump.*
*.dump
# Restart files are binary and large
restart.*
*.restart
# Binary files
*.bin

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
    no_validate: bool = typer.Option(
        False,
        "--no-validate",
        help="Skip file validation checks",
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
        
        # Run validators before staging (unless disabled)
        if not no_validate:
            console.print("[bold]Running validators...[/bold]\n")
            plugin_manager = PluginManager()
            plugin_manager.discover_plugins()
            plugin_manager.load_config(chemvcs_dir)
            validation_results = plugin_manager.run_validators(
                workspace_root, 
                [str(p) for p in path_objects],
                skip_disabled=True,
            )
            
            # Check if any validation failed
            has_errors = any(not result.passed for result in validation_results)
            
            if has_errors:
                console.print("\n[bold red]Staging aborted due to validation errors[/bold red]")
                console.print("  Use [bold]--no-validate[/bold] to skip validation\n")
                raise typer.Exit(1)
            
            if validation_results:
                console.print()  # Empty line after validation
        
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
            # Determine which files were tracked in the previous commit
            head_file_early = chemvcs_dir / "HEAD"
            early_parent = None
            if head_file_early.exists():
                early_content = head_file_early.read_text(encoding="utf-8").strip()
                if early_content:
                    early_parent = early_content

            if early_parent is None:
                console.print(
                    "[dim]No previous commit — --all has no tracked files to stage[/dim]"
                )
            else:
                _commit_builder_early = CommitBuilder(chemvcs_dir, object_store, metadata_db)
                parent_commit_data = _commit_builder_early.read_commit(early_parent)
                tracked_paths = [f["path"] for f in parent_commit_data.get("files", [])]
                auto_staged = 0
                auto_skipped = 0
                for tracked in tracked_paths:
                    disk_path = workspace_root / tracked
                    if disk_path.exists() and disk_path.is_file():
                        staging_manager.add([disk_path])
                        auto_staged += 1
                    else:
                        console.print(f"  [dim]skipped {tracked} (not found on disk)[/dim]")
                        auto_skipped += 1
                console.print(
                    f"  [dim]Auto-staged {auto_staged} tracked file(s)"
                    + (f", skipped {auto_skipped}" if auto_skipped else "")
                    + "[/dim]"
                )
        
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
        diff_engine = DiffEngine()
        
        # Prepare file list for commit and compute semantic diffs
        files_for_commit: List[Dict[str, Any]] = []
        semantic_changes: List[Dict[str, Any]] = []
        
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
            
            # Try to compute semantic diff if file type is supported
            if diff_engine.can_parse(rel_path) and parent_id:
                try:
                    # Get parent commit's version of this file
                    parent_commit = commit_builder.read_commit(parent_id)
                    if parent_commit:
                        parent_files = parent_commit.get("files", [])
                        parent_file = next(
                            (f for f in parent_files if f["path"] == rel_path),
                            None
                        )
                        
                        if parent_file:
                            # File exists in parent - compute diff
                            old_content = object_store.read_blob(parent_file["blob_hash"]).decode("utf-8")
                            new_content = content.decode("utf-8")
                            diff_entries = diff_engine.diff_files(old_content, new_content, rel_path)
                            
                            if diff_entries:
                                semantic_changes.append({
                                    "file": rel_path,
                                    "entries": diff_entries,
                                })
                except Exception:
                    # Silently ignore semantic diff errors
                    pass
            
            files_for_commit.append({
                "path": rel_path,
                "content": content,
                "file_type": file_type,
                "size_bytes": size_bytes,
                "is_reference": False,
            })
        
        # Build semantic_data dict from collected changes to persist in commit
        semantic_data: Optional[Dict[str, Any]] = None
        if semantic_changes:
            semantic_data = {
                "changes": [
                    {
                        "file": change["file"],
                        "diffs": [e.to_dict() for e in change["entries"]],
                    }
                    for change in semantic_changes
                ]
            }

        # Create commit
        commit_hash = commit_builder.create_commit(
            files=files_for_commit,
            message=message,
            author=author,
            parent_hash=parent_id,
            semantic_data=semantic_data,
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
        
        # Display semantic diff summary if available
        if semantic_changes:
            console.print("\n[bold]Semantic Changes:[/bold]")
            for change in semantic_changes:
                from chemvcs.parsers.base_parser import DiffEntry
                file_name = change["file"]
                entries: List[DiffEntry] = list(change["entries"])
                summary = diff_engine.summarize_diff(entries)
                
                # Show file and summary
                console.print(f"\n  [cyan]{file_name}[/cyan]:")
                
                # Show counts by significance
                sig = summary["by_significance"]
                if sig.get("critical", 0) > 0:
                    console.print(f"    ‼️  {sig['critical']} critical change(s)")
                if sig.get("major", 0) > 0:
                    console.print(f"    ⚠️  {sig['major']} major change(s)")
                if sig.get("minor", 0) > 0:
                    console.print(f"    ℹ️  {sig['minor']} minor change(s)")
                
                # Show top changes (limit to 3 per file)
                console.print("    [dim]Key changes:[/dim]")
                for entry in entries[:3]:
                    from chemvcs.parsers.base_parser import DiffEntry
                    if not isinstance(entry, DiffEntry):
                        continue
                    if entry.change_type == "added":
                        console.print(f"      [green]+[/green] {entry.path} = {entry.new_value}")
                    elif entry.change_type == "deleted":
                        console.print(f"      [red]-[/red] {entry.path} = {entry.old_value}")
                    else:
                        console.print(f"      [yellow]~[/yellow] {entry.path}: {entry.old_value} → {entry.new_value}")
                
                if len(entries) > 3:
                    console.print(f"    [dim]... and {len(entries) - 3} more change(s)[/dim]")
        
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
        # Initialize database
        metadata_db = MetadataDB(chemvcs_dir)
        metadata_db.open()
        
        try:
            # Get commit history
            commits = metadata_db.get_commit_history(limit=max_count)
            
            if not commits:
                console.print("[dim]No commits yet[/dim]")
                return
            
            # Determine output format
            use_oneline = oneline or format == "oneline"
            
            if format == "json":
                # JSON format output
                import json
                console.print(json.dumps(commits, indent=2))
            elif use_oneline:
                # One-line format: <short_hash> <message>
                for commit in commits:
                    commit_hash = commit["commit_hash"]
                    message = commit["message"].split("\n")[0]  # First line only
                    console.print(f"[yellow]{commit_hash[:7]}[/yellow] {message}")
            else:
                # Default format: detailed commit info
                for i, commit in enumerate(commits):
                    commit_hash = commit["commit_hash"]
                    author = commit["author"]
                    timestamp = commit["timestamp"]
                    message = commit["message"]
                    parent_hash = commit.get("parent_hash")
                    
                    # Format timestamp
                    from datetime import datetime
                    dt = datetime.fromisoformat(timestamp.replace("Z", "+00:00"))
                    date_str = dt.strftime("%Y-%m-%d %H:%M:%S")
                    
                    # Print commit header
                    console.print(f"[bold yellow]commit {commit_hash}[/bold yellow]")
                    
                    if parent_hash:
                        console.print(f"[dim]Parent: {parent_hash[:7]}[/dim]")
                    else:
                        console.print("[dim]Parent: (root commit)[/dim]")
                    
                    console.print(f"[bold]Author:[/bold] {author}")
                    console.print(f"[bold]Date:[/bold]   {date_str}")
                    console.print()
                    
                    # Print commit message (indented)
                    for line in message.split("\n"):
                        console.print(f"    {line}")
                    
                    # Add separator between commits (except after last one)
                    if i < len(commits) - 1:
                        console.print()
        finally:
            metadata_db.close()
    
    except Exception as e:
        console.print(
            f"[bold red]Error:[/bold red] {e}",
            style="red",
        )
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
    file: Optional[str] = typer.Option(
        None,
        "--file",
        help="Show diff for specific file only",
    ),
) -> None:
    """Compare two revisions or working directory."""
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
        metadata_db = MetadataDB(chemvcs_dir)
        metadata_db.open()
        
        try:
            commit_builder = CommitBuilder(chemvcs_dir, object_store, metadata_db)
            diff_engine = DiffEngine()
            
            # Determine source revision (rev1)
            if rev1 is None:
                # Default to HEAD
                head_file = chemvcs_dir / "HEAD"
                if not head_file.exists() or not head_file.read_text(encoding="utf-8").strip():
                    console.print(
                        "[bold yellow]Warning:[/bold yellow] No commits yet",
                        style="yellow",
                    )
                    raise typer.Exit(1)
                rev1 = head_file.read_text(encoding="utf-8").strip()
            
            # Read source commit
            commit1 = commit_builder.read_commit(rev1)
            if not commit1:
                console.print(
                    f"[bold red]Error:[/bold red] Commit not found: {rev1}",
                    style="red",
                )
                raise typer.Exit(1)
            
            # Determine target revision (rev2)
            compare_working_tree = rev2 is None
            if rev2 is None:
                commit2 = None
                files2 = {f["path"]: f for f in commit1.get("files", [])}
                files1 = {}

                # Compare HEAD-tracked files against the current working tree.
                for tracked_path in files2:
                    working_path = workspace_root / tracked_path
                    if not working_path.exists() or not working_path.is_file():
                        continue

                    files1[tracked_path] = {
                        "path": tracked_path,
                        "content": working_path.read_bytes(),
                    }

                # Allow explicit comparison for a working-tree file that is not
                # present in HEAD (reported as added).
                if file:
                    requested_path = workspace_root / file
                    if requested_path.exists() and requested_path.is_file() and file not in files1:
                        files1[file] = {
                            "path": file,
                            "content": requested_path.read_bytes(),
                        }
            else:
                # Compare two specific commits
                commit2 = commit_builder.read_commit(rev2)
                if not commit2:
                    console.print(
                        f"[bold red]Error:[/bold red] Commit not found: {rev2}",
                        style="red",
                    )
                    raise typer.Exit(1)
                
                files1 = {f["path"]: f for f in commit1.get("files", [])}
                files2 = {f["path"]: f for f in commit2.get("files", [])}
            
            # Get all file paths
            all_paths = set(files1.keys()) | set(files2.keys())
            if file:
                # Filter to specific file
                working_path = workspace_root / file
                if compare_working_tree and working_path.exists() and working_path.is_file():
                    all_paths = {file}
                elif file not in all_paths:
                    console.print(f"[yellow]File not found in either revision: {file}[/yellow]")
                    raise typer.Exit(1)
                else:
                    all_paths = {file}
            elif compare_working_tree:
                # Default working-tree diff only reports tracked files.
                all_paths = set(files2.keys())
            
            # Display header
            commit1_hash = commit1.get("hash", rev1)
            if compare_working_tree:
                console.print(f"[bold]Comparing:[/bold] {commit1_hash[:7]} → working tree\n")
            else:
                assert commit2 is not None
                commit2_hash = commit2.get("hash", rev2 or "")
                console.print(f"[bold]Comparing:[/bold] {commit2_hash[:7]} → {commit1_hash[:7]}\n")
            
            # Process each file
            has_changes = False
            for path in sorted(all_paths):
                file1_info = files1.get(path)
                file2_info = files2.get(path)

                diff_entries = None

                # Determine change type
                if file1_info and not file2_info:
                    change_type = "added"
                    has_changes = True
                elif file2_info and not file1_info:
                    change_type = "deleted"
                    has_changes = True
                else:
                    # Both versions exist — check for modifications.
                    # First read content bytes for both sides.
                    assert file2_info is not None
                    old_content_bytes = object_store.read_blob(file2_info["blob_hash"])
                    assert file1_info is not None
                    if compare_working_tree:
                        new_content_bytes = file1_info["content"]
                    else:
                        if file1_info["blob_hash"] == file2_info["blob_hash"]:
                            continue  # identical blobs → definitely unchanged
                        new_content_bytes = object_store.read_blob(file1_info["blob_hash"])

                    # Identical bytes → unchanged regardless of file type
                    if old_content_bytes == new_content_bytes:
                        continue

                    # For files with a semantic parser, semantic equality is the
                    # real test.  A cosmetic edit (whitespace, comments) that
                    # produces no parameter changes is NOT a modification in
                    # ChemVCS terms.
                    if diff_engine.can_parse(path):
                        try:
                            diff_entries = diff_engine.diff_files(
                                old_content_bytes.decode("utf-8"),
                                new_content_bytes.decode("utf-8"),
                                path,
                            )
                            if not diff_entries:
                                # Bytes differ but semantics are unchanged
                                continue
                        except Exception:
                            diff_entries = None  # parser error → fall through

                    change_type = "modified"
                    has_changes = True

                # Display file header
                if change_type == "added":
                    console.print(f"[bold green]+ {path}[/bold green] (added)")
                elif change_type == "deleted":
                    console.print(f"[bold red]- {path}[/bold red] (deleted)")
                else:
                    console.print(f"[bold yellow]~ {path}[/bold yellow] (modified)")

                # Show semantic diff detail for modified files
                if change_type == "modified" and diff_entries is not None:
                    if summary:
                        diff_summary = diff_engine.summarize_diff(diff_entries)
                        console.print(f"  Total: {diff_summary['total_changes']} change(s)")
                        sig = diff_summary["by_significance"]
                        if sig.get("critical", 0) > 0:
                            console.print(f"  ‼️  {sig['critical']} critical")
                        if sig.get("major", 0) > 0:
                            console.print(f"  ⚠️  {sig['major']} major")
                        if sig.get("minor", 0) > 0:
                            console.print(f"  ℹ️  {sig['minor']} minor")
                    else:
                        if format == "json":
                            import json
                            entries_dict = [e.to_dict() for e in diff_entries]
                            console.print(json.dumps(entries_dict, indent=2))
                        else:
                            style = "compact" if format == "unified" else "default"
                            formatted = diff_engine.format_diff(diff_entries, style=style)
                            for line in formatted.split("\n"):
                                console.print(f"  {line}")
                elif change_type == "modified":
                    # Non-parseable file (POSCAR, POTCAR, README, …) — contents differ
                    console.print("  [dim](binary or unrecognised format — contents differ)[/dim]")

                console.print()
            
            if not has_changes:
                console.print("[dim]No changes between revisions[/dim]")
        
        finally:
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
        metadata_db = MetadataDB(chemvcs_dir)
        metadata_db.open()
        
        try:
            commit_builder = CommitBuilder(
                chemvcs_dir=chemvcs_dir,
                object_store=object_store,
                metadata_db=metadata_db,
            )
            
            # Load commit
            commit_data = commit_builder.read_commit(revision)
            
            if commit_data is None:
                console.print(
                    f"[bold red]Error:[/bold red] Commit not found: {revision}",
                    style="red",
                )
                raise typer.Exit(1)
            
            # Determine output directory
            if output_dir:
                output_path = Path(output_dir)
            else:
                short_hash = revision[:7]
                output_path = workspace_root / f"reproduce_{short_hash}"
            
            # Create output directory
            output_path.mkdir(parents=True, exist_ok=True)
            
            console.print(f"[bold]Reproducing commit:[/bold] {revision[:7]}")
            console.print(f"[bold]Output directory:[/bold] {output_path}\n")
            
            # Extract files from commit
            files = commit_data.get("files", [])
            
            if not files:
                console.print("[yellow]No files in this commit[/yellow]")
                return
            
            console.print(f"Extracting {len(files)} file(s)...\n")
            
            # Restore each file
            for file_info in files:
                rel_path = file_info["path"]
                blob_hash = file_info["blob_hash"]
                
                try:
                    # Read blob content
                    blob_content = object_store.read_blob(blob_hash)
                    
                    # Write to output directory
                    output_file = output_path / rel_path
                    output_file.parent.mkdir(parents=True, exist_ok=True)
                    output_file.write_bytes(blob_content)
                    
                    # Get file size
                    size_bytes = len(blob_content)
                    if size_bytes < 1024:
                        size_str = f"{size_bytes} B"
                    elif size_bytes < 1024 * 1024:
                        size_str = f"{size_bytes / 1024:.1f} KB"
                    else:
                        size_str = f"{size_bytes / (1024 * 1024):.2f} MB"
                    
                    console.print(
                        f"  [green]✓[/green] {rel_path}  "
                        f"[dim]({size_str}, {blob_hash[:8]})[/dim]"
                    )
                    
                except Exception as e:
                    console.print(
                        f"  [red]✗[/red] {rel_path}  [red]Error: {e}[/red]"
                    )
            
            console.print(f"\n[bold green]✓[/bold green] Reproduced successfully to {output_path}")

            # Verification: --verify-potcar checks that restored file hashes match commit records
            if verify_potcar:
                import hashlib
                console.print("\n[bold]Verifying file integrity...[/bold]")
                all_ok = True
                for file_info in files:
                    rel_path = file_info["path"]
                    expected_hash = file_info["blob_hash"]
                    output_file = output_path / rel_path
                    if not output_file.exists():
                        console.print(f"  [red]✗[/red] {rel_path}  [red](file missing)[/red]")
                        all_ok = False
                        continue
                    actual_hash = hashlib.sha256(output_file.read_bytes()).hexdigest()
                    if actual_hash == expected_hash:
                        console.print(
                            f"  [green]✓[/green] {rel_path}  [dim]hash ok[/dim]"
                        )
                    else:
                        console.print(
                            f"  [red]✗[/red] {rel_path}  "
                            f"[red]hash mismatch (expected {expected_hash[:8]}, got {actual_hash[:8]})[/red]"
                        )
                        all_ok = False
                if not all_ok:
                    raise typer.Exit(1)

            # Verification: --verify-env compares stored environment against current system
            if verify_env:
                import sys as _sys
                stored_env = commit_data.get("environment", {})
                if not stored_env:
                    console.print(
                        "\n[dim]Note: no environment info recorded in this commit "
                        "(run with environment data to enable --verify-env)[/dim]"
                    )
                else:
                    console.print("\n[bold]Comparing environment...[/bold]")
                    stored_host = stored_env.get("hostname", "")
                    current_host = socket.gethostname()
                    if stored_host and stored_host != current_host:
                        console.print(
                            f"  [yellow]⚠[/yellow]  hostname: commit=[dim]{stored_host}[/dim]"
                            f", current=[dim]{current_host}[/dim]"
                        )
                    else:
                        console.print(f"  [green]✓[/green] hostname: {current_host}")

                    stored_py = stored_env.get("python_version", "")
                    current_py = f"{_sys.version_info.major}.{_sys.version_info.minor}.{_sys.version_info.micro}"
                    if stored_py and stored_py != current_py:
                        console.print(
                            f"  [yellow]⚠[/yellow]  python: commit=[dim]{stored_py}[/dim]"
                            f", current=[dim]{current_py}[/dim]"
                        )
                    else:
                        console.print(
                            f"  [green]✓[/green] python: {current_py}"
                        )
        
        finally:
            metadata_db.close()
    
    except Exception as e:
        console.print(
            f"[bold red]Error:[/bold red] {e}",
            style="red",
        )
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
        
        # Get HEAD commit
        head_file = chemvcs_dir / "HEAD"
        current_commit = None
        if head_file.exists():
            head_content = head_file.read_text(encoding="utf-8").strip()
            if head_content:
                current_commit = head_content
        
        # Get staged files
        staged_files = staging_manager.get_staged_files()
        
        if short:
            # Short format output (similar to git status --short)
            if staged_files:
                for rel_path in sorted(staged_files.keys()):
                    console.print(f"A  {rel_path}")
            else:
                # No output in short mode if nothing staged
                pass
        else:
            # Long format output
            console.print(f"[bold]On branch:[/bold] main  [dim](linear history)[/dim]")
            
            if current_commit:
                console.print(f"[bold]HEAD:[/bold] {current_commit[:7]}  [dim]({current_commit})[/dim]")
            else:
                console.print("[bold]HEAD:[/bold] [dim](no commits yet)[/dim]")
            
            console.print()
            
            if staged_files:
                console.print("[bold green]Changes to be committed:[/bold green]")
                console.print("  [dim](use \"chemvcs commit -m <message>\" to commit)[/dim]\n")
                
                for rel_path, file_info in sorted(staged_files.items()):
                    file_type = file_info.get("file_type", "unknown")
                    size_bytes = file_info.get("size_bytes", 0)
                    blob_hash = file_info["blob_hash"]
                    
                    # Format size
                    if size_bytes < 1024:
                        size_str = f"{size_bytes} B"
                    elif size_bytes < 1024 * 1024:
                        size_str = f"{size_bytes / 1024:.1f} KB"
                    else:
                        size_str = f"{size_bytes / (1024 * 1024):.2f} MB"
                    
                    console.print(
                        f"  [green]+[/green] {rel_path}  "
                        f"[dim]({size_str}, {file_type}, {blob_hash[:8]})[/dim]"
                    )
                
                console.print()
            else:
                if current_commit:
                    console.print("[dim]Nothing to commit (working tree clean)[/dim]\n")
                else:
                    console.print("[yellow]No files staged for commit[/yellow]")
                    console.print("  Use [bold]chemvcs add <file>[/bold] to stage files\n")
    
    except Exception as e:
        console.print(
            f"[bold red]Error:[/bold red] {e}",
            style="red",
        )
        raise typer.Exit(1)


@app.command()
def plugin(
    action: str = typer.Argument(
        ...,
        help="Action to perform: list, info, enable, disable",
    ),
    plugin_name: Optional[str] = typer.Argument(
        None,
        help="Plugin name (required for info, enable, disable)",
    ),
) -> None:
    """Manage ChemVCS plugins."""
    workspace_root = Path.cwd()
    
    # Only check repository for enable/disable actions
    # list and info can work without a repository
    if action in ["enable", "disable"]:
        chemvcs_dir = workspace_root / CHEMVCS_DIR
        if not chemvcs_dir.exists():
            console.print(
                "[bold red]Error:[/bold red] Not a ChemVCS repository",
                style="red",
            )
            console.print(
                f"  No .chemvcs/ directory found in {workspace_root}",
                style="dim",
            )
            raise typer.Exit(1)
    
    try:
        plugin_manager = PluginManager()
        plugin_manager.discover_plugins()
        
        # Load persisted config when a repo exists
        chemvcs_dir = workspace_root / CHEMVCS_DIR
        if chemvcs_dir.exists():
            plugin_manager.load_config(chemvcs_dir)

        if action == "list":
            # List all discovered plugins
            validators = plugin_manager.validators

            if not validators:
                console.print("[yellow]No plugins found[/yellow]")
                console.print("\nInstall plugins with: [bold]pip install chemvcs-validator[/bold]")
                return

            console.print(f"[bold]Discovered {len(validators)} validator plugin(s):[/bold]\n")

            for validator in validators.values():
                is_on = plugin_manager.is_validator_enabled(validator)
                enabled = "[green]✓ enabled[/green]" if is_on else "[dim]○ disabled[/dim]"
                console.print(f"  {enabled}  [cyan]{validator.name}[/cyan] v{validator.version}")
                console.print(f"           {validator.description}")
                console.print(f"           [dim]Priority: {validator.priority}[/dim]\n")
        
        elif action == "info":
            if not plugin_name:
                console.print("[bold red]Error:[/bold red] Plugin name required for 'info' action")
                raise typer.Exit(1)
            
            # Find the plugin
            validator = None
            for v in plugin_manager.validators.values():
                if v.name == plugin_name:
                    validator = v
                    break
            
            if not validator:
                console.print(f"[bold red]Error:[/bold red] Plugin '{plugin_name}' not found")
                console.print("\nAvailable plugins:")
                for v in plugin_manager.validators.values():
                    console.print(f"  - {v.name}")
                raise typer.Exit(1)
            
            # Display detailed info
            console.print(Panel(
                f"""[bold cyan]{validator.name}[/bold cyan] v{validator.version}

{validator.description}

[bold]Priority:[/bold] {validator.priority}
[bold]Enabled by default:[/bold] {'Yes' if validator.enabled_by_default else 'No'}
[bold]Type:[/bold] Validator Plugin
""",
                title=f"Plugin: {validator.name}",
                border_style="cyan",
            ))
        
        elif action == "enable":
            if not plugin_name:
                console.print("[bold red]Error:[/bold red] Plugin name required for 'enable' action")
                raise typer.Exit(1)
            if not chemvcs_dir.exists():
                console.print("[bold red]Error:[/bold red] Not a ChemVCS repository")
                raise typer.Exit(1)
            plugin_manager.load_config(chemvcs_dir)
            if plugin_manager.set_validator_enabled(plugin_name, True):
                console.print(f"[green]✓[/green] Plugin '[cyan]{plugin_name}[/cyan]' enabled")
            else:
                console.print(f"[bold red]Error:[/bold red] Plugin '{plugin_name}' not found")
                raise typer.Exit(1)

        elif action == "disable":
            if not plugin_name:
                console.print("[bold red]Error:[/bold red] Plugin name required for 'disable' action")
                raise typer.Exit(1)
            if not chemvcs_dir.exists():
                console.print("[bold red]Error:[/bold red] Not a ChemVCS repository")
                raise typer.Exit(1)
            plugin_manager.load_config(chemvcs_dir)
            if plugin_manager.set_validator_enabled(plugin_name, False):
                console.print(f"[dim]○[/dim] Plugin '[cyan]{plugin_name}[/cyan]' disabled")
            else:
                console.print(f"[bold red]Error:[/bold red] Plugin '{plugin_name}' not found")
                raise typer.Exit(1)
        
        else:
            console.print(f"[bold red]Error:[/bold red] Unknown action '{action}'")
            console.print("\nAvailable actions: list, info, enable, disable")
            raise typer.Exit(1)
    
    except Exception as e:
        console.print(
            f"[bold red]Error:[/bold red] {e}",
            style="red",
        )
        raise typer.Exit(1)


def main() -> None:
    """Entry point for the CLI."""
    app()


if __name__ == "__main__":
    main()
