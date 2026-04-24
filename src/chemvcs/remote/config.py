"""Remote repository configuration stored in .chemvcs/remotes.toml."""

from __future__ import annotations

from pathlib import Path

try:
    import tomllib  # type: ignore[import-not-found]  # Python 3.11+
except ImportError:
    import tomli as tomllib

from chemvcs.constants import REMOTES_FILE


class RemoteNotFoundError(ValueError):
    """Raised when a named remote is not found in remotes.toml."""


class RemoteAlreadyExistsError(ValueError):
    """Raised when trying to add a remote that already exists."""


def _remotes_path(chemvcs_dir: Path) -> Path:
    return chemvcs_dir / REMOTES_FILE


def list_remotes(chemvcs_dir: Path) -> dict[str, str]:
    """Return a mapping of ``{name: url}`` from the remotes config file.

    Returns an empty dict if the file does not exist.
    """
    path = _remotes_path(chemvcs_dir)
    if not path.exists():
        return {}
    data = tomllib.loads(path.read_text(encoding="utf-8"))
    result: dict[str, str] = {}
    for name, entry in data.items():
        result[name] = entry.get("url", "")
    return result


def add_remote(chemvcs_dir: Path, name: str, url: str) -> None:
    """Register a new remote in the remotes config file.

    Raises:
        RemoteAlreadyExistsError: If *name* is already registered.
    """
    remotes = list_remotes(chemvcs_dir)
    if name in remotes:
        raise RemoteAlreadyExistsError(f"Remote '{name}' already exists")

    remotes[name] = url
    _write_remotes(chemvcs_dir, remotes)


def remove_remote(chemvcs_dir: Path, name: str) -> None:
    """Remove a remote from the remotes config file.

    Raises:
        RemoteNotFoundError: If *name* is not registered.
    """
    remotes = list_remotes(chemvcs_dir)
    if name not in remotes:
        raise RemoteNotFoundError(f"Remote '{name}' not found")
    del remotes[name]
    _write_remotes(chemvcs_dir, remotes)


def get_remote_url(chemvcs_dir: Path, name: str) -> str:
    """Return the URL for a named remote.

    Raises:
        RemoteNotFoundError: If *name* is not registered.
    """
    remotes = list_remotes(chemvcs_dir)
    if name not in remotes:
        raise RemoteNotFoundError(f"Remote '{name}' not found")
    return remotes[name]


def _write_remotes(chemvcs_dir: Path, remotes: dict[str, str]) -> None:
    """Persist remotes dictionary to TOML format."""
    lines: list[str] = []
    for name, url in sorted(remotes.items()):
        lines.append(f"[{name}]")
        lines.append(f'url = "{url}"')
        lines.append("")
    path = _remotes_path(chemvcs_dir)
    path.write_text("\n".join(lines), encoding="utf-8")
