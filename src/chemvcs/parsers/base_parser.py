"""Base parser interface for VASP files.

This module defines the abstract base class for all VASP file parsers.
"""

from abc import ABC, abstractmethod
from typing import Any, Dict, List, Optional, Tuple


class ParserError(Exception):
    """Exception raised during parsing."""


class DiffEntry:
    """Represents a single change in semantic diff.
    
    Attributes:
        path: Dot-separated path to the changed field (e.g., "ENCUT", "positions.0.x")
        old_value: Previous value (None if added)
        new_value: New value (None if deleted)
        change_type: Type of change: "added", "deleted", "modified"
        significance: Importance level: "critical", "major", "minor"
    """
    
    def __init__(
        self,
        path: str,
        old_value: Any,
        new_value: Any,
        change_type: str,
        significance: str = "minor",
    ):
        self.path = path
        self.old_value = old_value
        self.new_value = new_value
        self.change_type = change_type
        self.significance = significance
    
    def __repr__(self) -> str:
        if self.change_type == "added":
            return f"DiffEntry({self.path}: +{self.new_value} [{self.significance}])"
        elif self.change_type == "deleted":
            return f"DiffEntry({self.path}: -{self.old_value} [{self.significance}])"
        else:
            return f"DiffEntry({self.path}: {self.old_value} → {self.new_value} [{self.significance}])"
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary representation."""
        return {
            "path": self.path,
            "old_value": self.old_value,
            "new_value": self.new_value,
            "change_type": self.change_type,
            "significance": self.significance,
        }


class BaseParser(ABC):
    """Abstract base class for VASP file parsers.
    
    Each parser implementation must provide:
    - parse(): Convert raw text to structured data
    - diff(): Compute semantic differences between two parsed structures
    - format_diff(): Human-readable diff output
    """
    
    @abstractmethod
    def parse(self, content: str) -> Dict[str, Any]:
        """Parse file content into structured data.
        
        Args:
            content: Raw file content as string
            
        Returns:
            Dictionary containing parsed structured data
            
        Raises:
            ParserError: If content is malformed or invalid
        """
        pass
    
    @abstractmethod
    def diff(
        self,
        old_data: Dict[str, Any],
        new_data: Dict[str, Any],
    ) -> List[DiffEntry]:
        """Compute semantic differences between two parsed structures.
        
        Args:
            old_data: Parsed data from old version
            new_data: Parsed data from new version
            
        Returns:
            List of DiffEntry objects describing changes
        """
        pass
    
    def format_diff(
        self,
        diff_entries: List[DiffEntry],
        style: str = "default",
    ) -> str:
        """Format diff entries for human-readable output.
        
        Args:
            diff_entries: List of diff entries
            style: Output style ("default", "compact", "detailed")
            
        Returns:
            Formatted diff string
        """
        if not diff_entries:
            return "No changes"
        
        lines = []
        
        if style == "compact":
            for entry in diff_entries:
                if entry.change_type == "added":
                    lines.append(f"+ {entry.path} = {entry.new_value}")
                elif entry.change_type == "deleted":
                    lines.append(f"- {entry.path} = {entry.old_value}")
                else:
                    lines.append(f"~ {entry.path}: {entry.old_value} → {entry.new_value}")
        
        elif style == "detailed":
            for entry in diff_entries:
                lines.append(f"Path: {entry.path}")
                lines.append(f"  Type: {entry.change_type}")
                lines.append(f"  Significance: {entry.significance}")
                if entry.old_value is not None:
                    lines.append(f"  Old: {entry.old_value}")
                if entry.new_value is not None:
                    lines.append(f"  New: {entry.new_value}")
                lines.append("")
        
        else:  # default
            for entry in diff_entries:
                sig_marker = {
                    "critical": "‼️",
                    "major": "⚠️",
                    "minor": "ℹ️",
                }.get(entry.significance, "ℹ️")
                
                if entry.change_type == "added":
                    lines.append(f"{sig_marker} Added {entry.path} = {entry.new_value}")
                elif entry.change_type == "deleted":
                    lines.append(f"{sig_marker} Deleted {entry.path} = {entry.old_value}")
                else:
                    lines.append(
                        f"{sig_marker} Modified {entry.path}: "
                        f"{entry.old_value} → {entry.new_value}"
                    )
        
        return "\n".join(lines)
    
    def validate(self, data: Dict[str, Any]) -> Tuple[bool, List[str]]:
        """Validate parsed data for correctness.
        
        Args:
            data: Parsed data dictionary
            
        Returns:
            Tuple of (is_valid, list_of_errors)
        """
        # Default implementation - subclasses can override
        return True, []
