"""Semantic diff engine for VASP files.

This module provides high-level diff functionality that automatically
selects and applies appropriate parsers based on file type.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional

from chemvcs.parsers.base_parser import BaseParser, DiffEntry, ParserError
from chemvcs.parsers.incar_parser import IncarParser
from chemvcs.parsers.kpoints_parser import KpointsParser


class DiffEngine:
    """Main diff engine for semantic comparison of VASP files.
    
    Automatically detects file types and applies appropriate parsers
    to compute semantic diffs rather than line-by-line text diffs.
    """
    
    def __init__(self):
        """Initialize the diff engine with available parsers."""
        self.parsers: Dict[str, BaseParser] = {
            "INCAR": IncarParser(),
            "KPOINTS": KpointsParser(),
            # POSCAR, POTCAR, OUTCAR parsers can be added later
        }
    
    def get_file_type(self, filename: str) -> Optional[str]:
        """Detect file type from filename.
        
        Args:
            filename: Name of the file
            
        Returns:
            File type string (e.g., "INCAR") or None if not recognized
        """
        # Get base filename without path
        base = Path(filename).name.upper()
        
        # Check for known VASP file types
        if base == "INCAR":
            return "INCAR"
        elif base == "KPOINTS":
            return "KPOINTS"
        elif base == "POSCAR" or base == "CONTCAR":
            return "POSCAR"
        elif base.startswith("POTCAR"):
            return "POTCAR"
        elif base == "OUTCAR":
            return "OUTCAR"
        
        return None
    
    def can_parse(self, filename: str) -> bool:
        """Check if file type is supported for parsing.
        
        Args:
            filename: Name of the file
            
        Returns:
            True if file can be parsed semantically
        """
        file_type = self.get_file_type(filename)
        return file_type is not None and file_type in self.parsers
    
    def diff_files(
        self,
        old_content: str,
        new_content: str,
        filename: str,
    ) -> Optional[List[DiffEntry]]:
        """Compute semantic diff between two file versions.
        
        Args:
            old_content: Content of old version
            new_content: Content of new version
            filename: Name of the file (for type detection)
            
        Returns:
            List of DiffEntry objects, or None if file type not supported
            
        Raises:
            ParserError: If parsing fails
        """
        file_type = self.get_file_type(filename)
        
        if file_type is None or file_type not in self.parsers:
            return None
        
        parser = self.parsers[file_type]
        
        # Parse both versions
        try:
            old_data = parser.parse(old_content)
            new_data = parser.parse(new_content)
        except ParserError:
            # If parsing fails, return None to fall back to text diff
            return None
        
        # Compute diff
        return parser.diff(old_data, new_data)
    
    def format_diff(
        self,
        diff_entries: List[DiffEntry],
        style: str = "default",
    ) -> str:
        """Format diff entries for display.
        
        Args:
            diff_entries: List of diff entries
            style: Format style ("default", "compact", "detailed")
            
        Returns:
            Formatted diff string
        """
        if not diff_entries:
            return "No changes"
        
        # Use the BaseParser's format_diff via any parser instance
        # (they all share the same base implementation)
        parser = next(iter(self.parsers.values()))
        return parser.format_diff(diff_entries, style=style)
    
    def summarize_diff(
        self,
        diff_entries: List[DiffEntry],
    ) -> Dict[str, Any]:
        """Generate summary statistics for a diff.
        
        Args:
            diff_entries: List of diff entries
            
        Returns:
            Dictionary with summary information
        """
        if not diff_entries:
            return {
                "total_changes": 0,
                "by_type": {},
                "by_significance": {},
            }
        
        by_type = {"added": 0, "deleted": 0, "modified": 0}
        by_significance = {"critical": 0, "major": 0, "minor": 0}
        
        for entry in diff_entries:
            by_type[entry.change_type] = by_type.get(entry.change_type, 0) + 1
            by_significance[entry.significance] = by_significance.get(entry.significance, 0) + 1
        
        return {
            "total_changes": len(diff_entries),
            "by_type": by_type,
            "by_significance": by_significance,
            "has_critical_changes": by_significance.get("critical", 0) > 0,
        }
    
    def validate_file(
        self,
        content: str,
        filename: str,
    ) -> tuple[bool, List[str]]:
        """Validate file content.
        
        Args:
            content: File content
            filename: Filename for type detection
            
        Returns:
            Tuple of (is_valid, list_of_errors)
        """
        file_type = self.get_file_type(filename)
        
        if file_type is None or file_type not in self.parsers:
            # Can't validate unknown file types
            return True, []
        
        parser = self.parsers[file_type]
        
        try:
            data = parser.parse(content)
            return parser.validate(data)
        except ParserError as e:
            return False, [str(e)]
