"""INCAR file parser for VASP input parameters.

INCAR files contain key-value pairs controlling VASP calculation parameters.
Format: KEY = VALUE [; KEY = VALUE ...]
"""

import re
from typing import Any, Dict, List, Optional

from chemvcs.parsers.base_parser import BaseParser, DiffEntry, ParserError


class IncarParser(BaseParser):
    """Parser for VASP INCAR files.
    
    INCAR files specify calculation parameters in key-value format:
    - One parameter per line (or semicolon-separated)
    - Comments start with # or !
    - Values can be integers, floats, strings, or boolean-like
    - Some tags accept multiple values (space-separated)
    
    Example:
        ENCUT = 520
        PREC = Accurate
        ISMEAR = 0; SIGMA = 0.05
        LDAU = .TRUE.
    """
    
    # Critical parameters that significantly affect results
    CRITICAL_TAGS = {
        "ENCUT", "PREC", "ALGO", "LREAL",
        "ISMEAR", "SIGMA", "EDIFF", "EDIFFG",
        "NSW", "IBRION", "POTIM", "ISIF",
        "NELM", "NELMIN",
    }
    
    # Major parameters that affect computational approach
    MAJOR_TAGS = {
        "LWAVE", "LCHARG", "LORBIT", "LVTOT", "LVHAR",
        "ISPIN", "MAGMOM", "NCORE", "KPAR",
        "LDAU", "LDAUTYPE", "LDAUL", "LDAUU", "LDAUJ",
        "IVDW", "LUSE_VDW", "AGGAC",
    }
    
    def parse(self, content: str) -> Dict[str, Any]:
        """Parse INCAR content into structured data.
        
        Args:
            content: Raw INCAR file content
            
        Returns:
            Dictionary mapping parameter names to values
            
        Raises:
            ParserError: If content is malformed
        """
        try:
            parameters = {}
            
            # Process line by line
            for line_num, line in enumerate(content.split("\n"), start=1):
                # Remove comments (# or !)
                line = re.split(r"[#!]", line)[0].strip()
                
                if not line:
                    continue
                
                # Split by semicolon for multiple assignments
                assignments = line.split(";")
                
                for assignment in assignments:
                    assignment = assignment.strip()
                    if not assignment:
                        continue
                    
                    # Parse KEY = VALUE
                    match = re.match(r"^([A-Z_][A-Z0-9_]*)\s*=\s*(.+)$", assignment, re.IGNORECASE)
                    
                    if not match:
                        # Not a valid assignment, skip
                        continue
                    
                    key = match.group(1).upper()
                    value_str = match.group(2).strip()
                    
                    # Parse value
                    value = self._parse_value(value_str)
                    parameters[key] = value
            
            return parameters
        
        except Exception as e:
            raise ParserError(f"Failed to parse INCAR: {e}") from e
    
    def _parse_value(self, value_str: str) -> Any:
        """Parse a value string to appropriate Python type.
        
        Args:
            value_str: Raw value string
            
        Returns:
            Parsed value (int, float, bool, str, or list)
        """
        value_str = value_str.strip()
        
        # Boolean-like values
        if value_str.upper() in (".TRUE.", "T", "TRUE"):
            return True
        if value_str.upper() in (".FALSE.", "F", "FALSE"):
            return False
        
        # Try parsing as number
        try:
            # Check if it's an integer
            if "." not in value_str and "e" not in value_str.lower():
                return int(value_str)
            else:
                return float(value_str)
        except ValueError:
            pass
        
        # Check for multiple values (space-separated)
        parts = value_str.split()
        if len(parts) > 1:
            # Try parsing as list of numbers (all parts must be numeric)
            all_numeric = True
            for p in parts:
                try:
                    float(p)
                except ValueError:
                    all_numeric = False
                    break
            
            if all_numeric:
                # Parse as list of numbers
                try:
                    # Try integers first
                    if all("." not in p and "e" not in p.lower() for p in parts):
                        return [int(p) for p in parts]
                    else:
                        return [float(p) for p in parts]
                except ValueError:
                    pass
            
            # Not all numeric, return as full string
            return value_str
        
        # Return as string
        return value_str
    
    def diff(
        self,
        old_data: Dict[str, Any],
        new_data: Dict[str, Any],
    ) -> List[DiffEntry]:
        """Compute semantic differences between two INCAR files.
        
        Args:
            old_data: Parsed INCAR from old version
            new_data: Parsed INCAR from new version
            
        Returns:
            List of DiffEntry objects describing changes
        """
        diff_entries = []
        
        # Get all keys
        all_keys = set(old_data.keys()) | set(new_data.keys())
        
        for key in sorted(all_keys):
            old_value = old_data.get(key)
            new_value = new_data.get(key)
            
            # Determine significance
            if key in self.CRITICAL_TAGS:
                significance = "critical"
            elif key in self.MAJOR_TAGS:
                significance = "major"
            else:
                significance = "minor"
            
            # Determine change type
            if old_value is None:
                # Added
                diff_entries.append(
                    DiffEntry(
                        path=key,
                        old_value=None,
                        new_value=new_value,
                        change_type="added",
                        significance=significance,
                    )
                )
            elif new_value is None:
                # Deleted
                diff_entries.append(
                    DiffEntry(
                        path=key,
                        old_value=old_value,
                        new_value=None,
                        change_type="deleted",
                        significance=significance,
                    )
                )
            elif old_value != new_value:
                # Modified
                diff_entries.append(
                    DiffEntry(
                        path=key,
                        old_value=old_value,
                        new_value=new_value,
                        change_type="modified",
                        significance=significance,
                    )
                )
        
        return diff_entries
    
    def validate(self, data: Dict[str, Any]) -> tuple[bool, List[str]]:
        """Validate INCAR parameters.
        
        Args:
            data: Parsed INCAR data
            
        Returns:
            Tuple of (is_valid, list_of_error_messages)
        """
        errors = []
        
        # Check ENCUT is positive
        if "ENCUT" in data:
            if not isinstance(data["ENCUT"], (int, float)) or data["ENCUT"] <= 0:
                errors.append("ENCUT must be a positive number")
        
        # Check ISMEAR range
        if "ISMEAR" in data:
            if isinstance(data["ISMEAR"], int):
                if data["ISMEAR"] < -5 or data["ISMEAR"] > 5:
                    errors.append("ISMEAR should be in range -5 to 5")
        
        # Check SIGMA is positive
        if "SIGMA" in data:
            if not isinstance(data["SIGMA"], (int, float)) or data["SIGMA"] <= 0:
                errors.append("SIGMA must be a positive number")
        
        # Check NSW is non-negative
        if "NSW" in data:
            if not isinstance(data["NSW"], int) or data["NSW"] < 0:
                errors.append("NSW must be a non-negative integer")
        
        # Warn if ISMEAR=0 but SIGMA is large
        if "ISMEAR" in data and "SIGMA" in data:
            if data["ISMEAR"] == 0 and isinstance(data["SIGMA"], (int, float)) and data["SIGMA"] > 0.2:
                errors.append("Warning: ISMEAR=0 (Gaussian) with large SIGMA (>0.2) may be inappropriate")
        
        is_valid = len(errors) == 0
        return is_valid, errors
