"""KPOINTS file parser for VASP k-point sampling.

KPOINTS files specify the k-point mesh for Brillouin zone sampling.
Multiple formats are supported: automatic, explicit, line-mode, etc.
"""

import re
from typing import Any, Dict, List, Optional, Tuple

from chemvcs.parsers.base_parser import BaseParser, DiffEntry, ParserError


class KpointsParser(BaseParser):
    """Parser for VASP KPOINTS files.
    
    KPOINTS files specify k-point sampling schemes:
    - Line 1: Comment/description
    - Line 2: Number of k-points (0 = automatic)
    - Line 3: Generation style (Automatic/Gamma/Monkhorst-Pack/etc.)
    - Line 4+: Grid or k-point specifications
    
    Formats supported:
    1. Automatic mesh (M/Monkhorst-Pack, G/Gamma)
    2. Explicit k-points
    3. Line-mode for band structure
    """
    
    def parse(self, content: str) -> Dict[str, Any]:
        """Parse KPOINTS content into structured data.
        
        Args:
            content: Raw KPOINTS file content
            
        Returns:
            Dictionary containing parsed k-point data
            
        Raises:
            ParserError: If content is malformed
        """
        try:
            lines = [l.strip() for l in content.split("\n") if l.strip()]
            
            if len(lines) < 3:
                raise ParserError("KPOINTS file too short (minimum 3 lines required)")
            
            # Line 1: Comment
            comment = lines[0]
            
            # Line 2: Number of k-points
            try:
                num_kpoints = int(lines[1])
            except ValueError:
                raise ParserError(f"Invalid number of k-points: {lines[1]}")
            
            # Line 3: Generation style
            style_line = lines[2].upper()
            
            # Determine generation style
            if num_kpoints == 0 or "AUTO" in style_line:
                # Automatic mesh generation
                return self._parse_automatic(lines, comment, style_line)
            elif "LINE" in style_line:
                # Line mode for band structure
                return self._parse_line_mode(lines, comment)
            else:
                # Explicit k-points
                return self._parse_explicit(lines, comment, num_kpoints)
        
        except ParserError:
            raise
        except Exception as e:
            raise ParserError(f"Failed to parse KPOINTS: {e}") from e
    
    def _parse_automatic(
        self,
        lines: List[str],
        comment: str,
        style_line: str,
    ) -> Dict[str, Any]:
        """Parse automatic k-point mesh generation.
        
        Format:
            Comment
            0
            Gamma (or Monkhorst-Pack)
            nx ny nz
            [shift_x shift_y shift_z]
        """
        if len(lines) < 4:
            raise ParserError("Automatic KPOINTS requires at least 4 lines")
        
        # Determine if Gamma-centered or Monkhorst-Pack
        gamma_centered = "G" in style_line
        
        # Parse grid dimensions
        grid_line = lines[3].split()
        if len(grid_line) < 3:
            raise ParserError(f"Invalid grid specification: {lines[3]}")
        
        try:
            grid = [int(grid_line[0]), int(grid_line[1]), int(grid_line[2])]
        except ValueError:
            raise ParserError(f"Invalid grid values: {lines[3]}")
        
        # Parse optional shift
        shift = [0.0, 0.0, 0.0]
        if len(lines) > 4:
            shift_line = lines[4].split()
            if len(shift_line) >= 3:
                try:
                    shift = [float(shift_line[0]), float(shift_line[1]), float(shift_line[2])]
                except ValueError:
                    pass  # Use default shift
        
        return {
            "type": "automatic",
            "comment": comment,
            "gamma_centered": gamma_centered,
            "grid": grid,
            "shift": shift,
        }
    
    def _parse_explicit(
        self,
        lines: List[str],
        comment: str,
        num_kpoints: int,
    ) -> Dict[str, Any]:
        """Parse explicit k-point list.
        
        Format:
            Comment
            N
            Cartesian (or Reciprocal)
            kx ky kz weight
            ...
        """
        if len(lines) < 3 + num_kpoints:
            raise ParserError(f"Expected {num_kpoints} k-points but file has fewer lines")
        
        # Determine coordinate system
        coord_line = lines[2].upper()
        cartesian = coord_line.startswith("C") or "CART" in coord_line
        
        # Parse k-points
        kpoints = []
        for i in range(num_kpoints):
            line = lines[3 + i].split()
            if len(line) < 3:
                raise ParserError(f"Invalid k-point line: {lines[3 + i]}")
            
            try:
                kx, ky, kz = float(line[0]), float(line[1]), float(line[2])
                weight = float(line[3]) if len(line) > 3 else 1.0
                kpoints.append({"coords": [kx, ky, kz], "weight": weight})
            except ValueError as e:
                raise ParserError(f"Invalid k-point values: {lines[3 + i]}") from e
        
        return {
            "type": "explicit",
            "comment": comment,
            "num_kpoints": num_kpoints,
            "cartesian": cartesian,
            "kpoints": kpoints,
        }
    
    def _parse_line_mode(
        self,
        lines: List[str],
        comment: str,
    ) -> Dict[str, Any]:
        """Parse line-mode for band structure calculations.
        
        Format:
            Comment
            N (points per segment)
            Line-mode
            Reciprocal
            kx ky kz  ! Label
            kx ky kz  ! Label
            (blank line separates segments)
        """
        # Find coordinate system line
        coord_line = lines[2].upper() if len(lines) > 2 else ""
        reciprocal = "REC" in coord_line or not ("CART" in coord_line)
        
        # Parse segments
        segments = []
        i = 3  # Start after header
        
        while i < len(lines):
            line = lines[i].strip()
            if not line:
                i += 1
                continue
            
            # Parse start point
            start_parts = line.split("!")
            start_coords = start_parts[0].split()
            if len(start_coords) < 3:
                i += 1
                continue
            
            start_label = start_parts[1].strip() if len(start_parts) > 1 else ""
            
            # Parse end point (next line)
            if i + 1 < len(lines):
                end_parts = lines[i + 1].split("!")
                end_coords = end_parts[0].split()
                if len(end_coords) >= 3:
                    end_label = end_parts[1].strip() if len(end_parts) > 1 else ""
                    
                    try:
                        segments.append({
                            "start": {
                                "coords": [float(start_coords[0]), float(start_coords[1]), float(start_coords[2])],
                                "label": start_label,
                            },
                            "end": {
                                "coords": [float(end_coords[0]), float(end_coords[1]), float(end_coords[2])],
                                "label": end_label,
                            },
                        })
                    except ValueError:
                        pass
            
            i += 2  # Move to next segment
        
        # Get number of points per segment
        num_points = 0
        try:
            num_points = int(lines[1])
        except (ValueError, IndexError):
            num_points = 40  # Default
        
        return {
            "type": "line",
            "comment": comment,
            "num_points": num_points,
            "reciprocal": reciprocal,
            "segments": segments,
        }
    
    def diff(
        self,
        old_data: Dict[str, Any],
        new_data: Dict[str, Any],
    ) -> List[DiffEntry]:
        """Compute semantic differences between two KPOINTS files.
        
        Args:
            old_data: Parsed KPOINTS from old version
            new_data: Parsed KPOINTS from new version
            
        Returns:
            List of DiffEntry objects describing changes
        """
        diff_entries = []
        
        # Check if type changed
        old_type = old_data.get("type")
        new_type = new_data.get("type")
        
        if old_type != new_type:
            diff_entries.append(
                DiffEntry(
                    path="type",
                    old_value=old_type,
                    new_value=new_type,
                    change_type="modified",
                    significance="critical",
                )
            )
            # Type changed, everything is different
            return diff_entries
        
        # Compare based on type
        if old_type == "automatic":
            diff_entries.extend(self._diff_automatic(old_data, new_data))
        elif old_type == "explicit":
            diff_entries.extend(self._diff_explicit(old_data, new_data))
        elif old_type == "line":
            diff_entries.extend(self._diff_line(old_data, new_data))
        
        return diff_entries
    
    def _diff_automatic(
        self,
        old_data: Dict[str, Any],
        new_data: Dict[str, Any],
    ) -> List[DiffEntry]:
        """Diff automatic mesh generation."""
        entries = []
        
        # Check gamma_centered
        if old_data.get("gamma_centered") != new_data.get("gamma_centered"):
            entries.append(
                DiffEntry(
                    path="gamma_centered",
                    old_value=old_data.get("gamma_centered"),
                    new_value=new_data.get("gamma_centered"),
                    change_type="modified",
                    significance="major",
                )
            )
        
        # Check grid
        if old_data.get("grid") != new_data.get("grid"):
            entries.append(
                DiffEntry(
                    path="grid",
                    old_value=old_data.get("grid"),
                    new_value=new_data.get("grid"),
                    change_type="modified",
                    significance="critical",
                )
            )
        
        # Check shift
        if old_data.get("shift") != new_data.get("shift"):
            entries.append(
                DiffEntry(
                    path="shift",
                    old_value=old_data.get("shift"),
                    new_value=new_data.get("shift"),
                    change_type="modified",
                    significance="major",
                )
            )
        
        return entries
    
    def _diff_explicit(
        self,
        old_data: Dict[str, Any],
        new_data: Dict[str, Any],
    ) -> List[DiffEntry]:
        """Diff explicit k-points."""
        entries = []
        
        # Check number of k-points
        if old_data.get("num_kpoints") != new_data.get("num_kpoints"):
            entries.append(
                DiffEntry(
                    path="num_kpoints",
                    old_value=old_data.get("num_kpoints"),
                    new_value=new_data.get("num_kpoints"),
                    change_type="modified",
                    significance="critical",
                )
            )
        
        # Check coordinate system
        if old_data.get("cartesian") != new_data.get("cartesian"):
            entries.append(
                DiffEntry(
                    path="cartesian",
                    old_value=old_data.get("cartesian"),
                    new_value=new_data.get("cartesian"),
                    change_type="modified",
                    significance="major",
                )
            )
        
        return entries
    
    def _diff_line(
        self,
        old_data: Dict[str, Any],
        new_data: Dict[str, Any],
    ) -> List[DiffEntry]:
        """Diff line-mode k-points."""
        entries = []
        
        # Check number of points
        if old_data.get("num_points") != new_data.get("num_points"):
            entries.append(
                DiffEntry(
                    path="num_points",
                    old_value=old_data.get("num_points"),
                    new_value=new_data.get("num_points"),
                    change_type="modified",
                    significance="minor",
                )
            )
        
        # Check number of segments
        old_segs = len(old_data.get("segments", []))
        new_segs = len(new_data.get("segments", []))
        if old_segs != new_segs:
            entries.append(
                DiffEntry(
                    path="num_segments",
                    old_value=old_segs,
                    new_value=new_segs,
                    change_type="modified",
                    significance="major",
                )
            )
        
        return entries
