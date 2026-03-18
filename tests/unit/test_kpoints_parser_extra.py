"""Extra coverage tests for kpoints_parser.py.

Targets uncovered lines:
    71-72  – diff: elif explicit branch + _diff_explicit call
   206    – _diff_explicit: cartesian field change
  221-247  – _diff_line body (comparison of line-mode KPOINTS)
"""

import pytest

from chemvcs.parsers.kpoints_parser import KpointsParser
from chemvcs.parsers.base_parser import DiffEntry


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

EXPLICIT_RECIPROCAL = """\
Explicit k-points
4
Reciprocal
0.0 0.0 0.0  1.0
0.5 0.0 0.0  1.0
0.5 0.5 0.0  1.0
0.5 0.5 0.5  1.0
"""

EXPLICIT_CARTESIAN = """\
Explicit k-points
4
Cartesian
0.0 0.0 0.0  1.0
0.5 0.0 0.0  1.0
0.5 0.5 0.0  1.0
0.5 0.5 0.5  1.0
"""

LINE_MODE_40 = """\
k-points for band structure
40
Line-mode
Reciprocal
0.0 0.0 0.0  ! Gamma
0.5 0.0 0.0  ! X

0.5 0.0 0.0  ! X
0.5 0.5 0.0  ! M
"""

LINE_MODE_50 = """\
k-points for band structure
50
Line-mode
Reciprocal
0.0 0.0 0.0  ! Gamma
0.5 0.0 0.0  ! X

0.5 0.0 0.0  ! X
0.5 0.5 0.0  ! M
"""

LINE_MODE_3SEG = """\
k-points for band structure
40
Line-mode
Reciprocal
0.0 0.0 0.0  ! Gamma
0.5 0.0 0.0  ! X

0.5 0.0 0.0  ! X
0.5 0.5 0.0  ! M

0.5 0.5 0.0  ! M
0.5 0.5 0.5  ! R
"""


# ---------------------------------------------------------------------------
# diff() – explicit KPOINTS branch (lines 71-72)
# ---------------------------------------------------------------------------

class TestDiffExplicit:
    def test_diff_explicit_no_change(self) -> None:
        """Two identical explicit KPOINTS: no diff entries."""
        parser = KpointsParser()
        data1 = parser.parse(EXPLICIT_RECIPROCAL)
        data2 = parser.parse(EXPLICIT_RECIPROCAL)
        entries = parser.diff(data1, data2)
        assert entries == []

    def test_diff_explicit_cartesian_change(self) -> None:
        """Switching from Reciprocal to Cartesian is a major change."""
        parser = KpointsParser()
        old = parser.parse(EXPLICIT_RECIPROCAL)
        new = parser.parse(EXPLICIT_CARTESIAN)
        entries = parser.diff(old, new)
        cartesian_entries = [e for e in entries if e.path == "cartesian"]
        assert len(cartesian_entries) == 1
        assert cartesian_entries[0].old_value is False
        assert cartesian_entries[0].new_value is True
        assert cartesian_entries[0].significance == "major"

    def test_diff_explicit_num_kpoints_change(self) -> None:
        """Changing number of explicit k-points is a critical change."""
        parser = KpointsParser()
        content_2pt = """\
Explicit k-points
2
Reciprocal
0.0 0.0 0.0  1.0
0.5 0.5 0.5  1.0
"""
        old = parser.parse(EXPLICIT_RECIPROCAL)
        new = parser.parse(content_2pt)
        entries = parser.diff(old, new)
        nk_entries = [e for e in entries if e.path == "num_kpoints"]
        assert len(nk_entries) == 1
        assert nk_entries[0].change_type == "modified"


# ---------------------------------------------------------------------------
# diff() – line-mode branch (lines 221-247)
# ---------------------------------------------------------------------------

class TestDiffLineMode:
    def test_diff_line_no_change(self) -> None:
        """Two identical line-mode KPOINTS: no diff."""
        parser = KpointsParser()
        d1 = parser.parse(LINE_MODE_40)
        d2 = parser.parse(LINE_MODE_40)
        entries = parser.diff(d1, d2)
        assert entries == []

    def test_diff_line_num_points_change(self) -> None:
        """Changing num_points is detected as modified (minor)."""
        parser = KpointsParser()
        old = parser.parse(LINE_MODE_40)
        new = parser.parse(LINE_MODE_50)
        entries = parser.diff(old, new)
        np_entries = [e for e in entries if e.path == "num_points"]
        assert len(np_entries) == 1
        assert np_entries[0].change_type == "modified"
        assert np_entries[0].old_value == 40
        assert np_entries[0].new_value == 50

    def test_diff_line_num_segments_change(self) -> None:
        """Adding a segment is detected as modified (major)."""
        parser = KpointsParser()
        old = parser.parse(LINE_MODE_40)   # 2 segments
        new = parser.parse(LINE_MODE_3SEG)  # 3 segments
        entries = parser.diff(old, new)
        seg_entries = [e for e in entries if e.path == "num_segments"]
        assert len(seg_entries) == 1
        assert seg_entries[0].change_type == "modified"
        assert seg_entries[0].old_value == 2
        assert seg_entries[0].new_value == 3
        assert seg_entries[0].significance == "major"

    def test_diff_line_both_changed(self) -> None:
        """Changing both num_points and num_segments in one diff."""
        parser = KpointsParser()
        old = parser.parse(LINE_MODE_40)
        new = parser.parse(LINE_MODE_3SEG.replace("40", "80"))
        entries = parser.diff(old, new)
        paths = {e.path for e in entries}
        assert "num_points" in paths
        assert "num_segments" in paths
