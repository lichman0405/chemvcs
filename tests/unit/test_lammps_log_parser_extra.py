"""Extra coverage tests for lammps_log_parser.py.

Targets uncovered lines (see coverage report):
    87-91   – regex fallback path in parse()
   122     – diff: added-field branch
   132     – diff: float-noise continue
   145     – validate: invalid atom count error
  179-180  – YAML ImportError guard  (yaml not installed → skip)
  194-195  – YAML: multiple '---' blocks (prev block flushed)
   211     – YAML: block open at EOF (no '...' terminator)
  223-224  – YAML: yaml.safe_load exception → continue
   227     – YAML: non-dict doc → continue
   229     – YAML: missing keywords/data → continue
   235     – YAML: empty kws/rows → continue
   247     – YAML: all blocks invalid → return None
   292     – columnar: row-length mismatch → break
   297     – columnar: no rows → continue
   307     – columnar: no valid header → return None
  315-329  – _try_regex_fallback body
"""

import pytest

from chemvcs.parsers.lammps_log_parser import LammpsLogParser
from chemvcs.parsers.base_parser import ParserError


# ---------------------------------------------------------------------------
# Helpers / fixtures content
# ---------------------------------------------------------------------------

# Content that has no YAML blocks AND no columnar "Step" header,
# but has a TotEng pattern for the regex fallback.
REGEX_FALLBACK_CONTENT = """\
LAMMPS (22 Jul 2025)
  500 atoms

TotEng -1234.5678
Temp 300.0
Total wall time: 0:00:05
"""

# Same but without wall-time (incomplete run)
REGEX_FALLBACK_NO_WALLTIME = """\
LAMMPS (22 Jul 2025)
  500 atoms
TotEng -1234.5678
"""

# Content with two consecutive YAML '---' blocks (exercises lines 194-195)
YAML_MULTI_BLOCK = """\
LAMMPS (22 Jul 2025)
  1000 atoms

---
keywords: ["Step", "Temp", "TotEng"]
data:
- [0, 300.0, -2000.0]
- [50, 299.0, -2001.0]
---
keywords: ["Step", "Temp", "TotEng"]
data:
- [100, 298.0, -2003.0]
- [200, 297.5, -2005.0]
...

Total wall time: 0:01:00
"""

# YAML block without a '...' terminator (block still open at EOF → line 211)
YAML_OPEN_AT_EOF = """\
LAMMPS (22 Jul 2025)
  1000 atoms

---
keywords: ["Step", "Temp", "TotEng"]
data:
- [0, 300.0, -2000.0]
- [100, 298.5, -2002.5]

Total wall time: 0:00:30
"""

# YAML block where yaml.safe_load raises (invalid YAML → exercises 223-224)
YAML_INVALID_SYNTAX = """\
LAMMPS (22 Jul 2025)
  100 atoms

---
keywords: unclosed_bracket: [
data: bad yaml
...

Step Temp TotEng
0 300.0 -500.0
100 299.0 -502.0
Total wall time: 0:00:01
"""

# YAML block where doc is not a dict (scalar value → exercises line 227)
YAML_NON_DICT = """\
LAMMPS (22 Jul 2025)
  100 atoms

---
just a plain string
...

Step Temp TotEng
0 300.0 -500.0
Total wall time: 0:00:01
"""

# YAML block without 'keywords' key (exercises line 229)
YAML_MISSING_KEYWORDS = """\
LAMMPS (22 Jul 2025)
  100 atoms

---
some_key: value
other_key: 123
...

Step Temp TotEng
0 300.0 -600.0
Total wall time: 0:00:01
"""

# YAML block with empty keywords list (exercises line 235)
YAML_EMPTY_KEYWORDS = """\
LAMMPS (22 Jul 2025)
  100 atoms

---
keywords: []
data: []
...

Step Temp TotEng
0 300.0 -700.0
Total wall time: 0:00:01
"""

# Content that looks like it has a YAML block but all blocks are invalid
# (so _try_parse_yaml_thermo returns None → line 247)
YAML_ALL_INVALID_BLOCKS = """\
LAMMPS (22 Jul 2025)
  100 atoms

---
not_keywords: something
not_data: []
...

Step Temp TotEng
0 300.0 -800.0
Total wall time: 0:00:01
"""

# Columnar content where the data row has wrong column count (exercises line 292)
# The second data row has only 2 values vs 3 in header → row length mismatch
COLUMNAR_ROW_LENGTH_MISMATCH = """\
LAMMPS
  200 atoms

Step Temp TotEng
       0   300.0  -1000.0
     100   299.0
Total wall time: 0:00:10
"""

# Columnar header with NO subsequent numeric rows (exercises line 297)
COLUMNAR_NO_ROWS = """\
LAMMPS
  200 atoms

Step Temp TotEng
This line is not numeric.
Neither is this.
Total wall time: 0:00:10
"""

# Content that has absolutely no thermo data at all → ParserError
CONTENT_NO_THERMO = """\
LAMMPS (22 Jul 2025)
  200 atoms
Some random output without any thermodynamic columns.
"""


# ---------------------------------------------------------------------------
# diff() edge cases
# ---------------------------------------------------------------------------

class TestLammpsLogDiff:
    def test_diff_added_field(self) -> None:
        """Field present in new but not old → added entry."""
        parser = LammpsLogParser()
        old: dict = {"final_etotal": -1000.0, "completed": True}
        new: dict = {"final_etotal": -1000.0, "completed": True, "wall_time": "0:01:00"}
        entries = parser.diff(old, new)
        added = [e for e in entries if e.path == "wall_time" and e.change_type == "added"]
        assert len(added) == 1
        assert added[0].new_value == "0:01:00"

    def test_diff_deleted_field(self) -> None:
        """Field present in old but not new → deleted entry."""
        parser = LammpsLogParser()
        old: dict = {"final_etotal": -1000.0, "num_atoms": 500, "completed": True}
        new: dict = {"final_etotal": -1000.0, "completed": True}
        entries = parser.diff(old, new)
        deleted = [e for e in entries if e.path == "num_atoms" and e.change_type == "deleted"]
        assert len(deleted) == 1
        assert deleted[0].old_value == 500

    def test_diff_float_noise_suppressed(self) -> None:
        """Energy difference below 1e-10 is suppressed."""
        parser = LammpsLogParser()
        old: dict = {"final_etotal": -1234.5678, "completed": True}
        new: dict = {"final_etotal": -1234.5678 + 1e-12, "completed": True}
        entries = parser.diff(old, new)
        energy_entries = [e for e in entries if e.path == "final_etotal"]
        assert len(energy_entries) == 0

    def test_diff_float_above_noise_threshold_kept(self) -> None:
        """Energy difference >= 1e-10 is NOT suppressed."""
        parser = LammpsLogParser()
        old: dict = {"final_etotal": -1234.0, "completed": True}
        new: dict = {"final_etotal": -1235.0, "completed": True}
        entries = parser.diff(old, new)
        energy_entries = [e for e in entries if e.path == "final_etotal"]
        assert len(energy_entries) == 1


# ---------------------------------------------------------------------------
# validate()
# ---------------------------------------------------------------------------

class TestLammpsLogValidate:
    def test_validate_incomplete_run_gives_error(self) -> None:
        parser = LammpsLogParser()
        data: dict = {"completed": False, "final_etotal": -1000.0}
        is_valid, errors = parser.validate(data)
        assert not is_valid
        assert any("complete" in e.lower() for e in errors)

    def test_validate_invalid_atom_count_gives_error(self) -> None:
        """num_atoms <= 0 should be flagged as invalid."""
        parser = LammpsLogParser()
        data: dict = {"completed": True, "num_atoms": 0}
        is_valid, errors = parser.validate(data)
        assert not is_valid
        assert any("atom" in e.lower() for e in errors)

    def test_validate_negative_atom_count_gives_error(self) -> None:
        parser = LammpsLogParser()
        data: dict = {"completed": True, "num_atoms": -5}
        is_valid, errors = parser.validate(data)
        assert not is_valid

    def test_validate_complete_valid_run_passes(self) -> None:
        parser = LammpsLogParser()
        data: dict = {"completed": True, "num_atoms": 500, "final_etotal": -1000.0}
        is_valid, errors = parser.validate(data)
        assert is_valid
        assert errors == []


# ---------------------------------------------------------------------------
# parse() – regex fallback path (lines 87-91, 315-329)
# ---------------------------------------------------------------------------

class TestLammpsLogParseRegexFallback:
    def test_regex_fallback_extracts_toteng(self) -> None:
        parser = LammpsLogParser()
        data = parser.parse(REGEX_FALLBACK_CONTENT)
        assert data.get("final_etotal") == pytest.approx(-1234.5678)
        assert data.get("num_atoms") == 500
        assert data.get("completed") is True

    def test_regex_fallback_incomplete_run(self) -> None:
        """Regex fallback without wall time → completed=False."""
        parser = LammpsLogParser()
        data = parser.parse(REGEX_FALLBACK_NO_WALLTIME)
        assert data.get("final_etotal") == pytest.approx(-1234.5678)
        assert data.get("completed") is False

    def test_no_thermo_data_raises_parser_error(self) -> None:
        """Content with no YAML, no columnar, no TotEng/Temp → ParserError."""
        parser = LammpsLogParser()
        with pytest.raises(ParserError):
            parser.parse(CONTENT_NO_THERMO)


# ---------------------------------------------------------------------------
# _try_parse_yaml_thermo edge cases
# ---------------------------------------------------------------------------

class TestYamlThermoEdgeCases:
    def test_yaml_multi_block(self) -> None:
        """Multiple '---' blocks (exercises lines 194-195)."""
        pytest.importorskip("yaml")
        parser = LammpsLogParser()
        data = parser.parse(YAML_MULTI_BLOCK)
        # Should pick up data from the last block
        assert data.get("completed") is True

    def test_yaml_block_open_at_eof(self) -> None:
        """YAML block still open at EOF (no '...') exercises line 211."""
        pytest.importorskip("yaml")
        parser = LammpsLogParser()
        data = parser.parse(YAML_OPEN_AT_EOF)
        assert data.get("completed") is True

    def test_yaml_invalid_syntax_falls_back_to_columnar(self) -> None:
        """Malformed YAML (safe_load raises) → lines 223-224; columnar wins."""
        pytest.importorskip("yaml")
        parser = LammpsLogParser()
        data = parser.parse(YAML_INVALID_SYNTAX)
        # Falls through to columnar parse
        assert data.get("final_etotal") is not None

    def test_yaml_non_dict_falls_back(self) -> None:
        """YAML doc is not a dict (line 227) → falls back to columnar."""
        pytest.importorskip("yaml")
        parser = LammpsLogParser()
        data = parser.parse(YAML_NON_DICT)
        assert data.get("final_etotal") is not None

    def test_yaml_missing_keywords_falls_back(self) -> None:
        """YAML doc missing 'keywords' key (line 229) → falls back to columnar."""
        pytest.importorskip("yaml")
        parser = LammpsLogParser()
        data = parser.parse(YAML_MISSING_KEYWORDS)
        assert data.get("final_etotal") is not None

    def test_yaml_empty_keywords_falls_back(self) -> None:
        """YAML block with empty kws (line 235) → falls back to columnar."""
        pytest.importorskip("yaml")
        parser = LammpsLogParser()
        data = parser.parse(YAML_EMPTY_KEYWORDS)
        assert data.get("final_etotal") is not None

    def test_yaml_all_invalid_blocks_returns_none_falls_back(self) -> None:
        """All YAML blocks missing keywords/data → _try_parse_yaml_thermo returns
        None (line 247) → columnar succeeds."""
        pytest.importorskip("yaml")
        parser = LammpsLogParser()
        data = parser.parse(YAML_ALL_INVALID_BLOCKS)
        # Falls back to columnar
        assert data.get("final_etotal") is not None


# ---------------------------------------------------------------------------
# _try_parse_columnar_thermo edge cases
# ---------------------------------------------------------------------------

class TestColumnarEdgeCases:
    def test_columnar_row_length_mismatch_still_parses_first_row(self) -> None:
        """Row with wrong column count stops data collection (line 292)."""
        parser = LammpsLogParser()
        data = parser.parse(COLUMNAR_ROW_LENGTH_MISMATCH)
        # Only the first row (0, 300.0, -1000.0) is good
        assert data.get("final_etotal") == pytest.approx(-1000.0)

    def test_columnar_no_rows_skips_header(self) -> None:
        """Header found but no following numeric rows (line 297).
        Parser falls through to None and then to regex/error."""
        parser = LammpsLogParser()
        # This content has no parseable data at all after the header
        with pytest.raises(ParserError):
            parser.parse(COLUMNAR_NO_ROWS)
