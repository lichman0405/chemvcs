"""Unit tests for LammpsLogParser."""

import pytest

from chemvcs.parsers.lammps_log_parser import LammpsLogParser
from chemvcs.parsers.base_parser import ParserError


# ---------------------------------------------------------------------------
# Sample log content fixtures
# ---------------------------------------------------------------------------

# Realistic columnar thermo output (NVT, real units)
SAMPLE_LOG_COLUMNAR = """\
LAMMPS (22 Jul 2025)
  using 1 OpenMP thread(s) per MPI task

Created orthogonal box = (0 0 0) to (20.0 20.0 20.0)
  1 by 1 by 1 MPI processor grid
reading atoms ...
  500 atoms
read_data   data.lammps
Setting up Verlet run ...

Step Temp E_pair E_mol TotEng Press
       0   300.0  -1234.5678       0  -1100.1234   1.23456
     100   301.2  -1238.9012       0  -1104.4568   1.19800
     200   299.8  -1236.4444       0  -1101.9800   1.21000
     300   300.5  -1237.1122       0  -1102.6478   1.20500
Loop time of 12.345 on 1 procs for 300 steps with 500 atoms

Performance: 12.3 ns/day, 1.95 hours/ns, 24.3 timesteps/s
Total wall time: 0:00:12
"""

# Partial log – no wall time (run did not complete)
SAMPLE_LOG_INCOMPLETE = """\
LAMMPS (22 Jul 2025)
  256 atoms

Step Temp TotEng Press
       0   300.0   -500.0   2.0
     100   302.5   -502.3   1.9
"""

# YAML thermo format (LAMMPS >= 24Mar2022)
SAMPLE_LOG_YAML = """\
LAMMPS (22 Jul 2025)
Running on 1 MPI tasks
  1000 atoms

---
keywords: ['Step', 'Temp', 'TotEng', 'PotEng', 'KinEng', 'Press', ]
data:
- [0, 300.0, -2000.5, -2300.5, 300.0, 1.5, ]
- [100, 298.5, -2005.2, -2305.2, 300.0, 1.4, ]
- [200, 299.1, -2003.8, -2303.8, 300.0, 1.45, ]
...

Total wall time: 0:01:00
"""

# Minimal log – only essential content
SAMPLE_LOG_MINIMAL = """\
100 atoms
Step Temp TotEng
0 300.0 -100.0
Total wall time: 0:00:01
"""


# ---------------------------------------------------------------------------
# Test cases
# ---------------------------------------------------------------------------

class TestLammpsLogParser:
    """Tests for LammpsLogParser."""

    def test_parse_columnar_returns_dict(self) -> None:
        parser = LammpsLogParser()
        data = parser.parse(SAMPLE_LOG_COLUMNAR)
        assert isinstance(data, dict)

    def test_parse_columnar_num_atoms(self) -> None:
        parser = LammpsLogParser()
        data = parser.parse(SAMPLE_LOG_COLUMNAR)
        assert data["num_atoms"] == 500

    def test_parse_columnar_completed(self) -> None:
        parser = LammpsLogParser()
        data = parser.parse(SAMPLE_LOG_COLUMNAR)
        assert data["completed"] is True
        assert data["wall_time"] == "0:00:12"

    def test_parse_columnar_final_toteng(self) -> None:
        parser = LammpsLogParser()
        data = parser.parse(SAMPLE_LOG_COLUMNAR)
        assert data["final_etotal"] == pytest.approx(-1102.6478)

    def test_parse_columnar_final_temp(self) -> None:
        parser = LammpsLogParser()
        data = parser.parse(SAMPLE_LOG_COLUMNAR)
        assert data["final_temp"] == pytest.approx(300.5)

    def test_parse_columnar_keywords(self) -> None:
        parser = LammpsLogParser()
        data = parser.parse(SAMPLE_LOG_COLUMNAR)
        assert "Step" in data.get("thermo_keywords", [])
        assert "TotEng" in data.get("thermo_keywords", [])

    def test_parse_incomplete_not_completed(self) -> None:
        parser = LammpsLogParser()
        data = parser.parse(SAMPLE_LOG_INCOMPLETE)
        assert data["completed"] is False

    def test_parse_incomplete_has_energy(self) -> None:
        parser = LammpsLogParser()
        data = parser.parse(SAMPLE_LOG_INCOMPLETE)
        assert data["final_etotal"] == pytest.approx(-502.3)

    def test_parse_yaml_format(self) -> None:
        pytest.importorskip("yaml")
        parser = LammpsLogParser()
        data = parser.parse(SAMPLE_LOG_YAML)
        assert data["completed"] is True
        assert data.get("num_atoms") == 1000

    def test_parse_yaml_final_values(self) -> None:
        pytest.importorskip("yaml")
        parser = LammpsLogParser()
        data = parser.parse(SAMPLE_LOG_YAML)
        assert data["final_etotal"] == pytest.approx(-2003.8)
        assert data["final_temp"] == pytest.approx(299.1)

    def test_parse_minimal(self) -> None:
        parser = LammpsLogParser()
        data = parser.parse(SAMPLE_LOG_MINIMAL)
        assert data["completed"] is True
        assert data["final_etotal"] == pytest.approx(-100.0)

    def test_parse_empty_raises(self) -> None:
        parser = LammpsLogParser()
        with pytest.raises(ParserError):
            parser.parse("")

    def test_parse_whitespace_only_raises(self) -> None:
        parser = LammpsLogParser()
        with pytest.raises(ParserError):
            parser.parse("   \n\n  ")

    # -- diff tests --

    def test_diff_detects_energy_change(self) -> None:
        parser = LammpsLogParser()
        old = parser.parse(SAMPLE_LOG_COLUMNAR)
        new = parser.parse(SAMPLE_LOG_INCOMPLETE)
        entries = parser.diff(old, new)
        paths = [e.path for e in entries]
        assert "final_etotal" in paths

    def test_diff_no_change(self) -> None:
        parser = LammpsLogParser()
        data = parser.parse(SAMPLE_LOG_COLUMNAR)
        entries = parser.diff(data, data)
        assert entries == []

    def test_diff_completion_status_change(self) -> None:
        parser = LammpsLogParser()
        old = parser.parse(SAMPLE_LOG_COLUMNAR)
        new = parser.parse(SAMPLE_LOG_INCOMPLETE)
        entries = parser.diff(old, new)
        completed_entries = [e for e in entries if e.path == "completed"]
        assert completed_entries
        assert completed_entries[0].old_value is True
        assert completed_entries[0].new_value is False

    def test_diff_critical_energy_significance(self) -> None:
        parser = LammpsLogParser()
        old = parser.parse(SAMPLE_LOG_COLUMNAR)
        new = parser.parse(SAMPLE_LOG_INCOMPLETE)
        entries = parser.diff(old, new)
        etotal_entry = next((e for e in entries if e.path == "final_etotal"), None)
        assert etotal_entry is not None
        assert etotal_entry.significance == "critical"

    # -- validate tests --

    def test_validate_complete_run(self) -> None:
        parser = LammpsLogParser()
        data = parser.parse(SAMPLE_LOG_COLUMNAR)
        is_valid, errors = parser.validate(data)
        assert is_valid
        assert errors == []

    def test_validate_incomplete_run(self) -> None:
        parser = LammpsLogParser()
        data = parser.parse(SAMPLE_LOG_INCOMPLETE)
        is_valid, errors = parser.validate(data)
        assert not is_valid
        assert any("complete" in e.lower() for e in errors)
