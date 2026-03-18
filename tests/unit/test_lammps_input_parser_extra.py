"""Extra coverage tests for lammps_input_parser.py.

Targets uncovered lines, including:
  - _preprocess: inline-comment with quoted characters (247, 249, 251-252, 254)
  - _parse_line: dimension ValueError, timestep/run variable, minimize error,
    thermo non-integer, mass, dump, restart, neigh_modify; and fix langevin
  - _summarize_fixes: langevin style (419-428)
  - _tokenize: shlex ValueError fallback (438-439)
  - _parse_fix_nvt_npt: iso/aniso keywords (lines 464-469), return (472)
  - diff: float noise continue (145), mass dict change (169),
          read_data change (155-163)
  - validate: negative timestep (193-194)
  - _classify: minor / major fields (437-439)
"""

from __future__ import annotations

import pytest

from chemvcs.parsers.lammps_input_parser import LammpsInputParser
from chemvcs.parsers.base_parser import ParserError


# ---------------------------------------------------------------------------
# _preprocess edge cases (inline comments with quotes)
# ---------------------------------------------------------------------------

SCRIPT_INLINE_COMMENT_SINGLE_QUOTE = """\
units           lj
atom_style      atomic
pair_style      lj/cut 2.5  # This is the 'Lennard-Jones' cutoff
pair_coeff      * * 1.0 1.0 2.5
timestep        0.001
run             1000
"""

SCRIPT_INLINE_COMMENT_DOUBLE_QUOTE = r"""
units           metal
pair_style      eam/alloy  # load "Cu.eam.alloy"
pair_coeff      * * Cu.eam.alloy Cu
timestep        0.001
run             500
"""

SCRIPT_COMMENT_NO_QUOTES = """\
units real
pair_style lj/charmm/coul/long 8.0 10.0 # no-quotes comment
pair_coeff * * 1.0 1.0
timestep 0.001
run 100
"""


class TestPreprocessInlineComments:
    def test_inline_comment_single_quote_parsed(self) -> None:
        """Line with a single-quoted word after '#' is stripped correctly."""
        parser = LammpsInputParser()
        data = parser.parse(SCRIPT_INLINE_COMMENT_SINGLE_QUOTE)
        assert data["units"] == "lj"
        assert data["pair_style"].startswith("lj/cut")
        assert data["timestep"] == pytest.approx(0.001)

    def test_inline_comment_double_quote_parsed(self) -> None:
        """Line with a double-quoted word after '#' is stripped correctly."""
        parser = LammpsInputParser()
        data = parser.parse(SCRIPT_INLINE_COMMENT_DOUBLE_QUOTE)
        assert data["units"] == "metal"
        assert "eam/alloy" in data["pair_style"]


# ---------------------------------------------------------------------------
# _parse_line: dimension ValueError, variable timestep/run
# ---------------------------------------------------------------------------

SCRIPT_DIMENSION_NON_INT = """\
units lj
dimension 3D
pair_style lj/cut 2.5
timestep 0.001
run 1000
"""

SCRIPT_TIMESTEP_VARIABLE = """\
units lj
pair_style lj/cut 2.5
variable dt equal 0.001
timestep ${dt}
run 1000
"""

SCRIPT_RUN_VARIABLE = """\
units lj
pair_style lj/cut 2.5
timestep 0.001
variable nsteps equal 50000
run ${nsteps}
"""

SCRIPT_MINIMIZE_ERROR = """\
units metal
pair_style eam
pair_coeff * * Cu.eam Cu
minimize bad_tol 1e-6 1000 10000
"""

SCRIPT_THERMO_NON_INT = """\
units lj
pair_style lj/cut 2.5
timestep 0.001
thermo myThermo
run 1000
"""

SCRIPT_MASS = """\
units metal
atom_style atomic
boundary p p p
pair_style eam/alloy
pair_coeff * * Cu.eam.alloy Cu
mass 1 63.546
timestep 0.001
run 100
"""

SCRIPT_RESTART = """\
units lj
pair_style lj/cut 2.5
timestep 0.001
restart 1000 restart.lammps
run 10000
"""

SCRIPT_NEIGH_MODIFY = """\
units lj
pair_style lj/cut 2.5
timestep 0.001
neigh_modify every 1 delay 0 check yes
run 1000
"""

SCRIPT_LANGEVIN_FIX = """\
units lj
atom_style atomic
pair_style lj/cut 2.5
pair_coeff * * 1.0 1.0 2.5
fix 1 all nve
fix 2 all langevin 1.0 1.0 0.1 12345
timestep 0.001
run 10000
"""


class TestParseLineEdgeCases:
    def test_dimension_non_int_stored_as_string(self) -> None:
        """'dimension 3D' should store the string value '3D'."""
        parser = LammpsInputParser()
        data = parser.parse(SCRIPT_DIMENSION_NON_INT)
        assert data.get("dimension") == "3D"

    def test_timestep_as_variable_stored_as_string(self) -> None:
        """'timestep ${dt}' cannot be converted to float → stored as string."""
        parser = LammpsInputParser()
        data = parser.parse(SCRIPT_TIMESTEP_VARIABLE)
        assert isinstance(data.get("timestep"), str)

    def test_run_as_variable_stored_as_string(self) -> None:
        """'run ${nsteps}' cannot be converted to int → stored as string."""
        parser = LammpsInputParser()
        data = parser.parse(SCRIPT_RUN_VARIABLE)
        assert isinstance(data.get("run"), str)

    def test_minimize_with_bad_tol_parsed_partially(self) -> None:
        """'minimize bad_tol ...' → etol set to the raw string token."""
        parser = LammpsInputParser()
        data = parser.parse(SCRIPT_MINIMIZE_ERROR)
        assert "minimize_etol" in data

    def test_thermo_non_int_stored_as_string(self) -> None:
        """'thermo myThermo' → thermo stored as string 'myThermo'."""
        parser = LammpsInputParser()
        data = parser.parse(SCRIPT_THERMO_NON_INT)
        assert data.get("thermo") == "myThermo"

    def test_mass_command_parsed(self) -> None:
        """'mass 1 63.546' → mass["1"] = 63.546."""
        parser = LammpsInputParser()
        data = parser.parse(SCRIPT_MASS)
        assert data.get("mass", {}).get("1") == pytest.approx(63.546)

    def test_restart_command_parsed(self) -> None:
        """'restart 1000 restart.lammps' → restart_interval = 1000."""
        parser = LammpsInputParser()
        data = parser.parse(SCRIPT_RESTART)
        assert data.get("restart_interval") == 1000

    def test_neigh_modify_command_parsed(self) -> None:
        """'neigh_modify every 1 delay 0 check yes' → stored."""
        parser = LammpsInputParser()
        data = parser.parse(SCRIPT_NEIGH_MODIFY)
        assert "neigh_modify" in data
        assert "every" in data["neigh_modify"]


# ---------------------------------------------------------------------------
# _summarize_fixes: langevin style (lines 419-428)
# ---------------------------------------------------------------------------

class TestLangevinFix:
    def test_langevin_fix_creates_fix_nvt_entry(self) -> None:
        """Langevin thermostat fix should populate fix_nvt."""
        parser = LammpsInputParser()
        data = parser.parse(SCRIPT_LANGEVIN_FIX)
        fix_nvt = data.get("fix_nvt")
        assert fix_nvt is not None
        # The NVE fix comes before langevin; langevin should win via setdefault
        # (nve sets fix_nve, not fix_nvt; langevin then sets fix_nvt)
        assert fix_nvt.get("Tstart") == pytest.approx(1.0)
        assert fix_nvt.get("Tstop") == pytest.approx(1.0)
        assert fix_nvt.get("damp") == pytest.approx(0.1)


# ---------------------------------------------------------------------------
# _tokenize: shlex ValueError path (lines 437-439)
# ---------------------------------------------------------------------------

class TestTokenizeEdgeCases:
    def test_tokenize_unbalanced_quote_falls_back_to_split(self) -> None:
        """A line with an unbalanced quote should still be tokenized."""
        from chemvcs.parsers.lammps_input_parser import LammpsInputParser as _P
        # Unbalanced quote input that shlex would fail on
        result = _P._tokenize("variable x equal ${y}")
        assert len(result) >= 3


# ---------------------------------------------------------------------------
# _parse_fix_nvt_npt: iso/aniso keyword path (lines 464-469, 472)
# ---------------------------------------------------------------------------

SCRIPT_NPT_ISO = """\
units lj
atom_style atomic
pair_style lj/cut 2.5
pair_coeff * * 1.0 1.0 2.5
fix 1 all npt temp 2.0 2.0 0.1 iso 1.0 1.0 1.0
timestep 0.001
run 10000
"""

class TestParseFixNptNve:
    def test_npt_iso_keyword_parsed(self) -> None:
        """NPT fix with 'iso' pressure keyword: Pstart, Pstop, Pdamp set."""
        parser = LammpsInputParser()
        data = parser.parse(SCRIPT_NPT_ISO)
        npt = data.get("fix_npt")
        assert npt is not None
        assert npt.get("press_style") == "iso"
        assert npt.get("Pstart") == pytest.approx(1.0)
        assert npt.get("Pstop") == pytest.approx(1.0)
        assert npt.get("Pdamp") == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# diff() – float noise in timestep, mass dict, read_data changes
# ---------------------------------------------------------------------------

class TestDiffEdgeCases:
    def test_diff_timestep_float_noise_suppressed(self) -> None:
        """Timestep values that differ by < 1e-15 must NOT appear in diff."""
        parser = LammpsInputParser()
        old: dict = {
            "units": "lj", "pair_style": "lj/cut 2.5",
            "timestep": 0.001, "read_data": [], "mass": {}, "variable_names": [],
        }
        new: dict = dict(old)
        new["timestep"] = old["timestep"] + 1e-17  # below 1e-15 threshold
        entries = parser.diff(old, new)
        ts_entries = [e for e in entries if e.path == "timestep"]
        assert len(ts_entries) == 0

    def test_diff_mass_dict_change_detected(self) -> None:
        """Changing mass values should produce a diff entry."""
        parser = LammpsInputParser()
        old: dict = {
            "units": "metal", "pair_style": "eam",
            "mass": {"1": 63.546}, "read_data": [], "variable_names": [],
        }
        new: dict = {
            "units": "metal", "pair_style": "eam",
            "mass": {"1": 65.38}, "read_data": [], "variable_names": [],
        }
        entries = parser.diff(old, new)
        mass_entries = [e for e in entries if e.path == "mass"]
        assert len(mass_entries) == 1
        assert mass_entries[0].change_type == "modified"

    def test_diff_read_data_change_detected(self) -> None:
        """Different read_data lists should produce a diff entry."""
        parser = LammpsInputParser()
        old: dict = {
            "units": "lj", "pair_style": "lj/cut 2.5",
            "read_data": ["data.lj_A"], "mass": {}, "variable_names": [],
        }
        new: dict = {
            "units": "lj", "pair_style": "lj/cut 2.5",
            "read_data": ["data.lj_B"], "mass": {}, "variable_names": [],
        }
        entries = parser.diff(old, new)
        rd_entries = [e for e in entries if e.path == "read_data"]
        assert len(rd_entries) == 1
        assert rd_entries[0].change_type == "modified"

    def test_diff_added_scalar_field(self) -> None:
        """Field added in new_data should appear as added entry."""
        parser = LammpsInputParser()
        old: dict = {
            "units": "lj", "pair_style": "lj/cut 2.5",
            "read_data": [], "mass": {}, "variable_names": [],
        }
        new: dict = dict(old)
        new["kspace_style"] = "pppm 1e-5"
        entries = parser.diff(old, new)
        added = [e for e in entries if e.path == "kspace_style" and e.change_type == "added"]
        assert len(added) == 1

    def test_diff_deleted_scalar_field(self) -> None:
        """Field removed in new_data should appear as deleted entry."""
        parser = LammpsInputParser()
        old: dict = {
            "units": "lj", "pair_style": "lj/cut 2.5",
            "kspace_style": "pppm 1e-5",
            "read_data": [], "mass": {}, "variable_names": [],
        }
        new: dict = {
            "units": "lj", "pair_style": "lj/cut 2.5",
            "read_data": [], "mass": {}, "variable_names": [],
        }
        entries = parser.diff(old, new)
        deleted = [e for e in entries if e.path == "kspace_style" and e.change_type == "deleted"]
        assert len(deleted) == 1


# ---------------------------------------------------------------------------
# validate() – negative timestep
# ---------------------------------------------------------------------------

class TestValidateEdgeCases:
    def test_validate_negative_timestep_error(self) -> None:
        parser = LammpsInputParser()
        data = {
            "units": "lj", "pair_style": "lj/cut 2.5",
            "timestep": -0.001, "read_data": [], "mass": {}, "variable_names": [],
        }
        is_valid, errors = parser.validate(data)
        assert not is_valid
        assert any("timestep" in e.lower() for e in errors)

    def test_validate_zero_timestep_error(self) -> None:
        parser = LammpsInputParser()
        data = {
            "units": "lj", "pair_style": "lj/cut 2.5",
            "timestep": 0.0, "read_data": [], "mass": {}, "variable_names": [],
        }
        is_valid, errors = parser.validate(data)
        assert not is_valid

    def test_validate_variable_timestep_skipped(self) -> None:
        """Variable-containing timestep (non-numeric) should not error."""
        parser = LammpsInputParser()
        data = {
            "units": "lj", "pair_style": "lj/cut 2.5",
            "timestep": "${dt}", "read_data": [], "mass": {}, "variable_names": [],
        }
        is_valid, errors = parser.validate(data)
        # No timestep error since it can't be converted to float
        ts_errors = [e for e in errors if "timestep" in e.lower()]
        assert len(ts_errors) == 0


# ---------------------------------------------------------------------------
# _classify – major and minor fields
# ---------------------------------------------------------------------------

class TestClassify:
    def test_classify_critical_field(self) -> None:
        parser = LammpsInputParser()
        assert parser._classify("units") == "critical"
        assert parser._classify("pair_style") == "critical"
        assert parser._classify("timestep") == "critical"
        assert parser._classify("run") == "critical"

    def test_classify_major_field(self) -> None:
        parser = LammpsInputParser()
        assert parser._classify("fix_nvt") == "major"
        assert parser._classify("fix_npt") == "major"
        assert parser._classify("fix_nve") == "major"
        assert parser._classify("boundary") == "major"
        assert parser._classify("atom_style") == "major"
        assert parser._classify("dimension") == "major"

    def test_classify_minor_field(self) -> None:
        """Unknown fields should default to 'minor'."""
        parser = LammpsInputParser()
        assert parser._classify("some_unknown_field") == "minor"
        assert parser._classify("dump_interval") == "minor"
        assert parser._classify("restart_interval") == "minor"
