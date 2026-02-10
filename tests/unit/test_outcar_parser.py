"""Tests for OUTCAR parser."""

import pytest

from chemvcs.parsers.outcar_parser import OutcarParser
from chemvcs.parsers.base_parser import ParserError


# A minimal but realistic OUTCAR excerpt for testing
SAMPLE_OUTCAR = """\
 running on    4 total cores
 distrk:  each k-point on    4 cores,    1 groups
 vasp.6.3.2 03Jul22

 POSCAR found type information on POSCAR: Si
 POSCAR found :  1 types and       2 ions

 E-fermi :   5.8921     XC(G=0):  -4.8951     alpha+bet : -3.8128

 Free energy of the ion-electron system (eV)
   ---------------------------------------------------
   alpha Z        PSCENC =         3.37178130
   Ewald energy   TEWEN  =      -255.20844012
   -1/2 Hartree   DENC   =       -35.06795802
   exchange EXHF  =         0.00000000
   -V(xc)+E(xc)   XCENC  =        38.04879640
   PAW double counting   =       247.45653230     -248.77183524
   entropy T*S    ECONTRIB =       -0.00026891
   eigenvalues    EBANDS =        -8.26149268
   atomic energy  EATOM  =       215.60634181
   ---------------------------------------------------
   free  energy   TOTEN  =       -10.82654316 eV

   energy  without entropy=     -10.82654312  energy(sigma->0) =     -10.82654312

 General timing and accounting informations for this job:

  Total CPU time used (sec):       42.358
  User time (sec):                 38.921
  System time (sec):                3.437
  Elapsed time (sec):              45.123

  Maximum memory used (kb):      125408
"""

SAMPLE_OUTCAR_V2 = """\
 running on    8 total cores
 vasp.6.3.2 03Jul22

 E-fermi :   5.9100     XC(G=0):  -4.9000     alpha+bet : -3.8200

 Free energy of the ion-electron system (eV)
   ---------------------------------------------------
   free  energy   TOTEN  =       -10.95000000 eV

   energy  without entropy=     -10.95000000  energy(sigma->0) =     -10.95000000

 General timing and accounting informations for this job:

  Total CPU time used (sec):       85.200
  Elapsed time (sec):              90.500

  Maximum memory used (kb):      256000
"""


class TestOutcarParser:
    """Test OUTCAR file parser."""

    def test_parse_extracts_energy(self) -> None:
        """Test that final energy is extracted."""
        parser = OutcarParser()
        data = parser.parse(SAMPLE_OUTCAR)

        assert data["final_energy"] == pytest.approx(-10.82654316)

    def test_parse_extracts_efermi(self) -> None:
        """Test that Fermi energy is extracted."""
        parser = OutcarParser()
        data = parser.parse(SAMPLE_OUTCAR)

        assert data["efermi"] == pytest.approx(5.8921)

    def test_parse_extracts_run_stats(self) -> None:
        """Test that run statistics are extracted."""
        parser = OutcarParser()
        data = parser.parse(SAMPLE_OUTCAR)

        assert "run_stats" in data
        run_stats = data["run_stats"]
        assert run_stats.get("Total CPU time used (sec)") == pytest.approx(42.358)
        assert run_stats.get("Elapsed time (sec)") == pytest.approx(45.123)

    def test_parse_extracts_memory(self) -> None:
        """Test that memory usage is extracted."""
        parser = OutcarParser()
        data = parser.parse(SAMPLE_OUTCAR)

        run_stats = data["run_stats"]
        assert run_stats.get("Maximum memory used (kb)") == 125408

    def test_parse_has_required_fields(self) -> None:
        """Test that all required fields are present."""
        parser = OutcarParser()
        data = parser.parse(SAMPLE_OUTCAR)

        for field in ["final_energy", "magnetization", "charge",
                       "run_stats", "is_stopped"]:
            assert field in data, f"Missing field: {field}"

    def test_parse_empty_content_fails(self) -> None:
        """Test that empty content raises ParserError."""
        parser = OutcarParser()

        with pytest.raises(ParserError):
            parser.parse("")

    def test_parse_garbage_content_fails(self) -> None:
        """Test that non-OUTCAR content raises ParserError."""
        parser = OutcarParser()

        with pytest.raises(ParserError):
            parser.parse("This is not an OUTCAR file at all.\nRandom text.")

    def test_diff_energy_change(self) -> None:
        """Test diff detects energy change."""
        parser = OutcarParser()
        old = parser.parse(SAMPLE_OUTCAR)
        new = parser.parse(SAMPLE_OUTCAR_V2)

        diff = parser.diff(old, new)

        energy_diffs = [d for d in diff if d.path == "final_energy"]
        assert len(energy_diffs) == 1
        assert energy_diffs[0].change_type == "modified"
        assert energy_diffs[0].significance == "critical"
        assert energy_diffs[0].old_value == pytest.approx(-10.82654316)
        assert energy_diffs[0].new_value == pytest.approx(-10.95)

    def test_diff_efermi_change(self) -> None:
        """Test diff detects Fermi energy change."""
        parser = OutcarParser()
        old = parser.parse(SAMPLE_OUTCAR)
        new = parser.parse(SAMPLE_OUTCAR_V2)

        diff = parser.diff(old, new)

        fermi_diffs = [d for d in diff if d.path == "efermi"]
        assert len(fermi_diffs) == 1
        assert fermi_diffs[0].significance == "critical"

    def test_diff_no_changes(self) -> None:
        """Test diff with identical data returns empty list."""
        parser = OutcarParser()
        data = parser.parse(SAMPLE_OUTCAR)

        diff = parser.diff(data, data)
        assert len(diff) == 0

    def test_diff_run_stats_change(self) -> None:
        """Test diff detects run statistics change."""
        parser = OutcarParser()
        old = parser.parse(SAMPLE_OUTCAR)
        new = parser.parse(SAMPLE_OUTCAR_V2)

        diff = parser.diff(old, new)

        stats_diffs = [d for d in diff if d.path == "run_stats"]
        assert len(stats_diffs) == 1
        assert stats_diffs[0].significance == "minor"

    def test_validate_normal_run(self) -> None:
        """Test validation of a normal completed run."""
        parser = OutcarParser()
        data = parser.parse(SAMPLE_OUTCAR)

        is_valid, errors = parser.validate(data)
        assert is_valid

    def test_validate_stopped_run(self) -> None:
        """Test validation warns about stopped run."""
        parser = OutcarParser()
        data = {"is_stopped": True, "final_energy": -10.0}

        is_valid, errors = parser.validate(data)
        assert not is_valid
        assert any("stopped" in err.lower() for err in errors)

    def test_validate_no_energy(self) -> None:
        """Test validation warns about missing energy."""
        parser = OutcarParser()
        data = {"is_stopped": False, "final_energy": None}

        is_valid, errors = parser.validate(data)
        assert not is_valid
        assert any("energy" in err.lower() for err in errors)

    def test_format_diff(self) -> None:
        """Test diff formatting output."""
        parser = OutcarParser()
        old = parser.parse(SAMPLE_OUTCAR)
        new = parser.parse(SAMPLE_OUTCAR_V2)

        diff = parser.diff(old, new)
        formatted = parser.format_diff(diff)

        assert "final_energy" in formatted
        assert "efermi" in formatted
