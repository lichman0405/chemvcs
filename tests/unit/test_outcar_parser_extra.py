"""Extra coverage tests for outcar_parser.py.

Targets uncovered lines:
    75-76  – except for regex fallback (needs regex to raise)
   104    – diff: added field branch
   108    – diff: deleted field branch
   152    – _parse_with_pymatgen: return _extract_pymatgen(outcar)
  156-157  – _parse_with_pymatgen: OSError in finally os.unlink
  162-184  – _extract_pymatgen body
   208    – _parse_essential_fields returns {} (returns None / no data)
"""

from unittest.mock import MagicMock, patch

import pytest

from chemvcs.parsers.base_parser import ParserError
from chemvcs.parsers.outcar_parser import OutcarParser

# ---------------------------------------------------------------------------
# validate()
# ---------------------------------------------------------------------------


class TestOutcarValidate:
    def test_validate_stopped_run_gives_error(self) -> None:
        parser = OutcarParser()
        data = {"is_stopped": True, "final_energy": -10.5}
        is_valid, errors = parser.validate(data)
        assert not is_valid
        assert any("stopped" in e.lower() for e in errors)

    def test_validate_no_energy_gives_error(self) -> None:
        parser = OutcarParser()
        data = {"is_stopped": False, "final_energy": None}
        is_valid, errors = parser.validate(data)
        assert not is_valid
        assert any("energy" in e.lower() for e in errors)

    def test_validate_stopped_and_no_energy_gives_two_errors(self) -> None:
        parser = OutcarParser()
        data = {"is_stopped": True, "final_energy": None}
        is_valid, errors = parser.validate(data)
        assert not is_valid
        assert len(errors) == 2

    def test_validate_clean_data_passes(self) -> None:
        parser = OutcarParser()
        data = {"is_stopped": False, "final_energy": -10.5}
        is_valid, errors = parser.validate(data)
        assert is_valid
        assert errors == []


# ---------------------------------------------------------------------------
# diff() – added / deleted / noise-filtered branches
# ---------------------------------------------------------------------------


class TestOutcarDiff:
    def test_diff_added_field(self) -> None:
        """Field appears only in new_data → added entry."""
        parser = OutcarParser()
        old: dict = {"final_energy": -10.0}
        new: dict = {"final_energy": -10.0, "total_magnetization": 0.5}
        entries = parser.diff(old, new)
        added = [e for e in entries if e.path == "total_magnetization"]
        assert len(added) == 1
        assert added[0].change_type == "added"
        assert added[0].new_value == 0.5

    def test_diff_deleted_field(self) -> None:
        """Field appears only in old_data → deleted entry."""
        parser = OutcarParser()
        old: dict = {"final_energy": -10.0, "total_magnetization": 1.5}
        new: dict = {"final_energy": -10.0}
        entries = parser.diff(old, new)
        deleted = [e for e in entries if e.path == "total_magnetization"]
        assert len(deleted) == 1
        assert deleted[0].change_type == "deleted"
        assert deleted[0].old_value == 1.5

    def test_diff_energy_noise_filtered_out(self) -> None:
        """Energy differences smaller than 1e-8 must be suppressed."""
        parser = OutcarParser()
        base_energy = -10.82654316
        old: dict = {"final_energy": base_energy}
        new: dict = {"final_energy": base_energy + 1e-10}  # below threshold
        entries = parser.diff(old, new)
        energy_entries = [e for e in entries if e.path == "final_energy"]
        assert len(energy_entries) == 0

    def test_diff_energy_above_noise_threshold_kept(self) -> None:
        """Energy differences >= 1e-8 must NOT be suppressed."""
        parser = OutcarParser()
        old: dict = {"final_energy": -10.0}
        new: dict = {"final_energy": -10.5}  # well above 1e-8
        entries = parser.diff(old, new)
        energy_entries = [e for e in entries if e.path == "final_energy"]
        assert len(energy_entries) == 1

    def test_diff_efermi_modification(self) -> None:
        parser = OutcarParser()
        old: dict = {"final_energy": -10.0, "efermi": 5.8}
        new: dict = {"final_energy": -10.0, "efermi": 6.1}
        entries = parser.diff(old, new)
        efermi = [e for e in entries if e.path == "efermi"]
        assert len(efermi) == 1
        assert efermi[0].change_type == "modified"


# ---------------------------------------------------------------------------
# _extract_pymatgen – called directly with a mock Outcar
# ---------------------------------------------------------------------------


class TestExtractPymatgen:
    """Call the static method directly so lines 162-184 are executed."""

    def test_basic_extraction(self) -> None:
        mock_outcar = MagicMock()
        mock_outcar.final_energy = -12.345
        mock_outcar.efermi = 5.9
        mock_outcar.magnetization = None
        mock_outcar.total_mag = 0.0
        mock_outcar.charge = None
        mock_outcar.run_stats = {"Total CPU time used (sec)": 42.0}
        mock_outcar.is_stopped = False

        data = OutcarParser._extract_pymatgen(mock_outcar)

        assert data["final_energy"] == -12.345
        assert data["efermi"] == 5.9
        assert data["is_stopped"] is False
        assert data["magnetization"] == []
        assert data["charge"] == []
        assert data["run_stats"]["Total CPU time used (sec)"] == 42.0

    def test_with_nonempty_magnetization(self) -> None:
        mock_outcar = MagicMock()
        mock_outcar.final_energy = -10.0
        mock_outcar.efermi = 5.5
        mock_outcar.magnetization = [{"s": 0.1, "p": 0.0, "d": 0.5, "tot": 0.6}]
        mock_outcar.total_mag = 0.6
        mock_outcar.charge = None
        mock_outcar.run_stats = {}
        mock_outcar.is_stopped = False

        data = OutcarParser._extract_pymatgen(mock_outcar)

        assert len(data["magnetization"]) == 1
        assert data["total_magnetization"] == 0.6

    def test_with_nonempty_charge(self) -> None:
        mock_outcar = MagicMock()
        mock_outcar.final_energy = -10.0
        mock_outcar.efermi = 5.5
        mock_outcar.magnetization = None
        mock_outcar.total_mag = None
        mock_outcar.charge = [{"s": 2.0, "p": 0.0, "tot": 2.0}]
        mock_outcar.run_stats = {}
        mock_outcar.is_stopped = False

        data = OutcarParser._extract_pymatgen(mock_outcar)

        assert len(data["charge"]) == 1

    def test_run_stats_none_becomes_empty_dict(self) -> None:
        """run_stats=None should be coerced to {}."""
        mock_outcar = MagicMock()
        mock_outcar.final_energy = -10.0
        mock_outcar.efermi = None
        mock_outcar.magnetization = None
        mock_outcar.total_mag = None
        mock_outcar.charge = None
        mock_outcar.run_stats = None
        mock_outcar.is_stopped = False

        data = OutcarParser._extract_pymatgen(mock_outcar)
        # None → {} due to `getattr(..., {}) or {}`
        assert data["run_stats"] == {}


# ---------------------------------------------------------------------------
# _parse_with_pymatgen – success path via mock PmgOutcar
# ---------------------------------------------------------------------------


class TestParseWithPymatgen:
    """Patch PmgOutcar so pymatgen 'succeeds', covering lines 146-157."""

    def test_parse_with_pymatgen_success_path(self) -> None:
        mock_outcar = MagicMock()
        mock_outcar.final_energy = -12.0
        mock_outcar.efermi = 5.5
        mock_outcar.magnetization = None
        mock_outcar.total_mag = 0.0
        mock_outcar.charge = None
        mock_outcar.run_stats = {}
        mock_outcar.is_stopped = False

        with patch("chemvcs.parsers.outcar_parser.PmgOutcar", return_value=mock_outcar):
            parser = OutcarParser()
            data = parser._parse_with_pymatgen("irrelevant content")

        assert data["final_energy"] == -12.0
        assert data["is_stopped"] is False

    def test_parse_with_pymatgen_unlink_oserror_silently_caught(self) -> None:
        """If os.unlink raises OSError in the finally block, it is silently
        swallowed (lines 156-157).  The parse result is unaffected."""
        mock_outcar = MagicMock()
        mock_outcar.final_energy = -7.0
        mock_outcar.efermi = 4.0
        mock_outcar.magnetization = None
        mock_outcar.total_mag = None
        mock_outcar.charge = None
        mock_outcar.run_stats = {}
        mock_outcar.is_stopped = False

        with (
            patch("chemvcs.parsers.outcar_parser.PmgOutcar", return_value=mock_outcar),
            patch("os.unlink", side_effect=OSError("simulated unlink error")),
        ):
            parser = OutcarParser()
            data = parser._parse_with_pymatgen("some content")

        assert data["final_energy"] == -7.0


# ---------------------------------------------------------------------------
# parse() – regex fallback (uses content with TOTEN pattern only)
# ---------------------------------------------------------------------------


class TestOutcarParseRegexFallback:
    def test_parse_regex_fallback_extracts_energy(self) -> None:
        """Content where pymatgen fails but regex succeeds."""
        content = " free  energy   TOTEN  =       -15.12345678 eV\n E-fermi :   6.1234\n"
        parser = OutcarParser()
        data = parser.parse(content)
        assert data.get("final_energy") == pytest.approx(-15.12345678)
        assert data.get("efermi") == pytest.approx(6.1234)

    def test_parse_regex_fallback_exception_raises_parser_error(self) -> None:
        """If _parse_essential_fields itself raises, ParserError is propagated."""
        parser = OutcarParser()
        with patch.object(
            OutcarParser,
            "_parse_essential_fields",
            side_effect=RuntimeError("boom"),
        ), pytest.raises(ParserError):
            parser.parse("some irrelevant content without patterns")
