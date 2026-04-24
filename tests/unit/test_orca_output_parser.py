"""Tests for OrcaOutputParser.

All tests rely on synthetic ORCA-like output strings that reproduce the
key marker lines.  No real ORCA installation is required.

cclib is called only when the content looks like valid ORCA output. For
unit tests we rely on the regex fallback path (synthetic snippets are too
short for cclib to fully parse).  The cclib primary path is exercised via
an integration fixture that builds a minimal but realistic header.
"""

from __future__ import annotations

import pytest

from chemvcs.parsers.base_parser import DiffEntry, ParserError
from chemvcs.parsers.orca_output_parser import OrcaOutputParser, _extract_wall_time

# ---------------------------------------------------------------------------
# Helpers – synthetic ORCA output fragments
# ---------------------------------------------------------------------------


def _make_orca_output(
    *,
    energy: float = -76.123456789,
    terminated: bool = True,
    scf_cycles: int = 12,
    scf_converged: bool = True,
    opt_converged: bool | None = None,
    charge: int = 0,
    mult: int = 1,
    version: str = "5.0.3",
    wall: str = "0 days 0 hours 1 minutes 23.456 seconds",
) -> str:
    lines = [
        f"       Program Version {version}",
        "",
        f"   Total Charge           Charge          ....  {charge:+d}",
        f"   Multiplicity           Mult            ....  {mult}",
        "",
    ]
    if scf_converged:
        lines.append(f"                  SUCCESS! SCF CONVERGED AFTER {scf_cycles} CYCLES")
    else:
        lines.append("                  WARNING: SCF NOT CONVERGED")

    lines.append(f"FINAL SINGLE POINT ENERGY        {energy:.10f}")

    if opt_converged is True:
        lines.append("              THE OPTIMIZATION HAS CONVERGED")
    elif opt_converged is False:
        lines.append("              THE OPTIMIZATION HAS NOT CONVERGED")

    if terminated:
        lines.append("                        ****ORCA TERMINATED NORMALLY****")
    else:
        lines.append("                        ERROR TERMINATION")

    lines.append(f"TOTAL RUN TIME: {wall}")
    return "\n".join(lines)


NORMAL_OUTPUT = _make_orca_output()
ERROR_OUTPUT = _make_orca_output(terminated=False)
OPT_CONVERGED = _make_orca_output(opt_converged=True)
OPT_NOT_CONVERGED = _make_orca_output(opt_converged=False)
HIGHER_ENERGY_OUTPUT = _make_orca_output(energy=-75.0)
RADICAL_OUTPUT = _make_orca_output(charge=0, mult=2)


# ---------------------------------------------------------------------------
# Parse tests (regex fallback path)
# ---------------------------------------------------------------------------


class TestOrcaOutputParserParse:
    def setup_method(self) -> None:
        self.parser = OrcaOutputParser()

    def test_parse_final_energy(self) -> None:
        data = self.parser.parse(NORMAL_OUTPUT)
        assert data["final_energy"] == pytest.approx(-76.123456789, rel=1e-8)

    def test_parse_termination_normal(self) -> None:
        data = self.parser.parse(NORMAL_OUTPUT)
        assert data["termination_status"] == "normal"

    def test_parse_termination_error(self) -> None:
        data = self.parser.parse(ERROR_OUTPUT)
        assert data["termination_status"] == "error"

    def test_parse_scf_converged(self) -> None:
        data = self.parser.parse(NORMAL_OUTPUT)
        assert data["scf_converged"] is True

    def test_parse_scf_cycles(self) -> None:
        data = self.parser.parse(NORMAL_OUTPUT)
        assert data["scf_cycles"] == 12

    def test_parse_scf_not_converged(self) -> None:
        data = self.parser.parse(_make_orca_output(scf_converged=False))
        assert data["scf_converged"] is False

    def test_parse_opt_converged(self) -> None:
        data = self.parser.parse(OPT_CONVERGED)
        assert data["optimization_converged"] is True

    def test_parse_opt_not_converged(self) -> None:
        data = self.parser.parse(OPT_NOT_CONVERGED)
        assert data["optimization_converged"] is False

    def test_parse_opt_absent(self) -> None:
        data = self.parser.parse(NORMAL_OUTPUT)
        # SP calculation — no opt markers
        assert data["optimization_converged"] is None

    def test_parse_orca_version(self) -> None:
        data = self.parser.parse(NORMAL_OUTPUT)
        assert data["orca_version"] == "5.0.3"

    def test_parse_wall_time(self) -> None:
        data = self.parser.parse(NORMAL_OUTPUT)
        assert data["wall_time"] is not None
        assert "seconds" in data["wall_time"].lower()

    def test_parse_charge(self) -> None:
        data = self.parser.parse(NORMAL_OUTPUT)
        assert data["charge"] == 0

    def test_parse_multiplicity(self) -> None:
        data = self.parser.parse(RADICAL_OUTPUT)
        assert data["multiplicity"] == 2

    def test_multiple_energies_last_wins(self) -> None:
        content = (
            "FINAL SINGLE POINT ENERGY  -76.100000000\n"
            "FINAL SINGLE POINT ENERGY  -76.123456789\n"
            "                        ****ORCA TERMINATED NORMALLY****\n"
        )
        data = self.parser.parse(content)
        assert data["final_energy"] == pytest.approx(-76.123456789, rel=1e-8)

    def test_empty_content_raises(self) -> None:
        with pytest.raises(ParserError):
            self.parser.parse("")

    def test_whitespace_only_raises(self) -> None:
        with pytest.raises(ParserError):
            self.parser.parse("   \n   \t   ")

    def test_unrelated_text_raises(self) -> None:
        with pytest.raises(ParserError):
            self.parser.parse("Hello, world!")


# ---------------------------------------------------------------------------
# Diff tests
# ---------------------------------------------------------------------------


class TestOrcaOutputParserDiff:
    def setup_method(self) -> None:
        self.parser = OrcaOutputParser()

    def test_no_changes(self) -> None:
        data = self.parser.parse(NORMAL_OUTPUT)
        entries = self.parser.diff(data, data)
        assert entries == []

    def test_energy_change_is_critical(self) -> None:
        old = self.parser.parse(NORMAL_OUTPUT)
        new = self.parser.parse(HIGHER_ENERGY_OUTPUT)
        entries = self.parser.diff(old, new)
        e_entry = next(e for e in entries if e.path == "final_energy")
        assert e_entry.significance == "critical"
        assert e_entry.old_value == pytest.approx(-76.123456789, rel=1e-8)
        assert e_entry.new_value == pytest.approx(-75.0, rel=1e-8)

    def test_termination_change_is_critical(self) -> None:
        old = self.parser.parse(NORMAL_OUTPUT)
        new = self.parser.parse(ERROR_OUTPUT)
        entries = self.parser.diff(old, new)
        term_e = next(e for e in entries if e.path == "termination_status")
        assert term_e.significance == "critical"
        assert term_e.old_value == "normal"
        assert term_e.new_value == "error"

    def test_opt_convergence_change_is_critical(self) -> None:
        old = self.parser.parse(OPT_CONVERGED)
        new = self.parser.parse(OPT_NOT_CONVERGED)
        entries = self.parser.diff(old, new)
        opt_e = next(e for e in entries if e.path == "optimization_converged")
        assert opt_e.significance == "critical"

    def test_scf_cycles_change_is_major(self) -> None:
        old = self.parser.parse(_make_orca_output(scf_cycles=12))
        new = self.parser.parse(_make_orca_output(scf_cycles=50))
        entries = self.parser.diff(old, new)
        cyc_e = next(e for e in entries if e.path == "scf_cycles")
        assert cyc_e.significance == "major"

    def test_version_change_is_minor(self) -> None:
        old = self.parser.parse(_make_orca_output(version="5.0.3"))
        new = self.parser.parse(_make_orca_output(version="6.0.0"))
        entries = self.parser.diff(old, new)
        ver_e = next(e for e in entries if e.path == "orca_version")
        assert ver_e.significance == "minor"

    def test_wall_time_change_is_minor(self) -> None:
        old = self.parser.parse(_make_orca_output(wall="0 days 0 hours 0 minutes 30.000 seconds"))
        new = self.parser.parse(_make_orca_output(wall="0 days 0 hours 1 minutes 20.000 seconds"))
        entries = self.parser.diff(old, new)
        wt_e = next((e for e in entries if e.path == "wall_time"), None)
        if wt_e is not None:
            assert wt_e.significance == "minor"

    def test_float_noise_suppressed(self) -> None:
        """Differences smaller than 1e-10 Hartree should not appear."""
        old = {"final_energy": -76.123456789123456}
        new = {"final_energy": -76.123456789123456 + 1e-15}
        entries = self.parser.diff(old, new)
        assert not any(e.path == "final_energy" for e in entries)

    def test_diff_returns_diff_entry_objects(self) -> None:
        old = self.parser.parse(NORMAL_OUTPUT)
        new = self.parser.parse(HIGHER_ENERGY_OUTPUT)
        entries = self.parser.diff(old, new)
        for e in entries:
            assert isinstance(e, DiffEntry)

    def test_change_type_added(self) -> None:
        old = {"termination_status": None, "final_energy": -76.0}
        new = {"termination_status": "normal", "final_energy": -76.0}
        entries = self.parser.diff(old, new)
        added = next(e for e in entries if e.path == "termination_status")
        assert added.change_type == "added"

    def test_change_type_deleted(self) -> None:
        old = {"termination_status": "normal", "final_energy": -76.0}
        new = {"termination_status": None, "final_energy": -76.0}
        entries = self.parser.diff(old, new)
        deleted = next(e for e in entries if e.path == "termination_status")
        assert deleted.change_type == "deleted"


# ---------------------------------------------------------------------------
# Validate tests
# ---------------------------------------------------------------------------


class TestOrcaOutputParserValidate:
    def setup_method(self) -> None:
        self.parser = OrcaOutputParser()

    def test_valid_normal_output(self) -> None:
        data = self.parser.parse(NORMAL_OUTPUT)
        ok, errors = self.parser.validate(data)
        assert ok
        assert errors == []

    def test_invalid_error_termination(self) -> None:
        data = self.parser.parse(ERROR_OUTPUT)
        ok, errors = self.parser.validate(data)
        assert not ok
        assert any("normally" in e.lower() or "terminat" in e.lower() for e in errors)

    def test_invalid_opt_not_converged(self) -> None:
        data = self.parser.parse(OPT_NOT_CONVERGED)
        ok, errors = self.parser.validate(data)
        assert not ok
        assert any("optimis" in e.lower() or "optim" in e.lower() for e in errors)

    def test_missing_energy(self) -> None:
        data = {"termination_status": "normal", "final_energy": None}
        ok, errors = self.parser.validate(data)
        assert not ok
        assert any("energy" in e.lower() for e in errors)


# ---------------------------------------------------------------------------
# _extract_wall_time helper
# ---------------------------------------------------------------------------


class TestExtractWallTime:
    def test_standard_format(self) -> None:
        content = "TOTAL RUN TIME: 0 days 0 hours 5 minutes 23.456 seconds"
        assert _extract_wall_time(content) is not None
        assert "5 minutes" in _extract_wall_time(content)

    def test_no_wall_time(self) -> None:
        assert _extract_wall_time("some random content") is None

    def test_case_insensitive(self) -> None:
        content = "TOTAL RUN TIME: 1 days 2 hours 3 minutes 4.5 seconds"
        assert _extract_wall_time(content.lower()) is not None
