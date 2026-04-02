"""Tests for OrcaInputParser."""

import pytest

from chemvcs.parsers.orca_input_parser import OrcaInputParser, _try_numeric
from chemvcs.parsers.base_parser import ParserError, DiffEntry


# ---------------------------------------------------------------------------
# Fixtures / shared content strings
# ---------------------------------------------------------------------------

SIMPLE_OPT = """\
! B3LYP def2-TZVP Opt TightSCF D3BJ RIJCOSX Grid4

%pal nprocs 8 end
%maxcore 4000

* xyz 0 1
  C    0.000   0.000   0.000
  H    1.089   0.000   0.000
  H   -0.363   1.027   0.000
  H   -0.363  -0.513   0.890
  H   -0.363  -0.513  -0.890
*
"""

SIMPLE_SP = """\
! CCSD(T) def2-TZVPP TightSCF RIJCOSX

%maxcore 8000

* xyz 0 1
  O    0.000   0.000   0.117
  H    0.000   0.757  -0.469
  H    0.000  -0.757  -0.469
*
"""

WITH_SCF_BLOCK = """\
! HF 6-31G* SP

%scf
   MaxIter 200
   Convergence Tight
end

* xyz 1 2
  N    0.000   0.000   0.000
*
"""

WITH_BLOCKS = """\
! PBE0 def2-SVP CPCM(water) Opt

%cpcm
   epsilon 80.4
   refrac  1.33
end

%geom
   MaxIter 50
end

%pal nprocs 4 end
%maxcore 2000

* xyz -1 1
  Cl   0.000   0.000   0.000
*
"""

XYZFILE = """\
! B3LYP def2-SVP SP

* xyzfile 0 1 molecule.xyz
"""

WITH_COMMENTS = """\
# This is an optimisation job
! B3LYP def2-TZVP Opt  # method line

%pal nprocs 4 end  # parallelism

* xyz 0 3  # triplet
  O    0.000   0.000   0.000
  O    1.200   0.000   0.000
*
"""

MULTI_STEP = """\
! HF STO-3G SP

* xyz 0 1
  H    0.000   0.000   0.000
  H    0.740   0.000   0.000
*

$new_job
! CCSD(T) cc-pVTZ SP

* xyz 0 1
  H    0.000   0.000   0.000
  H    0.740   0.000   0.000
*
"""


# ---------------------------------------------------------------------------
# Parse tests
# ---------------------------------------------------------------------------

class TestOrcaInputParserParse:
    def setup_method(self) -> None:
        self.parser = OrcaInputParser()

    # ---- keywords ----

    def test_parse_keyword_line(self) -> None:
        data = self.parser.parse(SIMPLE_OPT)
        assert "B3LYP" in data["keywords"]
        assert "DEF2-TZVP" in data["keywords"]
        assert "OPT" in data["keywords"]
        assert "TIGHTSCF" in data["keywords"]
        assert "D3BJ" in data["keywords"]
        assert "RIJCOSX" in data["keywords"]
        assert "GRID4" in data["keywords"]

    def test_keywords_are_upper_case(self) -> None:
        data = self.parser.parse(SIMPLE_OPT)
        for kw in data["keywords"]:
            assert kw == kw.upper(), f"keyword {kw!r} is not upper-case"

    def test_keywords_are_sorted(self) -> None:
        data = self.parser.parse(SIMPLE_OPT)
        assert data["keywords"] == sorted(data["keywords"])

    # ---- charge / multiplicity / coord format ----

    def test_parse_charge_and_mult(self) -> None:
        data = self.parser.parse(SIMPLE_OPT)
        assert data["charge"] == 0
        assert data["multiplicity"] == 1

    def test_parse_negative_charge(self) -> None:
        data = self.parser.parse(WITH_BLOCKS)
        assert data["charge"] == -1

    def test_parse_radical(self) -> None:
        data = self.parser.parse(WITH_SCF_BLOCK)
        assert data["charge"] == 1
        assert data["multiplicity"] == 2

    def test_parse_coord_format_xyz(self) -> None:
        data = self.parser.parse(SIMPLE_OPT)
        assert data["coord_format"] == "xyz"

    def test_parse_coord_format_xyzfile(self) -> None:
        data = self.parser.parse(XYZFILE)
        assert data["coord_format"] == "xyzfile"
        assert data["xyzfile"] == "molecule.xyz"
        assert data["charge"] == 0
        assert data["multiplicity"] == 1

    def test_parse_atoms(self) -> None:
        data = self.parser.parse(SIMPLE_OPT)
        assert len(data["atoms"]) == 5
        assert data["atoms"][0].startswith("C")

    # ---- parallel / memory ----

    def test_parse_maxcore(self) -> None:
        data = self.parser.parse(SIMPLE_OPT)
        assert data["maxcore"] == 4000

    def test_parse_nprocs_inline(self) -> None:
        data = self.parser.parse(SIMPLE_OPT)
        assert data["nprocs"] == 8

    def test_parse_nprocs_from_pal_block(self) -> None:
        content = """\
! HF STO-3G SP

%pal
   nprocs 16
end

* xyz 0 1
  H    0.000   0.000   0.000
  H    0.740   0.000   0.000
*
"""
        data = self.parser.parse(content)
        assert data["nprocs"] == 16

    # ---- blocks ----

    def test_parse_scf_block(self) -> None:
        data = self.parser.parse(WITH_SCF_BLOCK)
        assert "scf" in data["blocks"]
        scf = data["blocks"]["scf"]
        assert scf["maxiter"] == 200
        assert scf["convergence"] == "Tight"

    def test_parse_cpcm_block(self) -> None:
        data = self.parser.parse(WITH_BLOCKS)
        assert "cpcm" in data["blocks"]
        cpcm = data["blocks"]["cpcm"]
        assert cpcm["epsilon"] == pytest.approx(80.4)

    def test_parse_geom_block(self) -> None:
        data = self.parser.parse(WITH_BLOCKS)
        assert "geom" in data["blocks"]
        assert data["blocks"]["geom"]["maxiter"] == 50

    # ---- comments ----

    def test_strip_hash_comments(self) -> None:
        data = self.parser.parse(WITH_COMMENTS)
        # Inline comment on ! line should not appear as a keyword
        assert "This" not in data["keywords"]
        assert "B3LYP" in data["keywords"]

    def test_triplet_via_comments(self) -> None:
        data = self.parser.parse(WITH_COMMENTS)
        assert data["multiplicity"] == 3

    # ---- multi-step ----

    def test_multistep_parses_first_job_only(self) -> None:
        data = self.parser.parse(MULTI_STEP)
        # First job uses HF STO-3G; second is CCSD(T) cc-pVTZ
        assert "HF" in data["keywords"]
        assert "STO-3G" in data["keywords"]
        assert "CCSD(T)" not in data["keywords"]

    # ---- edge cases ----

    def test_empty_content_raises(self) -> None:
        with pytest.raises(ParserError):
            self.parser.parse("")

    def test_whitespace_only_raises(self) -> None:
        with pytest.raises(ParserError):
            self.parser.parse("   \n  \t  ")


# ---------------------------------------------------------------------------
# Diff tests
# ---------------------------------------------------------------------------

class TestOrcaInputParserDiff:
    def setup_method(self) -> None:
        self.parser = OrcaInputParser()

    def test_no_changes(self) -> None:
        data = self.parser.parse(SIMPLE_OPT)
        entries = self.parser.diff(data, data)
        assert entries == []

    def test_method_change_is_critical(self) -> None:
        old = self.parser.parse(SIMPLE_OPT)          # B3LYP
        new_content = SIMPLE_OPT.replace("B3LYP", "PBE0")
        new = self.parser.parse(new_content)
        entries = self.parser.diff(old, new)
        # B3LYP deleted (critical), PBE0 added (critical)
        paths = {e.path for e in entries}
        assert "keyword.B3LYP" in paths
        assert "keyword.PBE0" in paths
        b3lyp_entry = next(e for e in entries if e.path == "keyword.B3LYP")
        pbe0_entry = next(e for e in entries if e.path == "keyword.PBE0")
        assert b3lyp_entry.significance == "critical"
        assert pbe0_entry.significance == "critical"

    def test_basis_set_change_is_critical(self) -> None:
        old = self.parser.parse(SIMPLE_OPT)          # def2-TZVP
        new_content = SIMPLE_OPT.replace("def2-TZVP", "def2-SVP")
        new = self.parser.parse(new_content)
        entries = self.parser.diff(old, new)
        sigs = {e.path: e.significance for e in entries}
        assert sigs.get("keyword.DEF2-TZVP") == "critical"
        assert sigs.get("keyword.DEF2-SVP") == "critical"

    def test_dispersion_change_is_major(self) -> None:
        old = self.parser.parse(SIMPLE_OPT)          # D3BJ present
        new_content = SIMPLE_OPT.replace(" D3BJ", "")
        new = self.parser.parse(new_content)
        entries = self.parser.diff(old, new)
        d3bj_del = next((e for e in entries if e.path == "keyword.D3BJ"), None)
        assert d3bj_del is not None
        assert d3bj_del.significance == "major"

    def test_charge_change_is_critical(self) -> None:
        old = self.parser.parse(SIMPLE_OPT)          # charge 0
        new_content = SIMPLE_OPT.replace("xyz 0 1", "xyz 1 2")
        new = self.parser.parse(new_content)
        entries = self.parser.diff(old, new)
        charge_e = next(e for e in entries if e.path == "charge")
        mult_e = next(e for e in entries if e.path == "multiplicity")
        assert charge_e.significance == "critical"
        assert mult_e.significance == "critical"

    def test_maxcore_change_is_major(self) -> None:
        old = self.parser.parse(SIMPLE_OPT)          # maxcore 4000
        new_content = SIMPLE_OPT.replace("%maxcore 4000", "%maxcore 8000")
        new = self.parser.parse(new_content)
        entries = self.parser.diff(old, new)
        mc_e = next(e for e in entries if e.path == "maxcore")
        assert mc_e.significance == "major"
        assert mc_e.old_value == 4000
        assert mc_e.new_value == 8000

    def test_nprocs_change_is_major(self) -> None:
        old = self.parser.parse(SIMPLE_OPT)          # nprocs 8
        new_content = SIMPLE_OPT.replace("nprocs 8", "nprocs 16")
        new = self.parser.parse(new_content)
        entries = self.parser.diff(old, new)
        np_e = next(e for e in entries if e.path == "nprocs")
        assert np_e.significance == "major"

    def test_atom_count_change_is_critical(self) -> None:
        old = self.parser.parse(SIMPLE_OPT)
        new_content = SIMPLE_OPT.replace(
            "  H   -0.363  -0.513  -0.890\n", ""
        )
        new = self.parser.parse(new_content)
        entries = self.parser.diff(old, new)
        coord_e = next((e for e in entries if e.path == "coordinates"), None)
        assert coord_e is not None
        assert coord_e.significance == "critical"

    def test_block_key_change(self) -> None:
        old = self.parser.parse(WITH_SCF_BLOCK)
        new_content = WITH_SCF_BLOCK.replace("MaxIter 200", "MaxIter 500")
        new = self.parser.parse(new_content)
        entries = self.parser.diff(old, new)
        it_e = next(e for e in entries if e.path == "block.scf.maxiter")
        assert it_e.old_value == 200
        assert it_e.new_value == 500

    def test_keyword_added(self) -> None:
        old = self.parser.parse(SIMPLE_OPT)
        new_content = SIMPLE_OPT.replace("! B3LYP", "! B3LYP NMR")
        new = self.parser.parse(new_content)
        entries = self.parser.diff(old, new)
        nmr_e = next(e for e in entries if e.path == "keyword.NMR")
        assert nmr_e.change_type == "added"

    def test_diff_returns_diff_entry_objects(self) -> None:
        old = self.parser.parse(SIMPLE_OPT)
        new = self.parser.parse(SIMPLE_SP)
        entries = self.parser.diff(old, new)
        for e in entries:
            assert isinstance(e, DiffEntry)


# ---------------------------------------------------------------------------
# Format diff tests
# ---------------------------------------------------------------------------

class TestOrcaInputParserFormat:
    def setup_method(self) -> None:
        self.parser = OrcaInputParser()

    def test_format_no_changes(self) -> None:
        data = self.parser.parse(SIMPLE_OPT)
        result = self.parser.format_diff([])
        assert result == "No changes"

    def test_format_default_contains_marker(self) -> None:
        old = self.parser.parse(SIMPLE_OPT)
        new_content = SIMPLE_OPT.replace("B3LYP", "PBE0")
        new = self.parser.parse(new_content)
        entries = self.parser.diff(old, new)
        formatted = self.parser.format_diff(entries, style="default")
        assert "B3LYP" in formatted or "PBE0" in formatted

    def test_format_compact(self) -> None:
        old = self.parser.parse(SIMPLE_OPT)
        new_content = SIMPLE_OPT.replace("B3LYP", "PBE0")
        new = self.parser.parse(new_content)
        entries = self.parser.diff(old, new)
        formatted = self.parser.format_diff(entries, style="compact")
        assert "B3LYP" in formatted or "PBE0" in formatted

    def test_format_detailed(self) -> None:
        old = self.parser.parse(SIMPLE_OPT)
        new_content = SIMPLE_OPT.replace("B3LYP", "PBE0")
        new = self.parser.parse(new_content)
        entries = self.parser.diff(old, new)
        formatted = self.parser.format_diff(entries, style="detailed")
        assert "Significance" in formatted or "significance" in formatted.lower()


# ---------------------------------------------------------------------------
# Validate tests
# ---------------------------------------------------------------------------

class TestOrcaInputParserValidate:
    def setup_method(self) -> None:
        self.parser = OrcaInputParser()

    def test_valid_simple_opt(self) -> None:
        data = self.parser.parse(SIMPLE_OPT)
        ok, errors = self.parser.validate(data)
        assert ok
        assert errors == []

    def test_invalid_no_keywords(self) -> None:
        content = """\
* xyz 0 1
  H    0.000   0.000   0.000
  H    0.740   0.000   0.000
*
"""
        data = self.parser.parse(content)
        ok, errors = self.parser.validate(data)
        assert not ok
        assert any("keyword" in e.lower() for e in errors)

    def test_invalid_no_coordinates(self) -> None:
        content = "! B3LYP def2-SVP SP\n"
        data = self.parser.parse(content)
        ok, errors = self.parser.validate(data)
        assert not ok
        assert any("coordinate" in e.lower() for e in errors)

    def test_xyzfile_is_valid(self) -> None:
        data = self.parser.parse(XYZFILE)
        ok, errors = self.parser.validate(data)
        # xyzfile has coord_format set so the coordinate check passes
        assert ok


# ---------------------------------------------------------------------------
# _try_numeric helper
# ---------------------------------------------------------------------------

class TestTryNumeric:
    def test_int(self) -> None:
        assert _try_numeric("200") == 200

    def test_float(self) -> None:
        assert _try_numeric("80.4") == pytest.approx(80.4)

    def test_str(self) -> None:
        assert _try_numeric("Tight") == "Tight"

    def test_negative_int(self) -> None:
        assert _try_numeric("-5") == -5
