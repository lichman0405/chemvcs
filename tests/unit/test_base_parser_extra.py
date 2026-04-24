"""Extra coverage tests for base_parser.py.

Targets uncovered lines:
    49  – DiffEntry.__repr__ added branch
    80  – DiffEntry.__repr__ deleted branch
    97  – BaseParser.parse abstract body (pass)
   114  – BaseParser.diff abstract body (pass)
   121  – format_diff compact added branch
   123  – format_diff compact deleted branch
  128-136 – format_diff detailed style
   149  – BaseParser.validate default return
"""


from chemvcs.parsers.base_parser import BaseParser, DiffEntry
from chemvcs.parsers.incar_parser import IncarParser  # concrete implementation

# ---------------------------------------------------------------------------
# DiffEntry.__repr__
# ---------------------------------------------------------------------------


class TestDiffEntryRepr:
    def test_repr_added(self) -> None:
        e = DiffEntry("ENCUT", None, 500, "added", "critical")
        r = repr(e)
        assert "+" in r
        assert "ENCUT" in r
        assert "500" in r
        assert "critical" in r

    def test_repr_deleted(self) -> None:
        e = DiffEntry("ENCUT", 400, None, "deleted", "critical")
        r = repr(e)
        assert "-" in r
        assert "ENCUT" in r
        assert "400" in r
        assert "critical" in r

    def test_repr_modified(self) -> None:
        e = DiffEntry("ENCUT", 400, 500, "modified", "major")
        r = repr(e)
        assert "→" in r
        assert "400" in r
        assert "500" in r

    def test_to_dict_all_fields(self) -> None:
        e = DiffEntry("SIGMA", 0.1, 0.05, "modified", "minor")
        d = e.to_dict()
        assert d["path"] == "SIGMA"
        assert d["old_value"] == 0.1
        assert d["new_value"] == 0.05
        assert d["change_type"] == "modified"
        assert d["significance"] == "minor"

    def test_to_dict_added(self) -> None:
        e = DiffEntry("NSW", None, 100, "added", "major")
        d = e.to_dict()
        assert d["old_value"] is None
        assert d["new_value"] == 100
        assert d["change_type"] == "added"

    def test_to_dict_deleted(self) -> None:
        e = DiffEntry("MAGMOM", "1 1 1", None, "deleted", "minor")
        d = e.to_dict()
        assert d["old_value"] == "1 1 1"
        assert d["new_value"] is None
        assert d["change_type"] == "deleted"


# ---------------------------------------------------------------------------
# BaseParser.format_diff – compact style
# ---------------------------------------------------------------------------


class TestFormatDiffCompact:
    def setup_method(self) -> None:
        self.parser = IncarParser()

    def test_compact_added(self) -> None:
        e = DiffEntry("ENCUT", None, 500, "added", "critical")
        result = self.parser.format_diff([e], style="compact")
        assert result == "+ ENCUT = 500"

    def test_compact_deleted(self) -> None:
        e = DiffEntry("ENCUT", 400, None, "deleted", "major")
        result = self.parser.format_diff([e], style="compact")
        assert result == "- ENCUT = 400"

    def test_compact_modified(self) -> None:
        e = DiffEntry("ENCUT", 400, 500, "modified", "critical")
        result = self.parser.format_diff([e], style="compact")
        assert "~ ENCUT" in result
        assert "400" in result
        assert "500" in result

    def test_compact_multiple_entries(self) -> None:
        entries = [
            DiffEntry("ENCUT", None, 500, "added", "critical"),
            DiffEntry("SIGMA", 0.2, None, "deleted", "minor"),
            DiffEntry("NSW", 50, 100, "modified", "major"),
        ]
        result = self.parser.format_diff(entries, style="compact")
        lines = result.splitlines()
        assert len(lines) == 3
        assert lines[0].startswith("+")
        assert lines[1].startswith("-")
        assert lines[2].startswith("~")


# ---------------------------------------------------------------------------
# BaseParser.format_diff – detailed style
# ---------------------------------------------------------------------------


class TestFormatDiffDetailed:
    def setup_method(self) -> None:
        self.parser = IncarParser()

    def test_detailed_modified_has_all_fields(self) -> None:
        e = DiffEntry("ENCUT", 400, 500, "modified", "critical")
        result = self.parser.format_diff([e], style="detailed")
        assert "Path: ENCUT" in result
        assert "Type: modified" in result
        assert "Significance: critical" in result
        assert "Old: 400" in result
        assert "New: 500" in result
        # Empty-line separator between entries
        assert result.endswith("\n") or result.count("Path:") == 1

    def test_detailed_added_no_old_value(self) -> None:
        e = DiffEntry("NSW", None, 100, "added", "major")
        result = self.parser.format_diff([e], style="detailed")
        assert "Path: NSW" in result
        assert "New: 100" in result
        # old_value is None → "Old: ..." should NOT appear
        assert "Old:" not in result

    def test_detailed_deleted_no_new_value(self) -> None:
        e = DiffEntry("MAGMOM", "1 1 1", None, "deleted", "minor")
        result = self.parser.format_diff([e], style="detailed")
        assert "Old: 1 1 1" in result
        assert "New:" not in result

    def test_detailed_multiple_entries_separated(self) -> None:
        entries = [
            DiffEntry("ENCUT", 400, 500, "modified", "critical"),
            DiffEntry("SIGMA", 0.2, 0.05, "modified", "minor"),
        ]
        result = self.parser.format_diff(entries, style="detailed")
        assert result.count("Path:") == 2


# ---------------------------------------------------------------------------
# BaseParser.format_diff – default style (already mostly covered, extras)
# ---------------------------------------------------------------------------


class TestFormatDiffDefault:
    def setup_method(self) -> None:
        self.parser = IncarParser()

    def test_no_changes_returns_no_changes(self) -> None:
        assert self.parser.format_diff([]) == "No changes"

    def test_default_added_marker(self) -> None:
        e = DiffEntry("ENCUT", None, 500, "added", "minor")
        result = self.parser.format_diff([e])
        assert "Added ENCUT" in result

    def test_default_deleted_marker(self) -> None:
        e = DiffEntry("ENCUT", 400, None, "deleted", "minor")
        result = self.parser.format_diff([e])
        assert "Deleted ENCUT" in result


# ---------------------------------------------------------------------------
# BaseParser.validate – default implementation
# ---------------------------------------------------------------------------


class TestBaseParserValidate:
    def test_validate_returns_true_by_default(self) -> None:
        """IncarParser inherits BaseParser.validate which returns (True, [])
        unless overridden.  The INCAR parser's own validate may differ;
        we call the base implementation directly."""
        parser = IncarParser()
        # Call the base class method directly on a concrete instance
        is_valid, errors = BaseParser.validate(parser, {})
        assert is_valid is True
        assert errors == []


# ---------------------------------------------------------------------------
# BaseParser abstract method bodies (pass statements)
# ---------------------------------------------------------------------------


class TestAbstractMethodBodies:
    """Call abstract method bodies via class reference to exercise the
    ``pass`` statements (lines 97 and 114)."""

    def test_abstract_parse_body_returns_none(self) -> None:
        parser = IncarParser()
        result = BaseParser.parse(parser, "ISTART = 1")
        assert result is None

    def test_abstract_diff_body_returns_none(self) -> None:
        parser = IncarParser()
        result = BaseParser.diff(parser, {}, {})
        assert result is None
