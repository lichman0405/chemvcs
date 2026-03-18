"""Unit tests for LammpsDataParser."""

import pytest

from chemvcs.parsers.lammps_data_parser import LammpsDataParser
from chemvcs.parsers.base_parser import ParserError


# ---------------------------------------------------------------------------
# Sample data file content
# ---------------------------------------------------------------------------

# Atomic (no topology) – LJ fluid
SAMPLE_DATA_ATOMIC = """\
LJ fluid - 500 atoms - NVT test

500 atoms
1 atom types

0.0 20.0 xlo xhi
0.0 20.0 ylo yhi
0.0 20.0 zlo zhi

Masses

1 39.948

Atoms # atomic

1 1 3.56 7.12 10.5
2 1 4.20 1.50  2.3
"""

# Full atomistic with topology (molecular/polymer)
SAMPLE_DATA_MOLECULAR = """\
Water box - SPC/E model

2880 atoms
2 atom types
2880 bonds
1 bond types
1440 angles
1 angle types

0.0 40.0 xlo xhi
0.0 40.0 ylo yhi
0.0 60.0 zlo zhi

Masses

1  15.9994
2   1.0080

Atoms # full

1 1 1 -0.8476  1.0  1.0  1.0
2 1 2  0.4238  1.0  1.5  1.5
"""

# Triclinic box
SAMPLE_DATA_TRICLINIC = """\
Triclinic unit cell

100 atoms
1 atom types

0.0 10.0 xlo xhi
0.0 8.66 ylo yhi
0.0 14.14 zlo zhi
1.0 0.5 0.25 xy xz yz

Masses

1 26.98

Atoms # atomic

1 1 0.5 0.5 0.5
"""

# Missing atom count (invalid)
SAMPLE_DATA_NO_ATOMS = """\
Bad data file

1 atom types
0.0 10.0 xlo xhi
0.0 10.0 ylo yhi
0.0 10.0 zlo zhi
"""


# ---------------------------------------------------------------------------
# Test cases
# ---------------------------------------------------------------------------

class TestLammpsDataParser:

    def test_parse_atomic_returns_dict(self) -> None:
        parser = LammpsDataParser()
        data = parser.parse(SAMPLE_DATA_ATOMIC)
        assert isinstance(data, dict)

    def test_parse_atomic_num_atoms(self) -> None:
        parser = LammpsDataParser()
        data = parser.parse(SAMPLE_DATA_ATOMIC)
        assert data["num_atoms"] == 500

    def test_parse_atomic_num_types(self) -> None:
        parser = LammpsDataParser()
        data = parser.parse(SAMPLE_DATA_ATOMIC)
        assert data["num_atom_types"] == 1

    def test_parse_atomic_box_dimensions(self) -> None:
        parser = LammpsDataParser()
        data = parser.parse(SAMPLE_DATA_ATOMIC)
        assert data["lx"] == pytest.approx(20.0)
        assert data["ly"] == pytest.approx(20.0)
        assert data["lz"] == pytest.approx(20.0)

    def test_parse_atomic_volume(self) -> None:
        parser = LammpsDataParser()
        data = parser.parse(SAMPLE_DATA_ATOMIC)
        assert data["volume"] == pytest.approx(8000.0)

    def test_parse_atomic_no_bonds(self) -> None:
        parser = LammpsDataParser()
        data = parser.parse(SAMPLE_DATA_ATOMIC)
        assert data["has_bonds"] is False

    def test_parse_atomic_masses(self) -> None:
        parser = LammpsDataParser()
        data = parser.parse(SAMPLE_DATA_ATOMIC)
        masses = data.get("masses", {})
        assert "1" in masses
        assert masses["1"] == pytest.approx(39.948)

    def test_parse_atomic_sections_present(self) -> None:
        parser = LammpsDataParser()
        data = parser.parse(SAMPLE_DATA_ATOMIC)
        sections = data.get("sections", [])
        assert any("Masses" in s for s in sections)
        assert any("Atoms" in s for s in sections)

    def test_parse_molecular_has_bonds(self) -> None:
        parser = LammpsDataParser()
        data = parser.parse(SAMPLE_DATA_MOLECULAR)
        assert data["has_bonds"] is True
        assert data["num_bonds"] == 2880
        assert data["num_bond_types"] == 1

    def test_parse_molecular_has_angles(self) -> None:
        parser = LammpsDataParser()
        data = parser.parse(SAMPLE_DATA_MOLECULAR)
        assert data["num_angles"] == 1440
        assert data["num_angle_types"] == 1

    def test_parse_molecular_masses(self) -> None:
        parser = LammpsDataParser()
        data = parser.parse(SAMPLE_DATA_MOLECULAR)
        masses = data.get("masses", {})
        assert masses["1"] == pytest.approx(15.9994)
        assert masses["2"] == pytest.approx(1.0080)

    def test_parse_triclinic_tilt_factors(self) -> None:
        parser = LammpsDataParser()
        data = parser.parse(SAMPLE_DATA_TRICLINIC)
        assert data.get("xy") == pytest.approx(1.0)
        assert data.get("xz") == pytest.approx(0.5)
        assert data.get("yz") == pytest.approx(0.25)

    def test_parse_title_captured(self) -> None:
        parser = LammpsDataParser()
        data = parser.parse(SAMPLE_DATA_ATOMIC)
        assert "LJ fluid" in data.get("title", "")

    def test_parse_empty_raises(self) -> None:
        parser = LammpsDataParser()
        with pytest.raises(ParserError):
            parser.parse("")

    def test_parse_no_atom_count_raises(self) -> None:
        parser = LammpsDataParser()
        with pytest.raises(ParserError):
            parser.parse(SAMPLE_DATA_NO_ATOMS)

    # -- diff tests --

    def test_diff_no_change(self) -> None:
        parser = LammpsDataParser()
        data = parser.parse(SAMPLE_DATA_ATOMIC)
        entries = parser.diff(data, data)
        assert entries == []

    def test_diff_atom_count_change(self) -> None:
        parser = LammpsDataParser()
        old = parser.parse(SAMPLE_DATA_ATOMIC)
        new = parser.parse(SAMPLE_DATA_MOLECULAR)
        entries = parser.diff(old, new)
        paths = [e.path for e in entries]
        assert "num_atoms" in paths

    def test_diff_atom_count_is_critical(self) -> None:
        parser = LammpsDataParser()
        old = parser.parse(SAMPLE_DATA_ATOMIC)
        new = parser.parse(SAMPLE_DATA_MOLECULAR)
        entries = parser.diff(old, new)
        atom_entry = next(e for e in entries if e.path == "num_atoms")
        assert atom_entry.significance == "critical"

    def test_diff_bonds_change_is_major(self) -> None:
        parser = LammpsDataParser()
        old = parser.parse(SAMPLE_DATA_ATOMIC)
        new = parser.parse(SAMPLE_DATA_MOLECULAR)
        entries = parser.diff(old, new)
        bond_entry = next((e for e in entries if e.path == "has_bonds"), None)
        assert bond_entry is not None
        assert bond_entry.significance == "major"

    def test_diff_volume_change(self) -> None:
        parser = LammpsDataParser()
        old = parser.parse(SAMPLE_DATA_ATOMIC)
        new = parser.parse(SAMPLE_DATA_TRICLINIC)
        entries = parser.diff(old, new)
        paths = [e.path for e in entries]
        assert "volume" in paths

    # -- validate tests --

    def test_validate_good_data(self) -> None:
        parser = LammpsDataParser()
        data = parser.parse(SAMPLE_DATA_ATOMIC)
        is_valid, errors = parser.validate(data)
        assert is_valid
        assert errors == []

    def test_validate_zero_atoms_fails(self) -> None:
        parser = LammpsDataParser()
        data = {"num_atoms": 0, "num_atom_types": 1}
        is_valid, errors = parser.validate(data)
        assert not is_valid
        assert any("num_atoms" in e for e in errors)
