"""Unit tests for LammpsInputParser."""

import pytest

from chemvcs.parsers.lammps_input_parser import LammpsInputParser
from chemvcs.parsers.base_parser import ParserError


# ---------------------------------------------------------------------------
# Sample input scripts
# ---------------------------------------------------------------------------

# LJ fluid NVT – minimal but realistic
SAMPLE_INPUT_NVT = """\
# LJ fluid NVT simulation
# =======================

units           lj
atom_style      atomic
dimension       3
boundary        p p p

# --- Geometry ---
read_data       data.lj_fluid

# --- Interactions ---
pair_style      lj/cut 2.5
pair_coeff      * * 1.0 1.0 2.5

# --- Ensemble ---
velocity        all create 1.0 12345 mom yes rot yes
fix             1 all nvt temp 1.0 1.0 0.1

# --- Output ---
thermo          100
thermo_style    custom step temp etotal press
dump            1 all atom 1000 dump.nvt.lammpstrj

# --- Run ---
timestep        0.005
run             50000
"""

# Energy minimization
SAMPLE_INPUT_MINIMIZE = """\
units           metal
atom_style      charge
boundary        p p p

read_data       data.crystal

pair_style      reaxff NULL
pair_coeff      * * ffield.reax C H O

minimize        1.0e-4 1.0e-6 1000 10000

timestep        0.0005
"""

# NPT with variable
SAMPLE_INPUT_NPT = """\
units real
atom_style full
boundary p p p

read_data data.polymer
read_data data.solvent add append

pair_style      lj/charmm/coul/long 8.0 10.0
kspace_style    pppm 1.0e-5

variable        T equal 300.0
variable        P equal 1.0

fix NVT_EQ all nvt temp ${T} ${T} 100.0
run 10000

fix NPT_PROD all npt temp ${T} ${T} 100.0 iso ${P} ${P} 1000.0

thermo 500
dump 1 all custom 2000 dump.npt.lammpstrj id type x y z

timestep 1.0
run 500000
"""

# Continuation lines
SAMPLE_INPUT_CONTINUATION = """\
units           lj
atom_style      atomic
boundary        p p p
pair_style      lj/cut 2.5
pair_coeff      * * 1.0 1.0 2.5
fix             1 all nve
timestep        0.001
run             &
                10000
"""

# Empty input
SAMPLE_INPUT_EMPTY = ""


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestLammpsInputParser:

    def test_parse_nvt_returns_dict(self) -> None:
        parser = LammpsInputParser()
        data = parser.parse(SAMPLE_INPUT_NVT)
        assert isinstance(data, dict)

    def test_parse_nvt_units(self) -> None:
        parser = LammpsInputParser()
        data = parser.parse(SAMPLE_INPUT_NVT)
        assert data["units"] == "lj"

    def test_parse_nvt_atom_style(self) -> None:
        parser = LammpsInputParser()
        data = parser.parse(SAMPLE_INPUT_NVT)
        assert data["atom_style"] == "atomic"

    def test_parse_nvt_boundary(self) -> None:
        parser = LammpsInputParser()
        data = parser.parse(SAMPLE_INPUT_NVT)
        assert data["boundary"] == "p p p"

    def test_parse_nvt_pair_style(self) -> None:
        parser = LammpsInputParser()
        data = parser.parse(SAMPLE_INPUT_NVT)
        assert data["pair_style"].startswith("lj/cut")

    def test_parse_nvt_timestep(self) -> None:
        parser = LammpsInputParser()
        data = parser.parse(SAMPLE_INPUT_NVT)
        assert data["timestep"] == pytest.approx(0.005)

    def test_parse_nvt_run_steps(self) -> None:
        parser = LammpsInputParser()
        data = parser.parse(SAMPLE_INPUT_NVT)
        assert data["run"] == 50000

    def test_parse_nvt_thermo_interval(self) -> None:
        parser = LammpsInputParser()
        data = parser.parse(SAMPLE_INPUT_NVT)
        assert data["thermo"] == 100

    def test_parse_nvt_fix_nvt(self) -> None:
        parser = LammpsInputParser()
        data = parser.parse(SAMPLE_INPUT_NVT)
        nvt = data.get("fix_nvt")
        assert nvt is not None
        assert nvt["Tstart"] == pytest.approx(1.0)
        assert nvt["Tstop"] == pytest.approx(1.0)
        assert nvt["Tdamp"] == pytest.approx(0.1)

    def test_parse_nvt_read_data(self) -> None:
        parser = LammpsInputParser()
        data = parser.parse(SAMPLE_INPUT_NVT)
        assert "data.lj_fluid" in data.get("read_data", [])

    def test_parse_nvt_dump_interval(self) -> None:
        parser = LammpsInputParser()
        data = parser.parse(SAMPLE_INPUT_NVT)
        assert data.get("dump_interval") == 1000

    def test_parse_minimize_etol(self) -> None:
        parser = LammpsInputParser()
        data = parser.parse(SAMPLE_INPUT_MINIMIZE)
        assert data["minimize_etol"] == pytest.approx(1.0e-4)
        assert data["minimize_ftol"] == pytest.approx(1.0e-6)
        assert data["minimize_maxiter"] == 1000
        assert data["minimize_maxeval"] == 10000

    def test_parse_minimize_units(self) -> None:
        parser = LammpsInputParser()
        data = parser.parse(SAMPLE_INPUT_MINIMIZE)
        assert data["units"] == "metal"

    def test_parse_npt_kspace(self) -> None:
        parser = LammpsInputParser()
        data = parser.parse(SAMPLE_INPUT_NPT)
        assert "pppm" in data.get("kspace_style", "")

    def test_parse_npt_fix_npt(self) -> None:
        parser = LammpsInputParser()
        data = parser.parse(SAMPLE_INPUT_NPT)
        npt = data.get("fix_npt")
        assert npt is not None

    def test_parse_npt_multiple_read_data(self) -> None:
        parser = LammpsInputParser()
        data = parser.parse(SAMPLE_INPUT_NPT)
        read_data = data.get("read_data", [])
        assert len(read_data) == 2

    def test_parse_npt_variables(self) -> None:
        parser = LammpsInputParser()
        data = parser.parse(SAMPLE_INPUT_NPT)
        vars_ = data.get("variable_names", [])
        assert "T" in vars_
        assert "P" in vars_

    def test_parse_continuation_line(self) -> None:
        parser = LammpsInputParser()
        data = parser.parse(SAMPLE_INPUT_CONTINUATION)
        assert data["run"] == 10000

    def test_parse_nve_fix(self) -> None:
        parser = LammpsInputParser()
        data = parser.parse(SAMPLE_INPUT_CONTINUATION)
        assert data.get("fix_nve") is True

    def test_parse_empty_raises(self) -> None:
        parser = LammpsInputParser()
        with pytest.raises(ParserError):
            parser.parse(SAMPLE_INPUT_EMPTY)

    # -- diff tests --

    def test_diff_no_change(self) -> None:
        parser = LammpsInputParser()
        data = parser.parse(SAMPLE_INPUT_NVT)
        entries = parser.diff(data, data)
        assert entries == []

    def test_diff_units_change(self) -> None:
        parser = LammpsInputParser()
        old = parser.parse(SAMPLE_INPUT_NVT)
        new = parser.parse(SAMPLE_INPUT_MINIMIZE)
        entries = parser.diff(old, new)
        paths = [e.path for e in entries]
        assert "units" in paths

    def test_diff_units_is_critical(self) -> None:
        parser = LammpsInputParser()
        old = parser.parse(SAMPLE_INPUT_NVT)
        new = parser.parse(SAMPLE_INPUT_MINIMIZE)
        entries = parser.diff(old, new)
        u = next(e for e in entries if e.path == "units")
        assert u.significance == "critical"

    def test_diff_pair_style_change(self) -> None:
        parser = LammpsInputParser()
        old = parser.parse(SAMPLE_INPUT_NVT)
        new = parser.parse(SAMPLE_INPUT_MINIMIZE)
        entries = parser.diff(old, new)
        paths = [e.path for e in entries]
        assert "pair_style" in paths

    def test_diff_timestep_change(self) -> None:
        parser = LammpsInputParser()
        old = parser.parse(SAMPLE_INPUT_NVT)
        new = parser.parse(SAMPLE_INPUT_MINIMIZE)
        entries = parser.diff(old, new)
        paths = [e.path for e in entries]
        assert "timestep" in paths

    # -- validate tests --

    def test_validate_good_script(self) -> None:
        parser = LammpsInputParser()
        data = parser.parse(SAMPLE_INPUT_NVT)
        is_valid, errors = parser.validate(data)
        assert is_valid
        assert errors == []

    def test_validate_missing_units(self) -> None:
        parser = LammpsInputParser()
        # Manually remove units
        data = parser.parse(SAMPLE_INPUT_NVT)
        del data["units"]
        is_valid, errors = parser.validate(data)
        assert not is_valid
        assert any("units" in e for e in errors)

    def test_validate_missing_pair_style(self) -> None:
        parser = LammpsInputParser()
        data = parser.parse(SAMPLE_INPUT_NVT)
        del data["pair_style"]
        is_valid, errors = parser.validate(data)
        assert not is_valid
        assert any("pair_style" in e for e in errors)
