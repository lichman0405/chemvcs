"""Tests for DiffEngine LAMMPS file type detection and diff."""

from chemvcs.parsers.diff_engine import DiffEngine


class TestDiffEngineLammps:
    """Tests for LAMMPS file type detection in DiffEngine."""

    def setup_method(self) -> None:
        self.engine = DiffEngine()

    # -- File type detection -------------------------------------------

    def test_detect_in_lammps(self) -> None:
        assert self.engine.get_file_type("in.lammps") == "LAMMPS_INPUT"

    def test_detect_in_dot_prefix(self) -> None:
        assert self.engine.get_file_type("in.nvt") == "LAMMPS_INPUT"
        assert self.engine.get_file_type("in.minimize") == "LAMMPS_INPUT"

    def test_detect_lammps_suffix(self) -> None:
        assert self.engine.get_file_type("lj_fluid.lammps") == "LAMMPS_INPUT"

    def test_detect_lammps_in(self) -> None:
        assert self.engine.get_file_type("lammps.in") == "LAMMPS_INPUT"

    def test_detect_data_prefix(self) -> None:
        assert self.engine.get_file_type("data.lj_fluid") == "LAMMPS_DATA"
        assert self.engine.get_file_type("data.polymer") == "LAMMPS_DATA"

    def test_detect_lmp_suffix(self) -> None:
        assert self.engine.get_file_type("crystal.lmp") == "LAMMPS_DATA"

    def test_detect_data_suffix(self) -> None:
        assert self.engine.get_file_type("system.data") == "LAMMPS_DATA"

    def test_detect_log_lammps(self) -> None:
        assert self.engine.get_file_type("log.lammps") == "LAMMPS_LOG"

    def test_detect_log_prefix(self) -> None:
        assert self.engine.get_file_type("log.nvt") == "LAMMPS_LOG"

    def test_detect_lammps_log(self) -> None:
        assert self.engine.get_file_type("lammps.log") == "LAMMPS_LOG"

    def test_detect_vasp_still_works(self) -> None:
        assert self.engine.get_file_type("INCAR") == "INCAR"
        assert self.engine.get_file_type("OUTCAR") == "OUTCAR"
        assert self.engine.get_file_type("KPOINTS") == "KPOINTS"

    def test_can_parse_lammps_input(self) -> None:
        assert self.engine.can_parse("in.lammps") is True

    def test_can_parse_lammps_data(self) -> None:
        assert self.engine.can_parse("data.crystal") is True

    def test_can_parse_lammps_log(self) -> None:
        assert self.engine.can_parse("log.lammps") is True

    def test_dump_not_parseable(self) -> None:
        # dump files are tracked, not semantically parsed
        assert self.engine.can_parse("dump.nvt") is False

    # -- Diff via engine -----------------------------------------------

    LAMMPS_INPUT_V1 = """\
units lj
atom_style atomic
boundary p p p
pair_style lj/cut 2.5
pair_coeff * * 1.0 1.0 2.5
fix 1 all nvt temp 1.0 1.0 0.1
thermo 100
timestep 0.005
run 10000
"""

    LAMMPS_INPUT_V2 = """\
units lj
atom_style atomic
boundary p p p
pair_style lj/cut 3.0
pair_coeff * * 1.2 1.0 3.0
fix 1 all nvt temp 2.0 2.0 0.1
thermo 100
timestep 0.005
run 50000
"""

    def test_diff_files_lammps_input(self) -> None:
        entries = self.engine.diff_files(self.LAMMPS_INPUT_V1, self.LAMMPS_INPUT_V2, "in.lammps")
        assert entries is not None
        paths = [e.path for e in entries]
        assert "pair_style" in paths
        assert "run" in paths

    LAMMPS_LOG_V1 = """\
500 atoms
Step Temp TotEng
0 300.0 -1000.0
100 301.0 -1005.0
Total wall time: 0:00:10
"""

    LAMMPS_LOG_V2 = """\
500 atoms
Step Temp TotEng
0 300.0 -1000.0
100 305.0 -1050.0
Total wall time: 0:00:12
"""

    def test_diff_files_lammps_log(self) -> None:
        entries = self.engine.diff_files(self.LAMMPS_LOG_V1, self.LAMMPS_LOG_V2, "log.lammps")
        assert entries is not None
        paths = [e.path for e in entries]
        assert "final_etotal" in paths

    LAMMPS_DATA_V1 = """\
LJ box v1
500 atoms
1 atom types
0.0 20.0 xlo xhi
0.0 20.0 ylo yhi
0.0 20.0 zlo zhi
Masses
1 39.948
Atoms # atomic
1 1 1.0 2.0 3.0
"""

    LAMMPS_DATA_V2 = """\
LJ box v2
800 atoms
1 atom types
0.0 25.0 xlo xhi
0.0 25.0 ylo yhi
0.0 25.0 zlo zhi
Masses
1 39.948
Atoms # atomic
1 1 1.0 2.0 3.0
"""

    def test_diff_files_lammps_data(self) -> None:
        entries = self.engine.diff_files(self.LAMMPS_DATA_V1, self.LAMMPS_DATA_V2, "data.lammps")
        assert entries is not None
        paths = [e.path for e in entries]
        assert "num_atoms" in paths
        assert "volume" in paths
