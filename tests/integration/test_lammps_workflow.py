"""Integration tests for ChemVCS LAMMPS workflows."""

import json
import os
from pathlib import Path

from typer.testing import CliRunner

from chemvcs.cli.main import app
from chemvcs.constants import CHEMVCS_DIR

runner = CliRunner()


LAMMPS_INPUT_V1 = """\
units lj
atom_style atomic
boundary p p p
read_data data.lammps
pair_style lj/cut 2.5
pair_coeff 1 1 1.0 1.0 2.5
fix 1 all nvt temp 1.0 1.0 0.1
thermo 100
timestep 0.005
run 50000
"""


LAMMPS_INPUT_V2 = """\
units lj
atom_style atomic
boundary p p p
read_data data.lammps
pair_style lj/cut 2.5
pair_coeff 1 1 1.0 1.0 2.5
fix 1 all nvt temp 1.0 1.0 0.5
thermo 500
timestep 0.002
run 200000
"""


LAMMPS_DATA_V1 = """\
LJ fluid test

500 atoms
1 atom types

0.0 20.0 xlo xhi
0.0 20.0 ylo yhi
0.0 20.0 zlo zhi

Masses

1 39.948

Atoms # atomic

1 1 1.0 1.0 1.0
2 1 2.0 2.0 2.0
"""


LAMMPS_DATA_V2 = """\
LJ fluid test

512 atoms
1 atom types

0.0 21.0 xlo xhi
0.0 21.0 ylo yhi
0.0 21.0 zlo zhi

Masses

1 39.948

Atoms # atomic

1 1 1.0 1.0 1.0
2 1 2.0 2.0 2.0
"""


LAMMPS_LOG_V1 = """\
LAMMPS (22 Jul 2025)
  500 atoms

Step Temp TotEng PotEng KinEng Press
0 1.0000 -1004.1000 -1005.9000 1.8000 0.3300
100 0.9998 -1004.2800 -1005.4930 1.8130 0.3368
50000 1.0000 -1004.2850 -1005.4936 1.8130 0.3368
Loop time of 14.2 seconds with 1 procs for 50000 steps with 500 atoms

Total wall time: 0:00:14
"""


LAMMPS_LOG_V2 = """\
LAMMPS (22 Jul 2025)
  512 atoms

Step Temp TotEng PotEng KinEng Press
0 1.0000 -1004.2000 -1005.9500 1.7500 0.3100
500 1.0001 -1006.7000 -1007.9200 1.2200 0.3200
200000 1.0000 -1007.8214 -1009.0344 1.2130 0.3368
Loop time of 57.4 seconds with 1 procs for 200000 steps with 512 atoms

Total wall time: 0:00:57
"""


class TestLammpsWorkflow:
    """Integration coverage for end-to-end LAMMPS usage."""

    def test_add_detects_lammps_file_types(self, tmp_path: Path) -> None:
        """Staging should classify LAMMPS input, data, and log files correctly."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)

        try:
            runner.invoke(app, ["init", "--quiet"])

            (tmp_path / "in.lammps").write_text(LAMMPS_INPUT_V1)
            (tmp_path / "data.lammps").write_text(LAMMPS_DATA_V1)
            (tmp_path / "log.lammps").write_text(LAMMPS_LOG_V1)

            result = runner.invoke(
                app,
                ["add", "--no-validate", "in.lammps", "data.lammps", "log.lammps"],
            )

            assert result.exit_code == 0

            index_path = tmp_path / CHEMVCS_DIR / "index"
            index = json.loads(index_path.read_text())

            assert index["entries"]["in.lammps"]["file_type"] == "LAMMPS_INPUT"
            assert index["entries"]["data.lammps"]["file_type"] == "LAMMPS_DATA"
            assert index["entries"]["log.lammps"]["file_type"] == "LAMMPS_LOG"
        finally:
            os.chdir(original_cwd)

    def test_default_ignore_skips_lammps_dump_files(self, tmp_path: Path) -> None:
        """Default ignore rules should skip large LAMMPS dump files."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)

        try:
            runner.invoke(app, ["init", "--quiet"])

            (tmp_path / "in.lammps").write_text(LAMMPS_INPUT_V1)
            (tmp_path / "dump.lj_nvt.lammpstrj").write_text("trajectory data")

            result = runner.invoke(
                app,
                ["add", "--no-validate", "in.lammps", "dump.lj_nvt.lammpstrj"],
            )

            assert result.exit_code == 0
            assert "Ignored" in result.stdout
            assert "dump.lj_nvt.lammpstrj" in result.stdout

            index_path = tmp_path / CHEMVCS_DIR / "index"
            index = json.loads(index_path.read_text())
            assert "in.lammps" in index["entries"]
            assert "dump.lj_nvt.lammpstrj" not in index["entries"]
        finally:
            os.chdir(original_cwd)

    def test_commit_diff_log_and_reproduce_lammps_workflow(self, tmp_path: Path) -> None:
        """ChemVCS should support a realistic LAMMPS workflow across core commands."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)

        try:
            runner.invoke(app, ["init", "--quiet"])

            (tmp_path / "in.lammps").write_text(LAMMPS_INPUT_V1)
            (tmp_path / "data.lammps").write_text(LAMMPS_DATA_V1)
            (tmp_path / "log.lammps").write_text(LAMMPS_LOG_V1)

            add_v1 = runner.invoke(
                app,
                ["add", "--no-validate", "in.lammps", "data.lammps", "log.lammps"],
            )
            assert add_v1.exit_code == 0

            commit_v1 = runner.invoke(app, ["commit", "-m", "Initial LAMMPS equilibration"])
            assert commit_v1.exit_code == 0
            assert "Semantic Changes:" not in commit_v1.stdout

            (tmp_path / "in.lammps").write_text(LAMMPS_INPUT_V2)
            (tmp_path / "data.lammps").write_text(LAMMPS_DATA_V2)
            (tmp_path / "log.lammps").write_text(LAMMPS_LOG_V2)

            add_v2 = runner.invoke(
                app,
                ["add", "--no-validate", "in.lammps", "data.lammps", "log.lammps"],
            )
            assert add_v2.exit_code == 0

            commit_v2 = runner.invoke(app, ["commit", "-m", "Production LAMMPS run"])
            assert commit_v2.exit_code == 0
            assert "Semantic Changes:" in commit_v2.stdout
            assert "in.lammps" in commit_v2.stdout
            assert "data.lammps" in commit_v2.stdout
            assert "log.lammps" in commit_v2.stdout
            assert "timestep" in commit_v2.stdout
            assert "num_atoms" in commit_v2.stdout
            assert "final_etotal" in commit_v2.stdout

            # Modify the working tree after commit so `chemvcs diff` compares
            # HEAD against real current files instead of a clean tree.
            (tmp_path / "in.lammps").write_text(
                LAMMPS_INPUT_V2.replace("run 200000", "run 250000")
            )
            (tmp_path / "log.lammps").write_text(
                LAMMPS_LOG_V2.replace("-1007.8214", "-1008.5000")
            )

            diff_result = runner.invoke(app, ["diff", "--summary"])
            assert diff_result.exit_code == 0
            assert "Comparing:" in diff_result.stdout
            assert "in.lammps" in diff_result.stdout
            assert "log.lammps" in diff_result.stdout
            assert "data.lammps" not in diff_result.stdout
            assert "critical" in diff_result.stdout.lower()

            log_result = runner.invoke(app, ["log", "--oneline"])
            assert log_result.exit_code == 0
            assert "Production LAMMPS run" in log_result.stdout
            assert "Initial LAMMPS equilibration" in log_result.stdout
            assert log_result.stdout.index("Production LAMMPS run") < log_result.stdout.index(
                "Initial LAMMPS equilibration"
            )

            head_hash = (tmp_path / CHEMVCS_DIR / "HEAD").read_text().strip()
            reproduce_result = runner.invoke(app, ["reproduce", head_hash])
            assert reproduce_result.exit_code == 0

            reproduced_dir = tmp_path / f"reproduce_{head_hash[:7]}"
            assert (reproduced_dir / "in.lammps").read_text() == LAMMPS_INPUT_V2
            assert (reproduced_dir / "data.lammps").read_text() == LAMMPS_DATA_V2
            assert (reproduced_dir / "log.lammps").read_text() == LAMMPS_LOG_V2
        finally:
            os.chdir(original_cwd)

    def test_reproduce_keeps_unmodified_lammps_parent_files(self, tmp_path: Path) -> None:
        """Later commits should still reproduce unchanged files inherited from parents."""
        original_cwd = Path.cwd()
        os.chdir(tmp_path)

        try:
            runner.invoke(app, ["init", "--quiet"])

            (tmp_path / "in.lammps").write_text(LAMMPS_INPUT_V1)
            (tmp_path / "data.lammps").write_text(LAMMPS_DATA_V1)

            runner.invoke(app, ["add", "--no-validate", "in.lammps", "data.lammps"])
            commit_v1 = runner.invoke(app, ["commit", "-m", "Initial LAMMPS inputs"])
            assert commit_v1.exit_code == 0

            (tmp_path / "log.lammps").write_text(LAMMPS_LOG_V1)
            runner.invoke(app, ["add", "--no-validate", "log.lammps"])
            commit_v2 = runner.invoke(app, ["commit", "-m", "Finished NVT equilibration"])
            assert commit_v2.exit_code == 0

            (tmp_path / "in.lammps").write_text(LAMMPS_INPUT_V2)
            (tmp_path / "log.lammps").write_text(LAMMPS_LOG_V2)
            runner.invoke(app, ["add", "--no-validate", "in.lammps", "log.lammps"])
            commit_v3 = runner.invoke(app, ["commit", "-m", "Production LAMMPS run"])
            assert commit_v3.exit_code == 0

            head_hash = (tmp_path / CHEMVCS_DIR / "HEAD").read_text().strip()
            reproduce_result = runner.invoke(app, ["reproduce", head_hash])
            assert reproduce_result.exit_code == 0

            reproduced_dir = tmp_path / f"reproduce_{head_hash[:7]}"
            assert (reproduced_dir / "in.lammps").read_text() == LAMMPS_INPUT_V2
            assert (reproduced_dir / "log.lammps").read_text() == LAMMPS_LOG_V2
            assert (reproduced_dir / "data.lammps").read_text() == LAMMPS_DATA_V1
        finally:
            os.chdir(original_cwd)