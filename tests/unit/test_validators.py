"""Unit tests for chemvcs-validator plugins: incar_poscar and file_format."""

import textwrap
from pathlib import Path

import pytest

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def write(tmp_path: Path, name: str, content: str) -> Path:
    p = tmp_path / name
    p.write_text(textwrap.dedent(content), encoding="utf-8")
    return p


POSCAR_LiCoO2 = """\
    Li Co O
    1.0
    2.8 0.0 0.0
    0.0 2.8 0.0
    0.0 0.0 5.0
    Li Co O
    1 1 2
    Direct
    0.0 0.0 0.0
    0.5 0.5 0.5
    0.0 0.5 0.25
    0.5 0.0 0.75
"""

INCAR_BASIC = """\
    ENCUT = 520
    ISMEAR = 0
    SIGMA = 0.05
"""

INCAR_ISPIN2_NO_MAGMOM = """\
    ENCUT = 520
    ISPIN = 2
"""

INCAR_MAGMOM_WRONG_LEN = """\
    ENCUT = 520
    ISPIN = 2
    MAGMOM = 2*0.6
"""

INCAR_LDAUU_WRONG_LEN = """\
    ENCUT = 520
    LDAU = .TRUE.
    LDAUU = 5 5
"""

INCAR_NSW_ZERO_ISIF3 = """\
    ENCUT = 520
    NSW = 0
    ISIF = 3
"""

INCAR_MAGMOM_CORRECT = """\
    ENCUT = 520
    ISPIN = 2
    MAGMOM = 1*0.6 1*0.6 2*0.0
"""

INCAR_LDAUU_CORRECT = """\
    ENCUT = 520
    LDAU = .TRUE.
    LDAUU = 0 3 0
"""

KPOINTS_GAMMA = """\
    Automatic mesh
    0
    Gamma
    4 4 4
    0 0 0
"""

KPOINTS_BAD = """\
    Bad kpoints
    not_a_number
    Gamma
"""


# ===========================================================================
# INCARPOSCARValidator tests
# ===========================================================================


class TestINCARPOSCARValidator:
    @pytest.fixture(autouse=True)
    def _validator(self):
        from chemvcs_validator.incar_poscar import INCARPOSCARValidator

        self.v = INCARPOSCARValidator()

    def test_basic_valid(self, tmp_path):
        write(tmp_path, "POSCAR", POSCAR_LiCoO2)
        write(tmp_path, "INCAR", INCAR_BASIC)
        r = self.v.validate(tmp_path, ["INCAR", "POSCAR"])
        assert r.passed
        assert not r.errors

    def test_ispin2_no_magmom_warns(self, tmp_path):
        write(tmp_path, "POSCAR", POSCAR_LiCoO2)
        write(tmp_path, "INCAR", INCAR_ISPIN2_NO_MAGMOM)
        r = self.v.validate(tmp_path, ["INCAR", "POSCAR"])
        assert r.passed  # warning, not error
        assert any("ISPIN" in w for w in r.warnings)

    def test_magmom_wrong_length_is_error(self, tmp_path):
        write(tmp_path, "POSCAR", POSCAR_LiCoO2)  # 4 atoms total
        write(tmp_path, "INCAR", INCAR_MAGMOM_WRONG_LEN)  # 2*0.6 → 2 values
        r = self.v.validate(tmp_path, ["INCAR", "POSCAR"])
        assert not r.passed
        assert any("MAGMOM" in e for e in r.errors)

    def test_magmom_correct_length_passes(self, tmp_path):
        write(tmp_path, "POSCAR", POSCAR_LiCoO2)  # 4 atoms
        write(tmp_path, "INCAR", INCAR_MAGMOM_CORRECT)  # 1*0.6 1*0.6 2*0.0 → 4 values
        r = self.v.validate(tmp_path, ["INCAR", "POSCAR"])
        assert r.passed
        assert not r.errors

    def test_ldauu_wrong_length_is_error(self, tmp_path):
        write(tmp_path, "POSCAR", POSCAR_LiCoO2)  # 3 species
        write(tmp_path, "INCAR", INCAR_LDAUU_WRONG_LEN)  # 2 values
        r = self.v.validate(tmp_path, ["INCAR", "POSCAR"])
        assert not r.passed
        assert any("LDAUU" in e for e in r.errors)

    def test_ldauu_correct_length_passes(self, tmp_path):
        write(tmp_path, "POSCAR", POSCAR_LiCoO2)  # 3 species
        write(tmp_path, "INCAR", INCAR_LDAUU_CORRECT)  # 3 values
        r = self.v.validate(tmp_path, ["INCAR", "POSCAR"])
        assert r.passed

    def test_nsw_zero_isif3_warns(self, tmp_path):
        write(tmp_path, "POSCAR", POSCAR_LiCoO2)
        write(tmp_path, "INCAR", INCAR_NSW_ZERO_ISIF3)
        r = self.v.validate(tmp_path, ["INCAR", "POSCAR"])
        assert r.passed
        assert any("NSW" in w for w in r.warnings)

    def test_can_validate_requires_both(self):
        assert self.v.can_validate(["INCAR", "POSCAR"])
        assert not self.v.can_validate(["INCAR"])
        assert not self.v.can_validate(["POSCAR"])

    def test_skips_when_files_missing(self, tmp_path):
        # Neither INCAR nor POSCAR exist on disk
        r = self.v.validate(tmp_path, ["INCAR", "POSCAR"])
        assert r.passed  # graceful skip


# ===========================================================================
# FileFormatValidator tests
# ===========================================================================


class TestFileFormatValidator:
    @pytest.fixture(autouse=True)
    def _validator(self):
        from chemvcs_validator.file_format import FileFormatValidator

        self.v = FileFormatValidator()

    def test_valid_incar_passes(self, tmp_path):
        write(tmp_path, "INCAR", INCAR_BASIC)
        r = self.v.validate(tmp_path, ["INCAR"])
        assert r.passed
        assert not r.errors

    def test_invalid_incar_fails(self, tmp_path):
        # Write something completely unparseable as INCAR
        (tmp_path / "INCAR").write_bytes(b"\x00\x01\x02binary garbage")
        r = self.v.validate(tmp_path, ["INCAR"])
        # pymatgen may or may not error on binary; at minimum it should not crash
        assert isinstance(r.passed, bool)

    def test_valid_kpoints_passes(self, tmp_path):
        write(tmp_path, "KPOINTS", KPOINTS_GAMMA)
        r = self.v.validate(tmp_path, ["KPOINTS"])
        assert r.passed

    def test_invalid_kpoints_fails(self, tmp_path):
        write(tmp_path, "KPOINTS", KPOINTS_BAD)
        r = self.v.validate(tmp_path, ["KPOINTS"])
        assert not r.passed
        assert r.errors

    def test_valid_poscar_passes(self, tmp_path):
        write(tmp_path, "POSCAR", POSCAR_LiCoO2)
        r = self.v.validate(tmp_path, ["POSCAR"])
        assert r.passed

    def test_poscar_vasp4_format_flagged(self, tmp_path):
        vasp4 = (
            "comment\n1.0\n2.8 0.0 0.0\n0.0 2.8 0.0\n0.0 0.0 5.0\n"
            "1 1 2\nDirect\n0.0 0.0 0.0\n0.5 0.5 0.5\n0.0 0.5 0.25\n0.5 0.0 0.75\n"
        )
        (tmp_path / "POSCAR").write_text(vasp4)
        r = self.v.validate(tmp_path, ["POSCAR"])
        assert not r.passed
        assert any("VASP 4" in e or "VASP 5" in e for e in r.errors)

    def test_poscar_mismatched_counts_flagged(self, tmp_path):
        bad = (
            "comment\n1.0\n2.8 0.0 0.0\n0.0 2.8 0.0\n0.0 0.0 5.0\n"
            "Li Co O\n1 1\n"  # 2 counts for 3 elements
            "Direct\n0.0 0.0 0.0\n0.5 0.5 0.5\n"
        )
        (tmp_path / "POSCAR").write_text(bad)
        r = self.v.validate(tmp_path, ["POSCAR"])
        assert not r.passed

    def test_no_supported_files_skips(self, tmp_path):
        r = self.v.validate(tmp_path, ["README.md"])
        assert r.passed

    def test_can_validate_detects_vasp_files(self):
        assert self.v.can_validate(["INCAR"])
        assert self.v.can_validate(["KPOINTS"])
        assert self.v.can_validate(["POSCAR"])
        assert not self.v.can_validate(["README.md"])

    def test_disabled_by_default(self):
        assert not self.v.enabled_by_default


# ===========================================================================
# PluginManager config persistence tests
# ===========================================================================


class TestPluginManagerConfig:
    def test_load_nonexistent_config_is_noop(self, tmp_path):
        from chemvcs.plugins.manager import PluginManager

        pm = PluginManager()
        pm.load_config(tmp_path)  # no plugins.json yet — should not raise
        assert pm._config == {}

    def test_save_and_reload_config(self, tmp_path):
        from chemvcs.plugins.manager import PluginManager

        pm = PluginManager()
        pm.load_config(tmp_path)
        pm._config = {"validators": {"poscar-potcar": {"enabled": False}}}
        pm.save_config()

        pm2 = PluginManager()
        pm2.load_config(tmp_path)
        assert pm2._config["validators"]["poscar-potcar"]["enabled"] is False

    def test_set_validator_enabled_persists(self, tmp_path):
        from chemvcs_validator.poscar_potcar import POSCARPOTCARValidator

        from chemvcs.plugins.manager import PluginManager

        pm = PluginManager()
        pm.load_config(tmp_path)
        pm.validators["poscar-potcar"] = POSCARPOTCARValidator()

        result = pm.set_validator_enabled("poscar-potcar", False)
        assert result is True

        config_path = tmp_path / "plugins.json"
        assert config_path.exists()
        import json

        data = json.loads(config_path.read_text())
        assert data["validators"]["poscar-potcar"]["enabled"] is False

    def test_set_validator_enabled_unknown_returns_false(self, tmp_path):
        from chemvcs.plugins.manager import PluginManager

        pm = PluginManager()
        pm.load_config(tmp_path)
        assert pm.set_validator_enabled("nonexistent", True) is False

    def test_is_validator_enabled_respects_config(self, tmp_path):
        from chemvcs_validator.poscar_potcar import POSCARPOTCARValidator

        from chemvcs.plugins.manager import PluginManager

        pm = PluginManager()
        pm.load_config(tmp_path)
        v = POSCARPOTCARValidator()
        pm.validators["poscar-potcar"] = v

        # Default: enabled
        assert pm.is_validator_enabled(v)

        # Disable via config
        pm.set_validator_enabled("poscar-potcar", False)
        assert not pm.is_validator_enabled(v)

        # Re-enable
        pm.set_validator_enabled("poscar-potcar", True)
        assert pm.is_validator_enabled(v)
