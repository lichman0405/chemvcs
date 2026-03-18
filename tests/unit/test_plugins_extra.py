"""Extra coverage tests for plugins/manager.py.

Targets missing lines:
  60       – _discover_validators else-branch (Python < 3.10 fallback – unreachable on 3.10+;
             skipped via entry_points mock)
  72-75    – non-ValidatorPlugin warning path in _discover_validators
  89-90    – exception handler in _discover_validators
  136-141  – validator crash handler in run_validators
  158      – _is_validator_enabled explicit config lookup
  173-181  – _print_validation_result has_warnings path
  191-194  – load_config body (file exists)
  199-203  – save_config body
  218-223  – set_validator_enabled body
  227      – is_validator_enabled public wrapper
  235      – set_config
  246      – get_plugin
  254      – list_plugins
"""

import json
from pathlib import Path
from typing import Any
from unittest.mock import MagicMock, patch

import pytest

from chemvcs.plugins.base import Plugin, ValidationResult, ValidatorPlugin
from chemvcs.plugins.manager import PluginManager


# ---------------------------------------------------------------------------
# Minimal concrete plugin implementations for testing
# ---------------------------------------------------------------------------

class _TestValidator(ValidatorPlugin):
    @property
    def name(self) -> str:
        return "test-validator"

    @property
    def version(self) -> str:
        return "1.0.0"

    def validate(
        self, workspace_root: Path, files: list[str], **kwargs: Any
    ) -> ValidationResult:
        return ValidationResult(passed=True, message="OK")


class _CrashingValidator(ValidatorPlugin):
    @property
    def name(self) -> str:
        return "crashing-validator"

    @property
    def version(self) -> str:
        return "1.0.0"

    def validate(
        self, workspace_root: Path, files: list[str], **kwargs: Any
    ) -> ValidationResult:
        raise RuntimeError("Simulated crash")


class _WarnValidator(ValidatorPlugin):
    @property
    def name(self) -> str:
        return "warn-validator"

    @property
    def version(self) -> str:
        return "1.0.0"

    def validate(
        self, workspace_root: Path, files: list[str], **kwargs: Any
    ) -> ValidationResult:
        return ValidationResult(
            passed=True, message="OK with warning", warnings=["heads up"]
        )


class _FailValidator(ValidatorPlugin):
    @property
    def name(self) -> str:
        return "fail-validator"

    @property
    def version(self) -> str:
        return "1.0.0"

    def validate(
        self, workspace_root: Path, files: list[str], **kwargs: Any
    ) -> ValidationResult:
        return ValidationResult(passed=False, message="failed", errors=["err1"])


class _NonValidatorPlugin(Plugin):
    """A Plugin that is NOT a ValidatorPlugin."""

    @property
    def name(self) -> str:
        return "non-validator"

    @property
    def version(self) -> str:
        return "1.0.0"


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def chemvcs_dir(tmp_path: Path) -> Path:
    d = tmp_path / ".chemvcs"
    d.mkdir()
    return d


@pytest.fixture
def manager() -> PluginManager:
    return PluginManager()


# ---------------------------------------------------------------------------
# _discover_validators – entry point loading
# ---------------------------------------------------------------------------

class TestDiscoverPlugins:
    def test_discover_plugins_calls_discover_validators(
        self, manager: PluginManager
    ) -> None:
        """Covers the discover_plugins() top-level method (line 50)."""
        ep = MagicMock()
        ep.name = "ep"
        ep.load.return_value = _TestValidator
        with patch(
            "chemvcs.plugins.manager.metadata.entry_points",
            return_value=[ep],
        ):
            manager.discover_plugins()
        assert "test-validator" in manager.validators


class TestDiscoverValidators:
    @staticmethod
    def _fake_ep(plugin_class: type, name: str = "ep") -> MagicMock:
        ep = MagicMock()
        ep.name = name
        ep.load.return_value = plugin_class
        return ep

    def test_valid_validator_is_registered(self, manager: PluginManager) -> None:
        ep = self._fake_ep(_TestValidator)
        with patch(
            "chemvcs.plugins.manager.metadata.entry_points",
            return_value=[ep],
        ):
            manager._discover_validators()
        assert "test-validator" in manager.validators

    def test_non_validator_plugin_is_skipped(self, manager: PluginManager) -> None:
        """Non-ValidatorPlugin instance → warning printed, plugin not registered."""
        ep = self._fake_ep(_NonValidatorPlugin, "non-ep")
        with patch(
            "chemvcs.plugins.manager.metadata.entry_points",
            return_value=[ep],
        ):
            manager._discover_validators()
        assert "non-validator" not in manager.validators

    def test_entry_point_load_exception_is_caught(self, manager: PluginManager) -> None:
        """Entry point that raises on load() → warning printed, continues."""
        ep = MagicMock()
        ep.name = "bad-ep"
        ep.load.side_effect = ImportError("no module named foo")
        with patch(
            "chemvcs.plugins.manager.metadata.entry_points",
            return_value=[ep],
        ):
            manager._discover_validators()  # Must not raise
        assert len(manager.validators) == 0


# ---------------------------------------------------------------------------
# run_validators – crash handler (lines 136-141)
# ---------------------------------------------------------------------------

class TestRunValidatorsCrash:
    def test_crashing_validator_produces_failed_result(
        self, manager: PluginManager, tmp_path: Path
    ) -> None:
        manager.validators["crashing-validator"] = _CrashingValidator()
        results = manager.run_validators(tmp_path, ["INCAR"])
        assert len(results) == 1
        assert results[0].passed is False
        assert "crashed" in results[0].message.lower()

    def test_warn_validator_print_path(
        self, manager: PluginManager, tmp_path: Path
    ) -> None:
        """Covers the has_warnings branch in _print_validation_result (lines 173-181)."""
        manager.validators["warn-validator"] = _WarnValidator()
        results = manager.run_validators(tmp_path, ["INCAR"])
        assert len(results) == 1
        assert results[0].passed is True

    def test_fail_validator_print_path(
        self, manager: PluginManager, tmp_path: Path
    ) -> None:
        """Covers the passed=False branch in _print_validation_result."""
        manager.validators["fail-validator"] = _FailValidator()
        results = manager.run_validators(tmp_path, ["INCAR"])
        assert results[0].passed is False

    def test_skip_disabled_validator(
        self, manager: PluginManager, tmp_path: Path
    ) -> None:
        """Covers the skip_disabled=True continue (line 123)."""
        manager._config = {"validators": {"test-validator": {"enabled": False}}}
        manager.validators["test-validator"] = _TestValidator()
        results = manager.run_validators(tmp_path, ["INCAR"], skip_disabled=True)
        assert results == []

    def test_skip_cant_validate_files(
        self, manager: PluginManager, tmp_path: Path
    ) -> None:
        """Covers the can_validate() returning False continue (line 127)."""

        class _SelectiveValidator(_TestValidator):
            def can_validate(self, files: list[str]) -> bool:
                return False

        manager.validators["selective"] = _SelectiveValidator()
        results = manager.run_validators(tmp_path, ["INCAR"])
        assert results == []


# ---------------------------------------------------------------------------
# _is_validator_enabled – explicit config paths (line 158)
# ---------------------------------------------------------------------------

class TestIsValidatorEnabled:
    def test_explicitly_disabled(self, manager: PluginManager) -> None:
        manager._config = {"validators": {"test-validator": {"enabled": False}}}
        assert manager._is_validator_enabled(_TestValidator()) is False

    def test_explicitly_enabled_in_config(self, manager: PluginManager) -> None:
        manager._config = {"validators": {"test-validator": {"enabled": True}}}
        assert manager._is_validator_enabled(_TestValidator()) is True

    def test_not_in_config_uses_default(self, manager: PluginManager) -> None:
        assert manager._is_validator_enabled(_TestValidator()) is True


# ---------------------------------------------------------------------------
# _print_validation_result – standalone calls
# ---------------------------------------------------------------------------

class TestPrintValidationResult:
    def test_passed_no_warnings(self, manager: PluginManager) -> None:
        result = ValidationResult(passed=True, message="all good")
        manager._print_validation_result("v", result)

    def test_passed_with_warnings(self, manager: PluginManager) -> None:
        result = ValidationResult(passed=True, message="mostly OK", warnings=["note1"])
        manager._print_validation_result("v", result)

    def test_failed_with_errors(self, manager: PluginManager) -> None:
        result = ValidationResult(passed=False, message="bad", errors=["err1", "err2"])
        manager._print_validation_result("v", result)


# ---------------------------------------------------------------------------
# load_config / save_config (lines 191-203)
# ---------------------------------------------------------------------------

class TestLoadSaveConfig:
    def test_load_config_reads_existing_json(
        self, manager: PluginManager, chemvcs_dir: Path
    ) -> None:
        config = {"validators": {"x": {"enabled": False}}}
        (chemvcs_dir / "plugins.json").write_text(
            json.dumps(config), encoding="utf-8"
        )
        manager.load_config(chemvcs_dir)
        assert manager._config == config

    def test_load_config_missing_file_gives_empty(
        self, manager: PluginManager, chemvcs_dir: Path
    ) -> None:
        manager.load_config(chemvcs_dir)
        assert manager._config == {}

    def test_load_config_invalid_json_gives_empty(
        self, manager: PluginManager, chemvcs_dir: Path
    ) -> None:
        (chemvcs_dir / "plugins.json").write_text("INVALID{", encoding="utf-8")
        manager.load_config(chemvcs_dir)
        assert manager._config == {}

    def test_save_config_writes_json_file(
        self, manager: PluginManager, chemvcs_dir: Path
    ) -> None:
        manager.load_config(chemvcs_dir)
        manager._config = {"key": "value"}
        manager.save_config()
        data = json.loads(
            (chemvcs_dir / "plugins.json").read_text(encoding="utf-8")
        )
        assert data == {"key": "value"}

    def test_save_config_no_dir_is_noop(self, manager: PluginManager) -> None:
        """No _chemvcs_dir set → save_config must not raise."""
        manager.save_config()


# ---------------------------------------------------------------------------
# set_validator_enabled (lines 218-223)
# ---------------------------------------------------------------------------

class TestSetValidatorEnabled:
    def test_unknown_validator_returns_false(
        self, manager: PluginManager, chemvcs_dir: Path
    ) -> None:
        manager.load_config(chemvcs_dir)
        assert manager.set_validator_enabled("unknown", True) is False

    def test_known_validator_disables_in_config(
        self, manager: PluginManager, chemvcs_dir: Path
    ) -> None:
        manager.load_config(chemvcs_dir)
        manager.validators["test-validator"] = _TestValidator()
        result = manager.set_validator_enabled("test-validator", False)
        assert result is True
        assert manager._config["validators"]["test-validator"]["enabled"] is False


# ---------------------------------------------------------------------------
# is_validator_enabled public wrapper (line 227)
# ---------------------------------------------------------------------------

class TestIsValidatorEnabledPublic:
    def test_delegates_to_private_helper(self, manager: PluginManager) -> None:
        assert manager.is_validator_enabled(_TestValidator()) is True


# ---------------------------------------------------------------------------
# set_config (line 235)
# ---------------------------------------------------------------------------

class TestSetConfig:
    def test_set_config_replaces_config(self, manager: PluginManager) -> None:
        manager.set_config({"validators": {"x": {"enabled": False}}})
        assert manager._config["validators"]["x"]["enabled"] is False


# ---------------------------------------------------------------------------
# get_plugin / list_plugins (lines 246, 254)
# ---------------------------------------------------------------------------

class TestGetPluginListPlugins:
    def test_get_plugin_returns_registered_plugin(self, manager: PluginManager) -> None:
        v = _TestValidator()
        manager.plugins["test-validator"] = v
        assert manager.get_plugin("test-validator") is v

    def test_get_plugin_returns_none_for_unknown(self, manager: PluginManager) -> None:
        assert manager.get_plugin("nonexistent") is None

    def test_list_plugins_returns_all_info(self, manager: PluginManager) -> None:
        v = _TestValidator()
        manager.plugins["test-validator"] = v
        listing = manager.list_plugins()
        assert len(listing) == 1
        assert listing[0]["name"] == "test-validator"
        assert "version" in listing[0]
