"""Plugin manager for ChemVCS.

Discovers and loads plugins using Python entry points mechanism.
"""

import sys
from importlib import metadata
from pathlib import Path
from typing import Any, Optional

from rich.console import Console

from chemvcs.plugins.base import Plugin, ValidationResult, ValidatorPlugin

console = Console()


class PluginManager:
    """Manages ChemVCS plugins.
    
    Plugins are discovered via Python entry points:
    - chemvcs.validators: Validator plugins
    - chemvcs.parsers: Custom file parser plugins
    - chemvcs.hooks: Pre/post operation hooks
    
    Example:
        manager = PluginManager()
        manager.discover_plugins()
        
        # Run validators
        results = manager.run_validators(workspace_root, files)
    """
    
    ENTRY_POINT_GROUP_VALIDATORS = "chemvcs.validators"
    ENTRY_POINT_GROUP_PARSERS = "chemvcs.parsers"
    ENTRY_POINT_GROUP_HOOKS = "chemvcs.hooks"
    
    def __init__(self) -> None:
        """Initialize plugin manager."""
        self.plugins: dict[str, Plugin] = {}
        self.validators: dict[str, ValidatorPlugin] = {}
        self._config: dict[str, Any] = {}
    
    def discover_plugins(self) -> None:
        """Discover and load all installed plugins."""
        self._discover_validators()
        # Future: discover parsers, hooks, etc.
    
    def _discover_validators(self) -> None:
        """Discover validator plugins via entry points."""
        if sys.version_info >= (3, 10):
            # Python 3.10+ API
            entry_points = metadata.entry_points(group=self.ENTRY_POINT_GROUP_VALIDATORS)
        else:
            # Python 3.8-3.9 API
            entry_points = metadata.entry_points().get(self.ENTRY_POINT_GROUP_VALIDATORS, [])
        
        for entry_point in entry_points:
            try:
                # Load plugin class
                plugin_class = entry_point.load()
                
                # Instantiate plugin
                plugin = plugin_class()
                
                # Validate it's a ValidatorPlugin
                if not isinstance(plugin, ValidatorPlugin):
                    console.print(
                        f"[yellow]Warning: Plugin '{entry_point.name}' is not a ValidatorPlugin[/yellow]"
                    )
                    continue
                
                # Register plugin
                self.validators[plugin.name] = plugin
                self.plugins[plugin.name] = plugin
                
                # Activate plugin
                plugin.activate()
                
                console.print(
                    f"[dim]Loaded validator plugin: {plugin.name} v{plugin.version}[/dim]",
                    highlight=False,
                )
                
            except Exception as e:
                console.print(
                    f"[yellow]Warning: Failed to load plugin '{entry_point.name}': {e}[/yellow]"
                )
    
    def run_validators(
        self,
        workspace_root: Path,
        files: list[str],
        skip_disabled: bool = True,
        **kwargs: Any,
    ) -> list[ValidationResult]:
        """Run all applicable validators on specified files.
        
        Args:
            workspace_root: Root directory of the workspace
            files: List of files to validate (relative paths)
            skip_disabled: Skip validators disabled in config
            **kwargs: Additional context passed to validators
        
        Returns:
            List of ValidationResult objects
        """
        results: list[ValidationResult] = []
        
        # Sort validators by priority
        sorted_validators = sorted(
            self.validators.values(),
            key=lambda v: v.priority,
        )
        
        for validator in sorted_validators:
            # Check if validator is enabled
            if skip_disabled and not self._is_validator_enabled(validator):
                continue
            
            # Check if validator can handle these files
            if not validator.can_validate(files):
                continue
            
            try:
                result = validator.validate(workspace_root, files, **kwargs)
                results.append(result)
                
                # Print result
                self._print_validation_result(validator.name, result)
                
            except Exception as e:
                console.print(
                    f"[red]Validator '{validator.name}' failed: {e}[/red]"
                )
                # Create failed result
                results.append(
                    ValidationResult(
                        passed=False,
                        message=f"Validator crashed: {e}",
                        errors=[str(e)],
                    )
                )
        
        return results
    
    def _is_validator_enabled(self, validator: ValidatorPlugin) -> bool:
        """Check if a validator is enabled in configuration."""
        # Check global config
        validators_config = self._config.get("validators", {})
        
        # Check if explicitly disabled
        if validator.name in validators_config:
            return validators_config[validator.name].get("enabled", True)
        
        # Use default
        return validator.enabled_by_default
    
    def _print_validation_result(self, validator_name: str, result: ValidationResult) -> None:
        """Print validation result to console."""
        if result.passed:
            if result.has_warnings:
                console.print(
                    f"  [yellow]⚠[/yellow]  {validator_name}: {result.message}"
                )
                for warning in result.warnings:
                    console.print(f"      [dim]{warning}[/dim]")
            else:
                console.print(
                    f"  [green]✓[/green]  [dim]{validator_name}: passed[/dim]"
                )
        else:
            console.print(
                f"  [red]✗[/red]  {validator_name}: {result.message}"
            )
            for error in result.errors:
                console.print(f"      [red]{error}[/red]")
    
    def set_config(self, config: dict[str, Any]) -> None:
        """Set plugin configuration.
        
        Args:
            config: Configuration dictionary
        """
        self._config = config
    
    def get_plugin(self, name: str) -> Optional[Plugin]:
        """Get plugin by name.
        
        Args:
            name: Plugin name
            
        Returns:
            Plugin instance or None if not found
        """
        return self.plugins.get(name)
    
    def list_plugins(self) -> list[dict[str, Any]]:
        """List all loaded plugins.
        
        Returns:
            List of plugin info dictionaries
        """
        return [
            {
                "name": plugin.name,
                "version": plugin.version,
                "description": plugin.description,
                "type": plugin.__class__.__base__.__name__,
            }
            for plugin in self.plugins.values()
        ]
