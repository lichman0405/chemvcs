"""ChemVCS Plugin System.

Plugins can extend ChemVCS functionality through:
- Validators: Validate input file consistency
- Parsers: Add support for new file formats
- Hooks: Extend CLI commands behavior
"""

from chemvcs.plugins.base import Plugin, ValidatorPlugin
from chemvcs.plugins.manager import PluginManager

__all__ = ["Plugin", "ValidatorPlugin", "PluginManager"]
