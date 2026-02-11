"""Quick test of the plugin system with real POSCAR/POTCAR files."""

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).parent / "chemvcs-validator" / "src"))

from chemvcs.plugins.manager import PluginManager
from chemvcs_validator.poscar_potcar import POSCARPOTCARValidator

print("=" * 60)
print("ChemVCS Plugin System Demo")
print("=" * 60)

# Test 1: Direct validator usage
print("\n[Test 1] Direct Validator Usage")
print("-" * 60)

validator = POSCARPOTCARValidator()
print(f"✓ Loaded validator: {validator.name} v{validator.version}")
print(f"  Description: {validator.description}")
print(f"  Priority: {validator.priority}")

# Test with real files
workspace_root = Path(__file__).parent.parent / "sample-pospot"
files = ["POSCAR", "POTCAR"]

print(f"\nValidating files in: {workspace_root}")
print(f"Files: {files}")

result = validator.validate(workspace_root, files)

print(f"\nValidation Result:")
print(f"  Passed: {result.passed}")
print(f"  Message: {result.message}")
if result.errors:
    print(f"  Errors:")
    for error in result.errors:
        print(f"    - {error}")
if result.warnings:
    print(f"  Warnings:")
    for warning in result.warnings:
        print(f"    - {warning}")

# Test 2: Plugin Manager (manual loading)
print("\n\n[Test 2] Plugin Manager (Manual Registration)")
print("-" * 60)

manager = PluginManager()

# Manually register validator (simulating entry point discovery)
manager.validators[validator.name] = validator
manager.plugins[validator.name] = validator

print(f"✓ Registered {len(manager.validators)} validator(s)")

# List plugins
plugins = manager.list_plugins()
print("\nInstalled Plugins:")
for p in plugins:
    print(f"  • {p['name']} v{p['version']} - {p['description']}")

# Run all validators
print(f"\nRunning validators on: {workspace_root}")
results = manager.run_validators(workspace_root, files)

print(f"\nValidation Summary:")
print(f"  Total validators run: {len(results)}")
print(f"  Passed: {sum(1 for r in results if r.passed)}")
print(f"  Failed: {sum(1 for r in results if not r.passed)}")

# Test 3: Test with mismatched POSCAR/POTCAR (if available)
print("\n\n[Test 3] Test Mismatch Detection (Demo)")
print("-" * 60)
print("Creating a mock mismatch scenario...")

# Show what error message would look like
mock_poscar_elements = ["Li", "Co", "O"]
mock_potcar_elements = ["Li", "O", "Co"]

print(f"\nScenario:")
print(f"  POSCAR elements: {mock_poscar_elements}")
print(f"  POTCAR elements: {mock_potcar_elements}")

# Format error (using validator's private method)
error_msg = validator._format_mismatch_error(
    mock_poscar_elements,
    mock_potcar_elements,
)
print(f"\nError Message Would Be:")
print(error_msg)

print("\n" + "=" * 60)
print("Demo Complete!")
print("=" * 60)
