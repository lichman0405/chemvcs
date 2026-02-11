# ChemVCS Validator Plugin

Input file validators for [ChemVCS](https://github.com/lichman0405/chemvcs).

## Features

### ✅ POSCAR-POTCAR Validator
Validates that element order in POSCAR matches POTCAR concatenation order.

**Why it matters**: VASP does not validate this! If the order is wrong, your calculation will be silently incorrect (e.g., Co and O atoms swapped in LiCoO₂).

**Example**:
```bash
$ chemvcs add POSCAR POTCAR

  ✗  poscar-potcar: POSCAR-POTCAR element order mismatch

POSCAR elements (line 6): ['Li', 'Co', 'O']
POTCAR elements (TITEL):  ['Li', 'O', 'Co']

⚠️  Same elements but WRONG ORDER!
   Position 2: POSCAR=Co, POTCAR=O ← MISMATCH
   Position 3: POSCAR=O, POTCAR=Co ← MISMATCH

How to fix:
  1. Regenerate POTCAR with correct element order:
     cat Li/POTCAR Co/POTCAR O/POTCAR > POTCAR
  2. Or reorder elements in POSCAR line 6
```

### 🚧 INCAR-POSCAR Validator (Coming Soon)
Validates INCAR parameters match POSCAR structure (MAGMOM, LDAUU, etc.).

### 🚧 File Format Validator (Coming Soon)
Validates basic VASP file format correctness.

---

## Installation

### Option 1: Install from PyPI (after publishing)
```bash
pip install chemvcs-validator
```

### Option 2: Install from source
```bash
cd plugins/chemvcs-validator
pip install -e .
```

The plugins will be **automatically discovered** by ChemVCS via entry points.

---

## Usage

### Automatic Validation
Once installed, validators run automatically during `chemvcs add`:

```bash
$ chemvcs add POSCAR POTCAR INCAR

Validating input files...
  ✓  poscar-potcar: POSCAR-POTCAR match: ['Al', 'O', 'C', 'H']
  ✓  incar-poscar: passed

Adding 3 file(s) to staging area...
```

### Disable Validation
Skip validation for a single operation:
```bash
chemvcs add POSCAR POTCAR --no-validate
```

Disable specific validator globally:
```bash
chemvcs config validator.poscar-potcar.enabled false
```

### List Installed Validators
```bash
chemvcs plugin list

Installed plugins:
  • poscar-potcar v0.1.0 - Validates POSCAR-POTCAR element order
  • incar-poscar v0.1.0 - Validates INCAR-POSCAR consistency
  • file-format v0.1.0 - Validates file format correctness
```

---

## Development

### Creating a New Validator

1. **Inherit from `ValidatorPlugin`:**

```python
from chemvcs.plugins.base import ValidatorPlugin, ValidationResult

class MyValidator(ValidatorPlugin):
    @property
    def name(self) -> str:
        return "my-validator"
    
    @property
    def version(self) -> str:
        return "1.0.0"
    
    def validate(self, workspace_root, files, **kwargs):
        # Your validation logic
        if all_good:
            return ValidationResult(passed=True, message="All good!")
        else:
            return ValidationResult(
                passed=False,
                message="Found issues",
                errors=["Error 1", "Error 2"],
            )
```

2. **Register via entry point** in `pyproject.toml`:

```toml
[project.entry-points."chemvcs.validators"]
my-validator = "my_package.my_validator:MyValidator"
```

3. **Install and test:**
```bash
pip install -e .
chemvcs add POSCAR  # Your validator runs automatically!
```

---

## Configuration

Configure validators in ChemVCS config:

```bash
# Disable a validator
chemvcs config validator.poscar-potcar.enabled false

# Change priority (lower runs first)
chemvcs config validator.poscar-potcar.priority 5

# Validator-specific settings
chemvcs config validator.poscar-potcar.strict-mode true
```

---

## License

MIT License - see [LICENSE](../../LICENSE) file for details.

---

## Contributing

Contributions welcome! Please see [CONTRIBUTING.md](../../CONTRIBUTING.md) for guidelines.

New validator ideas:
- k-point density checker
- ENCUT convergence estimator
- Symmetry preservation validator
- Memory usage estimator
