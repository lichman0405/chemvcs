# ChemVCS Development Guidelines

> v1.0 | 2026-02-09  
> Scope: All contributors (core team + community)

---

## 1. Development Environment Setup

### 1.1 Prerequisites

- Python ≥3.8 (recommended 3.10 or 3.11)
- Git ≥2.30
- Operating System: Linux / macOS / Windows (WSL2 or PowerShell)

### 1.2 Initial Setup

```bash
# 1. Fork the repository (if you're an external contributor)
# Visit https://github.com/lichman0405/chemvcs and fork

# 2. Clone to local
git clone https://github.com/<your-username>/chemvcs.git
cd chemvcs

# 3. Add upstream remote
git remote add upstream https://github.com/lichman0405/chemvcs.git

# 4. Create virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate

# 5. Install development dependencies
pip install -e ".[dev]"

# 6. Install pre-commit hooks
pre-commit install

# 7. Verify installation
chemvcs --version
pytest --version
ruff --version
```

### 1.3 Development Dependencies

```toml
# dev dependencies in pyproject.toml
[project.optional-dependencies]
dev = [
    "pytest >=7.4.0",
    "pytest-cov >=4.1.0",
    "pytest-timeout >=2.1.0",
    "ruff >=0.1.9",
    "mypy >=1.8.0",
    "pre-commit >=3.5.0",
    "ipython >=8.12.0",  # REPL debugging
    "hypothesis >=6.92.0",  # property testing
]
```

---

## 2. Code Style

### 2.1 Automated Tool Configuration

#### Ruff (Linter + Formatter)

```toml
# pyproject.toml
[tool.ruff]
target-version = "py38"
line-length = 100
indent-width = 4

[tool.ruff.lint]
select = [
    "E",   # pycodestyle errors
    "W",   # pycodestyle warnings
    "F",   # pyflakes
    "I",   # isort
    "B",   # flake8-bugbear
    "C4",  # flake8-comprehensions
    "UP",  # pyupgrade
    "SIM", # flake8-simplify
]
ignore = [
    "E501",  # line too long (handled by formatter)
    "B008",  # do not perform function calls in argument defaults
]

[tool.ruff.format]
quote-style = "double"
indent-style = "space"
skip-magic-trailing-comma = false

[tool.ruff.lint.isort]
known-first-party = ["chemvcs"]
```

**Usage:**

```bash
# Check code
ruff check chemvcs/

# Auto-fix
ruff check --fix chemvcs/

# Format
ruff format chemvcs/
```

#### MyPy (Type Checking)

```toml
# pyproject.toml
[tool.mypy]
python_version = "3.8"
warn_return_any = true
warn_unused_configs = true
disallow_untyped_defs = true
disallow_any_generics = true
check_untyped_defs = true
no_implicit_optional = true
warn_redundant_casts = true
warn_unused_ignores = true
warn_no_return = true
strict_equality = true

[[tool.mypy.overrides]]
module = "pymatgen.*"
ignore_missing_imports = true
```

**Usage:**

```bash
mypy chemvcs/
```

### 2.2 Naming Conventions

| Type | Rule | Example |
|------|------|---------|
| **Modules** | lowercase + underscore | `object_store.py`, `incar_parser.py` |
| **Classes** | PascalCase | `ObjectStore`, `IncarParser` |
| **Functions/Methods** | lowercase + underscore | `write_blob()`, `parse_incar()` |
| **Constants** | uppercase + underscore | `DEFAULT_IGNORE_PATTERNS`, `MAX_FILE_SIZE` |
| **Private members** | prefix `_` | `_compute_hash()`, `_internal_state` |
| **Type variables** | PascalCase + `T` | `FileT`, `CommitT` |

### 2.3 Type Annotations

**All public APIs must include type annotations:**

```python
# ✅ Correct
def write_blob(content: bytes, objects_dir: Path) -> str:
    """Write blob and return SHA-256 hash.
    
    Args:
        content: File content
        objects_dir: Path to objects/ directory
    
    Returns:
        40-character SHA-256 hex string
    
    Raises:
        OSError: Write failed (permission/disk full)
    """
    ...

# ❌ Wrong (missing type annotations)
def write_blob(content, objects_dir):
    ...
```

**Use `typing` module:**

```python
from typing import Optional, Dict, List, Tuple, Union
from pathlib import Path

def parse_incar(path: Path) -> Dict[str, Union[int, float, bool, str, List[float]]]:
    ...

def get_commit(commit_id: str) -> Optional[Dict]:
    """Returns None when not found"""
    ...
```

### 2.4 Docstring Guidelines (Google Style)

```python
class IncarParser:
    """Semantic parser for INCAR files.
    
    Uses pymatgen as primary parser, falls back to regex on failure.
    
    Attributes:
        fallback_enabled: Whether to enable regex fallback mode
        param_meta: Parameter metadata (types, units, etc.)
    
    Example:
        >>> parser = IncarParser()
        >>> params = parser.parse(Path("INCAR"))
        >>> params["ENCUT"]
        520.0
    """
    
    def parse(self, path: Path) -> Dict[str, Any]:
        """Parse INCAR file.
        
        Args:
            path: Path to INCAR file
        
        Returns:
            Parameter dictionary with normalized values
        
        Raises:
            FileNotFoundError: File does not exist
            ValueError: File content cannot be parsed
        """
        ...
```

---

## 3. Branching Strategy

### 3.1 Branch Model

```
main (protected)          ← Stable version, tag on each release
  ↑
  merge ← dev             ← Development branch, integrates all features
           ↑
           merge ← feature/issue-123-incar-parser
                ← feature/issue-124-diff-poscar
                ← fix/issue-125-potcar-hash
```

### 3.2 Branch Naming Rules

| Type | Format | Example |
|------|--------|---------|
| Feature | `feature/<issue-number>-<short-description>` | `feature/42-semantic-diff` |
| Bug fix | `fix/<issue-number>-<short-description>` | `fix/55-potcar-hash-mismatch` |
| Documentation | `docs/<description>` | `docs/update-cli-spec` |
| Refactoring | `refactor/<description>` | `refactor/storage-layer` |
| Release prep | `release/v<version>` | `release/v0.1.0` |

### 3.3 Workflow

```bash
# 1. Update local dev branch
git checkout dev
git pull upstream dev

# 2. Create feature branch (from dev)
git checkout -b feature/42-semantic-diff

# 3. Develop + commit
# ... coding ...
git add chemvcs/parsers/incar_diff.py
git commit -m "feat: implement INCAR semantic diff

- Add IncarDiffer class
- Support numeric delta calculation
- Add unit tests (coverage: 92%)

Closes #42"

# 4. Local testing
pytest tests/unit/parsers/
ruff check chemvcs/
mypy chemvcs/

# 5. Push to fork
git push origin feature/42-semantic-diff

# 6. Create PR on GitHub (target: upstream/dev)
```

---

## 4. Commit Guidelines

### 4.1 Commit Message Format (Conventional Commits)

```
<type>(<scope>): <subject>

<body>

<footer>
```

**Type (required):**

| Type | Description | Example |
|------|-------------|---------|
| `feat` | New feature | `feat(cli): add chemvcs reproduce command` |
| `fix` | Bug fix | `fix(storage): handle disk full error in write_blob` |
| `docs` | Documentation | `docs: update CLI_SPEC with --format option` |
| `test` | Testing | `test(parsers): add edge case for INCAR comments` |
| `refactor` | Refactoring | `refactor(diff): extract common diff logic` |
| `perf` | Performance | `perf(commit): cache blob hash to avoid recalc` |
| `chore` | Build/tools | `chore: upgrade ruff to 0.2.0` |

**Scope (optional):** Affected module, e.g., `cli`, `storage`, `parsers`, `diff`

**Subject:** Short description (≤50 chars), use imperative mood (add, fix, update)

**Body (optional):** Detailed explanation (wrap at 72 chars)

**Footer (optional):** Issue references, Breaking Changes

**Example:**

```
feat(parsers): add KPOINTS parser with grid mode support

- Implement KpointsParser class using pymatgen
- Support Automatic (Gamma/Monkhorst-Pack) mode
- Add diff logic for grid size comparison
- Test coverage: 88%

Closes #42
```

### 4.2 Commit Atomicity

- Each commit should be **independently testable** (doesn't break CI)
- Single responsibility: one commit does one thing (implement feature vs refactor vs fix bug)
- Avoid "WIP", "fix typo" commits (squash locally before pushing)

---

## 5. Pull Request Process

### 5.1 PR Checklist

**Before submitting PR, verify:**

- [ ] All tests pass (`pytest tests/`)
- [ ] Code style check passes (`ruff check`)
- [ ] Type check passes (`mypy chemvcs/`)
- [ ] New code has unit tests (coverage ≥85%)
- [ ] Documentation updated (if API changed)
- [ ] PR description is clear and links to issue (`Closes #123`)

### 5.2 PR Template

```markdown
## Description
<!-- Brief description of what this PR does -->

Closes #<issue-number>

## Changes
<!-- List main changes -->

- Added X
- Fixed Y
- Refactored Z

## Testing
<!-- How were these changes tested? -->

- [ ] Unit tests added (coverage: __%)
- [ ] Integration test: <scenario>
- [ ] Manual testing on HPC: <cluster-name>

## Checklist
- [ ] Tests pass locally
- [ ] Ruff/mypy checks pass
- [ ] Documentation updated
- [ ] Breaking changes noted (if any)

## Screenshots/Output (if applicable)
```

### 5.3 Code Review Requirements

**All PRs must be reviewed by ≥1 person before merging.**

**Reviewer checklist:**

- [ ] Code logic is correct
- [ ] Edge cases handled (empty files, large files, concurrency, etc.)
- [ ] Error handling is comprehensive (meaningful error messages)
- [ ] No significant performance regression
- [ ] Tests cover critical paths
- [ ] Code readability (naming, comments)

**PR feedback response time:**

- Core team: First response within 24 hours
- Community contributors: First response within 48 hours

---

## 6. Testing Requirements

### 6.1 Unit Test Coverage

- **Minimum requirement**: New code coverage ≥85%
- **Target**: Core modules coverage ≥90% (`storage/`, `parsers/`, `cli/`)

### 6.2 Test File Organization

```
tests/
├── unit/
│   ├── storage/
│   │   ├── test_object_store.py
│   │   ├── test_metadata_db.py
│   │   └── test_commit_builder.py
│   ├── parsers/
│   │   ├── test_incar_parser.py
│   │   ├── test_poscar_parser.py
│   │   └── test_outcar_extractor.py
│   └── cli/
│       ├── test_init.py
│       ├── test_add.py
│       └── test_commit.py
├── integration/
│   ├── test_init_workflow.py
│   ├── test_commit_workflow.py
│   └── test_reproduce_workflow.py
├── fixtures/
│   └── ...
└── conftest.py
```

### 6.3 Test Naming Rules

```python
# Test function naming: test_<function>_<scenario>_<expected>
def test_write_blob_deduplication_skips_existing():
    """Test write_blob skips writing when content already exists"""
    ...

def test_parse_incar_with_comments_ignores_them():
    """Test parse_incar correctly ignores comments"""
    ...
```

### 6.4 Running Tests

```bash
# Run all tests
pytest

# Run specific module
pytest tests/unit/parsers/

# With coverage report
pytest --cov=chemvcs --cov-report=html

# Fast mode (skip slow tests)
pytest -m "not slow"

# Verbose mode (show full traceback on failure)
pytest -vv
```

---

## 7. Documentation Guidelines

### 7.1 Documentation Types

| Type | Location | Audience |
|------|----------|----------|
| User docs | `README.md`, `COMMANDS.md` | End users |
| API docs | Docstrings (auto-generated) | Developers |
| Architecture docs | `docs/TECH_ARCHITECTURE.md` | Contributors |
| Development guidelines | `CONTRIBUTING.md` (this file) | Contributors |

### 7.2 When to Update Documentation

**Must synchronize documentation when:**

- Adding CLI command → Update `COMMANDS.md`
- Modifying exit codes → Update exit code table in `COMMANDS.md`
- Adding storage fields → Update schema in `TECH_ARCHITECTURE.md`
- Breaking change → Mark specially in `CHANGELOG.md`

---

## 8. Release Process

### 8.1 Version Numbering (Semantic Versioning)

```
v<MAJOR>.<MINOR>.<PATCH>

v0.1.0 → Initial MVP release
v0.1.1 → Bug fixes
v0.2.0 → New HPC integration features
v1.0.0 → Stable version (production-ready)
```

### 8.2 Release Checklist

```markdown
- [ ] All CI passes (dev branch)
- [ ] Version number updated (`chemvcs/__init__.py`, `pyproject.toml`)
- [ ] CHANGELOG.md updated
- [ ] Create release/v<version> branch
- [ ] Manual testing on HPC environment (at least 2 filesystems)
- [ ] Merge release branch to main
- [ ] Tag on main: `git tag -a v0.1.0 -m "Release v0.1.0"`
- [ ] Push tag: `git push upstream v0.1.0`
- [ ] Build and upload to PyPI: `python -m build && twine upload dist/*`
- [ ] Create GitHub Release (with changelog)
- [ ] Announcement (mailing list / forum / Twitter)
```

---

## 9. Issue Tracking

### 9.1 Issue Templates

#### Bug Report

```markdown
**Environment:**
- ChemVCS version: 
- Python version: 
- OS: 
- File system (if HPC): 

**Steps to Reproduce:**
1. 
2. 

**Expected Behavior:**

**Actual Behavior:**

**Error Message/Traceback:**
```

#### Feature Request

```markdown
**Problem Statement:**
<!-- Describe the problem you want to solve -->

**Proposed Solution:**
<!-- Your suggested solution -->

**Alternatives Considered:**
<!-- Other approaches -->

**Additional Context:**
<!-- Any other information -->
```

### 9.2 Issue Labels

| Label | Description |
|-------|-------------|
| `bug` | Bug report |
| `enhancement` | New feature request |
| `documentation` | Documentation improvement |
| `good first issue` | Suitable for new contributors |
| `help wanted` | Community help needed |
| `priority: high` | High priority |
| `wontfix` | Will not be fixed |

---

## 10. Community Code of Conduct

### 10.1 Core Principles

- **Friendly & Respectful**: Welcome contributors of all skill levels
- **Constructive Feedback**: Critique code, not people; suggest solutions when raising issues
- **Patience**: Community members may be from different timezones and backgrounds
- **Inclusivity**: Avoid exclusionary language (e.g., "obviously", "everyone knows")

### 10.2 Unacceptable Behavior

The following behaviors are not tolerated:

- Personal attacks, insults, derogatory comments
- Publishing others' private information without permission
- Harassment (public or private)
- Other behavior violating professional ethics

**Reporting mechanism:** Email shadow.li981@gmail.com (handled confidentially)

---

## 11. FAQ for Contributors

### Q1: I'm a beginner, where should I start?

**A:** Look for issues labeled `good first issue`. These are typically independent and well-scoped tasks.

### Q2: Do I need to test in a real HPC environment?

**A:** For core functionality (blob storage, commit), local testing + unit tests are sufficient. HPC-specific features (e.g., SLURM integration) can be noted as "not HPC-tested" in PR, and maintainers will arrange testing.

### Q3: My PR was requested changes, but I disagree. What should I do?

**A:** Politely explain your reasoning in PR comments. If disagreement persists, tag `@chemvcs/maintainers` for a third-party opinion.

### Q4: How do I set up remote debugging (on HPC)?

**A:** Use `chemvcs --debug` to enable verbose logging. Or insert `import pdb; pdb.set_trace()` breakpoints in code (requires interactive shell).

### Q5: Tests fail but I can't reproduce locally. What now?

**A:** Check if it's OS-related (Linux vs Windows paths). Mark PR with `needs-investigation`. CI logs usually have detailed information.

---

## 12. Contact

- **GitHub Issues**: <https://github.com/lichman0405/chemvcs/issues>
- **Discussion Forum**: <https://github.com/lichman0405/chemvcs/discussions>
- **Email**: shadow.li981@gmail.com (prioritize GitHub for technical issues)

---

**Thank you for your contribution!** 🎉

All contributors will be listed in `CONTRIBUTORS.md`. We're committed to building a friendly and efficient open-source community for computational materials science.
