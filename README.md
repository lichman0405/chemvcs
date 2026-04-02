# ChemVCS

**Version Control for Computational Materials Science**

[![Python Version](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Tests](https://github.com/lichman0405/chemvcs/actions/workflows/test.yml/badge.svg)](https://github.com/lichman0405/chemvcs/actions/workflows/test.yml)
[![Lint](https://github.com/lichman0405/chemvcs/actions/workflows/lint.yml/badge.svg)](https://github.com/lichman0405/chemvcs/actions/workflows/lint.yml)
[![Code style: ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

ChemVCS provides Git-like version control specifically designed for computational materials science calculations, with semantic understanding of VASP, LAMMPS, and ORCA files.

## 🎯 Features (MVP)

- **Git-like Interface**: Familiar `init`, `add`, `commit`, `log`, `diff`, `reproduce` commands
- **VASP-Aware**: Semantic parsing of INCAR, KPOINTS, and OUTCAR files
- **LAMMPS-Aware**: Semantic parsing of input scripts, data files, and log files
- **ORCA-Aware**: Semantic parsing of ORCA input (`.inp`) and output (`.out`) files
- **Semantic Diff**: Human-readable comparison of calculation parameters — **semantic equality is the primary criterion**: cosmetic edits (whitespace, comments) do not trigger a "modified" report
  - **INCAR**: Automatic detection of critical (ENCUT, PREC, ISMEAR), major (LWAVE, LDAU), and minor parameter changes
  - **KPOINTS**: K-point grid changes, sampling type detection (Gamma/Monkhorst/Line-mode)
  - **LAMMPS input**: Thermostat/barostat, timestep, pair style, run/minimize settings
  - **LAMMPS data**: Atom counts, topology counts, box dimensions, masses
  - **LAMMPS log**: Final energies, temperature, pressure, completion status, wall time
  - **ORCA input**: Method, basis set, run type, charge, multiplicity, dispersion, RI, block settings
  - **ORCA output**: Final energy (Hartrees), termination status, SCF convergence, optimisation result
  - Change significance indicators: ‼️ critical, ⚠️ major, ℹ️ minor
- **Plugin System**: Extensible validation framework
  - **POSCAR-POTCAR validator**: Automatic element order consistency checking
  - **Pluggable architecture**: Easy to add custom validators for domain-specific checks
  - **Optional validation**: Use `--no-validate` flag to skip checks when needed
- **Reproducibility**: Exact restoration of calculation inputs from any commit
- **HPC-Native**: Zero infrastructure dependencies, works on any filesystem
- **Lightweight**: Content-addressable storage with automatic deduplication and compression

## 🚀 Quick Start

### Installation

```bash
pip install chemvcs
```

### Basic Usage

```bash
# Initialize a repository
cd /scratch/username/LCO_doping
chemvcs init

# Add files to staging
chemvcs add INCAR POSCAR KPOINTS POTCAR  # Automatic POSCAR-POTCAR validation

# Commit snapshot
chemvcs commit -m "PBE+U, U=3.5eV, kpoints 4x4x4"

# View plugin status
chemvcs plugin list

# Compare working tree against HEAD (default)
chemvcs diff

# Compare two specific commits
chemvcs diff <hash1> <hash2>

# Reproduce previous calculation
chemvcs reproduce abc1234
cd reproduce_abc1234
mpirun vasp_std  # Ready to run
```

### LAMMPS Workflow

```bash
# Initialize a repository for a LAMMPS run directory
cd /scratch/username/lj_fluid
chemvcs init

# Stage LAMMPS inputs
chemvcs add in.lammps data.lammps
chemvcs commit -m "Initial LJ fluid equilibration"

# After the run, track the primary thermo log
chemvcs add log.lammps
chemvcs commit -m "Completed NVT run at T*=1.0"

# Compare parameter or thermodynamic changes
chemvcs diff --summary

# Reproduce the latest tracked snapshot
chemvcs reproduce HEAD
```

## 📖 Documentation

- **[Command Reference](COMMANDS.md)** - Complete guide to all commands and options
- **[Demo Guide](demo/GUIDE.md)** - Step-by-step tutorial with Si convergence test example
- **[Advanced Demo](demo_advanced/GUIDE_ADVANCED.md)** - MOF-acetylene adsorption workflow with plugin validation
- **[LAMMPS Demo](demo_lammps/GUIDE_LAMMPS.md)** - Step-by-step LAMMPS workflow with semantic diffs
- **[ORCA Demo](demo_orca/GUIDE_ORCA.md)** - ORCA quantum chemistry workflow (DFT opt → CCSD(T) single point)
- **[Plugin README](plugins/chemvcs-validator/README.md)** - Built-in validator plugin overview
- [Contributing](CONTRIBUTING.md) - Development guide and coding standards

## 🛠️ Development

### Setup

```bash
# Clone repository
git clone https://github.com/lichman0405/chemvcs.git
cd chemvcs

# Create virtual environment
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate

# Install in development mode
pip install -e ".[dev]"

# Install pre-commit hooks
pre-commit install
```

### Running Tests

```bash
# Run all tests
pytest

# With coverage
pytest --cov=chemvcs --cov-report=html

# Run specific test module
pytest tests/unit/parsers/

# Run linting
ruff check chemvcs/

# Type checking
mypy chemvcs/
```

## 🗺️ Roadmap

### Phase 1: Core Version Control ✅ **COMPLETE**
- ✅ Project setup and scaffolding
- ✅ Storage layer (content-addressable blob store, SQLite metadata DB)
- ✅ Semantic diff engine for VASP, LAMMPS, and ORCA file types
- ✅ ORCA support: `OrcaInputParser` (pure regex, `.inp`) and `OrcaOutputParser` (cclib + regex fallback, `.out`)
- ✅ Plugin system for extensible validation (POSCAR-POTCAR validator included)
- ✅ Comprehensive test suite (520+ tests, 87%+ coverage)
  - ✅ Interactive demo and documentation
  - ✅ Real working-tree semantic diff (cosmetic-only changes correctly ignored)

### Phase 1.5: Plugin Ecosystem ✅ **COMPLETE**
- ✅ POSCAR-POTCAR element order validator
- ✅ INCAR-POSCAR consistency validator (MAGMOM/LDAUU length, ISPIN vs magnetism, NSW vs ISIF)
- ✅ File format validator (POSCAR structure, INCAR/KPOINTS syntax)
- ✅ Configuration system for plugin enable/disable (`chemvcs plugin enable/disable <name>`)

### Phase 2: HPC Integration 🔲 PLANNED
- SLURM job tracking
- Automatic commit on job completion
- Resource usage statistics

### Phase 3: Collaboration 🔲 PLANNED
- Remote repositories
- Push/pull/clone
- Conflict resolution
- Community plugin repository

## 🤝 Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for:

- Code style guidelines
- Development workflow
- Testing requirements
- Pull request process

## 📄 License

MIT License - see [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

This project builds on:

- [pymatgen](https://pymatgen.org/) - Materials analysis framework (VASP output parsing)
- [cclib](https://cclib.github.io/) - Computational chemistry logfile parser (ORCA output parsing)
- [Typer](https://typer.tiangolo.com/) - CLI framework
- [Rich](https://rich.readthedocs.io/) - Terminal formatting

Inspired by version control needs in computational materials science community.

## 📧 Contact

- **Issues**: [GitHub Issues](https://github.com/lichman0405/chemvcs/issues)
- **Discussions**: [GitHub Discussions](https://github.com/lichman0405/chemvcs/discussions)

---

**Note**: This project is in active development (v0.1.0). APIs may change before v1.0.0.
