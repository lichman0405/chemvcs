# ChemVCS Command Reference

Complete reference for all ChemVCS commands, options, and usage patterns.

---

## Table of Contents

- [Global Options](#global-options)
- [chemvcs init](#chemvcs-init)
- [chemvcs add](#chemvcs-add)
- [chemvcs status](#chemvcs-status)
- [chemvcs commit](#chemvcs-commit)
- [chemvcs log](#chemvcs-log)
- [chemvcs diff](#chemvcs-diff)
- [chemvcs reproduce](#chemvcs-reproduce)
- [chemvcs version](#chemvcs-version)

---

## Global Options

These options work with any ChemVCS command:

```bash
chemvcs --help          # Show help for all commands
chemvcs <command> --help  # Show help for specific command
```

---

## chemvcs init

**Purpose**: Initialize a new ChemVCS repository in the current directory.

### Synopsis

```bash
chemvcs init [OPTIONS]
```

### Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--force` | flag | `False` | Overwrite existing `.chemvcs/` directory (**dangerous!**) |
| `--quiet` | flag | `False` | Suppress success message, only show errors |

### Examples

```bash
# Initialize repository in current directory
chemvcs init

# Reinitialize repository (deletes all history!)
chemvcs init --force

# Silent initialization (for scripts)
chemvcs init --quiet
```

### What It Creates

```
.chemvcs/
├── objects/        # Content-addressable blob store
├── commits/        # Commit objects
├── metadata.db     # SQLite database for indexing
├── index           # Staging area (JSON)
└── HEAD            # Current commit hash
```

### Notes

- Creates `.chemvcs/` directory to store all version control data
- Safe to run in existing directories (won't touch your files)
- Use `--force` carefully—it deletes all commit history

---

## chemvcs add

**Purpose**: Add files to the staging area for the next commit.

### Synopsis

```bash
chemvcs add <paths>... [OPTIONS]
```

### Arguments

| Argument | Required | Description |
|----------|----------|-------------|
| `<paths>` | Yes | One or more files or directories to add |

### Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--force` | flag | `False` | Override `.chemvcsignore` rules |
| `--dry-run` | flag | `False` | Show what would be added without modifying staging area |

### Examples

```bash
# Add single file
chemvcs add INCAR

# Add multiple files
chemvcs add INCAR POSCAR KPOINTS POTCAR

# Add all VASP input files
chemvcs add INCAR POSCAR KPOINTS POTCAR

# Add directory recursively
chemvcs add calculations/

# Preview what would be added
chemvcs add OUTCAR vasprun.xml --dry-run

# Force add ignored file
chemvcs add debug.log --force
```

### File Type Detection

ChemVCS automatically detects VASP file types:

| Filename Pattern | Detected Type |
|-----------------|---------------|
| `INCAR*` | INCAR (calculation parameters) |
| `POSCAR*`, `CONTCAR*` | POSCAR (structure) |
| `KPOINTS*` | KPOINTS (k-point sampling) |
| `POTCAR*` | POTCAR (pseudopotentials) |
| `OUTCAR*` | OUTCAR (calculation output) |
| Others | Unknown (stored as binary blob) |

### Notes

- Files are hashed and stored in `.chemvcs/objects/`
- Deduplication: identical files are stored only once
- Large files are automatically compressed
- Use `.chemvcsignore` to exclude files (similar to `.gitignore`)

---

## chemvcs status

**Purpose**: Show the current state of the working directory and staging area.

### Synopsis

```bash
chemvcs status [OPTIONS]
```

### Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--short` / `-s` | flag | `False` | Show condensed output (1 line per file) |

### Examples

```bash
# Full status display
chemvcs status

# Short format
chemvcs status --short
chemvcs status -s
```

### Output Format

#### Default (Verbose)

```
On branch: main  (linear history)
HEAD: abc1234 "SCF calculation completed"

Changes to be committed:
  (use "chemvcs commit -m <message>" to commit)

  + INCAR  (209 B, INCAR, a9b8d5a1)
  + KPOINTS  (40 B, KPOINTS, 0a7cbf0a)
```

#### Short Format (`--short`)

```
M INCAR
A KPOINTS
```

**Short format symbols**:
- `A` = Added (new file)
- `M` = Modified (file changed)

### Notes

- Shows current HEAD commit if one exists
- Lists all staged files with their sizes and hashes
- Helps verify what will be included in next commit

---

## chemvcs commit

**Purpose**: Create a new commit snapshot with all staged files.

### Synopsis

```bash
chemvcs commit [OPTIONS]
```

### Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--message` / `-m` | text | **Required** | Commit message describing this snapshot |
| `--author` | text | `$USER@$HOSTNAME` | Override default author |
| `--allow-empty` | flag | `False` | Create commit even if staging area is empty |

### Examples

```bash
# Basic commit
chemvcs commit -m "Initial VASP setup for LCO"

# Multi-line commit message
chemvcs commit -m "ENCUT convergence test

Increased ENCUT from 400 to 520 eV to test
energy convergence. K-points kept at 4x4x4."

# Custom author
chemvcs commit -m "Fixed ISMEAR" --author "jane@university.edu"

# Allow empty commit (useful for milestones)
chemvcs commit -m "Checkpoint before optimization" --allow-empty
```

### Semantic Diff Display

When committing VASP files that have changed since the parent commit, ChemVCS shows **semantic changes**:

```
Semantic Changes:

  INCAR:
    ‼️  2 critical change(s)
    ⚠️  1 major change(s)
    Key changes:
      ~ ENCUT: 400 → 520
      ~ ISMEAR: 0 → 1
      ~ LWAVE: True → False
```

**Significance Levels**:
- `‼️ critical` - Parameters that significantly affect calculation results (ENCUT, PREC, ISMEAR, SIGMA, k-point grid)
- `⚠️ major` - Important settings (LWAVE, LDAU, ALGO)
- `ℹ️ minor` - Less critical parameters (NSW, NCORE)

### Output

Each commit receives:
- SHA-256 hash (64 characters)
- Timestamp (ISO 8601)
- Author information
- Parent commit hash (if not root)

### Notes

- Commit message is **required** (no default)
- All staged files are included in commit
- Creates immutable snapshot—files cannot be modified after commit
- Staging area is cleared after successful commit

---

## chemvcs log

**Purpose**: Display commit history.

### Synopsis

```bash
chemvcs log [OPTIONS]
```

### Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--max-count` / `-n` | int | None | Limit number of commits to show |
| `--oneline` | flag | `False` | Show each commit on a single line (short hash + message) |
| `--format` | choice | `default` | Output format: `default`, `oneline`, `json` |

### Examples

```bash
# Show all commits (most recent first)
chemvcs log

# Show last 5 commits
chemvcs log -n 5
chemvcs log --max-count 5

# One-line format (short hash + message)
chemvcs log --oneline

# JSON output (for scripts/tools)
chemvcs log --format json

# Combine options
chemvcs log --oneline -n 10
```

### Output Formats

#### Default (Verbose)

```
commit d77f4cbbff22472fe29a0af7639e9d6f1f7c7e9215ca02a7212840a43b7e71c9
Parent: 3fcec3e
Author: user@hostname
Date:   2026-02-10 07:11:49

    K-point convergence: 8x8x8, ISMEAR=1, E=-10.8489 eV

commit 3fcec3e48312bcc9c8a7ce1b80b70403041d3a19ac2ba7b83c7b341060562948
Parent: 9efe155
Author: user@hostname
Date:   2026-02-10 07:11:49

    ENCUT convergence: 400->520 eV, E=-10.8452 eV
```

#### Oneline Format

```
d77f4cb K-point convergence: 8x8x8, ISMEAR=1, E=-10.8489 eV
3fcec3e ENCUT convergence: 400->520 eV, E=-10.8452 eV
9efe155 SCF completed: E=-10.8265 eV (ENCUT=400)
41c41d0 Initial setup: Si bulk SCF (ENCUT=400, K=4x4x4)
```

#### JSON Format

```json
[
  {
    "hash": "d77f4cb...",
    "parent": "3fcec3e...",
    "author": "user@hostname",
    "timestamp": "2026-02-10T07:11:49",
    "message": "K-point convergence: 8x8x8, ISMEAR=1, E=-10.8489 eV",
    "files": [...]
  }
]
```

### Notes

- Commits are shown in reverse chronological order (newest first)
- Short hashes (7 characters) are sufficient for most operations
- Root commits show `Parent: (root commit)`

---

## chemvcs diff

**Purpose**: Show semantic differences between commits or with working directory.

### Synopsis

```bash
chemvcs diff [rev1] [rev2] [OPTIONS]
```

### Arguments

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `rev1` | No | `HEAD` | Source commit hash |
| `rev2` | No | `parent of rev1` | Target commit hash |

### Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--format` | choice | `semantic` | Output format: `semantic`, `unified`, `json` |
| `--summary` | flag | `False` | Show only change summary (counts by significance) |
| `--file` | text | None | Show diff for specific file only |

### Examples

```bash
# Compare HEAD with its parent (most common)
chemvcs diff

# Compare two specific commits
chemvcs diff abc1234 def5678

# Compare with older commit
chemvcs diff HEAD~2 HEAD

# Show only summary statistics
chemvcs diff --summary

# JSON output for scripting
chemvcs diff --format json

# Diff specific file only
chemvcs diff --file INCAR

# Compact unified diff format
chemvcs diff --format unified
```

### Output Formats

#### Semantic (Default)

```
Comparing: 3fcec3e → d77f4cb

~ INCAR (modified)
  ‼️ Modified ENCUT: 400 → 520
  ‼️ Modified ISMEAR: 0 → 1
  ⚠️ Modified LWAVE: True → False
  ‼️ Modified NELM: 100 → 200

~ KPOINTS (modified)
  ‼️ Modified grid: [4, 4, 4] → [8, 8, 8]
```

**Change indicators**:
- `~` = Modified file
- `+` = Added file
- `-` = Deleted file

#### Summary Format

```
~ INCAR (modified)
  Total: 4 change(s)
  ‼️  3 critical
  ⚠️  1 major

~ KPOINTS (modified)
  Total: 1 change(s)
  ‼️  1 critical
```

#### JSON Format

```json
{
  "INCAR": {
    "status": "modified",
    "changes": [
      {
        "path": "ENCUT",
        "old_value": 400,
        "new_value": 520,
        "change_type": "modified",
        "significance": "critical"
      }
    ]
  }
}
```

### Semantic Analysis

ChemVCS performs **parameter-level analysis** for supported file types:

| File Type | Semantic Features |
|-----------|-------------------|
| **INCAR** | Parameter changes with significance levels, value comparisons |
| **KPOINTS** | Grid size changes, k-point type changes (Gamma/Monkhorst), line-mode detection |
| **POSCAR** | *(Future)* Structure changes, lattice parameter shifts |
| **OUTCAR** | *(Future)* Energy comparisons, convergence metrics |

### Notes

- Unlike text diffs, semantic diffs **understand** parameter meanings
- Critical parameters are automatically flagged (ENCUT, PREC, ISMEAR, k-grid, etc.)
- Unsupported files fall back to checksum comparison
- Working directory comparison not yet implemented (use commits only)

---

## chemvcs reproduce

**Purpose**: Extract files from a commit to reproduce a calculation exactly.

### Synopsis

```bash
chemvcs reproduce <revision> [OPTIONS]
```

### Arguments

| Argument | Required | Description |
|----------|----------|-------------|
| `<revision>` | Yes | Commit hash (full or short, min 7 chars) |

### Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--output-dir` / `-o` | text | `reproduce_<rev>/` | Output directory for extracted files |
| `--verify-potcar` | flag | `True` | Verify POTCAR file hashes (future feature) |
| `--verify-env` | flag | `True` | Compare environment info (future feature) |

### Examples

```bash
# Reproduce with default output directory
chemvcs reproduce abc1234
# → Creates reproduce_abc1234/

# Custom output directory
chemvcs reproduce abc1234 -o old_calculation

# Short hash is sufficient (min 7 chars)
chemvcs reproduce abc1234

# Use commit from log
chemvcs log --oneline
# abc1234 Initial setup
chemvcs reproduce abc1234
```

### What Gets Extracted

All files from the commit are restored to the output directory:

```
reproduce_abc1234/
├── INCAR
├── POSCAR
├── KPOINTS
├── POTCAR
└── OUTCAR  (if included in commit)
```

### Workflow Example

```bash
# Find old commit
chemvcs log --oneline
# abc1234 Converged parameters for LCO

# Reproduce it
chemvcs reproduce abc1234

# Run calculation
cd reproduce_abc1234
mpirun -np 16 vasp_std
```

### Notes

- Files are **exact byte-for-byte copies** from the commit
- Output directory must not exist (prevents accidental overwrite)
- Useful for:
  - Verifying published results
  - Rerunning calculations with different resources
  - Comparing with newer parameter sets
  - Sharing exact inputs with collaborators

---

## chemvcs version

**Purpose**: Show ChemVCS version information.

### Synopsis

```bash
chemvcs version
```

### Example

```bash
chemvcs version
# Output: ChemVCS version 0.1.0
```

---

## File Patterns & Ignore Rules

### .chemvcsignore

Similar to `.gitignore`, create a `.chemvcsignore` file to exclude files:

```
# Ignore temporary files
*.tmp
*.bak

# Ignore large output files by default
vasprun.xml
CHG
CHGCAR
WAVECAR

# Ignore job scripts
*.sh
*.pbs
*.slurm

# Ignore logs
slurm-*.out
```

**Default exclusions** (always ignored):
- `.chemvcs/`
- Hidden files (`.*)` unless explicitly added with `--force`

---

## Common Workflows

### Initial Setup

```bash
cd /path/to/calculation
chemvcs init
chemvcs add INCAR POSCAR KPOINTS POTCAR
chemvcs commit -m "Initial setup: PBE+U for LiCoO2"
```

### After Calculation Completes

```bash
chemvcs add OUTCAR
chemvcs commit -m "SCF converged: E=-123.45 eV"
```

### Parameter Tuning

```bash
# Modify INCAR (e.g., change ENCUT)
nano INCAR

chemvcs add INCAR
chemvcs commit -m "ENCUT convergence test: 400 -> 520 eV"
# → Shows semantic diff automatically
```

### Review History

```bash
chemvcs log --oneline
chemvcs diff HEAD~1 HEAD
```

### Reproduce Old Calculation

```bash
chemvcs log --oneline
# Find commit: abc1234 "Optimized structure"

chemvcs reproduce abc1234
cd reproduce_abc1234
# Ready to run
```

---

## Exit Codes

| Code | Meaning |
|------|---------|
| `0` | Success |
| `1` | Error (not initialized, file not found, etc.) |
| `2` | Invalid arguments |

---

## Environment Variables

Currently, ChemVCS does not use environment variables. Configuration is workspace-local.

---

## Tips & Best Practices

### Commit Messages

✅ **Good**:
```bash
chemvcs commit -m "ENCUT convergence: 400->520 eV, E=-10.845 eV"
chemvcs commit -m "Optimized geometry: a=3.82Å, b=3.82Å, c=13.15Å"
```

❌ **Avoid**:
```bash
chemvcs commit -m "update"
chemvcs commit -m "fix"
```

### What to Commit

**Always commit**:
- Input files (INCAR, POSCAR, KPOINTS, POTCAR)
- Key outputs (OUTCAR with final energy)

**Consider committing**:
- CONTCAR (relaxed structure)
- OSZICAR (convergence history)

**Don't commit** (unless necessary):
- vasprun.xml (large, redundant)
- WAVECAR, CHGCAR (huge files)
- Temporary files

### Commit Frequency

- **After each significant change** (parameter adjustment, structure modification)
- **After calculation completion** (capture final energy/structure)
- **Before major changes** (create checkpoint)

---

## Troubleshooting

### "Not a ChemVCS repository"

**Problem**: Running commands outside initialized directory.

**Solution**:
```bash
chemvcs init
```

### "No files staged for commit"

**Problem**: Forgot to run `chemvcs add`.

**Solution**:
```bash
chemvcs add INCAR POSCAR KPOINTS
chemvcs commit -m "message"
```

### Semantic diff not showing

**Problem**: Files not present in parent commit.

**Solution**: Ensure each commit includes all files:
```bash
chemvcs add INCAR POSCAR KPOINTS OUTCAR  # Include all files
chemvcs commit -m "message"
```

---

## See Also

- [README.md](README.md) - Project overview
- [Demo Guide](demo/GUIDE.md) - Step-by-step tutorial
- [CONTRIBUTING.md](CONTRIBUTING.md) - Development guide

---

**Note**: ChemVCS is under active development. Some features (marked "future feature") are planned but not yet implemented.
