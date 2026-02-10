# ChemVCS Demo - Manual Walkthrough Guide

> **Scenario**: Si (Silicon) bulk DFT convergence test  
> **Objective**: Step-by-step demonstration of all ChemVCS MVP features

---

## Prerequisites

```powershell
# Create a temporary demo directory
mkdir $env:TEMP\si_convergence
cd $env:TEMP\si_convergence
```

> 💡 All commands below are executed in this directory.

---

## Step 1: Initialize Repository

**Narration**: ChemVCS provides version control for computational chemistry projects. First step is to initialize a repository.

```powershell
chemvcs init
```

**Expected Output**: Shows repository initialized successfully, including `.chemvcs/` directory path and usage hints.

---

## Step 2: Create Initial VASP Input Files

**Narration**: We'll perform a DFT calculation on Si bulk. Prepare three standard VASP input files:
- **POSCAR** — Silicon diamond structure (lattice constant 5.43 Å)
- **INCAR** — Calculation parameters (ENCUT=400 eV, ISMEAR=0)
- **KPOINTS** — 4×4×4 Gamma k-point grid

```powershell
# Copy from pre-made files
# ⚠️ Replace the path below with your actual chemvcs project path
$DEMO = "C:\Users\lishi\code\chemvcs\demo\vasp_files"

Copy-Item "$DEMO\step1_initial\POSCAR" .
Copy-Item "$DEMO\step1_initial\INCAR"  .
Copy-Item "$DEMO\step1_initial\KPOINTS" .
```

Check file contents:

```powershell
type INCAR
```

**Expected Output**:
```
# Si bulk - initial SCF calculation
SYSTEM  = Si-diamond
PREC    = Accurate
ENCUT   = 400
ISMEAR  = 0
SIGMA   = 0.05
...
```

---

## Step 3: Add Files to Staging Area

**Narration**: `chemvcs add` computes file content hashes and identifies VASP file types.

```powershell
chemvcs add POSCAR INCAR KPOINTS
```

**Expected Output**: Shows 3 files added, each displaying size and type (POSCAR / INCAR / KPOINTS).

---

## Step 4: Check Repository Status

**Narration**: Similar to `git status`, displays which files are staged and ready to commit.

```powershell
chemvcs status
```

**Expected Output**: Shows 3 files to be committed with their hash values.

---

## Step 5: Create Initial Commit

**Narration**: Save the current file snapshot to version history.

```powershell
chemvcs commit -m "Initial setup: Si bulk SCF (ENCUT=400, K=4x4x4)"
```

**Expected Output**: Shows commit hash, author, timestamp. Note this is a root commit (no parent).

---

## Step 6: Simulate Calculation Completion — Add OUTCAR

**Narration**: Assume VASP calculation finished, producing OUTCAR. Total energy E = **-10.8265 eV**. We'll commit all files together to save the complete calculation snapshot.

```powershell
Copy-Item "$DEMO\step2_scf_done\OUTCAR" .

# Check energy result
Select-String "energy.*sigma" OUTCAR
```

```powershell
chemvcs add POSCAR INCAR KPOINTS OUTCAR
chemvcs commit -m "SCF completed: E=-10.8265 eV (ENCUT=400)"
```

**Expected Output**: 4 files committed successfully.

---

## Step 7: ENCUT Convergence Test (400 → 520 eV) ⭐

**Narration**: Energy cutoff ENCUT is the most critical parameter in DFT. We'll increase ENCUT from 400 to 520 eV to check energy convergence.

```powershell
# Update INCAR and OUTCAR
Copy-Item "$DEMO\step3_encut_conv\INCAR"  . -Force
Copy-Item "$DEMO\step3_encut_conv\OUTCAR" . -Force

# Check INCAR changes
type INCAR
```

```powershell
chemvcs add POSCAR INCAR KPOINTS OUTCAR
```

```powershell
chemvcs commit -m "ENCUT convergence: 400->520 eV, E=-10.8452 eV"
```

**⭐ Key Observation**: The commit output will show a **Semantic Changes** section:
```
Semantic Changes:

  INCAR:
    ‼️  1 critical change(s)
    Key changes:
      ~ ENCUT: 400 → 520
```

**Key Talking Point**: ChemVCS automatically detected that ENCUT is a **critical**-level parameter change! This is far more meaningful than plain text diffs from `git diff`.

---

## Step 8: K-point Convergence Test (4×4×4 → 8×8×8) ⭐⭐

**Narration**: Next, test k-point grid density. Also adjust ISMEAR (0→1) and SIGMA (0.05→0.1).

```powershell
# Update INCAR, KPOINTS, OUTCAR
Copy-Item "$DEMO\step4_kpoint_conv\INCAR"   . -Force
Copy-Item "$DEMO\step4_kpoint_conv\KPOINTS" . -Force
Copy-Item "$DEMO\step4_kpoint_conv\OUTCAR"  . -Force
```

```powershell
chemvcs add POSCAR INCAR KPOINTS OUTCAR
```

```powershell
chemvcs commit -m "K-point convergence: 8x8x8, ISMEAR=1, E=-10.8489 eV"
```

**⭐⭐ Key Observation**: Detects semantic changes in both INCAR and KPOINTS:
```
Semantic Changes:

  INCAR:
    ‼️  5 critical change(s)
    ⚠️  1 major change(s)
    Key changes:
      ~ ISMEAR: 0 → 1
      ~ SIGMA: 0.05 → 0.1
      ~ EDIFF: 1e-06 → 1e-07
    ... and 3 more change(s)

  KPOINTS:
    ‼️  1 critical change(s)
    Key changes:
      ~ grid: [4, 4, 4] → [8, 8, 8]
```

**Key Talking Points**:
- INCAR shows 5 critical + 1 major changes
- KPOINTS k-point grid change automatically marked as critical
- Helps researchers quickly assess: this modification will significantly affect calculation results

---

## Step 9: View Commit History

**Narration**: View complete calculation version history, documenting each parameter adjustment.

```powershell
chemvcs log
```

Compact mode:

```powershell
chemvcs log --oneline
```

**Expected Output** (similar to):
```
d77f4cb K-point convergence: 8x8x8, ISMEAR=1, E=-10.8489 eV
3fcec3e ENCUT convergence: 400->520 eV, E=-10.8452 eV
9efe155 SCF completed: E=-10.8265 eV (ENCUT=400)
41c41d0 Initial setup: Si bulk SCF (ENCUT=400, K=4x4x4)
```

**Key Talking Point**: Energy converges progressively from -10.8265 → -10.8452 → -10.8489 eV.

---

## Step 10: Semantic Diff ⭐⭐⭐

**Narration**: This is ChemVCS's core feature — semantic-level comparison of VASP files.

### 10a. Detailed diff (default)

```powershell
chemvcs diff
```

**Expected Output**: Shows complete semantic differences between HEAD and previous commit, including ‼️ critical / ⚠️ major markers for each parameter.

### 10b. Summary mode

```powershell
chemvcs diff --summary
```

**Expected Output**: Only displays change count statistics.

### 10c. JSON output (machine-readable)

```powershell
chemvcs diff --format json
```

**Key Talking Point**: JSON format can be parsed by other tools, facilitating automated analysis of parameter changes in convergence tests.

---

## Step 11: Reproduce Historical Calculation ⭐⭐⭐

**Narration**: One of the most critical features — given a commit hash, completely reproduce that version's calculation inputs.

```powershell
# First find the initial commit hash
chemvcs log --oneline
```

```powershell
# Reproduce using the initial commit hash (replace with your actual hash)
chemvcs reproduce <initial-commit-hash> -o reproduce_initial
```

```powershell
# View reproduced files
Get-ChildItem reproduce_initial

# Verify INCAR is the initial version (ENCUT=400)
type reproduce_initial\INCAR
```

**Key Talking Points**:
- Completely exports all input files from that version to a new directory
- Ensures 100% calculation reproducibility
- Critical for paper review and result verification

---

## Step 12: Summary

| Feature | Command | Problem Solved |
|---------|---------|----------------|
| Initialize | `chemvcs init` | Establish version control for computational projects |
| Track files | `chemvcs add` | Automatically identify VASP file types |
| Check status | `chemvcs status` | Understand workspace changes |
| Commit snapshot | `chemvcs commit` | Save complete snapshot of parameters + results |
| Semantic Diff | `chemvcs diff` | **Understand** parameter changes, not just text diffs |
| Version history | `chemvcs log` | Track calculation evolution process |
| Reproduce | `chemvcs reproduce` | Ensure calculation reproducibility |

### Energy Convergence Process

```
ENCUT=400, K=4×4×4  →  E = -10.8265 eV
ENCUT=520, K=4×4×4  →  E = -10.8452 eV  (ΔE = 18.7 meV)
ENCUT=520, K=8×8×8  →  E = -10.8489 eV  (ΔE =  3.7 meV → Converged!)
```

---

## Cleanup

```powershell
cd ~
Remove-Item $env:TEMP\si_convergence -Recurse -Force
```

---

## File Structure Reference

```
demo/
├── GUIDE.md                      ← This file (manual walkthrough)
├── run_demo.py                   ← Automated demo script (alternative)
└── vasp_files/                   ← Pre-made VASP files for each stage
    ├── step1_initial/            ← Initial inputs
    │   ├── POSCAR                   Si diamond structure
    │   ├── INCAR                    ENCUT=400, ISMEAR=0
    │   └── KPOINTS                  4×4×4 Gamma
    ├── step2_scf_done/           ← After calculation completion
    │   └── OUTCAR                   E=-10.8265 eV
    ├── step3_encut_conv/         ← ENCUT convergence test
    │   ├── INCAR                    ENCUT=520 (changed!)
    │   └── OUTCAR                   E=-10.8452 eV
    └── step4_kpoint_conv/        ← K-point convergence test
        ├── INCAR                    ISMEAR=1, SIGMA=0.1 (multiple changes!)
        ├── KPOINTS                  8×8×8 (changed!)
        └── OUTCAR                   E=-10.8489 eV
```
