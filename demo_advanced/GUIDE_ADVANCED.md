# Advanced Demo: MOF-Acetylene Adsorption Workflow

This demo demonstrates tracking a realistic computational chemistry workflow using ChemVCS. The scenario simulates studying acetylene (C₂H₂) adsorption in HKUST-1 (Cu-BTC), a prototypical metal-organic framework.

## Directory Structure

```
demo_advanced/
├── vasp_files/          # Working directory (starts with step1 files)
│   ├── INCAR
│   ├── KPOINTS
│   ├── POSCAR
│   └── OUTCAR
├── snapshots/           # Reference files for each step
│   ├── step1_mof_opt/
│   ├── step2_add_c2h2/
│   ├── step3_vdw_opt/
│   └── step4_high_precision/
└── GUIDE_ADVANCED.md    # This file
```

**Workflow**: Copy files from `snapshots/step*/` to `vasp_files/`, then commit. ChemVCS tracks parameter evolution across the same files.

## Scientific Context

**Metal-Organic Frameworks (MOFs)** are porous materials with applications in gas storage, separation, and catalysis. HKUST-1 features copper paddle-wheel units connected by benzene-tricarboxylate linkers, creating a 3D porous structure ideal for small molecule adsorption.

**Acetylene adsorption** in MOFs is important for separating C₂H₂ from CO₂ and C₂H₄ in petrochemical processes.

## Workflow Steps

### Step 1: MOF Geometry Optimization (step1_mof_opt)
**Purpose**: Optimize the pristine MOF structure before adding guest molecules.

**Key Parameters**:
- `ENCUT = 450` eV: Moderate plane-wave cutoff for organic systems
- `ISIF = 3`: Relax both cell shape and ionic positions
- `EDIFFG = -0.02`: Forces converged to 0.02 eV/Å
- `ISPIN = 2`: Spin-polarized (Cu²⁺ has unpaired electrons)
- `MAGMOM = 4*5.0 ...`: Initial magnetic moments on Cu centers

**Output**: Final energy = -200.53 eV

---

### Step 2: Add Acetylene Molecule (step2_add_c2h2)
**Purpose**: Introduce C₂H₂ into the MOF pore and relax ionic positions.

**Key Changes**:
- POSCAR: Added 2 C atoms and 2 H atoms (C₂H₂) at center of unit cell
- `ISIF = 2`: Fix cell, relax ions only
- `MAGMOM`: Updated to include C₂H₂ atoms (non-magnetic)

**Why this matters**: 
- Pure PBE functional underestimates weak dispersion interactions
- This step captures electrostatic and orbital overlap effects only

**Output**: Final energy = -211.79 eV

---

### Step 3: Van der Waals Optimization (step3_vdw_opt)
**Purpose**: Include dispersion corrections critical for physisorption.

**Key Changes**:
- `IVDW = 12`: Grimme's D3 method with BJ damping
- `LVDW = .TRUE.`: Enable vdW corrections
- `EDIFFG = -0.01`: Tighter convergence for weak interactions
- `KPOINTS`: Increased to 3×3×3 for better accuracy

**Scientific Impact**:
- D3-BJ adds ~16 eV of attractive dispersion energy
- Captures π-π stacking between C₂H₂ and aromatic linkers
- Crucial for predicting realistic adsorption energies

**Output**: Final energy = -228.65 eV (note dispersion contribution)

---

### Step 4: High-Precision Single Point (step4_high_precision)
**Purpose**: Benchmark with hybrid functional for publication-quality results.

**Key Changes**:
- `LHFCALC = .TRUE.`: Enable hybrid functional (HSE06)
- `HFSCREEN = 0.2`: Screening parameter for HSE06
- `AEXX = 0.25`: 25% exact exchange mixing
- `ENCUT = 520` eV: Higher cutoff for precision
- `EDIFF = 1E-7`: Electronic convergence to 10⁻⁷ eV
- `KPOINTS`: 4×4×4 mesh
- `NSW = 0`: No relaxation, just single-point energy
- `LCHARG = .TRUE.`: Save charge density for analysis

**Why hybrid functionals?**:
- PBE underestimates bandgaps and Cu-C₂H₂ orbital mixing
- HSE06 provides better description of metal d-orbitals
- Computationally expensive (~30× slower than PBE)

**Output**: Final energy = -231.57 eV

---

## Tracking with ChemVCS

This demo uses a **single-directory workflow** where you modify files in `vasp_files/` and commit changes at each step. ChemVCS tracks the evolution and shows semantic diffs.

### Step 1: Initialize and Track Initial MOF Structure

```bash
cd demo_advanced
chemvcs init

# vasp_files/ already contains step1 files (pristine MOF)
# ChemVCS will automatically validate POSCAR-POTCAR element order
chemvcs add vasp_files/INCAR vasp_files/KPOINTS vasp_files/POSCAR vasp_files/POTCAR vasp_files/OUTCAR
chemvcs commit -m "Step 1: MOF geometry optimization (PBE, ISIF=3)"
```

**Note on Validation**: ChemVCS automatically validates that POSCAR and POTCAR element orders match. If there's a mismatch, staging will be aborted with a clear error message. Use `--no-validate` to skip this check if needed.

### Step 2: Add C₂H₂ Guest Molecule

```bash
# Replace files with step2 versions (snapshots provided in snapshots/step2_add_c2h2/)
cp snapshots/step2_add_c2h2/INCAR vasp_files/
cp snapshots/step2_add_c2h2/POSCAR vasp_files/
cp snapshots/step2_add_c2h2/POTCAR vasp_files/
cp snapshots/step2_add_c2h2/OUTCAR vasp_files/

# Track changes (validation runs automatically)
chemvcs add vasp_files/INCAR vasp_files/POSCAR vasp_files/POTCAR vasp_files/OUTCAR
chemvcs commit -m "Step 2: Add C2H2 guest molecule (ISIF=2, extended MAGMOM)"
```

**View semantic diff from Step 1:**
```bash
chemvcs log  # Get commit hashes
chemvcs diff <step1_hash> <step2_hash>
```

### Step 3: Enable vdW Dispersion Corrections

```bash
# Replace files with step3 versions
cp snapshots/step3_vdw_opt/INCAR vasp_files/
cp snapshots/step3_vdw_opt/POTCAR vasp_files/
cp snapshots/step3_vdw_opt/OUTCAR vasp_files/

# Track changes
chemvcs add vasp_files/INCAR vasp_files/KPOINTS vasp_files/POTCAR
chemvcs add vasp_files/INCAR vasp_files/KPOINTS vasp_files/OUTCAR
chemvcs commit -m "Step 3: Enable DFT-D3 dispersion (IVDW=12, tighter convergence)"
```

**View semantic diff from Step 2:**
```bash
chemvcs diff <step2_hash> <step3_hash>
```

Expected output shows:
- **INCAR changes**: `+IVDW=12`, `+LVDW=.TRUE.`, `EDIFF` 1E-5→1E-6, `EDIFFG` -0.02→-0.01
- **KPOINTS changes**: Grid 2×2×2 → 3×3×3
- **OUTCAR changes**: Energy -211.79 → -228.65 eV (**~17 eV dispersion contribution!**)

### Step 4: High-Precision HSE06 Single Point

```bash
# Replace files with step4 versions
cp snapshots/step4_high_precision/INCAR vasp_files/
cp snapshots/step4_high_precision/POTCAR vasp_files/
cp snapshots/step4_high_precision/OUTCAR vasp_files/

# Track changes
chemvcs add vasp_files/INCAR vasp_files/KPOINTS vasp_files/POTCAR
chemvcs add vasp_files/INCAR vasp_files/KPOINTS vasp_files/OUTCAR
chemvcs commit -m "Step 4: HSE06 hybrid functional (LHFCALC, ENCUT=520, 4x4x4 grid)"
```

**View semantic diff from Step 3:**
```bash
chemvcs diff <step3_hash> <step4_hash>
```

Expected output shows:
- **INCAR changes**: `+LHFCALC=.TRUE.`, `+HFSCREEN=0.2`, `+AEXX=0.25`, `ENCUT` 450→520, `EDIFF` 1E-6→1E-7, `NSW` 200→0
- **KPOINTS changes**: Grid 3×3×3 → 4×4×4
- **OUTCAR changes**: Energy -228.65 → -231.57 eV, Runtime 2234 → 8734 sec

### Compare Any Two Steps

```bash
# Compare step1 vs step4 (full evolution)
chemvcs diff <step1_hash> <step4_hash>
```

Shows cumulative changes across entire workflow!

---Plugin System: Automatic Validation

ChemVCS includes a plugin system for extensible validation. By default, the **POSCAR-POTCAR validator** ensures element order consistency.

### Viewing Available Plugins

```bash
chemvcs plugin list
```

Output:
```
Discovered 3 validator plugin(s):

  ○ disabled  file-format v0.1.0
           Validates basic VASP file format correctness
           Priority: 5

  ✓ enabled  incar-poscar v0.1.0
           Validates INCAR-POSCAR parameter consistency
           Priority: 20

  ✓ enabled  poscar-potcar v0.1.0
           Validates POSCAR-POTCAR element order consistency
           Priority: 10
```

### Plugin Information

```bash
chemvcs plugin info poscar-potcar
```

### Validation in Action

When you run `chemvcs add`, validators automatically check files:

```bash
chemvcs add vasp_files/POSCAR vasp_files/POTCAR
```

Output (success):
```
Running validators...
  ✓  poscar-potcar: passed

Staging files...
  + vasp_files/POSCAR  (1.5 KB, POSCAR)
  + vasp_files/POTCAR  (777 B, POTCAR)
```

Output (failure - element mismatch):
```
Running validators...
  ✗  poscar-potcar: Element order mismatch

Staging aborted due to validation errors
  Use --no-validate to skip validation
```

### Skipping Validation

If you need to bypass validators (e.g., working with non-standard files):

```bash
chemvcs add vasp_files/* --no-validate
```

---

## 

## Key Insights from ChemVCS Tracking

1. **Energy Progression**:
   - MOF alone: -200.53 eV
   - +C₂H₂ (PBE): -211.79 eV → Δ = -11.26 eV
   - +vdW (PBE-D3): -228.65 eV → Δ = -16.86 eV dispersion
   - Hybrid (HSE06): -231.57 eV → Δ = -2.92 eV electronic correction

2. **Parameter Evolution**:
   - Convergence criteria tightened as workflow progresses
   - K-point mesh density increases for final calculations
   - Dispersion corrections essential (contributes ~60% of total binding)

3. **Computational Cost**:
   - Step 1: 1.3k sec
   - Step 2: 1.5k sec  
   - Step 3: 2.2k sec (D3 overhead ~40%)
   - Step 4: 8.7k sec (HSE06 ~4× slower)

4. **Reproducibility**:
   - Each commit captures complete calculation state
   - Semantic diffs highlight **what changed** and **why it matters**
   - Easy to reproduce any step or branch alternative methods

---

## Next Steps

- Calculate adsorption energy: E(MOF+C₂H₂) - E(MOF) - E(C₂H₂)
- Compare with experimental values (~30-40 kJ/mol)
- Test different loadings (multiple C₂H₂ molecules)
- Screen other MOFs (UiO-66, MIL-101, etc.)

ChemVCS enables tracking all these variations systematically!
