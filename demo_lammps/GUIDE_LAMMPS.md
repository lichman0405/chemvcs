# LAMMPS Demo: LJ Fluid NVT Convergence Study

This demo shows how to use ChemVCS to track a series of LAMMPS simulations
for a Lennard-Jones (LJ) fluid: first finding a good thermostat damping
parameter, then running production.

## Files

```
demo_lammps/
├── GUIDE_LAMMPS.md          # this file
└── lammps_files/
    ├── step1_initial/        # initial equilibration (NVT, T=1.0)
    │   ├── in.lammps
    │   ├── data.lj_fluid
    │   └── log.lammps
    ├── step2_production/     # production run (longer, larger damping)
    │   ├── in.lammps
    │   └── log.lammps
    └── step3_npt/            # NPT production (pressure coupling)
        ├── in.lammps
        └── log.lammps
```

## Step-by-step walkthrough

```bash
# Navigate to your LAMMPS working directory
cd /scratch/username/lj_fluid

# Initialize ChemVCS
chemvcs init

# --- Step 1: Stage initial equilibration setup ---
cp /path/to/demo_lammps/lammps_files/step1_initial/* .
chemvcs add in.lammps data.lj_fluid
chemvcs commit -m "Initial NVT equilibration: T=1.0, damp=0.1, 50k steps"

# Run LAMMPS
lmp -in in.lammps

# Stage and commit results
chemvcs add log.lammps
chemvcs commit -m "Equilibration complete: TotEng=-1005.2, T=1.003"

# --- Step 2: Production run with larger damping ---
cp /path/to/demo_lammps/lammps_files/step2_production/in.lammps .
chemvcs add in.lammps
chemvcs commit -m "Production NVT: damp=0.5, 200k steps"

lmp -in in.lammps

chemvcs add log.lammps
chemvcs commit -m "Production complete: TotEng=-1007.8"

# --- View history ---
chemvcs log

# --- Compare equilibration vs production input scripts ---
chemvcs diff HEAD~3 HEAD~1

# --- Reproduce exact equilibration inputs ---
chemvcs reproduce HEAD~3
```

## Expected diff output

When comparing the two input scripts, ChemVCS shows:

```
Semantic Diff: in.lammps
  ‼️  Modified timestep: 0.005 → 0.002
  ⚠️  Modified fix_nvt.Tdamp: 0.1 → 0.5
  ℹ️  Modified run: 50000 → 200000
  ℹ️  Modified thermo: 100 → 500
```
