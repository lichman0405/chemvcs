# Snapshots Directory

This directory contains reference VASP files for each step of the MOF-C₂H₂ adsorption workflow. These files are used to populate the `vasp_files/` working directory at each stage of the demo.

## Usage

Copy files from the appropriate step directory to `vasp_files/` before committing:

```bash
# For step 2:
cp snapshots/step2_add_c2h2/* vasp_files/

# For step 3:
cp snapshots/step3_vdw_opt/* vasp_files/

# For step 4:
cp snapshots/step4_high_precision/* vasp_files/
```

## Step Contents

- **step1_mof_opt**: Pristine Cu-BTC MOF optimization (already in vasp_files/)
- **step2_add_c2h2**: MOF + C₂H₂ guest molecule (added atoms in POSCAR, modified INCAR)
- **step3_vdw_opt**: With DFT-D3 dispersion corrections (IVDW=12, denser k-mesh)
- **step4_high_precision**: HSE06 hybrid functional (LHFCALC, higher ENCUT)

Each directory contains: INCAR, KPOINTS, POSCAR (step2-4 only), OUTCAR
