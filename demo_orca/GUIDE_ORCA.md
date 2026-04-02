# ORCA Demo: Water Molecule — From Optimization to High-Level Single Point

This demo walks through a typical quantum chemistry workflow with ORCA:

1. **Step 1** – Geometry optimisation of water at B3LYP/def2-TZVP level (DFT)
2. **Step 2** – High-level single point energy at CCSD(T)/def2-TZVPP on the optimised geometry

All ORCA input/output files in this directory are synthetic examples generated
to illustrate ChemVCS features without needing a real ORCA installation.

## Workflow

```bash
# Initialise a new ChemVCS repository
chemvcs init

# --- Step 1: DFT geometry optimisation ---
chemvcs add step1_opt/water.inp
chemvcs commit -m "DFT geometry optimisation B3LYP/def2-TZVP"

# View workflow
chemvcs log

# --- Step 2: High-level single point ---
chemvcs add step2_sp/water.inp
chemvcs commit -m "CCSD(T)/def2-TZVPP single point on optimised geometry"

# Compare the two input files
chemvcs diff step2_sp/water.inp
```

The semantic diff will highlight:
- **‼️ CRITICAL**: Method changed from `B3LYP` → `CCSD(T)`, basis set from `DEF2-TZVP` → `DEF2-TZVPP`
- **⚠️ MAJOR**: RI approximation `RIJCOSX` removed (not available for CCSD(T))

## Files

| File | Description |
|------|-------------|
| `step1_opt/water.inp` | B3LYP/def2-TZVP geometry optimisation input |
| `step1_opt/water.out` | Synthetic DFT optimisation output |
| `step2_sp/water.inp`  | CCSD(T)/def2-TZVPP single point input (optimised geometry) |
| `step2_sp/water.out`  | Synthetic CCSD(T) single point output |
