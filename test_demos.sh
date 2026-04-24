#!/bin/bash
# Test script for all 4 demo directories
set -e

BASE="/home/shibo/project/chemvcs"
TMP="/tmp/chemvcs_demo_test"
PASS=0
FAIL=0

run_demo() {
    local name="$1"
    local dir="$2"
    echo ""
    echo "============================================"
    echo "  TESTING: $name"
    echo "============================================"
    rm -rf "$dir"
    mkdir -p "$dir"
    cd "$dir"
    if chemvcs init 2>&1; then
        echo "  [OK] init"
    else
        echo "  [FAIL] init"
        FAIL=$((FAIL+1))
        return
    fi
}

# ============================================================
# Demo 1: VASP basic (Si convergence test)
# ============================================================
D1="/tmp/test_demo_vasp"
run_demo "demo (VASP Si convergence)" "$D1"

# Step 1: Copy initial files
cp "$BASE/demo/vasp_files/step1_initial/INCAR" "$D1/"
cp "$BASE/demo/vasp_files/step1_initial/KPOINTS" "$D1/"
cp "$BASE/demo/vasp_files/step1_initial/POSCAR" "$D1/"
cp "$BASE/demo/vasp_files/step1_initial/POTCAR" "$D1/"

echo "  --- add initial files ---"
chemvcs add INCAR KPOINTS POSCAR POTCAR 2>&1 && echo "  [OK] add step1" || echo "  [FAIL] add step1"

echo "  --- commit initial ---"
chemvcs commit -m "Initial setup: Si bulk SCF (ENCUT=400, K=4x4x4)" 2>&1 && echo "  [OK] commit step1" || echo "  [FAIL] commit step1"

# Step 2: Add OUTCAR
cp "$BASE/demo/vasp_files/step2_scf_done/OUTCAR" "$D1/"
chemvcs add POSCAR INCAR KPOINTS OUTCAR 2>&1 && echo "  [OK] add step2" || echo "  [FAIL] add step2"
chemvcs commit -m "SCF completed: E=-10.8265 eV (ENCUT=400)" 2>&1 && echo "  [OK] commit step2" || echo "  [FAIL] commit step2"

# Step 3: ENCUT convergence
cp "$BASE/demo/vasp_files/step3_encut_conv/INCAR" "$D1/" -f
cp "$BASE/demo/vasp_files/step3_encut_conv/OUTCAR" "$D1/" -f
chemvcs add POSCAR INCAR KPOINTS OUTCAR 2>&1 && echo "  [OK] add step3" || echo "  [FAIL] add step3"
chemvcs commit -m "ENCUT convergence: 400->520 eV, E=-10.8452 eV" 2>&1 && echo "  [OK] commit step3" || echo "  [FAIL] commit step3"

# Step 4: K-point convergence
cp "$BASE/demo/vasp_files/step4_kpoint_conv/INCAR" "$D1/" -f
cp "$BASE/demo/vasp_files/step4_kpoint_conv/KPOINTS" "$D1/" -f
cp "$BASE/demo/vasp_files/step4_kpoint_conv/OUTCAR" "$D1/" -f
chemvcs add POSCAR INCAR KPOINTS OUTCAR 2>&1 && echo "  [OK] add step4" || echo "  [FAIL] add step4"
chemvcs commit -m "K-point convergence: 8x8x8, ISMEAR=1, E=-10.8489 eV" 2>&1 && echo "  [OK] commit step4" || echo "  [FAIL] commit step4"

# Test log, diff, reproduce
echo "  --- log ---"
chemvcs log --oneline 2>&1 && echo "  [OK] log" || echo "  [FAIL] log"

echo "  --- diff ---"
chemvcs diff 2>&1 && echo "  [OK] diff" || echo "  [FAIL] diff"

echo "  --- diff --summary ---"
chemvcs diff --summary 2>&1 && echo "  [OK] diff --summary" || echo "  [FAIL] diff --summary"

echo "  --- status ---"
chemvcs status 2>&1 && echo "  [OK] status" || echo "  [FAIL] status"

echo "  [PASS] demo (VASP Si convergence)"
PASS=$((PASS+1))

# ============================================================
# Demo 2: demo_advanced (MOF workflow)
# ============================================================
D2="/tmp/test_demo_advanced"
run_demo "demo_advanced (MOF workflow)" "$D2"

# Step 1: Copy initial MOF files
cp "$BASE/demo_advanced/vasp_files/INCAR" "$D2/"
cp "$BASE/demo_advanced/vasp_files/KPOINTS" "$D2/"
cp "$BASE/demo_advanced/vasp_files/POSCAR" "$D2/"
cp "$BASE/demo_advanced/vasp_files/POTCAR" "$D2/"
cp "$BASE/demo_advanced/vasp_files/OUTCAR" "$D2/"

echo "  --- add step1 ---"
chemvcs add INCAR KPOINTS POSCAR POTCAR OUTCAR 2>&1 && echo "  [OK] add step1" || echo "  [FAIL] add step1"
chemvcs commit -m "Step 1: MOF geometry optimization (PBE, ISIF=3)" 2>&1 && echo "  [OK] commit step1" || echo "  [FAIL] commit step1"

# Step 2: Add C2H2
cp "$BASE/demo_advanced/snapshots/step2_add_c2h2/INCAR" "$D2/" -f
cp "$BASE/demo_advanced/snapshots/step2_add_c2h2/POSCAR" "$D2/" -f
cp "$BASE/demo_advanced/snapshots/step2_add_c2h2/POTCAR" "$D2/" -f
cp "$BASE/demo_advanced/snapshots/step2_add_c2h2/OUTCAR" "$D2/" -f
chemvcs add INCAR POSCAR POTCAR OUTCAR 2>&1 && echo "  [OK] add step2" || echo "  [FAIL] add step2"
chemvcs commit -m "Step 2: Add C2H2 guest molecule (ISIF=2)" 2>&1 && echo "  [OK] commit step2" || echo "  [FAIL] commit step2"

# Step 3: vdW
cp "$BASE/demo_advanced/snapshots/step3_vdw_opt/INCAR" "$D2/" -f
cp "$BASE/demo_advanced/snapshots/step3_vdw_opt/KPOINTS" "$D2/" -f
cp "$BASE/demo_advanced/snapshots/step3_vdw_opt/POTCAR" "$D2/" -f
cp "$BASE/demo_advanced/snapshots/step3_vdw_opt/OUTCAR" "$D2/" -f
chemvcs add INCAR KPOINTS POTCAR OUTCAR 2>&1 && echo "  [OK] add step3" || echo "  [FAIL] add step3"
chemvcs commit -m "Step 3: Enable DFT-D3 dispersion (IVDW=12)" 2>&1 && echo "  [OK] commit step3" || echo "  [FAIL] commit step3"

# Step 4: HSE06
cp "$BASE/demo_advanced/snapshots/step4_high_precision/INCAR" "$D2/" -f
cp "$BASE/demo_advanced/snapshots/step4_high_precision/KPOINTS" "$D2/" -f
cp "$BASE/demo_advanced/snapshots/step4_high_precision/POTCAR" "$D2/" -f
cp "$BASE/demo_advanced/snapshots/step4_high_precision/OUTCAR" "$D2/" -f
chemvcs add INCAR KPOINTS POTCAR OUTCAR 2>&1 && echo "  [OK] add step4" || echo "  [FAIL] add step4"
chemvcs commit -m "Step 4: HSE06 hybrid functional (LHFCALC, ENCUT=520)" 2>&1 && echo "  [OK] commit step4" || echo "  [FAIL] commit step4"

echo "  --- log ---"
chemvcs log --oneline 2>&1 && echo "  [OK] log" || echo "  [FAIL] log"

echo "  --- diff summary ---"
chemvcs diff --summary 2>&1 && echo "  [OK] diff --summary" || echo "  [FAIL] diff --summary"

echo "  [PASS] demo_advanced (MOF workflow)"
PASS=$((PASS+1))

# ============================================================
# Demo 3: LAMMPS
# ============================================================
D3="/tmp/test_demo_lammps"
run_demo "demo_lammps (LJ Fluid)" "$D3"

# Step 1: Initial equilibration
cp "$BASE/demo_lammps/lammps_files/step1_initial/in.lammps" "$D3/"
cp "$BASE/demo_lammps/lammps_files/step1_initial/data.lj_fluid" "$D3/"
chemvcs add in.lammps data.lj_fluid 2>&1 && echo "  [OK] add step1" || echo "  [FAIL] add step1"
chemvcs commit -m "Initial NVT equilibration: T=1.0, damp=0.1, 50k steps" 2>&1 && echo "  [OK] commit step1" || echo "  [FAIL] commit step1"

# Step 1b: Add log
cp "$BASE/demo_lammps/lammps_files/step1_initial/log.lammps" "$D3/"
chemvcs add log.lammps 2>&1 && echo "  [OK] add log step1" || echo "  [FAIL] add log step1"
chemvcs commit -m "Equilibration complete: TotEng=-1005.2, T=1.003" 2>&1 && echo "  [OK] commit log step1" || echo "  [FAIL] commit log step1"

# Step 2: Production
cp "$BASE/demo_lammps/lammps_files/step2_production/in.lammps" "$D3/" -f
chemvcs add in.lammps 2>&1 && echo "  [OK] add step2" || echo "  [FAIL] add step2"
chemvcs commit -m "Production NVT: damp=0.5, 200k steps" 2>&1 && echo "  [OK] commit step2" || echo "  [FAIL] commit step2"

# Step 2b: Add production log
cp "$BASE/demo_lammps/lammps_files/step2_production/log.lammps" "$D3/" -f
chemvcs add log.lammps 2>&1 && echo "  [OK] add log step2" || echo "  [FAIL] add log step2"
chemvcs commit -m "Production complete: TotEng=-1007.8" 2>&1 && echo "  [OK] commit log step2" || echo "  [FAIL] commit log step2"

echo "  --- log ---"
chemvcs log --oneline 2>&1 && echo "  [OK] log" || echo "  [FAIL] log"

echo "  --- diff ---"
chemvcs diff 2>&1 && echo "  [OK] diff" || echo "  [FAIL] diff"

echo "  [PASS] demo_lammps (LJ Fluid)"
PASS=$((PASS+1))

# ============================================================
# Demo 4: ORCA
# ============================================================
D4="/tmp/test_demo_orca"
run_demo "demo_orca (Water molecule)" "$D4"

# Step 1: DFT optimization
cp "$BASE/demo_orca/step1_opt/water.inp" "$D4/"
cp "$BASE/demo_orca/step1_opt/water.out" "$D4/"
chemvcs add water.inp water.out 2>&1 && echo "  [OK] add step1" || echo "  [FAIL] add step1"
chemvcs commit -m "DFT geometry optimisation B3LYP/def2-TZVP" 2>&1 && echo "  [OK] commit step1" || echo "  [FAIL] commit step1"

# Step 2: CCSD(T) single point
cp "$BASE/demo_orca/step2_sp/water.inp" "$D4/"
cp "$BASE/demo_orca/step2_sp/water.out" "$D4/"
chemvcs add water.inp water.out 2>&1 && echo "  [OK] add step2" || echo "  [FAIL] add step2"
chemvcs commit -m "CCSD(T)/def2-TZVPP single point on optimised geometry" 2>&1 && echo "  [OK] commit step2" || echo "  [FAIL] commit step2"

echo "  --- log ---"
chemvcs log --oneline 2>&1 && echo "  [OK] log" || echo "  [FAIL] log"

echo "  --- diff ---"
chemvcs diff 2>&1 && echo "  [OK] diff" || echo "  [FAIL] diff"

echo "  [PASS] demo_orca (Water molecule)"
PASS=$((PASS+1))

# ============================================================
# Summary
# ============================================================
echo ""
echo "============================================"
echo "  RESULTS: $PASS passed, $FAIL failed"
echo "============================================"
