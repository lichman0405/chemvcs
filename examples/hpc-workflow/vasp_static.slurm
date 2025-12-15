#!/bin/bash
#SBATCH --job-name=water-static
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --time=12:00:00
#SBATCH --partition=batch
#SBATCH --output=vasp_static.out
#SBATCH --error=vasp_static.err

# Load required modules
module load vasp/6.3.0
module load intel/2021.4
module load mkl/2021.4

# Set threading
export OMP_NUM_THREADS=1

# Run VASP single-point calculation
mpirun -np 28 vasp_std

echo "Static calculation complete"
