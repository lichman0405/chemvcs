#!/bin/bash
#SBATCH --job-name=water-relax
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=28
#SBATCH --time=24:00:00
#SBATCH --partition=batch
#SBATCH --output=vasp_relax.out
#SBATCH --error=vasp_relax.err

# Load required modules
module load vasp/6.3.0
module load intel/2021.4
module load mkl/2021.4

# Set threading
export OMP_NUM_THREADS=1

# Run VASP geometry optimization
mpirun -np 56 vasp_std

echo "Relaxation complete"
