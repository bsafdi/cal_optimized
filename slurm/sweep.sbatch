#!/bin/bash
#SBATCH --job-name=cal_sweep
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr7
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --output=slurm/sweep-%j.out
#SBATCH --error=slurm/sweep-%j.err

set -euo pipefail
cd /global/scratch/projects/pc_heptheory/brsafdi/github/cal_optimized

module load miniforge3/25.9.1
module load texlive/2024

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

echo "host: $(hostname)"
echo "nproc: $OMP_NUM_THREADS"
date
python3 code/run_sweep.py --out outputs/production_sweep.csv --nt 15 --wire ${WIRE_TYPE:-production}
date
