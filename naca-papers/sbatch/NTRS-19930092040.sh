#!/bin/bash
#SBATCH --job-name=naca-92040
#SBATCH --account=def-rahimk
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=0:20:00
#SBATCH --output=/scratch/rahimk/naca-sims/logs/%x_%j.out

# NTRS-19930092040 — Lessen (1950), "On Stability of Free Laminar Boundary
# Layer Between Parallel Streams": direct Orr-Sommerfeld solution of the free
# shear layer. CPU-only. Pattern: one sbatch per paper, named by NTRS id.
# Solver: ../sim/NTRS-19930092040.py (this script cd's into sim/ and runs it).

set -euo pipefail
mkdir -p /scratch/rahimk/naca-sims/logs

module --force purge
module load StdEnv/2026 scipy-stack/2026b
which python && python -V

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MPLCONFIGDIR="${SLURM_TMPDIR:-/tmp}/mpl"

REPO=/scratch/rahimk/repos/example-using-agenticai
cd "$REPO/naca-papers/sim"

python -m py_compile NTRS-19930092040.py
python NTRS-19930092040.py --results-dir "$REPO/naca-papers/results/NTRS-19930092040"

echo "=== outputs ==="
ls -la "$REPO/naca-papers/results/NTRS-19930092040"
