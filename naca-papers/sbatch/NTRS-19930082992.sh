#!/bin/bash
#SBATCH --job-name=naca-82992
#SBATCH --account=def-rahimk
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=0:10:00
#SBATCH --output=/scratch/rahimk/naca-sims/logs/%x_%j.out

# NTRS-19930082992 — Pekeris (1951), "Effect of Quadratic Terms in
# Differential Equations of Atmospheric Oscillations": exact verification of
# the quadratic-term estimate + re-application to the US76 atmosphere.
# CPU-only, < 1 min. Pattern: one sbatch per paper, named by NTRS id.
# Solver: ../sim/NTRS-19930082992.py (this script cd's into sim/ and runs it).
#
# PORTABILITY (other clusters): the two "module load" lines are Vulcan-
# specific. Run `module spider scipy-stack` (or python/numpy equivalents)
# on the new cluster and adjust; the python script itself is pure
# numpy/scipy/matplotlib. See naca-papers/PROCESS.md step 5 + gotcha 1.

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

python -m py_compile NTRS-19930082992.py
python NTRS-19930082992.py --results-dir "$REPO/naca-papers/results/NTRS-19930082992"

echo "=== outputs ==="
ls -la "$REPO/naca-papers/results/NTRS-19930082992"
