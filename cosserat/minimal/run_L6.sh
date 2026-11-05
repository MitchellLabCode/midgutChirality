#!/bin/bash
#SBATCH --job-name=S14_L6        # create a short name for your job
#SBATCH -o ./logs/run_L6.sh.log-%A_%a
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --partition=caslake
#SBATCH -c 7                     # one core per run: 7 α
#SBATCH --time=6:00:00           # ~2.3 h per task
#SBATCH --account=pi-npmitchell
#SBATCH --array=0-7

# S14 family grown from L/Le = 1.05 to 6 in the hard-walled tube (R = 0.3 Le) and without
# confinement, η τE = 0, 10, 100, 1000, 7 values of α. The rod remeshes from 40 to ~160 edges;
# saved shapes use 260 edges to keep that resolution.
#   sbatch run_L6.sh              -> 20261020_tube_S14_{hard,free}_L6.0/
#   sbatch run_L6.sh nonoise      -> 20261020_tube_S14_{hard,free}_L6.0_nonoise/

module load julia
module load hdf5

TASKS=("hard 0" "hard 10" "hard 100" "hard 1000" "free 0" "free 10" "free 100" "free 1000")
read CONF ETA <<< "${TASKS[$SLURM_ARRAY_TASK_ID]}"
SIGMA=1e-3; [ "$1" = "nonoise" ] && SIGMA=0

TUBE_SIGMA=$SIGMA TUBE_NSAVE=260 TUBE_LEND=6 TUBE_ONLY_ETA=$ETA julia -p 7 tube_S14.jl $CONF
