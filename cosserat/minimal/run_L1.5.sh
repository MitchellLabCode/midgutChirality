#!/bin/bash
#SBATCH --job-name=S14_L1.5      # create a short name for your job
#SBATCH -o ./logs/run_L1.5.sh.log-%A_%a
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --partition=caslake
#SBATCH -c 21                    # one core per run: 21 α
#SBATCH --time=1:00:00           # ~15 min per task
#SBATCH --account=pi-npmitchell
#SBATCH --array=0-7

# S14 family grown from L/Le = 1.05 to 1.5 in the hard-walled tube (R = 0.3 Le) and without
# confinement, η τE = 0, 10, 100, 1000, 21 values of α (-1:0.1:1).
#   sbatch run_L1.5.sh            -> 20261020_tube_S14_{hard,free}_L1.5_na21/
#   sbatch run_L1.5.sh nonoise    -> 20261020_tube_S14_{hard,free}_L1.5_na21_nonoise/

module load julia
module load hdf5

TASKS=("hard 0" "hard 10" "hard 100" "hard 1000" "free 0" "free 10" "free 100" "free 1000")
read CONF ETA <<< "${TASKS[$SLURM_ARRAY_TASK_ID]}"
SIGMA=1e-3; [ "$1" = "nonoise" ] && SIGMA=0

TUBE_SIGMA=$SIGMA TUBE_NALPHA=21 TUBE_LEND=1.5 TUBE_ONLY_ETA=$ETA julia -p 21 tube_S14.jl $CONF
