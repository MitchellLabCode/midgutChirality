#!/bin/bash
#SBATCH --job-name=prodV6 # create a short name for your job
#SBATCH -o logs/prod_V6_batch.sh.log-%j
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --partition=caslake
#SBATCH -c 36                     # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --time=8:00:00
#SBATCH --account=pi-npmitchell

# usage: sbatch run_V6_batch.sh            (all samples)
#        sbatch run_V6_batch.sh 4 5 6      (selected samples)

module load julia
module load hdf5

julia -p 36 QTpTw_batch_samples_V6.jl "$@"
