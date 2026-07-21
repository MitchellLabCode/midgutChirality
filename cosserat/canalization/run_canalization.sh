#!/bin/bash
#SBATCH --job-name=canalization  # create a short name for your job
#SBATCH -o ./logs/canalization.sh.log-%j
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --partition=caslake
#SBATCH -c 12                     # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --time=8:00:00
#SBATCH --account=pi-npmitchell

module load julia
module load hdf5

julia -p 12 canalization.jl