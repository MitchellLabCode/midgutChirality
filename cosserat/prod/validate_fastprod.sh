#!/bin/bash
#SBATCH --job-name=validate_fastprod  # create a short name for your job
#SBATCH -o ./logs/validate_fastprod.sh.log-%j
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --partition=caslake
#SBATCH -c 1                      # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --time=2:00:00
#SBATCH --account=pi-npmitchell

# fastprod.jl vs run_simulation_plastic_remeshed on V6 data (short T), with timing

module load julia
module load hdf5

julia validate_fastprod.jl
