#!/bin/bash
#SBATCH --job-name=contspl  # create a short name for your job
#SBATCH -o logs/contactSpline.sh.log-%j
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --partition=caslake
#SBATCH -c 12                     # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --time=4:00:00
#SBATCH --account=pi-npmitchell

module load julia
module load hdf5

julia -p 12 Qtwp_contact_spline.jl