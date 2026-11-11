#!/bin/bash
#SBATCH --job-name=canalization # create a short name for your job
#SBATCH -o ./logs/canalization.sh.log-%j
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --partition=caslake
#SBATCH -c 12                     # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --time=4:00:00           # fast solver: ~15-20 min per sweep; ref: ~2.5-3 h
#SBATCH --account=pi-npmitchell

module load julia
module load hdf5

# CANAL_GROWTH = common (default) | <final L/L0> | none (control), e.g. sbatch --export=ALL,CANAL_GROWTH=none run_canalization.sh
export CANAL_GROWTH=${CANAL_GROWTH:-common}
export CANAL_INPUT=${CANAL_INPUT:-SRVbox3f_s6s3} # sample 6 <-> sample 3
export CANAL_START=${CANAL_START:-xi}            # xi: first frame of the series; x0: the x0 field
export CANAL_ALPHA=${CANAL_ALPHA:-300}           # self-contact stiffness
export CANAL_SOLVER=${CANAL_SOLVER:-fast}         # fast (fastprod.jl, ~10x) | ref (rod_contact2.jl)
julia -p 12 canalization.jl