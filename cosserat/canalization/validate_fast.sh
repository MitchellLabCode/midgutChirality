#!/bin/bash
#SBATCH --job-name=canalfastval
#SBATCH -o ./logs/validate_fast.sh.log-%j
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=caslake
#SBATCH -c 4
#SBATCH --time=01:00:00
#SBATCH --account=pi-npmitchell

# Fast solver (../prod/fastprod.jl) against the reference one (../src/rod_contact2.jl) on three runs of
# the sample 6 <-> 3 canalization sweep, short time T = 0.1, 10 saves:
#   1. growth normalized, start from x0 (tests total_growth and the x0 start)
#   2. control (no normalization)
# Prints "VALIDATE run n: max |dx|, max |dQ|, time ref / fast" per run. Output goes to validate_* folders.

module load julia
module load hdf5
cd "$SLURM_SUBMIT_DIR" || exit 1

export CANAL_SOLVER=both CANAL_T=0.1 CANAL_NSAVES=10 CANAL_ONLY=1,25,50 CANAL_DATE=validate CANAL_INPUT=SRVbox3f_s6s3
CANAL_GROWTH=common CANAL_START=x0 julia -p 3 canalization.jl
CANAL_GROWTH=none   CANAL_START=xi julia -p 3 canalization.jl
grep -h "VALIDATE" logs/validate_fast.sh.log-$SLURM_JOB_ID | sed 's/.*VALIDATE/VALIDATE/'
