#!/bin/bash
#SBATCH --job-name=canalframes   # create a short name for your job
#SBATCH -o ./logs/canalization_frames.sh.log-%j
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --partition=caslake
#SBATCH -c 4                     # MATLAB rendering is essentially single-threaded
#SBATCH --mem=16G
#SBATCH --time=02:00:00          # 50 frames
#SBATCH --account=pi-npmitchell

module load matlab/2023a

# plot_canalization_frames.m uses relative paths (the run folder)
cd "$SLURM_SUBMIT_DIR" || exit 1

# the writhe panel reads this; produce it with precompute_writhe.py first
# usage: sbatch run_canalization_frames.sh [sweep folder]   (default 20261110_canalization_box3f_s6s3_Lnorm)
FOLDER=${1:-20261110_canalization_box3f_s6s3_Lnorm}
WR=$FOLDER/wr_vs_len.mat
if [ ! -f "$WR" ]; then
    echo "missing $WR -- run: python precompute_writhe.py $FOLDER"
    exit 1
fi

# -batch: no desktop, exits with a nonzero code if the script errors.
# -softwareopengl: the nodes have no display; gouraud lighting needs OpenGL.
matlab -softwareopengl -batch "folder_override = \"$FOLDER\"; plot_canalization_frames"
