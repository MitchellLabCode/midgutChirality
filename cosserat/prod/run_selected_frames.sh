#!/bin/bash
#SBATCH --job-name=selframes     # create a short name for your job
#SBATCH -o logs/selected_frames.sh.log-%j
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --partition=caslake
#SBATCH -c 4                     # MATLAB rendering is essentially single-threaded
#SBATCH --mem=16G
#SBATCH --time=04:00:00          # ~360 frames; generous margin
#SBATCH --account=pi-npmitchell

module load matlab/2023a

# plot_selected_frames.m uses relative paths (selected/, ../data/)
cd "$SLURM_SUBMIT_DIR" || exit 1

# -batch: no desktop, exits with a nonzero code if the script errors.
# -softwareopengl: the nodes have no display, and the FaceAlpha transparency and
#   gouraud lighting need OpenGL, so use MATLAB's bundled software renderer.
matlab -softwareopengl -batch "plot_selected_frames"
