#!/bin/bash
#SBATCH --job-name=diagplots     # create a short name for your job
#SBATCH -o logs/diagnostic_plots.sh.log-%j
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --partition=caslake
#SBATCH -c 4                     # MATLAB rendering is essentially single-threaded
#SBATCH --mem=32G                # meshes + all sweep files of a sample
#SBATCH --time=04:00:00
#SBATCH --account=pi-npmitchell

# Tube-grid diagnostic figures (rows: T x tw, cols: p x view) for every output folder
# matching folder_glob in diagnostic_plots_tube_QTpTw_quiver_batch.m; one PNG per sample
# and view, saved next to the folders: <folder>_tubegridq_allT_<view>.png

module load matlab/2023a

# the script uses relative paths (output folders, ../data/, vik.mat)
cd "$SLURM_SUBMIT_DIR" || exit 1

# -batch: no desktop, exits with a nonzero code if the script errors.
# -softwareopengl: the nodes have no display, and the FaceAlpha transparency and
#   gouraud lighting need OpenGL, so use MATLAB's bundled software renderer.
matlab -softwareopengl -batch "diagnostic_plots_tube_QTpTw_quiver_batch"
