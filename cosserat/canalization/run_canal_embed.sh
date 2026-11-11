#!/bin/bash
#SBATCH --job-name=canalembed
#SBATCH -o ./logs/canal_embed.sh.log-%j
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=caslake
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --account=pi-npmitchell

# PCA, metric MDS and t-SNE of the curvature / Q-measure profiles of a canalization sweep (canal_pca.py)
# usage: sbatch run_canal_embed.sh [options of canal_pca.py]   (default: --features kq --mds --tsne)

module load python/anaconda-2022.05
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

cd "$SLURM_SUBMIT_DIR" || exit 1
python canal_pca.py 20261110_canalization_box3f_s6s3_Lnorm ${@:---features kq --mds --tsne}
