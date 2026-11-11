#!/bin/bash

#SBATCH --job-name=job_name      
#SBATCH --partition=caslake
#SBATCH --time=36:00:00
#SBATCH --mem=0G
#SBATCH --cpus-per-task=40
#SBATCH --account=account_name
#SBATCH --array=0-90 # range of timepoints for job array

bash Deconvolve_timepoints.sh
