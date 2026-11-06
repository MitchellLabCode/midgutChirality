#!/usr/bin/env bash
# #####################################################################################################################
# This script calls a Fiji macro and passes the timepoints to deconvolve.
# Specify start and end timepoints.
#
# $ bash deconvolve_timepoints.sh
#
# NPMitchell 2019
#
# Edited and adapted to Midway3 at UChicago RCC, Chris Anto 2024, canto@uchicago.edu
# Edited and adapted for parallel processing, Wenjie Shao 2024 wjsh@uchicago.edu
# #####################################################################################################################
source ~/.bashrc
module load fiji



echo "executing deconvolution"
xvfb-run /project/npmitchell/canto/fijiattempt/fiji-linux64/Fiji.app/ImageJ-linux64 -Dimage.updater.disableAutocheck=true --ij2 --run batch_ch0.ijm

# BACK_PID=$!
# wait $BACK_PID
echo "done with deconvolution ch0"

sleep 20s

xvfb-run /project/npmitchell/canto/fijiattempt/fiji-linux64/Fiji.app/ImageJ-linux64 -Dimage.updater.disableAutocheck=true --ij2 --run  batch_ch1.ijm

echo "done with deconvolution ch1"

xvfb-run /project/npmitchell/canto/fijiattempt/fiji-linux64/Fiji.app/ImageJ-linux64 -Dimage.updater.disableAutocheck=true --ij2 --run batch_ch2.ijm

echo "done with deconvolution ch2"

sleep 20s
