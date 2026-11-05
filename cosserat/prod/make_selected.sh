#!/bin/bash
# Build selected/ for plot_selected_frames.m from the V6 batch outputs, with the same
# choice as the V5 selection: T = 3, Bt/B = 1, ητE = 10, and the same start times t0
# (s7 from t0 = 20, s10 from t0 = 35). Even samples are WT, odd samples myo.
#   usage (from prod/):  bash make_selected.sh [date_str]
date_str=${1:-20261019}
run=cutoff01contact_T3.0_tw1.0_p10.0.mat
declare -A t0=([1]=38 [2]=25 [3]=45 [4]=36 [5]=25 [6]=54 [7]=20 [8]=37 [9]=45 [10]=35)
mkdir -p selected
for s in 1 2 3 4 5 6 7 8 9 10; do
    if (( s % 2 == 0 )); then lab=WT; else lab=myo; fi
    src=${date_str}_QTpTw_contact2_t${t0[$s]}_w3_a3e2_s${s}V6/$run
    dst=selected/sample${s}_${lab}_T3.0_tw1.0_p10.0.mat
    if [ -f "$src" ]; then cp -p "$src" "$dst" && echo "$dst <- $src"; else echo "MISSING $src"; fi
done
