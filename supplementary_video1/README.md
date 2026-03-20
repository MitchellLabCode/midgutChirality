# Gut Chirality Supplementary Video Pipeline

This folder contains the segmentation and masking workflow used to generate embryo and midgut RGB volumes for movies.

## Overview

1. Prepare an ilastik-ready base H5 from TIFF (`mask_prepare_ilastik.py`).
2. Train ilastik model for midgut-only segmentation.
3. Merge ilastik probabilities into the base data H5.
4. Run whole-embryo masking (`mask_embryo.py`).
5. Run midgut-only masking (`mask_midgut.py`).

## Step 0: Prepare ilastik training input (TIFF -> H5)

Run:

```bash
.venv/bin/python mask_prepare_ilastik.py \
  --tiff-path /path/to/Time_000005.tif
```

Default behavior:

- Writes `/path/to/Time_000005.h5` (same folder, same stem).
- Writes dataset name `inputData`.
- Downsamples by 4 in `xy` only (keeps `z` unchanged).
- Writes chunked H5 without compression by default (ilastik-friendly).

Useful options:

- `--h5-path /custom/output/Time_000005.h5`
- `--dataset-name inputData`
- `--downsample-factor 4`
- `--downsample-axes xy` (default) or `--downsample-axes xyz`
- `--compression none` (default) or e.g. `--compression gzip --compression-opts 4`

This H5 is the base file consumed by the merge step below.

## Step 1: Train ilastik segmentation

Train and save the ilastik project:

- `20240531_masking_midgutOnly_step1.ilp`

Then run ilastik batch export to produce probability H5 files (for each timepoint), typically named like:

- `Time_000005-inputData_Probabilities.h5`

Use dataset name `exported_data` in the probability file (this is what the scripts expect).

## Step 2: Merge ilastik probabilities with base input data

Run:

```bash
.venv/bin/python merge_ilastik_probabilities.py \
  --datadir /Users/npmitchell/Dropbox/Soft_Matter_Biophysics/PAPER/gut_chirality/supplementary_videos/20240531_142410_22C_HandGFPbynGAL4klar_UASmChCAAXH2AviRFP \
  --timepoint Time_000005 \
  --write-mode copy
```

This writes:

- `Time_000005_merged.h5`

with a 7-channel `inputData` dataset:

- channels 0-2: original image channels
- channels 3-6: ilastik probability channels

## Step 3: Whole-embryo masking

Run:

```bash
.venv/bin/python mask_embryo.py
```

Default key inputs:

- probability file: `Time_000005-inputData_Probabilities_embryo_vs_bg.h5`
- data file: `Time_000005_2x.tif`

Default key outputs:

- masked data: `Time_000005_WT_masked_embryo.tif`
- RGB data: `Time_000005_WT_masked_embryo_rgb.tif`

Notes:

- Mask is built from thresholding + largest connected component + hole fill (+ optional dilation).
- RGB brightness uses per-channel scaling and per-channel background subtraction.

## Step 4: Midgut-only masking

Run:

```bash
.venv/bin/python mask_midgut.py
```

Default key inputs:

- probability file: `Time_000005_merged-inputData_Probabilities_midgutOnly.h5`
- data file: `Time_000005_2x.tif`

Default key outputs:

- masked data: `Time_000005_WT_masked_midgutOnly.tif`
- RGB data: `Time_000005_WT_masked_midgutOnly_rgb.tif`

Notes:

- Uses the same masking logic as embryo masking.
- Default RGB channels for composition are `[0, 2]`.

## Optional next step: rotate RGB outputs (Fiji-equivalent pipeline)

Run:

```bash
.venv/bin/python rotate_anatomical_rgb.py
```

This script applies the same ordered operations to both embryo and midgut RGB stacks:

1. Rotate in-plane by `A1` (clockwise-positive, Fiji convention)
2. Reslice from Top
3. Rotate in-plane by `A2`
4. Reslice from Left
5. Rotate in-plane by `A3`
6. Final Reslice from Left
7. Final Rotate 90 Degrees Right

The script then crops to the tight nonzero bounding box after **each** step.

Default key options:

- `--fiji-angles A1 A2 A3` (default from script: `[13, -8, 140]`)
- `--interp-order` (default from script: `3`, cubic spline)
- `--enlarge/--no-enlarge` (default from script)

Default key outputs:

- rotated embryo RGB: `Time_000005_WT_masked_embryo_rgb_anatomical.tif`
- rotated midgut RGB: `Time_000005_WT_masked_midgutOnly_rgb_anatomical.tif`
- rotation report: `Time_000005_anatomical_rotation_report.txt`
- embryo crop history JSON: `Time_000005_WT_masked_embryo_rgb_anatomical_crop_history.json`
- midgut crop history JSON: `Time_000005_WT_masked_midgutOnly_rgb_anatomical_crop_history.json`

The crop-history JSON files contain per-step crop bounds and trim values so the same transforms can be applied to meshes.

To inspect settings and all CLI flags:

```bash
.venv/bin/python rotate_anatomical_rgb.py --help
```

## Practical reminders

- If you change filenames or timepoints, pass overrides via CLI flags.
- Keep `prob-dataset=exported_data` and `data-dataset=inputData` unless your files differ.
- Use `--help` on each script to see all options.
