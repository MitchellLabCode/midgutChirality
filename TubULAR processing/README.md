# Midgut surface pipeline

MATLAB scripts used to process the light-sheet movies in *Mechanical canalization of chiral gut morphogenesis*. Starting from fused, deconvolved multiview light-sheet stacks of the embryonic *Drosophila* midgut, they produce:

- segmented surface meshes in an anatomical (APDV) frame;
- cylindrical parameterizations and pullback images;
- surface kinematics.

The scripts are built on [TubULAR](https://github.com/npmitchell/tubular) (Mitchell & Cislo, *Nat. Methods* 2023) and follow its example scripts closely. See the [TubULAR documentation](https://npmitchell.github.io/tubular/) for detailed explanations.

Each sample gets its own copy of both scripts (e.g. `midgut_tubular_pipeline_WT1.m`). Run them in this order:

1. `midgut_tubular_pipeline.m`: pre-processing, surface segmentation, APDV frame, cylindrical cut meshes, pullbacks.
2. `midgut_tubular_kinematics.m`: time-smoothed meshes and pullbacks, PIV, 3D velocities and metric kinematics, then pathlines and strain from a reference time set by the cross-sample time alignment.

## Requirements

- **MATLAB.** TubULAR has been tested with releases 2018–2025.
- **TubULAR**, including its `DECLab` and `TexturePatch` submodules.
  - A plain clone does not pull the submodules; use `git clone --recursive`.
  - CGAL is needed for full functionality; see TubULAR's installation guide.
- **gptoolbox**: the copy bundled in TubULAR's `external/` folder, or a separate install.
- **ilastik**, for pixel classification (Autocontext, 2-stage workflow).
- **PIVlab**, which TubULAR's `measurePIV2d` calls by default.
- the file of per-sample time shifts from the cross-sample alignment (Supplement Sec. IV.3);
- **Input data**: fused, deconvolved 32-bit TIFFs, one file per timepoint and channel. Ours were acquired as six views at 60° increments on a Luxendo MuVi SPIM, then registered, fused and deconvolved as described in Methods.

## Usage

1. Copy both scripts for your sample and edit their settings:
   - pipeline script, **USER SETTINGS**: code paths, data directory, reference timepoints;
   - pipeline script, **MASTER SETTINGS**: resolution, timepoints, channels, intensity scaling, filename pattern of the input TIFFs;
   - kinematics script: the same paths, plus `shiftsFile`.
2. Run each script section by section in the MATLAB editor.
   - Most steps skip outputs that already exist unless `overwrite = true`.
   - After restarting MATLAB, rerun the settings sections first.
3. Do the pipeline script's two ilastik steps where it indicates them:
   - segment the gut on the h5 files written by `prepareIlastik`;
   - mark anterior, posterior and dorsal landmarks for the APDV frame.

Outputs follow TubULAR's standard layout under `<16-bit data directory>/tubular_output/`:

- meshes at the top level;
- (s, φ) meshes, pullbacks, PIV and kinematics in `gridCoords_nU0100_nV0100/` (for nU = nV = 100).

Our channel order:

- channel 0: HandGFP (circular visceral muscle nuclei);
- channel 1: membrane or RFP marker (mCherry-CAAX, Myo1C-RFP or HLH54F-RFP, depending on genotype);
- channel 2: H2Av-iRFP (endoderm nuclei).

Filenames use 0-based channel indices, while TubULAR uses 1-based ones; the scripts convert between them as needed.

## Differences from the standard TubULAR workflow

- **Pre-processing.** The pipeline script starts from 32-bit stacks and, before TubULAR sees the data:
  - converts them to 16-bit;
  - collates the channels into one TIFF per timepoint.
- **ilastik.**
  - The h5 files for ilastik are written in `cxyz` order so that ilastik's three orthogonal views display properly.
  - `getMeshes` then permutes the exported probabilities back to `xyzc` by matching axis lengths. This requires x, y and z to differ in length after subsampling.
  - The segmentation is trained with anterior and posterior endcap labels that close the tube for the level-set step (Supplement Sec. IV.1).
- **Surface and pullback settings.**
  - Meshes are smoothed in MATLAB (`smooth_with_matlab > 0`). It is also smoothed in Meshlab as detailed in this paper's repository in midgutChirality/MeshSmoothing.
- **Reference times.**
  - The phi0 reference of the (s, φ) maps (`t0_for_phi0`) is set in the pipeline script.
  - The Lagrangian reference time t0 for pathlines and strain comes from the cross-sample time alignment. The kinematics script writes it to `tubular_output/t0_referenceTime.txt` and uses a conformal (Ricci-flow) mesh at t0 as the pathline reference.
