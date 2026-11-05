# Midgut surface pipeline

MATLAB pipeline used to process the light-sheet movies in *Mechanical canalization of chiral gut morphogenesis*. Starting from deconvolved multiview light-sheet stacks of the embryonic *Drosophila* midgut, it produces:

- segmented surface meshes in an anatomical (APDV) frame,
- cylindrical parameterizations and pullback images,
- surface kinematics.

The pipeline is built on [TubULAR](https://github.com/npmitchell/tubular) (Mitchell & Cislo, *Nat. Methods* 2023) and follows its example scripts closely. See the [TubULAR documentation](https://npmitchell.github.io/tubular/) for anything not covered here.

`midgut_tubular_pipeline.m` is a template: each sample gets its own copy (e.g. `midgut_tubular_pipeline_WT1.m`).

## Requirements

- **MATLAB.** TubULAR has been tested with releases 2018–2021.
- **TubULAR**, including its `DECLab` and `TexturePatch` submodules.
  - A plain clone does not pull the submodules: use `git clone --recursive https://github.com/npmitchell/tubular`.
  - CGAL is needed for full functionality; see TubULAR's installation guide.
  - TubULAR also provides the pre-processing helpers used here (`convert32to16bit`, `makeMips`, `stabilizeImages`, `collateColors`) and `exportOBJ`.
  - We used commit `<commit>`.
- **gptoolbox**, either the copy bundled in TubULAR's `external/` folder or a separate install.
- **ilastik**, for pixel classification (Autocontext, 2-stage workflow).
- **PIVlab**, which TubULAR's `measurePIV2d` calls by default. A note in TubULAR's source reports that v3.02 works and v3.09 does not.
- **Input data**: fused, deconvolved 32-bit TIFFs, one file per timepoint and channel. Ours were acquired as six views at 60° increments on a Luxendo MuVi SPIM, then registered, fused and deconvolved as described in Methods.

## Usage

1. Copy the script for your sample and edit two blocks:
   - **USER SETTINGS**: code paths, data directory, reference timepoints.
   - **MASTER SETTINGS**: resolution, timepoints, channels, intensity scaling, filename pattern of the input TIFFs.
2. Run the script section by section in the MATLAB editor.
   - Most steps skip outputs that already exist unless `overwrite = true`.
   - After restarting MATLAB, rerun the two settings sections first.
3. Do the two ilastik steps where the script indicates them:
   - segment the gut on the h5 files written by `prepareIlastik`;
   - mark anterior, posterior and dorsal landmarks for the APDV frame.

Outputs follow TubULAR's standard layout under `<16-bit data directory>/tubular_output/`:

- meshes and cut meshes at the top level;
- (s, φ) meshes, pullbacks, PIV and kinematics in `gridCoords_nU0100_nV0100/` (for nU = nV = 100).

The channel-dependent defaults (`pbOptions.channels`, `falseColors`) assume our channel order:

- channel 0: HandGFP (circular visceral muscle nuclei);
- channel 1: membrane or RFP marker (mCherry-CAAX, Myo1C-RFP or HLH54F-RFP, depending on genotype);
- channel 2: H2Av-iRFP (endoderm nuclei).

Filenames use 0-based channel indices, while TubULAR uses 1-based ones; the script converts between them as needed.

## Differences from the standard TubULAR workflow

- **Pre-processing.** The script starts from 32-bit stacks, and then:
  - converts them to 16-bit;
  - collates the channels into one TIFF per timepoint.
- **ilastik.**
  - The h5 files for ilastik are written in `cxyz` order so that ilastik's three orthogonal views display properly.
  - `getMeshes` then permutes the exported probabilities back to `xyzc` by matching axis lengths (`detectOptions.useAxesPixels4Permutation = true`; off by default in TubULAR). This requires x, y and z to differ in length after subsampling, which the script checks.
  - The segmentation is trained with anterior and posterior endcap labels that close the tube for the level-set step (Supplement Sec. IV.1).
- **Parameter choices.**
  - Meshes are smoothed in MATLAB (`smooth_with_matlab > 0`). After this step, see https://github.com/MitchellLabCode/midgutChirality/MeshSmoothing for how they are further refined using Meshlab.
  - Pullbacks are made only for selected channels, as maximum intensity projections over 20 + 20 onion layers.
- **Additions.**
  - Surface area and enclosed volume over time (`measureSurfaceAreaVolume`, Fig. S7).
  - OBJ export of the (s, φ) meshes for visualization.
