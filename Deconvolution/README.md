# Multiview fusion and deconvolution

Fiji macros and SLURM scripts used to fuse and deconvolve the multiview light-sheet datasets in *Mechanical canalization of chiral gut morphogenesis* (Methods, "Light-sheet imaging and processing").

- These ran on shared nodes at the University of Chicago Research Computing Center's Midway3 cluster, as SLURM job arrays with one task per timepoint.
- Each task deconvolves the three channels of its timepoint in turn and writes one 32-bit TIFF stack per channel.
- These outputs are the input to `midgut_tubular_pipeline.m`.

## Files

- `CPU_deconvolution.sh`: SLURM submission script. It requests one array task per timepoint (`--array=0-90` for example) and runs `Deconvolve_timepoints.sh`.
- `Deconvolve_timepoints.sh`: starts Fiji under `xvfb-run` once for each channel macro, in sequence.
- `batch_ch0.ijm`, `batch_ch1.ijm`, `batch_ch2.ijm`: each fuses and deconvolves one channel of the timepoint given by `SLURM_ARRAY_TASK_ID`, then quits Fiji.

## Requirements

- **A SLURM cluster.** We used RCC's Midway3, partition `caslake` or `amd`.
- **Fiji with the legacy SPIM_Registration plugin** ("Multiview Reconstruction"), which provides the `Fuse/Deconvolve Dataset` command.
  - BigStitcher replaces this command with separate Image Fusion and Deconvolution commands, so the macros need adapting for a Fiji that only has BigStitcher.
  - We used `<Fiji / plugin version>`.
- **`xvfb-run`**, which gives Fiji the virtual display that its dialog-based plugins need on a compute node.
- **A registered dataset (`dataset.xml`).** The views were registered beforehand, using the fluorescent microspheres as interest points (Methods).
- **Measured PSFs**, one per excitation wavelength: `Mean_of_PSFs_488.tif`, `Mean_of_PSFs_561.tif` and `Mean_of_PSFs_642.tif`.

## Channels

| Channel | Excitation | Marker | Iterations |
|---|---|---|---|
| 0 | 488 nm | HandGFP: circular visceral muscle nuclei (gold in the renderings) | 20 |
| 1 | 561 nm | mCherry-CAAX, Myo1C-RFP or HLH54F-RFP, depending on genotype | 9 |
| 2 | 642 nm | H2Av-iRFP: endodermal nuclei (blue in the renderings) | 16 |

Each macro sets three things at the top:

- its channel;
- `excitation`, which selects the PSF file;
- the iteration count (`itenum`).

## Deconvolution settings

All channels share these settings:

- **Method:** multi-view deconvolution of all angles and illuminations, using "Efficient Bayesian – Optimization I" (Preibisch et al. 2014), Tikhonov regularization (λ = 0.006) and OSEM acceleration 1 (balanced).
- **Computation:** virtual weights; CellImg containers for images and FFTs; save-memory mode; the entire image is computed at once on the CPU (Java).
- **PSF:** one per wavelength, used for all angles and illuminations and transformed into each view.
- **Bounding box:** set by hand for each sample, in registered coordinates (`minx` to `maxz` at the top of each macro).
- **Output:** 32-bit TIFF stacks in `<dataDir>/deconvolved32bit/`, named `TP<t>_Ch<c>_Ill0_Ang<angles>.tif`. This is the pattern `file32Base` expects in `midgut_tubular_pipeline.m`.

## Usage

1. Edit the paths and sample settings:
   - in each macro: `dataDir` (the folder holding `dataset.xml`), the bounding box and the PSF folder;
   - in `Deconvolve_timepoints.sh`: the path to the Fiji launcher (`ImageJ-linux64`);
   - in `CPU_deconvolution.sh`: `--account`, and `--array` set to your timepoints. Array task IDs are used directly as timepoint IDs, so they must match the IDs in `dataset.xml`.
2. Submit from the folder containing all five files, since the scripts call each other by relative path:

   ```
   sbatch CPU_deconvolution.sh
   ```

Each task requested the full memory of a node (`--mem=0`) and up to 36 h.

- Fiji is started with `-Dimage.updater.disableAutocheck=true`, so it skips the update check at start-up. 
- Each macro ends with `System.exit(0)`, so Fiji quits and the next channel can start.

## References

- S. Preibisch, S. Saalfeld, J. Schindelin and P. Tomancak, Software for bead-based registration of selective plane illumination microscopy data, *Nat. Methods* **7**, 418 (2010).
- S. Preibisch, F. Amat, E. Stamataki, M. Sarov, R. H. Singer, E. Myers and P. Tomancak, Efficient Bayesian-based multiview deconvolution, *Nat. Methods* **11**, 645 (2014).

Computing resources were provided by the University of Chicago Research Computing Center.
