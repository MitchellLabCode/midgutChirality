# processing

MATLAB code that turns the TubULAR reconstructions of the 10 imaged embryos into
the framed-curve dataset, aligns the embryos in time, and makes the Fig. 1 and
Fig. 2 panels, movies and figshare data. The rod simulations are in
[`../cosserat`](../cosserat).

## Setup

- MATLAB R2023b or newer, with the Image Processing Toolbox (`medfilt2`,
  `imgaussfilt`) and the Statistics and Machine Learning Toolbox (`pca`).
- [TubULAR](https://github.com/npmitchell/tubular). Only needed for
  `data_generation/` and `timeAlign/`.
- Python 3 with numpy and Pillow, plus ffmpeg, for `Fig1/make_movie.py`.

All paths and per-sample settings are in [`config_paths.m`](config_paths.m),
which every script calls. For machine-specific overrides, add a
`config_paths_local.m` (gitignored) containing `function P = config_paths_local(P)`
instead of editing `config_paths.m`. Every script is self-contained: open it and
run it. Outputs go to `Fig1/output/` and `Fig2/output/` (gitignored).

**Samples.** Samples 1 to 10 are numbered as in `config_paths.m`. Odd samples are
byn>Myo1C and even samples are WT. WT #k and byn>Myo1C #k are the k-th embryo of
each genotype, so sample 4 is WT #2 and sample 1 is byn>Myo1C #1.

**Offscreen rendering.** In `matlab -batch` without a display, the TeX labels in
exported PDFs can come out garbled. Use a (virtual) X display, e.g.
`DISPLAY=:1`, and `set(groot,'defaultFigureVisible','off')`.

## Pipeline

| step | script | reads | writes |
|---|---|---|---|
| 1 | `data_generation/extractNucleiStripes_noTracks.m` | TubULAR pullbacks (one sample, `idx_sample`) | hand-traced myocyte stripes, `tubular_output/cellTracking/muscleStripes/muscleNucleiStripes.mat` |
| 2 | `data_generation/tubi_to_meshrs.m` | TubULAR `spcutMeshSmRS` meshes | `P.meshes/meshes_sample<N>.mat` |
| 3 | `data_generation/process_all_V6.m` | centerlines, stripes, folds, meshes | `P.framesV6` = `all_framesInterpV6_10samples.mat` |
| 4 | `timeAlign/timelag_fromRadii.m` | TubULAR radii and folds | `timeAlign/timeshifts_fromRadiusKymos_10samples.mat` (tracked) |
| 5 | `bounding_contact/generate_bounding_contact.m` | `P.framesV6` | walls and contact radius for the simulations, `../cosserat/data/refIV6_sample<N>_walls3_cball.mat`, plus heatmaps |
| 6 | `Fig1/*`, `Fig2/*` | `P.framesV6`, time shifts | figures, movies, figshare data |

Steps 1 to 4 need the raw data and TubULAR. Step 1 is interactive (drawing on
the pullbacks), but already-traced stripes are reused, so rerunning it only
re-embeds them. The dataset (`all_framesInterpV6_10samples.mat`) and the
figshare files are distributed separately; the dataset's structure is described
in [`Fig1/figshare_readme.txt`](Fig1/figshare_readme.txt).

### Fig. 1
- `fig1_chirality_10samples.m`: for all embryos, 3D renders coloured by the
  chirality q (dots, and a tube version), the chirality kymographs q(s/L, t), and
  a PCA of the kymographs.
- `fig1_wr.m`: writhe and twist vs time, Tw vs Wr and Tw + Wr (data only), and
  `wr_tw_10samples_V6.csv`.
- `illustration_qmeasure.m`: the q illustration (curve and tantrix with the
  signed area).
- Movie frames (PNG per frame), assembled with `make_movie.py`:
  `fig1_tantrix_movie.m` (tangent sweeping the unit sphere),
  `fig1_qmeasure_frame.m` (title card), `fig1_chiral_movie.m` (q on two embryos,
  ventral and side views), `fig1_coordinates_movie.m` (longitudinal and
  azimuthal surface coordinates), and `fig1_frame_movie.m` (material frame d1, d2).
  The combined Fig. 1 movie command is in the `make_movie.py` docstring.
- `export_figshare_fig1.m`: `chiral_kymographs.h5`, the CSV and `readme.txt`.

### Fig. 2
- `fig2_deformation_measures.m`: all time-aligned curvature, rotation-rate and
  curvature-change summaries. These are E (bending and twisting energy densities,
  integrals and total angles), ω and μ = |∂ₜκ|, each split into d1/d2/twist
  components and into anterior/posterior/whole gut. The computation is in
  `utils/compute_deformation_measures.m`; the plotting is in
  `utils/plot_timealigned.m`.
- `fig2_radius.m`: mean radius vs time, and the radius kymograph of sample 4.
- `fig2_growthrate.m`: kymographs of edge growth.
- `asym_constrictions.m`: the three constrictions on the surface at given
  aligned times, ventral and side views.
- `fig2_wr_twist_sim.m`: Tw/Wr of data and simulations (`P.sims`).
- `export_figshare_fig2.m`: `fig2_kymographs.h5`, three CSVs and `readme.txt`.
- The contact-radius heatmaps come from step 5.

### utils
Shared functions:
- `qmeasure.m`: the signed spherical area of a set of tangents (formerly `stare.m`).
- `coiling_profile.m`, `coiling_kymograph.m`: q along one curve, and over time.
- `makeFramedCurve.m`, `logrot.m`: the material frame, the Darboux vector κ and
  the spin ω.
- `compute_writhe.m`, `compute_twist.m`, `compute_wr_tw.m`.
- `get_time_lags.m`, `aligned_time.m`: the time alignment.
- `compute_deformation_measures.m`, `contact_radius.m`.
- Plotting helpers (`tubeplot`, `arrow3d`, `spherical_patch`, `colorbar_extend`,
  `robust_ylim`, `sample_color`, `sample_label`, `save_frame`) and colormaps
  (`managua_br`, `coolwarm`, `blues`).

Third-party files and their licenses are in `utils/third_party/`: `curvspace`,
`distance2curve` and `savgol`.

## Time alignment and time axes

Each embryo i has a lag in frames (`get_time_lags.m`). It is the sum of two
shifts: the shift relative to its genotype's reference embryo, and for
byn>Myo1C the cross-genotype shift. The sum is then offset so that t = 0 is the
onset of the 2nd fold of sample 4. The lags are −16, −5, −21, −9, +1, −26, +9,
−6, −18 and −32 frames for samples 1 to 10, with 2 min between frames.

The published figures place frame f at different times
(`aligned_time.m`, `P.time_convention = 'published'`):

| kind | time [min] | used by |
|---|---|---|
| `frame` | 2f + 2·lag | Wr/Tw, radius, growth-rate y-limits, contact radius, the movies' timestamps |
| `frame0` | 2(f−1) + 2·lag | chirality kymographs (and the h5), E measures |
| `interval0` | 2(f−1) + 2·lag | ω (interval f→f+1) |
| `interval0mid` | 2(f−½) + 2·lag | μ (interval f→f+1) |
| `interval` | 2(f+½) + 2·lag | edge growth rate |

Only `frame` puts the reference fold onset at exactly t = 0. The other kinds
are shifted by up to 2 min. Setting `P.time_convention = 'unified'` places every
frame at 2f + 2·lag and every interval at its midpoint, 2(f+½) + 2·lag. This
moves the `frame0`, `interval0` and `interval0mid` axes by +1 to +2 min, and
changes the frame choice in `asym_constrictions.m` (see below). Note that
`figshare_readme.txt` describes the published convention.

## Known issues (kept as published)

- `asym_constrictions.m`: the published frame lookup also added
  `fold_onset(2)` of sample 4 (9 frames) on top of the lags, which already
  contain it. The panel labelled T min is therefore at T + 18 min on the axis
  of the other figures. With `'unified'`, the frame at T is used.
- `fig2_wr_twist_sim.m`: the byn>Myo1C simulation (started from sample 1) is
  drawn on sample 3's frames and lag.
- Sample 6, last frame: Tw drops from 0.94 to 0.08. This is the discrete twist
  wrapping by one turn at a very tight bend, not a real change.
- Bending rates ω_b and μ_b use the norm √(x₁² + x₂²), not |x₁| + |x₂| as in
  the Supplement. The "total energies" in Fig. S8 are integrals [µm⁻¹], not
  means per unit length.
- `fig2_growthrate.m` now places each interval row at its midpoint time. The
  original spread the n−1 rows over the n frame times, which differs by less
  than 1 min.

## Verification

With the published inputs, the following outputs were regenerated and are
identical to the published files:
- `process_all_V6.m`: all fields of all 10 samples.
- `tubi_to_meshrs.m`: the meshes and radii.
- `timelag_fromRadii.m`: all shifts and correlations.
- The stripe embedding of sample 4.
- `generate_bounding_contact.m`: the reference files in `cosserat/data`.
- The figshare h5 and CSV datasets.

The figures and movie frames match the earlier outputs. The tantrix and frame
movies are pixel-identical; the other figures differ only in layout details.
