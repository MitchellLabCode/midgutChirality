# cosserat

Cosserat-rod growth simulations of the midgut. The folder contains:
- the source and the reproducible Julia/conda environment spec;
- the V6 reference data;
- three sets of runs: the production runs on the measured samples (`prod/`), the
  canalization runs (`canalization/`), and a minimal version of the tube runs of the coil rods
  (`minimal/`).

## Layout

```
cosserat/
├── src/                       # simulation source
│   ├── rod_contact2.jl        #   main solver (rest curvature in the material frame, growth of
│   │                          #   remeshed edges weighted by w_seg, optional prescribed total growth);
│   │                          #   include()s walls.jl + contact.jl
│   ├── fastkernel.jl          #   non-allocating rod kernel (same arithmetic as rod_contact2.jl)
│   ├── walls.jl               #   confinement / boundary polyhedron
│   ├── contact.jl             #   self-contact forces
│   ├── measurements.py        #   curvature, torsion, writhe (python)
│   └── rods/                  # environment spec (see "Environment" below)
│       ├── Project.toml       #   Julia deps
│       ├── Manifest.toml      #   pinned Julia versions (generated on Julia 1.9)
│       └── CondaPkg.toml      #   python deps (matplotlib, numpy) via CondaPkg.jl
├── data/                      # reference inputs: refIV6_sample<1-10>_walls3_cball.mat
├── prod/                      # production runs on the V6 data
│   ├── QTpTw_batch_samples_V6.jl  # batch over samples × (T, twist ratio, plastic rate)
│   ├── prod_common.jl         #   sample registry (reference file, x0, t0) and run setup
│   ├── fastprod.jl            #   fast solver (~10x) with walls, self-contact, remeshing; also used by
│   │                          #   canalization/
│   ├── validate_fastprod.jl   #   check of fastprod.jl against rod_contact2.jl
│   ├── x0_s*V5_t*.mat         #   initial centerlines
│   ├── run_V6_batch.sh, validate_fastprod.sh, run_diagnostic_plots.sh, run_selected_frames.sh
│   ├── make_selected.sh       #   collects the selected runs (T = 3, Bt/B = 1, ητE = 10) into selected/
│   ├── *.m                    #   MATLAB plots: tube grids, selected frames, dynamics
│   └── logs/
├── canalization/              # rods grown from shapes interpolated between sample 6 (WT) and sample 3
│   ├── canalization.jl        #   driver: one run per input file (options as environment variables)
│   ├── run_canalization.sh    #   SLURM launcher (12 workers); validate_fast.sh: fast vs reference solver
│   ├── SRVbox3f_s6s3/         #   inputs: 50 interpolated time series (tau = -0.1 .. 1.1), start shapes,
│   │                          #   circumradii, and walls containing every interpolated shape
│   │                          #   (made by ../processing/canalization/prep_interp_data.m)
│   ├── precompute_writhe.py   #   Wr vs L/L0 of a sweep -> <folder>/wr_vs_len.mat (flags self-passages)
│   ├── plot_wr_vs_len.py      #   Wr vs length of every run
│   ├── canal_pca.py           #   PCA / MDS / t-SNE of curvature and Q-measure profiles; run_canal_embed.sh
│   ├── plot_canalization_frames.m, run_canalization_frames.sh   # movie frames: tubes + Wr vs L/L0
│   ├── plot_tube_canal.m, plot_initial_tube.m, diagnostic_canalization.m, plot_tube_canal_overlay.m
│   │                          #   final / initial shapes as tubes, d1 quivers, overlays of ±alpha
│   └── logs/
└── minimal/                   # coil rods grown in a tube, hard / free, L = 1.5 and 6, ± noise
    ├── tube_coil.jl            #   driver; run_L1.5.sh, run_L6.sh
    ├── tubeplots.py           #   shared plotting helpers (κ²τ, Q measure, 3D tubes)
    ├── plot_wr_vs_L.py, plot_3d.py, plot_perversions.py
    └── README.md              #   parameters, measures and the commands for every figure
```

## Running

Each simulation script activates the shared environment and pulls source from `../src/`, so
scripts must be run from their own folder (`prod/`, `canalization/` or `minimal/`):

```bash
cd prod
sbatch run_V6_batch.sh               # all samples; or: sbatch run_V6_batch.sh 4 5 6
cd ../canalization
sbatch run_canalization.sh           # or e.g. sbatch --export=ALL,CANAL_GROWTH=none run_canalization.sh
cd ../minimal
sbatch run_L6.sh                     # or: sbatch run_L6.sh nonoise; see minimal/README.md
```

### Canalization

`canalization.jl` runs the 50 inputs of `SRVbox3f_s6s3/` (tau = 0: sample 6, WT; tau = 1: sample 3,
byn>Myo1C; the figures use alpha = 2 tau - 1, so alpha = -1 is the WT end). Each run starts from an
interpolated shape, follows the interpolated series at its clamped ends, grows with the growth field of
that series and is confined by the walls of the input file. Options (environment variables):

| variable | values |
|---|---|
| `CANAL_GROWTH` | `common` (default): every run grows to the same final length, R = mean input growth of the two measured ends; a number: that R; `none`: the growth of the input series |
| `CANAL_START` | `xi` (default): first frame of the series; `x0`: the `x0` shape of the input file |
| `CANAL_ALPHA` | self-contact stiffness (default 300) |
| `CANAL_SOLVER` | `fast` (default, `../prod/fastprod.jl`); `ref` (`run_simulation_plastic_remeshed`); `both` (compares them) |
| `CANAL_INPUT`, `CANAL_DATE`, `CANAL_ONLY`, `CANAL_T`, `CANAL_NSAVES` | input folder, date stamp of the output folder, runs to do, simulation time, number of saves |

Output: `<date>_canalization_box3f_s6s3_<Lnorm|ctrl>[_x0][_a<alpha>]/cutoff01contact_n<k>_tau<tau>.mat`
(xsaves, Qsaves, parameters). A sweep takes about 20 min on 12 cores with the fast solver. Then:

```bash
python precompute_writhe.py <folder>                 # Wr vs L/L0 -> <folder>/wr_vs_len.mat
python plot_wr_vs_len.py <folder> [--exclude-jump]   # Wr vs length
python canal_pca.py <folder> [--mds] [--tsne]        # shape-space embeddings
sbatch run_canalization_frames.sh <folder>           # movie frames
matlab -softwareopengl -batch "folder_override = \"<folder>\"; plot_tube_canal"   # also plot_initial_tube, diagnostic_canalization
```

On first run the scripts call `Pkg.activate("../src/rods")` + `Pkg.instantiate()`, which
rebuilds the Julia packages from `Manifest.toml` and the python/conda packages from
`CondaPkg.toml`. Output is written to a new dated folder created by the script itself.

## Environment

Only the environment *spec* files are included (`Project.toml`, `Manifest.toml`,
`CondaPkg.toml`). The packages are installed automatically by `Pkg.instantiate()` +
CondaPkg.jl on first run (into `src/rods/.CondaPkg/` and the Julia depot). The Manifest was
generated on Julia 1.9.

## Not included

- Simulation **output** folders and figures: they are produced by running the scripts.
- The midgut surface meshes (`data/meshes_spRS/meshes_sample<N>.mat`, ~209 MB). Only the MATLAB
  plots of `prod/` and `canalization/diagnostic_canalization.m`, `plot_initial_tube.m` and
  `plot_tube_canal.m` draw them; the simulations and the other plots do not need them.
