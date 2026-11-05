# cosserat

Cosserat-rod growth simulations of the midgut, extracted from `noah/rods`. The folder contains:
- the source and the reproducible Julia/conda environment spec;
- the V6 reference data;
- three sets of runs: the production runs on the measured samples (`prod/`), the
  canalization runs (`canalization/`), and a minimal version of the tube runs of the S14 rods
  (`minimal/`).

## Layout

```
cosserat/
├── src/                       # simulation source
│   ├── rod_contact2.jl        #   main solver (rest curvature in the material frame, growth of
│   │                          #   remeshed edges weighted by w_seg); include()s walls.jl + contact.jl
│   ├── rod_contact.jl         #   previous solver, used by canalization/canalization_notail.jl
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
│   ├── fastprod.jl            #   fast kernel with walls, self-contact, remeshing
│   ├── validate_fastprod.jl   #   check of fastprod.jl against rod_contact2.jl
│   ├── x0_s*V5_t*.mat         #   initial centerlines
│   ├── run_V6_batch.sh, validate_fastprod.sh, run_diagnostic_plots.sh, run_selected_frames.sh
│   ├── make_selected.sh       #   collects the selected runs (T = 3, Bt/B = 1, ητE = 10) into selected/
│   ├── *.m                    #   MATLAB plots: tube grids, selected frames, dynamics
│   └── logs/
├── canalization/              # canalization runs
│   ├── canalization.jl, canalization_spline.jl, canalization_notail.jl
│   ├── run_*.sh               #   SLURM launchers
│   ├── SRVbox*/               #   .mat input ensembles (looped over via readdir)
│   ├── data_interp/           #   interpolation inputs
│   ├── precompute_writhe.py   #   Wr vs L/L0 of a sweep -> <folder>/wr_vs_len.mat
│   ├── plot_canalization_frames.m, plot_tube_canal.m, plot_tube_canal_overlay.m, ...   # MATLAB plots
│   └── logs/
└── minimal/                   # S14 rods grown in a tube, hard / free, L = 1.5 and 6, ± noise
    ├── tube_S14.jl            #   driver; run_L1.5.sh, run_L6.sh
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
sbatch run_canalization.sh
cd ../minimal
sbatch run_L6.sh                     # or: sbatch run_L6.sh nonoise; see minimal/README.md
```

On first run the scripts call `Pkg.activate("../src/rods")` + `Pkg.instantiate()`, which
rebuilds the Julia packages from `Manifest.toml` and the python/conda packages from
`CondaPkg.toml`. Output is written to a new dated folder created by the script itself.

## Environment

Only the environment *spec* files are included (`Project.toml`, `Manifest.toml`,
`CondaPkg.toml`). The actual installed binaries in `noah/rods` — `src/rods/.CondaPkg/`
(~1.2 GB) and `noah/rods/env/` (~1.9 GB) — are **not** copied: they are large and contain
machine-specific absolute paths. They are regenerated automatically by
`Pkg.instantiate()` + CondaPkg.jl on first run. The Manifest was generated on Julia 1.9.

## Not included (results / unused inputs)

To keep the bundle lean, the following from `noah/rods` were intentionally omitted:
- Dated simulation **output** folders, `.png`/`.pdf` figures, and result `.zip` archives —
  these are produced by running the sims.
- `data/meshes_spRS/` (~209 MB). The MATLAB plots of `prod/` and some of `canalization/` draw the
  measured midgut meshes from it; copy it to `data/meshes_spRS/` to run them.
- The redundant `SRVbox*.zip` archives (the unzipped `SRVbox*/` dirs are included instead).

If you need `meshes_spRS` or any output folder, copy it from `noah/rods/`.
