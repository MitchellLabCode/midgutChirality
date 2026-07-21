# cosserat

Cosserat-rod growth simulations of the midgut, extracted from `noah/rods`. Contains
the source, the reproducible Julia/conda environment spec, the input data, and the
two simulation campaigns (`contact/` and `canalization/`).

## Layout

```
cosserat/
├── src/                     # simulation source (Julia)
│   ├── rod_contact.jl       #   main solver; include()s walls.jl + contact.jl
│   ├── walls.jl             #   confinement / boundary polyhedron
│   ├── contact.jl           #   self-contact forces
│   └── rods/                # environment spec (see "Environment" below)
│       ├── Project.toml     #   Julia deps
│       ├── Manifest.toml    #   pinned Julia versions (generated on Julia 1.9)
│       └── CondaPkg.toml    #   python deps (matplotlib, numpy) via CondaPkg.jl
├── data/                    # reference .mat inputs (refIV5_sample*_walls*.mat, ...)
├── Q_vi/                    # extra x0 initial conditions referenced as ../Q_vi/ by
│                            #   contact/Qtwp_contact_spline.jl (samples 7 & 8)
├── contact/                 # self-contact + twist-plasticity campaign
│   ├── *.jl                 #   Qtwp_contact.jl, Qtwp_contact_spline.jl, QTp_contact_spline.jl
│   ├── run_*.sh             #   SLURM launchers (julia -p 12 <script>.jl)
│   ├── x0_*.mat             #   initial conditions
│   ├── *.m                  #   MATLAB plotting / diagnostics
│   └── logs/                #   SLURM log output dir (referenced by run_*.sh)
└── canalization/            # canalization campaign
    ├── *.jl                 #   canalization.jl, canalization_spline.jl, canalization_notail.jl
    ├── run_*.sh             #   SLURM launchers
    ├── SRVbox*/             #   .mat input ensembles (looped over via readdir)
    ├── data_interp/         #   interpolation inputs
    ├── *.m                  #   MATLAB plotting / diagnostics
    └── logs/                #   SLURM log output dir
```

## Running

Each simulation script activates the shared environment and pulls source from `../src/`,
so scripts must be run from their own directory (`contact/` or `canalization/`):

```bash
cd contact          # or: cd canalization
sbatch run_cspline.sh          # SLURM
# or interactively:
julia -p 12 Qtwp_contact_spline.jl
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
- Dated simulation **output** folders (e.g. `20260304_Q_twp_contact_*`), `.png`/`.pdf`
  figures, and result `.zip` archives — these are produced by running the sims.
- `data/meshes_spRS/` (~182 MB) — not referenced by the contact or canalization scripts.
- The redundant `SRVbox*.zip` archives (the unzipped `SRVbox*/` dirs are included instead).

If you need `meshes_spRS` or any output folder, copy it from `noah/rods/`.
