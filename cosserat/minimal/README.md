# minimal: S14 rods grown in a tube

Pared-down code for the S14 tube runs: rods grown in a hard-walled tube (R = 0.3 Le) and without
confinement, to L/Le = 1.5 and 6, with and without rotational noise. The folder also has the
plots of those runs.

Each run grows an overdamped Cosserat rod, clamped at both ends, from L0 = 1.05 Le with:
- Maxwell remodeling of the rest curvature at fluidity η τE = 0, 10, 100, 1000;
- self-contact;
- remeshing.

The initial shape is the lowest clamped pair (S14) at initial asymmetry α ∈ [−1, 1], where
α = 0 is planar and α = ±1 are mirror images.

`tube_S14.jl` is the S14 / hard / free part of `noah/rods/reduced/tube_full.jl`. It uses the same
arithmetic, seeds and folder names, and reproduces that driver's output exactly: `xsaves`,
`Qsaves`, Tw and Wr were checked identical on short runs with and without noise.

## Files

| file | what |
|---|---|
| `tube_S14.jl` | simulation driver (Julia). It uses `../src/rod_contact2.jl`, `../src/fastkernel.jl` and `../prod/fastprod.jl`, and the environment in `../src/rods`. |
| `run_L1.5.sh`, `run_L6.sh` | SLURM launchers (array 0–7 = hard/free × 4 η). Add the argument `nonoise` for runs without noise. |
| `tubeplots.py` | shared helpers: loading, colormaps, 3D tube rendering, κ²τ and the Q measure |
| `plot_wr_vs_L.py` | Writhe during growth, rows η, columns hard / free, one line per α |
| `plot_3d.py` | 3D renders: perspective view, α < 0 / α > 0 overlay, and renders coloured by κ²τ or the Q measure |
| `plot_perversions.py` | number of perversions vs L/L0 at one η (default η τE = 10), with the wall contact marked |
| `vik.mat`, `managua_br.m` | colormaps (vik for α; managua_br for the chirality measures) |

## Running the simulations

From this folder:

```bash
sbatch run_L1.5.sh            # -> 20261020_tube_S14_{hard,free}_L1.5_na21/          (21 α, ~15 min/task)
sbatch run_L1.5.sh nonoise    # -> 20261020_tube_S14_{hard,free}_L1.5_na21_nonoise/
sbatch run_L6.sh              # -> 20261020_tube_S14_{hard,free}_L6.0/               (7 α, ~2.3 h/task)
sbatch run_L6.sh nonoise      # -> 20261020_tube_S14_{hard,free}_L6.0_nonoise/
```

Each output file is `eta<η>_a<i>.mat` (MAT v7.3). It contains:
- 200 snapshots, equally spaced in time: `xsaves`, `Qsaves`;
- `Lrest` (rest length L/Le), `Tw` (Wr = −Tw), `Wr_fuller`, `rho_max`;
- the parameters.

Runs are resumable: existing files are skipped. Options can be set by environment variable;
see the header of `tube_S14.jl`.

## Parameters

| quantity | value |
|---|---|
| initial / final length | 1.05 Le → 1.5 Le or 6 Le |
| growth rate | g τE = 1 (uniform, L = L0 e^{gt}) |
| fluidity | η τE = 0, 10, 100, 1000 |
| bending, twisting stiffness | B1 = B3 = 1 |
| stretch/shear stiffness | S = 6×10³ |
| drag | γL = γR = 1 |
| rotational noise | σ = 10⁻³, or 0 (`nonoise`) |
| tube radius, wall stiffness | R = 0.3 Le, W = 10⁴ per unit length |
| self-contact radius, stiffness | 0.14 Le, 10⁴ per unit length |
| edges | 40 initially, split at twice the initial length |
| saved shapes | 200, resampled on 120 (L = 1.5) or 260 (L = 6) edges |

## Plots

The plots read the run folders directly. Here `H` and `F` stand for the hard-tube and the free
folder of one set, e.g. `H=20261020_tube_S14_hard_L6.0_nonoise` and
`F=20261020_tube_S14_free_L6.0_nonoise`. `--at 3` takes the snapshot at L/L0 = 3 as the final
shape. Output goes to the current folder unless `--out` is given.

```bash
python plot_wr_vs_L.py $H $F                                      # Wr vs L/L0 (--x Le for L/Le)
python plot_3d.py $H $F --view persp                              # perspective view, coloured by α
python plot_3d.py $H $F --view above --overlay --at 3             # α<0 and α>0 overlaid, az −15°, el 82°
python plot_3d.py $H $F --view above --colour k2tau --at 3        # coloured by κ²τ
python plot_3d.py $H $F --view above --colour qmeasure --pct --at 3 --etas 10   # coloured by the Q measure
python plot_perversions.py $H --xlim 1,4 --ylim=-0.3,6.3          # perversion count, η τE = 10
python plot_perversions.py $F --thr 2.08e-4 --xlim 1,4 --ylim=-0.3,6.3   # on the threshold of the noisy hard runs
```

### Chirality measures

- **κ²τ** is the curvature-weighted torsion (r′×r″)·r‴/|r′|⁶. It is computed after resampling
  and Gaussian smoothing over 4 nodes.
- **Q measure** is the signed solid angle (sr) enclosed by the unit tangents over sliding windows
  of 0.1 Le. Each window is fanned from its own first tangent, after resampling at 0.01 Le and
  Gaussian smoothing over 0.08 Le. For a short window it is ≈ κ²τ ℓ³/12, but it needs only first
  derivatives, so it is much less noisy. Both measures are > 0 for right-handed stretches.

### Perversion count

The count is made on the Q-measure profile of every saved shape:
1. Threshold θ = 5% of the 95th percentile of |Q| over all shapes of all runs of the folder. For the
   noisy hard L = 6 runs this gives θ = 2.08×10⁻⁴ sr (1.77×10⁻⁴ sr without noise). `--thr` sets θ
   directly, e.g. `--thr 2.08e-4` to compare folders on one threshold.
2. Points with |Q| > θ are labelled by their sign; runs of equal sign form stretches.
3. Stretches shorter than 0.1 Le, and sub-threshold stretches, are dropped.
4. The count is the number of sign changes between the stretches that remain.

The black line is the median over α. The α = 0 runs are left out of the lines and the median
unless `--with-alpha0`, because without noise they stay planar; they still count in the
threshold. For hard-tube folders, the dashed line marks the first wall contact
(max ρ ≥ R, median over α; L/L0 ≈ 1.16).

Python requirements: numpy, scipy, h5py, matplotlib.
