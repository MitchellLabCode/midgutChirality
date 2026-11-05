"""Number of perversions along the rod during growth, for one fluidity (default η τE = 10): one line
per α (vik colours), the median over α in black, and for hard-tube runs a dashed line at the first
wall contact (first saved length with max ρ ≥ R, median over α).

Perversions are counted on the Q measure profile (tubeplots.qmeasure_profile: windows of
ell_w = 0.1 Le, smoothing σ_s = 0.08 Le):
  1. threshold θ = 5% of the 95th percentile of |Q| over all saved shapes of all runs of the folder
     (or --thr, e.g. the hard-tube value 2.08e-4 sr to compare folders on one threshold);
  2. each point is +1 (Q > θ), -1 (Q < -θ) or 0; runs of equal sign form stretches;
  3. stretches of one handedness shorter than ell_min = 0.1 Le, and the sub-threshold ones, are dropped;
  4. perversions = sign changes between consecutive kept stretches.
The α = 0 runs (planar without noise) are left out of the lines and the median unless --with-alpha0;
they still count in the threshold, which is over all runs of the folder.

usage: python plot_perversions.py <folder> [--eta 10] [--thr 2.08e-4] [--xlim 1,4] [--ylim=-0.3,6.3] [--out f.pdf]
  e.g. python plot_perversions.py 20261020_tube_coil_hard_L6.0 --xlim 1,4 --ylim=-0.3,6.3
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import tubeplots as tp  # noqa: E402

p = argparse.ArgumentParser()
p.add_argument("folder")
p.add_argument("--eta", type=float, default=10.0)
p.add_argument("--thr", type=float, default=None, help="fixed threshold (sr) instead of 5%% of the folder's 95th percentile")
p.add_argument("--ell-w", type=float, default=0.1, help="window length (Le)")
p.add_argument("--sigma-s", type=float, default=0.08, help="smoothing (Le)")
p.add_argument("--ell-min", type=float, default=0.1, help="minimum stretch length (Le)")
p.add_argument("--with-alpha0", action="store_true")
p.add_argument("--xlim", default=None, help="L/L0 range, e.g. 1,4 ('max' = final value)")
p.add_argument("--ylim", default=None)
p.add_argument("--out", default=None)
a = p.parse_args()


def perversions(s, v, thr, lmin):
    """sign changes between stretches with |v| > thr of length ≥ lmin"""
    sg = np.where(np.abs(v) > thr, np.sign(v), 0.0)
    ds = np.gradient(s)
    runs = []                                # [sign, length]
    for x, d in zip(sg, ds):
        if runs and runs[-1][0] == x:
            runs[-1][1] += d
        else:
            runs.append([x, d])
    kept = [r[0] for r in runs if r[0] != 0 and r[1] >= lmin]
    return int(sum(1 for p_, q in zip(kept[:-1], kept[1:]) if p_ != q))


runs = tp.load_runs(a.folder)
tp.orient(runs)
hard = "_hard" in os.path.basename(os.path.normpath(a.folder))
for r in runs:
    r["prof"] = [tp.qmeasure_profile(xs, a.ell_w, a.sigma_s) for xs in r["xsaves"]]
    rho = np.hypot(r["xsaves"][..., 0], r["xsaves"][..., 1]).max(axis=1)
    hit = np.nonzero(rho >= r["tube_radius"])[0] if hard else []
    r["i_wall"] = int(hit[0]) if len(hit) else None
thr = a.thr if a.thr is not None else 0.05 * np.percentile(np.abs(np.concatenate([v for r in runs for _, v in r["prof"]])), 95)
sub = sorted([r for r in runs if r["eta"] == a.eta and (a.with_alpha0 or abs(r["alpha"]) > 1e-9)], key=lambda r: r["a"])
if not sub:
    raise SystemExit(f"no runs at η τE = {a.eta:g}")
for r in sub:
    r["perv"] = np.array([perversions(s, v, thr, a.ell_min) for s, v in r["prof"]])
X = lambda r: r["Lrest"] / r["Lrest"][0]
print(f"threshold |Q| > {thr:.3g} sr, stretches ≥ {a.ell_min:g} Le")
print(f"final perversions (α = {sub[0]['a']:+.2f} .. {sub[-1]['a']:+.2f}): {' '.join(str(r['perv'][-1]) for r in sub)}")

cmap = tp.vik()
F = plt.figure(figsize=(3.4, 2.9))
gx = F.add_axes([0.16, 0.16, 0.62, 0.62 * 3.4 / 2.9])   # square axes
for r in sub:
    gx.plot(X(r), r["perv"], color=tp.alpha_colour(cmap, r["a"]), lw=0.9, alpha=0.85)
wl = [X(r)[r["i_wall"]] for r in sub if r["i_wall"] is not None]
if wl:
    gx.axvline(float(np.median(wl)), color="0.4", lw=1.0, ls="--", label="wall contact", zorder=0)
    print(f"wall contact at L/L0 = {np.median(wl):.3f} (median over α)")
Lg = X(sub[0])
gx.plot(Lg, np.median([np.interp(Lg, X(r), r["perv"]) for r in sub], axis=0), color="k", lw=2.2, label="median over α")
gx.set_xlabel(rf"$L/L_0$ ($L_0 = {sub[0]['Lrest'][0]:.3g}\,L_e$)", fontsize=9)
gx.set_ylabel("number of perversions", fontsize=9)
gx.set_title(tp.eta_label(a.eta), fontsize=9)
gx.tick_params(labelsize=8, direction="in", top=True, right=True)
if a.xlim:
    x0, x1 = (max(X(r)[-1] for r in sub) if v.strip() == "max" else float(v) for v in a.xlim.split(","))
    gx.set_xlim(x0, x1)
gx.set_ylim(*([float(v) for v in a.ylim.split(",")] if a.ylim else (-0.3, max(r["perv"].max() for r in sub) + 0.3)))
gx.legend(fontsize=7, frameon=False, loc="upper left")
cax = F.add_axes([0.81, 0.16, 0.03, 0.62 * 3.4 / 2.9])
F.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(-1, 1)), cax=cax)
cax.set_ylabel(r"initial asymmetry $\alpha$", fontsize=8); cax.tick_params(labelsize=7, direction="in")
out = a.out or os.path.basename(os.path.normpath(a.folder)) + f"_perversions_eta{a.eta:g}" + (f"_thr{a.thr:.3g}" if a.thr else "") + ".pdf"
tp.savefig(F, out)
