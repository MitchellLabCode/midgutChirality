"""Writhe during growth for the tube_coil.jl runs: rows η τE, columns hard-walled tube | no
confinement; one line per α (vik colours, α > 0 initially right-handed), length increasing
downward. Wr = -Tw, relative to the tube axis.

usage: python plot_wr_vs_L.py <hard_folder> <free_folder> [--x L0|Le] [--at L/L0] [--out file.pdf]
  e.g. python plot_wr_vs_L.py 20261020_tube_coil_hard_L6.0 20261020_tube_coil_free_L6.0 --x L0
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import tubeplots as tp  # noqa: E402

p = argparse.ArgumentParser()
p.add_argument("hard"); p.add_argument("free")
p.add_argument("--x", default="L0", choices=["L0", "Le"], help="length axis L/L0 (default) or L/Le")
p.add_argument("--at", type=float, default=None, help="stop at L/L0 = AT")
p.add_argument("--out", default=None)
a = p.parse_args()

confs = [("hard-walled tube", tp.load_runs(a.hard, a.at)), ("no confinement", tp.load_runs(a.free, a.at))]
tp.orient(*[rs for _, rs in confs])
cmap = tp.vik()
Lx = (lambda r: r["Lrest"] / r["Lrest"][0]) if a.x == "L0" else (lambda r: r["Lrest"])
Lname = r"$L/L_0$" if a.x == "L0" else r"$L/L_e$"
etas = [e for e in tp.ETAS if any(r["eta"] == e for _, rs in confs for r in rs)]
wr_lim = 1.08 * max(np.abs(r["Wr"]).max() for _, rs in confs for r in rs)

F, ax = plt.subplots(len(etas), 2, figsize=(5.2, 2.3 * len(etas)), sharex=True, sharey=True, squeeze=False)
for i, eta in enumerate(etas):
    for j, (cname, rs) in enumerate(confs):
        g = ax[i, j]
        for r in sorted([r for r in rs if r["eta"] == eta], key=lambda r: r["a"]):
            g.plot(r["Wr"], Lx(r), color=tp.alpha_colour(cmap, r["a"]), lw=0.8)
        g.axvline(0, color="k", lw=0.4)
        g.set_title(f"{tp.eta_label(eta)}, {cname}", fontsize=8)
        g.tick_params(labelsize=7, direction="in", top=True, right=True)
        g.set_xlim(-wr_lim, wr_lim)
    ax[i, 0].set_ylabel(Lname, fontsize=8)
r0 = confs[0][1][0]
ax[0, 0].set_ylim(Lx(r0)[-1], Lx(r0)[0])     # length increasing downward
for j in range(2):
    ax[-1, j].set_xlabel("Writhe Wr", fontsize=8)
cax = F.add_axes([0.93, 0.80, 0.015, 0.12])
F.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(-1, 1)), cax=cax)
cax.set_title(r"$\alpha$", fontsize=8); cax.tick_params(labelsize=7)
out = a.out or os.path.basename(os.path.normpath(a.hard)) + "_free" + (f"_at{a.at:g}" if a.at else "") + f"_Wr_vs_L{a.x}.pdf"
tp.savefig(F, out)
