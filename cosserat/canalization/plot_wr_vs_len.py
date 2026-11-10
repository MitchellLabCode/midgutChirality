"""Writhe vs length of every run of a canalization sweep, from <folder>/wr_vs_len.mat
(precompute_writhe.py: writhe of every 4th snapshot; length = summed edge length, initial length 1).

alpha = 2*tau - 1 (alpha = -1: tau = 0, the WT end, blue; +1: tau = 1, red), on vik over [-1.2, 1.2].

usage: python plot_wr_vs_len.py [folder] [--exclude-jump] [--zorder abs]
  --exclude-jump  leave out the run flagged in wr_vs_len.mat (writhe discontinuity), as in the movie
  --zorder abs    stack the curves by |alpha|: alpha = 0 at the bottom, the largest |alpha| on top
                  (default: drawn in order of tau, so the alpha = +1 end is on top)
writes <folder>_wr_vs_len[_nojump].pdf and the same with _xlimcrop (xlim -0.5..0.7, as in cell 40)
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402
from scipy.io import loadmat  # noqa: E402

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"

args = [a for a in sys.argv[1:] if not a.startswith("--")]
folder = os.path.normpath(args[0] if args else "20261110_canalization_box3f_s6s3_Lnorm")
exclude = "--exclude-jump" in sys.argv
by_abs = "--zorder" in sys.argv and sys.argv[sys.argv.index("--zorder") + 1] == "abs"
args = [x for x in args if x != "abs"]

W = loadmat(os.path.join(folder, "wr_vs_len.mat"))
tau = W["tau"].ravel()
alpha = 2 * tau - 1
wr, length = W["wr"], W["len_ratio"]
excluded = W["excluded"].ravel().astype(bool)
here = os.path.dirname(os.path.abspath(__file__))
vik = ListedColormap(loadmat(os.path.join(here, "vik.mat"))["vik"])
mapnorm = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(-1.2, 1.2), cmap=vik)

for crop in (False, True):
    plt.figure(figsize=(5, 5), layout="constrained")
    for k in (np.argsort(np.abs(alpha), kind="stable") if by_abs else np.argsort(tau)):
        if exclude and excluded[k]:
            continue
        plt.plot(wr[k], length[k], label=rf"$\alpha={alpha[k]:.3g}$", color=mapnorm.to_rgba(alpha[k]),
                 zorder=2 + abs(alpha[k]) if by_abs else None)
    plt.ylabel("Length")
    plt.xlabel(r"Writhe $Wr$")
    ax = plt.gca()
    ax.set_box_aspect(1)
    plt.gcf().colorbar(mapnorm, ax=ax, orientation="vertical", label=r"$\alpha$")
    if crop:
        plt.xticks([-.5, -0.25, 0., 0.25, 0.5])
        plt.xlim(-0.5, 0.7)
    out = f"{os.path.basename(folder)}_wr_vs_len{'_nojump' if exclude else ''}{'_zabs' if by_abs else ''}{'_xlimcrop' if crop else ''}.pdf"
    plt.savefig(out)
    plt.close()
    print("saved", out)
if exclude:
    print("left out: alpha =", ", ".join(f"{a:+.3f}" for a in alpha[excluded]))
