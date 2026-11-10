"""Canalization in shape space: PCA of curvature and Q-measure profiles of the smoothed rods.

Every shape (every STRIDE-th snapshot of every run of a canalization folder) is described by
  κ(s)  curvature                                  (1/L0)
  Q(s)  Q measure: signed solid angle swept by the unit tangents over windows of ELL_W (sr; > 0
        right-handed); chiral, unlike κ
on NU points of normalised arclength s/L, after resampling the centerline at H and Gaussian smoothing
over SIGMA_S (lengths in units of the initial length L0). The two blocks are scaled to the same total
variance and concatenated; PCA is fitted on all shapes of all runs minus the runs flagged for a writhe
discontinuity in wr_vs_len.mat, which are also left out of the plots and the MDS / t-SNE (×
in panel f).

Figure (<folder>_pca_<features>.pdf):
  a  explained variance
  b  PC1 vs length, one line per run, coloured by α  (the analogue of Wr vs length)
  c  PC2 vs length
  d  PC1-PC2 trajectories (o = start, ● = end)
  e  PC1 and PC2 loadings along s/L, κ block and Q block
  f  outcome vs initial α at a common length L* (--at; default the shortest final length, since the
     runs end at different lengths): PC1 and Wr at L*; the shapes at L* clustered (Ward, --nclust
     clusters, on the PCs that explain 90% of the variance), cluster as marker shape

α = 2τ − 1 (α = −1: τ = 0, the WT end, blue), vik over [−1.2, 1.2] as in plot_wr_vs_len.py.

usage: python canal_pca.py [folder] [--features kq|k|q] [--stride 4] [--at L] [--nclust 2] [--mds] [--tsne]
  --mds  also <folder>_mds_<features>.pdf: metric MDS (SMACOF, sklearn) of the same embedding
         (Euclidean distances between the feature vectors; classical MDS would equal the PCA),
         rotated/reflected onto PC1-PC2 (orthogonal Procrustes):
           a  all shapes, trajectories (o start, ◆ L*, ● end)   b  MDS 1 vs length
           c  MDS of the shapes at L* alone, cluster as marker shape
  --tsne same three panels with t-SNE (sklearn; --perplexity 30 for all shapes, --perplexity-lstar
         10 for the 50 shapes at L*), PCA-initialised and rotated onto PC1-PC2 likewise
"""
import argparse
import glob
import os

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402
from scipy.cluster.hierarchy import fcluster, linkage  # noqa: E402
from scipy.io import loadmat  # noqa: E402
from scipy.ndimage import gaussian_filter1d  # noqa: E402

p = argparse.ArgumentParser()
p.add_argument("folder", nargs="?", default="20261110_canalization_box3f_s6s3_Lnorm")
p.add_argument("--features", default="kq", choices=["kq", "k", "q"])
p.add_argument("--stride", type=int, default=4)
p.add_argument("--sigma-s", type=float, default=0.04, help="smoothing (L0)")
p.add_argument("--ell-w", type=float, default=0.15, help="Q-measure window (L0)")
p.add_argument("--h", type=float, default=0.005, help="resampling step (L0)")
p.add_argument("--nu", type=int, default=100, help="points along s/L")
p.add_argument("--at", type=float, default=None, help="common length L* of the outcome (L0)")
p.add_argument("--nclust", type=int, default=2)
p.add_argument("--mds", action="store_true")
p.add_argument("--tsne", action="store_true")
p.add_argument("--perplexity", type=float, default=30.0, help="t-SNE, all shapes")
p.add_argument("--perplexity-lstar", type=float, default=10.0, help="t-SNE, the shapes at L*")
a = p.parse_args()
folder = os.path.normpath(a.folder)

plt.rcParams.update({"pdf.fonttype": 42, "xtick.direction": "in", "ytick.direction": "in"})
here = os.path.dirname(os.path.abspath(__file__))
vik = ListedColormap(loadmat(os.path.join(here, "vik.mat"))["vik"])
norm = mpl.colors.Normalize(-1.2, 1.2)
col = lambda al: vik(norm(al))


# ---- geometry of one shape -------------------------------------------------------
def smooth(xs):
    """resampled (step ≈ H) and smoothed centerline, its arclength"""
    s = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(xs, axis=0), axis=1))])
    m = max(int(round(s[-1] / a.h)), 4)
    su = np.linspace(0, s[-1], m + 1)
    xu = np.stack([np.interp(su, s, xs[:, j]) for j in range(3)], axis=1)
    return su, gaussian_filter1d(xu, a.sigma_s / su[1], axis=0, mode="nearest")


def profiles(xs):
    """κ and Q on NU points of s/L"""
    su, xu = smooth(xs)
    hh, L = su[1], su[-1]
    t = np.diff(xu, axis=0)
    t /= np.linalg.norm(t, axis=1)[:, None]
    # curvature: turning angle between consecutive tangents / step, at the interior nodes
    kap = np.arccos(np.clip(np.einsum("ij,ij->i", t[:-1], t[1:]), -1, 1)) / hh
    # Q measure: each window of n tangents fanned from its first tangent
    n = max(int(round(a.ell_w / hh)), 3)
    k = len(t) - n + 1
    ap = t[:k]
    Q = np.zeros(k)
    for j in range(1, n - 1):
        b, c = t[j:j + k], t[j + 1:j + 1 + k]
        num = np.einsum("ij,ij->i", np.cross(ap, b), c)
        den = 1 + np.einsum("ij,ij->i", ap, b) + np.einsum("ij,ij->i", b, c) + np.einsum("ij,ij->i", ap, c)
        Q += 2 * np.arctan2(num, den)
    u = np.linspace(0, 1, a.nu)
    return (np.interp(u * L, su[1:-1], kap),
            np.interp(u * L, (np.arange(k) + n / 2.0) * hh, Q), L)


# ---- all shapes -----------------------------------------------------------------
files = sorted(glob.glob(os.path.join(folder, "cutoff*.mat")))
tau = np.array([float(os.path.basename(f).split("_")[2][3:-4]) for f in files])
order = np.argsort(tau)
files, tau = [files[i] for i in order], tau[order]
alpha = 2 * tau - 1
K, Qm, Ls, corrupt = [], [], [], []
for f in files:
    with h5py.File(f, "r") as h:
        xs = h["xsaves"][()].T[::a.stride]
    # corrupt snapshots (non-finite, or a length far above the final one) cannot be resampled: such
    # runs get the profiles of their first shape and are left out (as in precompute_writhe.py)
    lens = np.linalg.norm(np.diff(xs, axis=1), axis=2).sum(axis=1)
    bad = (not np.isfinite(xs).all()) or lens.max() > 1.5 * lens[-1]
    corrupt.append(bad)
    pr = [profiles(xs[0] if bad else x) for x in xs]
    K.append([q[0] for q in pr]); Qm.append([q[1] for q in pr]); Ls.append([q[2] for q in pr])
K, Qm, Ls = np.array(K), np.array(Qm), np.array(Ls)          # runs × shapes × NU
nr, ns = K.shape[:2]
W = loadmat(os.path.join(folder, "wr_vs_len.mat"))
assert np.allclose(W["tau"].ravel(), tau) and W["wr"].shape[1] == ns, "wr_vs_len.mat does not match"
wr, excl = W["wr"], W["excluded"].ravel().astype(bool) | np.array(corrupt)

blocks = {"k": K, "q": Qm}
fit = ~excl
X = []
for key in a.features:
    B = blocks[key].reshape(nr * ns, -1)
    mu = B.reshape(nr, ns, -1)[fit].reshape(-1, B.shape[1]).mean(0)
    sc = np.sqrt(((B.reshape(nr, ns, -1)[fit].reshape(-1, B.shape[1]) - mu) ** 2).sum(1).mean())
    X.append((B - mu) / sc)
X = np.concatenate(X, axis=1)
Xf = X.reshape(nr, ns, -1)[fit].reshape(-1, X.shape[1])
mu = Xf.mean(0)
U, S, Vt = np.linalg.svd(Xf - mu, full_matrices=False)
var = S ** 2 / (S ** 2).sum()
P = ((X - mu) @ Vt.T).reshape(nr, ns, -1)                    # scores of every shape
# sign: PC1 > 0 for the α = +1 end at the last snapshot
for c in range(3):
    if P[np.argmax(alpha), -1, c] < 0:
        P[..., c] *= -1; Vt[c] *= -1
npc = int(np.searchsorted(np.cumsum(var), 0.9) + 1)
Lstar = a.at if a.at is not None else Ls[fit, -1].min()   # over the runs kept in the fit
iL = np.array([int(np.argmin(np.abs(Ls[r] - Lstar))) for r in range(nr)])   # snapshot at L* per run
Pst = P[np.arange(nr), iL]
wst = wr[np.arange(nr), iL]
Z = linkage(Pst[fit, :npc], "ward")
cl = fcluster(Z, a.nclust, "maxclust")
cl_all = np.zeros(nr, int); cl_all[fit] = cl
print(f"{folder}: {nr} runs × {ns} shapes, features {a.features}; "
      f"PC1..5 explain {', '.join(f'{v:.0%}' for v in var[:5])}; {npc} PCs for 90%")
print(f"  outcome at L* = {Lstar:.3f} L0")
for nc in sorted({2, 3, a.nclust}):
    cc = fcluster(Z, nc, "maxclust")
    print(f"  {nc} clusters: " + ";  ".join(
        f"{(cc == c).sum()} runs, α {alpha[fit][cc == c].min():+.2f}..{alpha[fit][cc == c].max():+.2f}"
        for c in sorted(set(cc), key=lambda c: alpha[fit][cc == c].mean())))

# ---- figure -----------------------------------------------------------------------
F = plt.figure(figsize=(15, 9.5), layout="constrained")
gs = F.add_gridspec(2, 4, width_ratios=[1, 1, 1, 0.06])
ax = {k: F.add_subplot(gs[i, j]) for k, (i, j) in zip("abcdef", [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)])}
draw = lambda r: dict(color=col(alpha[r]), lw=1.0 if not excl[r] else 0.8, ls="-" if not excl[r] else "--")

g = ax["a"]
g.bar(np.arange(1, 11), var[:10], color="0.4")
g.plot(np.arange(1, 11), np.cumsum(var[:10]), "k.-", lw=1)
g.set_xlabel("PC"); g.set_ylabel("explained variance (bars), cumulative (line)")
g.set_title(f"a  features: {'κ and Q' if a.features == 'kq' else ('κ' if a.features == 'k' else 'Q')}")
for key, c in (("b", 0), ("c", 1)):
    g = ax[key]
    for r in np.where(fit)[0]:
        g.plot(P[r, :, c], Ls[r], **draw(r))
    g.axvline(0, color="0.5", lw=0.6)
    g.set_xlabel(f"PC{c + 1}"); g.set_ylabel("Length")
    g.set_title(f"{key}  PC{c + 1} vs length")
g = ax["d"]
for r in np.where(fit)[0]:
    g.plot(P[r, :, 0], P[r, :, 1], **draw(r))
    g.plot(P[r, 0, 0], P[r, 0, 1], "o", mfc="none", mec=col(alpha[r]), ms=4)
    g.plot(P[r, -1, 0], P[r, -1, 1], "o", color=col(alpha[r]), ms=5, mec="k", mew=0.4)
    g.plot(Pst[r, 0], Pst[r, 1], "D", color=col(alpha[r]), ms=3.5, mec="0.2", mew=0.3)
g.set_xlabel("PC1"); g.set_ylabel("PC2"); g.set_title(f"d  trajectories (o start, ◆ L = {Lstar:.2f}, ● end)")
g = ax["e"]
u = np.linspace(0, 1, a.nu)
for c, ls in ((0, "-"), (1, "--")):
    for ib, key in enumerate(a.features):
        g.plot(u, Vt[c, ib * a.nu:(ib + 1) * a.nu], ls=ls, color="C0" if key == "k" else "C3",
               label=f"PC{c + 1}, {'κ' if key == 'k' else 'Q'}")
g.axhline(0, color="0.5", lw=0.6)
g.set_xlabel("s / L"); g.set_ylabel("loading"); g.legend(fontsize=8, frameon=False)
g.set_title("e  loadings")
g = ax["f"]
g2 = g.twinx()
marks = "osDv^<>p"
for r in range(nr):
    mk = marks[(cl_all[r] - 1) % len(marks)] if cl_all[r] else "x"
    g.plot(alpha[r], Pst[r, 0], mk, color=col(alpha[r]), mec="k", mew=0.4, ms=6)
    g2.plot(alpha[r], wst[r], "^", color="none", mec="0.35", mew=0.8, ms=5)
g.set_xlabel(r"initial asymmetry $\alpha$"); g.set_ylabel(f"PC1 at L = {Lstar:.2f} (marker = cluster)")
g2.set_ylabel(f"Wr at L = {Lstar:.2f} (△)", color="0.35"); g2.tick_params(axis="y", colors="0.35")
g.set_title(f"f  outcome at L = {Lstar:.2f} vs initial shape ({a.nclust} clusters)")
cax = F.add_subplot(gs[:, 3])
F.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=vik), cax=cax, label=r"$\alpha$ (−1: WT end, +1: byn>Myo1C end)")
F.suptitle(f"{os.path.basename(folder)}: PCA of smoothed rods (σ = {a.sigma_s:g} L0, Q window {a.ell_w:g} L0)", fontsize=11)
out = f"{os.path.basename(folder)}_pca_{a.features}.pdf"
F.savefig(out)
F.savefig(out.replace(".pdf", ".png"), dpi=110)
print("saved", out)


# ---- 2D embeddings of the same feature space: metric MDS, t-SNE -----------------------------
def embed(Y, ref, method, perplexity=30.0):
    """2D embedding of the rows of Y, started from ref (n × 2, PC1-PC2) and rotated/reflected onto it
    (orthogonal Procrustes; the axes of both methods are arbitrary)"""
    from scipy.linalg import orthogonal_procrustes
    if method == "mds":
        from sklearn.manifold import MDS
        E = MDS(n_components=2, metric=True, n_init=1, max_iter=500, random_state=0,
                dissimilarity="euclidean").fit_transform(Y, init=ref.copy())
    else:
        from sklearn.manifold import TSNE
        init = (ref - ref.mean(0)) / ref[:, 0].std() * 1e-4      # PCA initialisation, sklearn scaling
        E = TSNE(n_components=2, perplexity=min(perplexity, (len(Y) - 1) / 3), init=init,
                 learning_rate="auto", random_state=0).fit_transform(Y)
    E = E - E.mean(0)
    R, _ = orthogonal_procrustes(E, ref - ref.mean(0))
    return E @ R


def plot_embedding(method):
    name = {"mds": "metric MDS", "tsne": f"t-SNE (perplexity {a.perplexity:g})"}[method]
    lab = {"mds": "MDS", "tsne": "t-SNE"}[method]
    keep = np.where(fit)[0]                     # runs flagged for a writhe jump are left out
    E = np.full((nr, ns, 2), np.nan)
    E[keep] = embed(X.reshape(nr, ns, -1)[keep].reshape(-1, X.shape[1]), P[keep, :, :2].reshape(-1, 2),
                    method, a.perplexity).reshape(len(keep), ns, 2)
    Xst = X.reshape(nr, ns, -1)[np.arange(nr), iL]
    Est = np.full((nr, 2), np.nan)
    Est[keep] = embed(Xst[keep], Pst[keep, :2], method, a.perplexity_lstar)
    F = plt.figure(figsize=(15, 4.9), layout="constrained")
    gs = F.add_gridspec(1, 4, width_ratios=[1, 1, 1, 0.05])
    g = F.add_subplot(gs[0])
    for r in keep:
        g.plot(E[r, :, 0], E[r, :, 1], **draw(r))
        g.plot(E[r, 0, 0], E[r, 0, 1], "o", mfc="none", mec=col(alpha[r]), ms=4)
        g.plot(E[r, iL[r], 0], E[r, iL[r], 1], "D", color=col(alpha[r]), ms=3.5, mec="0.2", mew=0.3)
        g.plot(E[r, -1, 0], E[r, -1, 1], "o", color=col(alpha[r]), ms=5, mec="k", mew=0.4)
    g.set_xlabel(f"{lab} 1"); g.set_ylabel(f"{lab} 2"); g.set_aspect("equal", "datalim")
    g.set_title(f"a  all shapes (o start, ◆ L = {Lstar:.2f}, ● end)")
    g = F.add_subplot(gs[1])
    for r in keep:
        g.plot(E[r, :, 0], Ls[r], **draw(r))
    g.set_xlabel(f"{lab} 1"); g.set_ylabel("Length"); g.set_title(f"b  {lab} 1 vs length")
    g = F.add_subplot(gs[2])
    for r in keep:
        mk = marks[(cl_all[r] - 1) % len(marks)] if cl_all[r] else "x"
        g.plot(Est[r, 0], Est[r, 1], mk, color=col(alpha[r]), mec="k", mew=0.4, ms=7)
    for r in (0, int(np.argmin(np.abs(alpha))), nr - 1):
        g.annotate(rf"$\alpha$ = {alpha[r]:+.2f}", Est[r], textcoords="offset points", xytext=(6, 4), fontsize=8)
    g.set_xlabel(f"{lab} 1"); g.set_ylabel(f"{lab} 2"); g.set_aspect("equal", "datalim")
    g.set_title(f"c  shapes at L = {Lstar:.2f} alone ({a.nclust} clusters"
                + (f", perplexity {min(a.perplexity_lstar, (nr - 1) / 3):g})" if method == "tsne" else ")"))
    F.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=vik), cax=F.add_subplot(gs[3]),
               label=r"$\alpha$ (−1: WT end, +1: byn>Myo1C end)")
    F.suptitle(f"{os.path.basename(folder)}: {name} of the κ/Q embedding (features {a.features})", fontsize=11)
    out = f"{os.path.basename(folder)}_{method}_{a.features}.pdf"
    F.savefig(out)
    F.savefig(out.replace(".pdf", ".png"), dpi=110)
    print("saved", out)


for m in (["mds"] if a.mds else []) + (["tsne"] if a.tsne else []):
    plot_embedding(m)
