"""Shared helpers for the plots of the tube_coil.jl runs: loading, colormaps, 3D tube rendering and
the two local chirality measures (κ²τ and the windowed Q measure).

Conventions
  α is oriented so that α > 0 is an initially right-handed shape (Wr > 0).
  Writhe relative to the tube axis: Wr = -Tw (the link of the clamped rod stays 0).
  Lengths in units of Le; L0 = initial length (1.05 Le).
"""
import glob
import os
import re

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LightSource, ListedColormap
from scipy.io import loadmat
from scipy.ndimage import gaussian_filter1d

here = os.path.dirname(os.path.abspath(__file__))
ETAS = [0.0, 10.0, 100.0, 1000.0]
R_DRAW = 0.035          # drawing radius of the centerline tube (Le)
# 3D views: (elev, azim) and the lab axes shown on the plot (x, y, z) axes
#   persp: oblique perspective, rod axis horizontal
#   above: from above in perspective, MATLAB view(az = -15°, el = 82°) of (lab z, lab y, lab x)
VIEWS = {"persp": ((24, -58), [2, 0, 1]), "above": ((82, -15 - 90), [2, 1, 0])}


# ---- data ------------------------------------------------------------------------
def load_runs(folder, at_LL0=None):
    """runs of one output folder of tube_coil.jl, as dicts with alpha, a (oriented α), eta, Lrest,
    Wr (= -Tw), xsaves (n_saves, nodes, 3), tube_radius; at_LL0: stop at the snapshot closest to
    L/L0 = at_LL0 (treated as the final shape)"""
    runs = []
    for f in sorted(glob.glob(os.path.join(folder, "eta*_a*.mat"))):
        with h5py.File(f, "r") as h:
            d = {k: np.squeeze(np.array(h[k]).T) for k in ["alpha", "eta", "Lrest", "Tw", "xsaves", "tube_radius"]}
        d["alpha"], d["eta"], d["tube_radius"] = float(d["alpha"]), float(d["eta"]), float(d["tube_radius"])
        d["Wr"] = -d["Tw"]
        if at_LL0 is not None:
            n = int(np.argmin(np.abs(d["Lrest"] / d["Lrest"][0] - at_LL0))) + 1
            for k in ("xsaves", "Lrest", "Wr", "Tw"):
                d[k] = d[k][:n]
        runs.append(d)
    if not runs:
        raise SystemExit(f"no runs in {folder}")
    return runs


def orient(*run_lists):
    """sets r["a"] = ±α so that α > 0 is initially right-handed (Wr(t = 0) > 0)"""
    allruns = [r for rs in run_lists for r in rs]
    o = sum(r["alpha"] * r["Wr"][0] for r in allruns)
    sgn = np.sign(o) if abs(o) > 1e-3 else 1.0
    for r in allruns:
        r["a"] = round(sgn * r["alpha"], 6)


def nearest_alphas(runs, targets=(-1, -1 / 3, 1 / 3, 1)):
    alphas = sorted({r["a"] for r in runs})
    return list(dict.fromkeys(min(alphas, key=lambda x: abs(x - t)) for t in targets))


def eta_label(eta):
    return "elastic" if eta == 0 else rf"$\eta\tau_E = {eta:g}$"


# ---- colormaps -------------------------------------------------------------------
def vik():
    """vik (Crameri), for α"""
    m = loadmat(os.path.join(here, "vik.mat"))
    return ListedColormap(next(v for k, v in m.items() if not k.startswith("__")))


def managua_br():
    """managua_br (managua_br.m): < 0 blue (left-handed), > 0 orange (right-handed), black at 0"""
    txt = open(os.path.join(here, "managua_br.m")).read()
    body = txt[txt.index("[") + 1: txt.rindex("]")]
    vals = np.array([float(v) for v in re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", body)])
    return ListedColormap(vals.reshape(-1, 3))


def alpha_colour(cmap, a):
    return cmap(0.5 * (a + 1))


# ---- local chirality ---------------------------------------------------------------
def arclength(xs):
    return np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(xs, axis=0), axis=1))])


def k2tau(xs, sigma=4.0):
    """curvature-weighted torsion κ²τ = (r'×r'')·r'''/|r'|⁶ (Le⁻³; > 0 right-handed) on the nodes
    of xs (n, 3): uniform arclength grid, Gaussian smoothing over sigma nodes (removes the kinks
    of the saved, linearly interpolated shapes from the third derivative)"""
    s = arclength(xs)
    su = np.linspace(0, s[-1], len(xs))
    xu = np.stack([np.interp(su, s, xs[:, j]) for j in range(3)], axis=1)
    xu = gaussian_filter1d(xu, sigma, axis=0, mode="nearest")
    d1 = np.gradient(xu, su, axis=0); d2 = np.gradient(d1, su, axis=0); d3 = np.gradient(d2, su, axis=0)
    val = np.einsum("ij,ij->i", np.cross(d1, d2), d3) / np.linalg.norm(d1, axis=1) ** 6
    return np.interp(s, su, val)


def qmeasure_profile(xs, ell_w=0.1, sigma_s=0.08, h=0.01):
    """Q measure: signed solid angle (sr; > 0 right-handed) enclosed by the unit tangents over
    sliding windows of fixed arclength ell_w. The centerline is resampled at spacing ≈ h and
    Gaussian-smoothed over sigma_s (all in Le, so the measure does not depend on the node spacing
    of the saved shapes); each window of n = ell_w/h tangents t_k .. t_k+n-1 is fanned from its
    own first tangent, Ω_k = Σ_j 2 atan2(t_k·(t_j×t_j+1), 1 + t_k·t_j + t_j·t_j+1 + t_j+1·t_k).
    For a short window Ω ≈ κ²τ ell_w³/12. Returns (s, Ω) on the resampled grid (value at each
    window centre, constant within ell_w/2 of the ends)."""
    xs = np.asarray(xs, float)
    s = arclength(xs)
    m = max(int(round(s[-1] / h)), 2)
    su = np.linspace(0, s[-1], m + 1)
    hh = su[1]
    xu = np.stack([np.interp(su, s, xs[:, j]) for j in range(3)], axis=1)
    xu = gaussian_filter1d(xu, sigma_s / hh, axis=0, mode="nearest")
    t = np.diff(xu, axis=0)
    t /= np.linalg.norm(t, axis=1)[:, None]
    n = max(int(round(ell_w / hh)), 3)
    k = len(t) - n + 1                       # number of windows
    a = t[:k]                                # apex of each window
    A = np.zeros(k)
    for j in range(1, n - 1):                # fan triangles (t_k, t_k+j, t_k+j+1), vectorized over k
        b, c = t[j:j + k], t[j + 1:j + 1 + k]
        num = np.einsum("ij,ij->i", np.cross(a, b), c)
        den = 1 + np.einsum("ij,ij->i", a, b) + np.einsum("ij,ij->i", b, c) + np.einsum("ij,ij->i", a, c)
        A += 2 * np.arctan2(num, den)
    centre = (np.arange(k) + n / 2.0) * hh   # arclength of the window centres
    return su, np.interp(su, centre, A)


def qmeasure(xs, ell_w=0.1, sigma_s=0.08):
    """Q measure (see qmeasure_profile) on the nodes of xs"""
    su, v = qmeasure_profile(xs, ell_w, sigma_s)
    return np.interp(arclength(xs), su, v)


MEASURES = {"k2tau": (k2tau, r"$\kappa^2\tau$ ($L_e^{-3}$)"),
            "qmeasure": (qmeasure, "Q measure,\n0.1 $L_e$ windows (sr)")}


# ---- 3D rendering ------------------------------------------------------------------
def tube_surface(P, radius, nphi=14):
    """surface of a tube of given radius around the polyline P (n, 3), parallel-transported frame"""
    t = np.gradient(P, axis=0)
    t /= np.linalg.norm(t, axis=1)[:, None]
    n = np.cross(t[0], [1.0, 0.0, 0.0] if abs(t[0][0]) < 0.9 else [0.0, 1.0, 0.0])
    n /= np.linalg.norm(n)
    N = [n]
    for i in range(1, len(P)):
        n = N[-1] - np.dot(N[-1], t[i]) * t[i]
        N.append(n / np.linalg.norm(n))
    N = np.array(N)
    Bn = np.cross(t, N)
    phi = np.linspace(0, 2 * np.pi, nphi)
    return P[:, None, :] + radius * (np.cos(phi)[None, :, None] * N[:, None, :] + np.sin(phi)[None, :, None] * Bn[:, None, :])


def draw_rod(ax, xs, view, zlim, extent, color=None, values=None, vnorm=None, vcmap=None, R_wall=None, opacity=1.0):
    """the centerline xs as a shaded tube, in one colour or coloured along its length by values
    (vnorm, vcmap), with the hard wall (radius R_wall) as a translucent cylinder"""
    (elev, azim), axes = VIEWS[view]
    P = xs[:, axes]
    S = tube_surface(P, R_DRAW)
    if values is None:
        ax.plot_surface(S[..., 0], S[..., 1], S[..., 2], color=color, linewidth=0, antialiased=True, alpha=opacity,
                        shade=True, lightsource=LightSource(azdeg=225, altdeg=60), rstride=1, cstride=1)
    else:                                    # colours along the rod, Lambert shading from the normals
        nrm = (S - P[:, None, :]) / R_DRAW
        light = np.array([-0.35, -0.55, 0.76]); light /= np.linalg.norm(light)
        shade = 0.5 + 0.5 * np.clip(nrm @ light, 0, 1)
        rgb = vcmap(vnorm(values))[:, None, :3] * shade[..., None]
        fc = np.concatenate([rgb, np.full(rgb.shape[:2] + (1,), opacity)], axis=2)
        ax.plot_surface(S[..., 0], S[..., 1], S[..., 2], facecolors=fc, linewidth=0, antialiased=True,
                        shade=False, rstride=1, cstride=1)
    if R_wall is not None:
        Z, PH = np.meshgrid(np.array(zlim), np.linspace(0, 2 * np.pi, 40))
        ax.plot_surface(Z, R_wall * np.cos(PH), R_wall * np.sin(PH), color="0.6", alpha=0.12, linewidth=0, shade=False)
    ax.set_xlim(*zlim); ax.set_ylim(-extent, extent); ax.set_zlim(-extent, extent)
    ax.set_box_aspect((zlim[1] - zlim[0], 2 * extent, 2 * extent))
    ax.view_init(elev=elev, azim=azim)
    ax.set_proj_type("persp", focal_length=0.45)
    ax.set_axis_off()


def savefig(fig, out):
    fig.savefig(out, bbox_inches="tight", dpi=200)
    print("saved", out)
