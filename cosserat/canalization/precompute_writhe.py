"""Writhe vs L/L0 for every run of a canalization folder, for plot_canalization_frames.m.

MATLAB has no writhe routine, so this computes it with the formula of
../src/measurements.py (measurements._writhe: Klenin & Langowski 2000,
method 1a, with NaN pair contributions set to 0), vectorized over segment pairs.
--verify checks the vectorized version against measurements.compute_wr.

Conventions:
    tau     parsed from the filename (cutoff01contact_n<k>_tau<tau>.mat)
    alpha   = 2*tau - 1   (alpha = -1: tau = 0, the WT end; +1: tau = 1)
    Wr      every `stride`-th snapshot (default 4)
    L/L0    summed edge length; positions are in units of the initial length

Also flags the run with an outlying jump in Wr(t) (max |dWr| between
consecutive samples, compared with the median over runs) so it can be left out
of the writhe panel.

    python precompute_writhe.py 20261110_canalization_box3f_s6s3_Lnorm [--verify]
"""

import os
import sys
import numpy as np
import h5py
import scipy.io as sio


def writhe(x):
    """Writhe of a polyline x (N, 3). Same as measurements._writhe, vectorized."""
    N = x.shape[0]
    i, j = np.triu_indices(N - 1, k=1)          # segment pairs i < j, as in _writhe
    i, j = i[i <= N - 3], j[i <= N - 3]
    p1, p2, p3, p4 = x[i], x[i + 1], x[j], x[j + 1]
    r12, r34 = p2 - p1, p4 - p3
    r13, r14, r23, r24 = p3 - p1, p4 - p1, p3 - p2, p4 - p2

    def unit(v):
        with np.errstate(invalid="ignore", divide="ignore"):
            return v / np.linalg.norm(v, axis=1, keepdims=True)

    n1 = unit(np.cross(r13, r14))
    n2 = unit(np.cross(r14, r24))
    n3 = unit(np.cross(r24, r23))
    n4 = unit(np.cross(r23, r13))
    dot = lambda a, b: np.einsum("ij,ij->i", a, b)
    with np.errstate(invalid="ignore"):
        om = (np.arcsin(dot(n1, n2)) + np.arcsin(dot(n2, n3))
              + np.arcsin(dot(n3, n4)) + np.arcsin(dot(n4, n1)))
    om = np.where(np.isnan(om), 0.0, om)        # _writhe: NaN -> 0
    sgn = np.sign(dot(np.cross(r34, r12), r13))
    return float(np.sum(om * sgn) / (2 * np.pi))


def main(folder, stride=4, verify=False, jump_factor=5.0):
    files = sorted(f for f in os.listdir(folder)
                   if f.startswith("cutoff") and f.endswith(".mat"))
    tau, wr, ln = [], [], []
    for fn in files:
        tau.append(float(fn.split("_")[2][3:-4]))          # tau from the file name
        with h5py.File(os.path.join(folder, fn), "r") as f:
            xs = np.array(f["xsaves"]).T                    # (nt, N+1, 3)
        xs = xs[::stride]
        wr.append([writhe(xs[t]) for t in range(xs.shape[0])])
        ln.append(np.linalg.norm(xs[:, 1:] - xs[:, :-1], axis=2).sum(axis=1))
        if verify and fn == files[0]:
            sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
            from measurements import compute_wr
            ref = compute_wr(xs[[0, xs.shape[0] // 2, -1]])
            got = np.array(wr[-1])[[0, xs.shape[0] // 2, -1]]
            print(f"verify vs measurements.compute_wr: max|diff| = "
                  f"{np.abs(ref - got).max():.2e}")

    order = np.argsort(tau)
    tau = np.array(tau)[order]
    wr = np.array(wr)[order]
    ln = np.array(ln)[order]
    files = [files[k] for k in order]
    alpha = 2 * tau - 1

    jump = np.abs(np.diff(wr, axis=1)).max(axis=1)
    med = np.median(jump)
    # a step of ~2 in Wr between samples is a self-passage (the rod crossing itself); the relative test
    # alone misses them when most runs of a sweep have one, so any step > 1 is flagged as well
    excluded = (jump > jump_factor * med) | (jump > 1.0)
    # also runs with corrupt snapshots: non-finite values, or a length well above the final one
    badlen = ~np.isfinite(wr).all(axis=1) | ~np.isfinite(ln).all(axis=1) | (np.nanmax(ln, axis=1) > 1.5 * ln[:, -1])
    if badlen.any():
        print(f"  {badlen.sum()} runs with corrupt snapshots (length > 1.5x final or non-finite), excluded: alpha = "
              + ", ".join(f"{v:+.3f}" for v in alpha[badlen]))
    excluded = excluded | badlen
    snap_idx = np.arange(0, wr.shape[1]) * stride + 1        # 1-based, for MATLAB

    print(f"{folder}: {len(files)} runs, {wr.shape[1]} samples each "
          f"(snapshots 1:{stride}:{snap_idx[-1]})")
    print(f"  L/L0 {ln[:, 0].mean():.3f} -> {ln[:, -1].min():.2f}..{ln[:, -1].max():.2f}")
    print(f"  max |dWr| per run: median {med:.4f}; flagged > {jump_factor:g}x median:")
    for k in np.argsort(-jump)[:4]:
        mark = "  <- EXCLUDED" if excluded[k] else ""
        print(f"    alpha={alpha[k]:+.3f} (tau={tau[k]:.4f})  max|dWr|={jump[k]:.4f}"
              f" ({jump[k] / med:4.1f}x){mark}")

    out = os.path.join(folder, "wr_vs_len.mat")
    sio.savemat(out, dict(tau=tau, alpha=alpha, wr=wr, len_ratio=ln,
                          snap_idx=snap_idx, excluded=excluded,
                          files=np.array(files, dtype=object), stride=stride))
    print(f"  wrote {out}")


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    main(args[0] if args else "20261110_canalization_box3f_s6s3_Lnorm",
         verify="--verify" in sys.argv)
