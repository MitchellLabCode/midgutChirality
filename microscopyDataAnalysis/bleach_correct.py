#!/usr/bin/env python3
"""
Bleach correction for sparse-pulse dI/I time series, prior to temporal color coding.

For each .tif/.tiff (T, Y, X) stack in the input folder:
  1. Compute the P-th percentile of pixel values in every frame (default P = 99.9).
     Pulses are sparse on a dark background, so a high percentile tracks pulse
     brightness rather than background.
  2. Fit a single exponential  P(t) = A * exp(-k t)  as a straight line to
     log P(t)  (least squares on the log values).
  3. Rescale every frame:  I_corr(t) = I(t) * F(0) / F(t),  F = fitted curve.
     This is a pure multiplicative correction (dI/I background ~0 stays ~0).
     The smooth fit is used rather than the raw per-frame percentile, so a frame
     that happens to contain more or brighter pulses is NOT flattened away.

Outputs (in <input>/bleachCorrected/):
  <name>_bleachCorr.tif   float32 ImageJ hyperstack, same frames / calibration
  <name>_bleachFit.png    log-percentile vs. time with fit + correction factor
  <name>_bleachFit.csv    per-frame percentile, fit, and correction factor

Usage:
  python3 bleach_correct.py                       # default folder below
  python3 bleach_correct.py /path/to/folder --percentile 99.5

Needs: numpy, tifffile, matplotlib
"""

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
import tifffile
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DEFAULT_DIR = Path(__file__).resolve().parent / "30independentFrames"


def read_stack(path):
    """Return (data[T,Y,X] float32, imagej metadata dict, (xres, yres) or None)."""
    with tifffile.TiffFile(path) as tif:
        series = tif.series[0]
        data = series.asarray()
        axes = series.axes
        ij = dict(tif.imagej_metadata or {})
        page = tif.pages[0]
        res = None
        if "XResolution" in page.tags and "YResolution" in page.tags:
            xr = page.tags["XResolution"].value
            yr = page.tags["YResolution"].value
            res = (xr[0] / xr[1], yr[0] / yr[1])

    data = np.squeeze(data)
    if data.ndim != 3:
        raise ValueError(f"expected a single-channel time series (T,Y,X); got axes "
                         f"'{axes}' with shape {series.shape}")
    return data.astype(np.float32, copy=False), ij, res


def fit_log_exponential(t, p):
    """Least-squares line through log(p). Returns (A, k, mask of frames used)."""
    used = np.isfinite(p) & (p > 0)
    if used.sum() < 2:
        raise ValueError("fewer than 2 frames have a positive percentile value; "
                         "try a higher --percentile")
    slope, intercept = np.polyfit(t[used], np.log(p[used]), 1)
    return np.exp(intercept), -slope, used


def process(path, out_dir, pct):
    print(f"Processing {path.name}")
    data, ij, res = read_stack(path)
    n = data.shape[0]
    t = np.arange(n, dtype=float)                       # frame index (0-based)

    p = np.array([np.percentile(frame, pct) for frame in data])
    A, k, used = fit_log_exponential(t, p)
    fit = A * np.exp(-k * t)
    factor = fit[0] / fit                               # 1.0 at first frame

    if (~used).any():
        print(f"  note: {(~used).sum()} frame(s) had percentile <= 0 and were "
              f"left out of the fit")
    if k <= 0:
        print(f"  WARNING: fitted decay rate k = {k:.4g} per frame is not positive "
              f"(no bleaching detected); correction will dim later frames")

    finterval = ij.get("finterval")
    unit_t = ij.get("tunit", "s")
    half_life = np.log(2) / k if k > 0 else np.inf
    msg = f"  k = {k:.4g} /frame, half-life = {half_life:.3g} frames"
    if finterval:
        msg += f" (= {half_life * finterval:.3g} {unit_t})"
    print(msg + f"; last-frame correction factor = {factor[-1]:.3g}x")

    corrected = data * factor[:, None, None].astype(np.float32)

    # --- write corrected stack as an ImageJ hyperstack ---
    meta = {"axes": "TYX"}
    for key in ("unit", "finterval", "tunit"):
        if key in ij:
            meta[key] = ij[key]
    out_tif = out_dir / f"{path.stem}_bleachCorr.tif"
    kwargs = {"imagej": True, "metadata": meta}
    if res:
        kwargs["resolution"] = res
    tifffile.imwrite(out_tif, corrected.astype(np.float32), **kwargs)
    print(f"  saved {out_tif}")

    # --- CSV ---
    out_csv = out_dir / f"{path.stem}_bleachFit.csv"
    with open(out_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["frame", f"p{pct}", "fit", "correction_factor",
                    f"p{pct}_after_correction"])
        for i in range(n):
            w.writerow([i + 1, p[i], fit[i], factor[i], p[i] * factor[i]])

    # --- plot ---
    x = t * finterval if finterval else t + 1
    xlabel = f"time ({unit_t})" if finterval else "frame"
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 7), sharex=True)
    ax1.semilogy(x[used], p[used], "o", color="0.3", label=f"{pct}th percentile")
    ax1.semilogy(x, fit, "-", color="tab:red",
                 label=f"fit: A·exp(−k·t), k = {k:.3g}/frame")
    ax1.semilogy(x[used], (p * factor)[used], "s", ms=4, mfc="none",
                 color="tab:blue", label="after correction")
    ax1.set_ylabel(f"{pct}th percentile (log scale)")
    ax1.set_title(path.name, fontsize=9)
    ax1.legend(fontsize=8)
    ax2.plot(x, factor, "-o", ms=3, color="tab:red")
    ax2.set_ylabel("correction factor F(0)/F(t)")
    ax2.set_xlabel(xlabel)
    fig.tight_layout()
    out_png = out_dir / f"{path.stem}_bleachFit.png"
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    print(f"  saved {out_png}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folder", nargs="?", type=Path, default=DEFAULT_DIR,
                    help=f"folder of .tif stacks (default: {DEFAULT_DIR})")
    ap.add_argument("--percentile", type=float, default=99.9,
                    help="per-frame percentile to track (default 99.9)")
    ap.add_argument("--outdir", type=Path, default=None,
                    help="output folder (default: <folder>/bleachCorrected)")
    args = ap.parse_args()

    folder = args.folder.expanduser().resolve()
    files = sorted(f for f in folder.iterdir()
                   if f.is_file() and f.suffix.lower() in (".tif", ".tiff")
                   and not f.name.startswith("."))
    if not files:
        sys.exit(f"No .tif/.tiff files in {folder}")

    out_dir = args.outdir or folder / "bleachCorrected"
    out_dir.mkdir(parents=True, exist_ok=True)

    for f in files:
        try:
            process(f, out_dir, args.percentile)
        except Exception as e:  # keep going with the other files
            print(f"  ERROR on {f.name}: {e}")


if __name__ == "__main__":
    main()
