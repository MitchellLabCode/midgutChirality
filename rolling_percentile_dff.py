#!/usr/bin/env python3
"""Compute rolling-percentile dF/F from a TIFF stack.

For each pixel (x, y), computes a rolling 10th-percentile baseline
I_10(x, y, t) from a symmetric window [t-dT, t+dT], then returns:

    dF/F(x, y, t) = (I(x, y, t) - I_10(x, y, t)) / I_10(x, y, t)

Usage
-----
    python rolling_percentile_dff.py --in MIPS_dFF_right_early.tif --dT 20
    python rolling_percentile_dff.py --in foo.tif --dT 20 --percentile 10 --out foo_dff.tif

Outputs a float32 TIFF stack with the same T×Y×X layout as the input.

NPMitchell 2026
"""
from pathlib import Path
import argparse

import numpy as np
import tifffile
from scipy.ndimage import percentile_filter
from PIL import Image, ImageDraw, ImageFont


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def load_tyx(path: Path) -> np.ndarray:
    """Load TIFF as float32 array with shape (T, Y, X)."""
    data = tifffile.imread(str(path)).astype(np.float32)
    if data.ndim == 2:
        data = data[np.newaxis, ...]          # single frame → (1, Y, X)
    elif data.ndim != 3:
        raise ValueError(f"Expected 2-D or 3-D TIFF, got shape {data.shape}")
    return data  # (T, Y, X)


def rolling_percentile_baseline(stack: np.ndarray, dT: int, pct: float) -> np.ndarray:
    """Return rolling percentile baseline with the same shape as *stack*.

    Parameters
    ----------
    stack : (T, Y, X) float32
    dT    : half-window size in frames (full window = 2*dT + 1)
    pct   : percentile value (0–100)

    Uses scipy percentile_filter with a footprint that spans only the
    time axis, so each pixel's baseline is computed independently.
    """
    T, Y, X = stack.shape
    window = 2 * dT + 1
    # footprint: (window, 1, 1) → filter only along the time axis
    footprint = np.zeros((window, 1, 1), dtype=bool)
    footprint[:, 0, 0] = True
    baseline = percentile_filter(stack, percentile=pct, footprint=footprint, mode='nearest')
    return baseline


def compute_dff(stack: np.ndarray, baseline: np.ndarray, eps: float = 1e-6) -> np.ndarray:
    """Return (stack - baseline) / baseline, guarded against near-zero baseline."""
    safe_baseline = np.where(np.abs(baseline) < eps, eps, baseline)
    return (stack - safe_baseline) / safe_baseline


_FONT_CANDIDATES = [
    "/System/Library/Fonts/Helvetica.ttc",
    "/System/Library/Fonts/Monaco.ttf",
    "/Library/Fonts/Arial.ttf",
    "/usr/share/fonts/truetype/liberation/LiberationMono-Regular.ttf",
    "/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf",
]


def _load_font(size: int) -> ImageFont.ImageFont:
    for path in _FONT_CANDIDATES:
        if Path(path).exists():
            try:
                return ImageFont.truetype(path, size=size)
            except Exception:
                pass
    try:
        return ImageFont.load_default(size=size)   # Pillow >= 10
    except TypeError:
        return ImageFont.load_default()


def burn_timestamps(stack: np.ndarray, dt_seconds: float,
                    margin_px: int = 12,
                    font_size: int = 20) -> np.ndarray:
    """Burn a mm:ss timestamp into the top-left corner of every frame."""
    out = stack.copy()
    bar_value = float(np.nanpercentile(stack, 99))
    T, Y, X = out.shape
    font = _load_font(font_size)

    for t in range(T):
        total_sec = t * dt_seconds
        mm = int(total_sec // 60)
        ss = int(total_sec % 60)
        text = f"{mm:02d}:{ss:02d}"

        img = Image.new("L", (X, Y), 0)
        ImageDraw.Draw(img).text((margin_px, margin_px), text, fill=255, font=font)
        mask = np.array(img) > 0
        out[t][mask] = bar_value

    return out


def burn_scale_bar(stack: np.ndarray, um_per_pixel: float,
                   scale_bar_um: float = 100.0,
                   bar_height_px: int = 5,
                   margin_px: int = 12) -> np.ndarray:
    """Burn a filled scale bar into the bottom-right corner of every frame.

    The bar value is set to the 99th-percentile of the stack so it appears
    bright regardless of the dF/F range.
    """
    out = stack.copy()
    bar_value = float(np.nanpercentile(stack, 99))
    bar_px = max(1, round(scale_bar_um / um_per_pixel))
    T, Y, X = out.shape

    y0 = Y - margin_px - bar_height_px
    y1 = Y - margin_px
    x1 = X - margin_px
    x0 = x1 - bar_px

    if x0 < 0:
        raise ValueError(
            f"Scale bar ({bar_px} px = {scale_bar_um} µm) wider than image ({X} px). "
            "Reduce --scale_bar_um."
        )

    print(f"Scale bar: {scale_bar_um:.0f} µm = {bar_px} px  "
          f"(value={bar_value:.4f}, corner [{y0}:{y1}, {x0}:{x1}])")
    out[:, y0:y1, x0:x1] = bar_value
    return out


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

DEFAULT_IN = "/Users/npmitchell/Dropbox/Soft_Matter_Biophysics/PAPER/gut_chirality/calcium/MIPS_dFF_right_early.tif"


def parse_args():
    p = argparse.ArgumentParser(description="Rolling-percentile dF/F for a TIFF stack.")
    p.add_argument("--in", dest="input", default=DEFAULT_IN,
                   help="Path to input TIFF stack (T×Y×X). Default: %(default)s")
    p.add_argument("--out", dest="output", default=None,
                   help="Path for output TIFF. Default: <input stem>_dff.tif")
    p.add_argument("--dT", type=int, default=10,
                   help="Half-window size in frames (full window = 2*dT+1). Default: %(default)s")
    p.add_argument("--percentile", type=float, default=10.0,
                   help="Percentile for baseline (0–100). Default: %(default)s")
    p.add_argument("--eps", type=float, default=1e-6,
                   help="Small value to avoid division by zero. Default: %(default)s")
    p.add_argument("--um_per_pixel", type=float, default=0.195,
                   help="Microns per pixel, used to size the scale bar. Default: %(default)s")
    p.add_argument("--scale_bar_um", type=float, default=100.0,
                   help="Scale bar length in microns. Default: %(default)s")
    p.add_argument("--no_scale_bar", action="store_true",
                   help="Skip burning in the scale bar.")
    p.add_argument("--dt_seconds", type=float, default=5.0,
                   help="Time between frames in seconds, used for the timestamp. Default: %(default)s")
    p.add_argument("--font_size", type=int, default=20,
                   help="Font size in pixels for the timestamp. Default: %(default)s")
    p.add_argument("--no_timestamp", action="store_true",
                   help="Skip burning in the timestamp.")
    return p.parse_args()


def main():
    args = parse_args()

    in_path = Path(args.input)
    out_path = Path(args.output) if args.output else in_path.with_name(in_path.stem + f"_dff_dT{args.dT}.tif")

    print(f"Loading {in_path} ...")
    stack = load_tyx(in_path)
    T, Y, X = stack.shape
    print(f"  shape: T={T}, Y={Y}, X={X}  (dtype={stack.dtype})")

    print(f"Computing rolling {args.percentile:.0f}th-percentile baseline (dT={args.dT} frames) ...")
    baseline = rolling_percentile_baseline(stack, dT=args.dT, pct=args.percentile)

    print("Computing dF/F ...")
    dff = compute_dff(stack, baseline, eps=args.eps)

    print(f"dF/F range: [{dff.min():.4f}, {dff.max():.4f}]")

    if not args.no_scale_bar:
        dff = burn_scale_bar(dff, um_per_pixel=args.um_per_pixel,
                             scale_bar_um=args.scale_bar_um)

    if not args.no_timestamp:
        print(f"Burning timestamps (dt={args.dt_seconds} s/frame, font size={args.font_size} px) ...")
        dff = burn_timestamps(dff, dt_seconds=args.dt_seconds, font_size=args.font_size)

    print(f"Saving to {out_path} ...")
    tifffile.imwrite(str(out_path), dff, photometric='minisblack')
    print("Done.")


if __name__ == "__main__":
    main()
