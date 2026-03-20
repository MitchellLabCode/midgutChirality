#!/usr/bin/env python3
"""Colorize a masked TIFF volume into an RGB uint16 multi-page TIFF.

Reads a masked TIFF (ZYX or ZYXC/CZYX-like), detects channel axis, applies
per-channel RGB colors and writes an RGB multi-page TIFF.

Usage: python colorize_masked_tiff.py --in Time_000005_WT_masked.tif

Minimal, standalone utility by NPMitchell 2026
"""
from pathlib import Path
import argparse
from typing import Tuple

import numpy as np
import tifffile


DEFAULT_IN = "Time_000005_WT_masked.tif"
DEFAULT_OUT = "Time_000005_WT_masked_rgb_py.tif"

# color choices (same as MATLAB script)
ENDODERM = np.array([100, 200, 255]) / 255.0
CVM = np.array([255, 200, 100]) / 255.0
LVM = np.array([1, 0, 0])
DEFAULT_CHANNEL_COLORS = (CVM, LVM, ENDODERM)


def detect_and_move_channel_axis(arr: np.ndarray) -> Tuple[np.ndarray, int]:
    """Return array with channel axis moved to last position and the original axis index.

    If input is 3D, returns same array and channel axis = -1 (interpreted as single channel).
    For 4D, chooses the axis with small size (<=8) as channel axis if available, else assumes last axis.
    """
    if arr.ndim == 3:
        return arr, -1

    if arr.ndim == 4:
        # candidate channel axes: small dimension (<=8)
        candidates = [i for i, s in enumerate(arr.shape) if s <= 8]
        ch_axis = candidates[0] if candidates else (arr.ndim - 1)
        if ch_axis == arr.ndim - 1:
            return arr, ch_axis
        moved = np.moveaxis(arr, ch_axis, -1)
        return moved, ch_axis

    raise ValueError(f"Unsupported array ndim {arr.ndim}; expected 3 or 4")


def colorize_channels_to_rgb_u16(masked_data_zyxc_u16: np.ndarray, channel_colors) -> np.ndarray:
    """Colorize channels into RGB uint16 volume (ZYX -> ZYX3 or ZYXC -> ZYX3).

    Accepts 3D (ZYX) single-channel or 4D (ZYXC) multi-channel uint16 input.
    """
    if masked_data_zyxc_u16.ndim == 3:
        gray = masked_data_zyxc_u16[..., None]
        return np.repeat(gray, 3, axis=-1)

    if masked_data_zyxc_u16.ndim != 4:
        raise ValueError(f"Expected ZYXC or ZYX input for RGB colorization, got {masked_data_zyxc_u16.shape}")

    num_channels = masked_data_zyxc_u16.shape[-1]
    colors = np.asarray(channel_colors, dtype=np.float32)
    if colors.ndim != 2 or colors.shape[1] != 3:
        raise ValueError("channel_colors must be a sequence of (R,G,B) triplets")
    if colors.shape[0] < num_channels:
        raise ValueError(f"Need at least {num_channels} channel colors, got {colors.shape[0]}")

    data_norm = masked_data_zyxc_u16.astype(np.float32) / 65535.0
    rgb = np.zeros(masked_data_zyxc_u16.shape[:3] + (3,), dtype=np.float32)
    for channel_index in range(num_channels):
        rgb += data_norm[..., channel_index, None] * colors[channel_index][None, None, None, :]

    rgb = np.clip(rgb, 0.0, 1.0)
    return np.rint(rgb * 65535.0).astype(np.uint16)


def main() -> None:
    p = argparse.ArgumentParser(description="Colorize masked TIFF to RGB uint16 TIFF")
    p.add_argument("--in", dest="infile", type=Path, default=Path(DEFAULT_IN))
    p.add_argument("--out", dest="outfile", type=Path, default=Path(DEFAULT_OUT))
    p.add_argument("--colors", nargs="*", help="Optional override RGB triplets as floats 0-1", default=None)
    args = p.parse_args()

    infile: Path = args.infile
    outfile: Path = args.outfile

    if not infile.exists():
        raise FileNotFoundError(f"Input file not found: {infile}")

    arr = tifffile.imread(infile)
    print(f"Read {infile}: shape={arr.shape} dtype={arr.dtype}")

    arr_moved, orig_ch_axis = detect_and_move_channel_axis(arr)
    if orig_ch_axis == -1:
        print("Detected single-channel volume (ZYX)")
    else:
        print(f"Detected channel axis at original index {orig_ch_axis}; moved to last axis")

    if arr_moved.dtype != np.uint16:
        # try to convert sensibly
        print(f"Converting input dtype {arr_moved.dtype} -> uint16")
        arr_moved = np.clip(np.rint(arr_moved), 0, 65535).astype(np.uint16)

    colors = DEFAULT_CHANNEL_COLORS
    if args.colors:
        # parse triples
        vals = [float(x) for x in args.colors]
        if len(vals) % 3 != 0:
            raise ValueError("--colors must be multiples of 3 floats (r g b ...)")
        colors = [tuple(vals[i : i + 3]) for i in range(0, len(vals), 3)]

    rgb_u16 = colorize_channels_to_rgb_u16(arr_moved, colors)

    # tifffile can write a stack of RGB frames; ensure shape is (Z, Y, X, 3)
    print(f"Writing RGB uint16 TIFF to {outfile}: shape={rgb_u16.shape} dtype={rgb_u16.dtype}")
    tifffile.imwrite(outfile, rgb_u16, photometric="rgb")
    print("Done.")


if __name__ == "__main__":
    main()
