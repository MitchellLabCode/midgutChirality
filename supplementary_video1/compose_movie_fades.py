#!/usr/bin/env python3
"""Compose supplementary movie parts with manual intensity scaling and fades.

This high-level script stitches three precomputed RGB TIFF stacks into one movie stack:

1. Embryo flythrough stack:
   - Input: ``Time_000005_WT_masked_embryo_rgb_anatomical_labeled.tif``
   - Every frame is scaled by ``scale_embryo_flythrough``.

2. Embryo rotation stack:
   - Input: ``Time_000005_WT_masked_embryo_rgb_anatomical_labeled_xrot_mip_stack_fade_s220_e280.tif``
   - Every frame is scaled by ``scale_embryo_rotation``.
   - Before this stack, a fade is generated from black to its first frame.

3. Midgut rotation stack:
   - Input: ``Time_000005_WT_masked_midgutOnly_rgb_anatomical_labeled_xrot_mip_stack_fade_s220_e280.tif``
   - Every frame is scaled by ``scale_midgut_rotation``.
   - A cross-fade is generated from the last embryo-rotation frame to the first
     midgut-rotation frame.

Output is a single RGB uint16 TIFF stack.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import tifffile


DEFAULT_DATADIR = Path(
    "/Users/npmitchell/Dropbox/Soft_Matter_Biophysics/PAPER/gut_chirality/supplementary_videos/"
    "20240531_142410_22C_HandGFPbynGAL4klar_UASmChCAAXH2AviRFP"
)

DEFAULT_EMBRYO_FLYTHROUGH_FILE = "Time_000005_WT_masked_embryo_rgb_anatomical_labeled.tif"
DEFAULT_EMBRYO_ROTATION_FILE = "Time_000005_WT_masked_embryo_rgb_anatomical_labeled_xrot_mip_stack_fade_s220_e280.tif"
DEFAULT_MIDGUT_ROTATION_FILE = "Time_000005_WT_masked_midgutOnly_rgb_anatomical_labeled_xrot_mip_stack_fade_s220_e280.tif"
DEFAULT_OUTPUT_FILE = "Time_000005_supplementary_movie_with_fades.tif"

DEFAULT_SCALE_EMBRYO_FLYTHROUGH = 2.0
DEFAULT_SCALE_EMBRYO_ROTATION = 1.0
DEFAULT_SCALE_MIDGUT_ROTATION = 1.0

DEFAULT_FADE_BLACK_TO_EMBRYO_ROT_FRAMES = 24
DEFAULT_FADE_EMBRYO_ROT_TO_MIDGUT_FRAMES = 24
DEFAULT_FLYTHROUGH_STEP = 4


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compose a supplementary RGB movie from embryo flythrough, embryo rotation, and "
            "midgut rotation stacks with manual per-part intensity scales and fades."
        )
    )
    parser.add_argument("--datadir", type=Path, default=DEFAULT_DATADIR)
    parser.add_argument("--embryo-flythrough-file", default=DEFAULT_EMBRYO_FLYTHROUGH_FILE)
    parser.add_argument("--embryo-rotation-file", default=DEFAULT_EMBRYO_ROTATION_FILE)
    parser.add_argument("--midgut-rotation-file", default=DEFAULT_MIDGUT_ROTATION_FILE)
    parser.add_argument("--output-file", default=DEFAULT_OUTPUT_FILE)

    parser.add_argument(
        "--scale-embryo-flythrough",
        type=float,
        default=DEFAULT_SCALE_EMBRYO_FLYTHROUGH,
        help="Multiply all flythrough frames by this factor.",
    )
    parser.add_argument(
        "--scale-embryo-rotation",
        type=float,
        default=DEFAULT_SCALE_EMBRYO_ROTATION,
        help="Multiply all embryo rotation frames by this factor.",
    )
    parser.add_argument(
        "--scale-midgut-rotation",
        type=float,
        default=DEFAULT_SCALE_MIDGUT_ROTATION,
        help="Multiply all midgut rotation frames by this factor.",
    )

    parser.add_argument(
        "--flythrough-step",
        type=int,
        default=DEFAULT_FLYTHROUGH_STEP,
        help="Use every Nth frame from the embryo flythrough stack (1 = all frames).",
    )
    parser.add_argument(
        "--fade-black-to-embryo-rot-frames",
        type=int,
        default=DEFAULT_FADE_BLACK_TO_EMBRYO_ROT_FRAMES,
        help="Number of frames in the fade from black to first embryo-rotation frame.",
    )
    parser.add_argument(
        "--fade-embryo-rot-to-midgut-frames",
        type=int,
        default=DEFAULT_FADE_EMBRYO_ROT_TO_MIDGUT_FRAMES,
        help="Number of frames in the cross-fade from last embryo-rotation frame to first midgut frame.",
    )

    return parser.parse_args()


def load_rgb_stack(path: Path) -> np.ndarray:
    """Load a TIFF stack and return data as (T, Y, X, 3)."""
    arr = tifffile.imread(path)
    if arr.ndim != 4:
        raise ValueError(f"Expected 4D RGB stack, got {arr.shape} from {path}")

    if arr.shape[-1] == 3:
        return arr
    if arr.shape[1] == 3:
        return np.transpose(arr, (0, 2, 3, 1))
    if arr.shape[0] == 3:
        return np.transpose(arr, (1, 2, 3, 0))

    raise ValueError(
        "Could not identify RGB channel axis (length 3); "
        f"got shape {arr.shape} from {path}"
    )


def scale_stack(stack_tyxc: np.ndarray, scale: float) -> np.ndarray:
    """Apply manual intensity scaling and clip to uint16 range."""
    scaled = stack_tyxc.astype(np.float32, copy=False) * float(scale)
    return np.clip(np.rint(scaled), 0, 65535).astype(np.uint16)


def fade_black_to_frame(first_frame_yxc: np.ndarray, n_frames: int) -> np.ndarray:
    """Generate fade frames from black to a target frame.

    The target frame itself is not included, so it can be appended once as part
    of its source stack without duplication.
    """
    if n_frames <= 0:
        return np.zeros((0,) + first_frame_yxc.shape, dtype=np.uint16)

    target = first_frame_yxc.astype(np.float32, copy=False)
    alphas = np.linspace(0.0, 1.0, n_frames + 1, dtype=np.float32)[1:-1]
    if alphas.size == 0:
        return np.zeros((0,) + first_frame_yxc.shape, dtype=np.uint16)

    frames = [np.clip(np.rint(alpha * target), 0, 65535).astype(np.uint16) for alpha in alphas]
    return np.stack(frames, axis=0)


def crossfade_frames(from_frame_yxc: np.ndarray, to_frame_yxc: np.ndarray, n_frames: int) -> np.ndarray:
    """Generate cross-fade frames between two RGB frames.

    Endpoints are excluded so caller can keep source stack endpoints unchanged.
    """
    if n_frames <= 0:
        return np.zeros((0,) + from_frame_yxc.shape, dtype=np.uint16)

    frm = from_frame_yxc.astype(np.float32, copy=False)
    to = to_frame_yxc.astype(np.float32, copy=False)
    alphas = np.linspace(0.0, 1.0, n_frames + 2, dtype=np.float32)[1:-1]
    if alphas.size == 0:
        return np.zeros((0,) + from_frame_yxc.shape, dtype=np.uint16)

    frames = [
        np.clip(np.rint((1.0 - alpha) * frm + alpha * to), 0, 65535).astype(np.uint16)
        for alpha in alphas
    ]
    return np.stack(frames, axis=0)


def ensure_same_frame_shape(*stacks: np.ndarray) -> None:
    """Require that all stacks share identical frame dimensions/channels."""
    shapes = [arr.shape[1:] for arr in stacks]
    if any(shape != shapes[0] for shape in shapes[1:]):
        raise ValueError(
            "All stacks must share identical frame shape (Y, X, C). "
            f"Got: {shapes}"
        )


def main() -> None:
    args = parse_args()

    embryo_flythrough_path = args.datadir / args.embryo_flythrough_file
    embryo_rotation_path = args.datadir / args.embryo_rotation_file
    midgut_rotation_path = args.datadir / args.midgut_rotation_file
    out_path = args.datadir / args.output_file

    print(f"Loading flythrough stack: {embryo_flythrough_path}")
    embryo_flythrough = load_rgb_stack(embryo_flythrough_path)
    print(f"Loading embryo-rotation stack: {embryo_rotation_path}")
    embryo_rotation = load_rgb_stack(embryo_rotation_path)
    print(f"Loading midgut-rotation stack: {midgut_rotation_path}")
    midgut_rotation = load_rgb_stack(midgut_rotation_path)

    ensure_same_frame_shape(embryo_flythrough, embryo_rotation, midgut_rotation)

    if args.flythrough_step > 1:
        embryo_flythrough = embryo_flythrough[:: args.flythrough_step]
    embryo_flythrough_scaled = scale_stack(embryo_flythrough, args.scale_embryo_flythrough)
    embryo_rotation_scaled = scale_stack(embryo_rotation, args.scale_embryo_rotation)
    midgut_rotation_scaled = scale_stack(midgut_rotation, args.scale_midgut_rotation)

    black_to_embryo_rot = fade_black_to_frame(
        embryo_rotation_scaled[0],
        n_frames=args.fade_black_to_embryo_rot_frames,
    )
    embryo_rot_to_midgut = crossfade_frames(
        embryo_rotation_scaled[-1],
        midgut_rotation_scaled[0],
        n_frames=args.fade_embryo_rot_to_midgut_frames,
    )

    movie_stack = np.concatenate(
        [
            embryo_flythrough_scaled,
            black_to_embryo_rot,
            embryo_rotation_scaled,
            embryo_rot_to_midgut,
            midgut_rotation_scaled,
        ],
        axis=0,
    )

    print(
        "Composed stack shape: "
        f"{movie_stack.shape} (flythrough={len(embryo_flythrough_scaled)}, "
        f"black_to_embryo_rot={len(black_to_embryo_rot)}, "
        f"embryo_rotation={len(embryo_rotation_scaled)}, "
        f"embryo_to_midgut={len(embryo_rot_to_midgut)}, "
        f"midgut_rotation={len(midgut_rotation_scaled)})"
    )
    print(f"Writing: {out_path}")
    tifffile.imwrite(out_path, movie_stack, photometric="rgb")
    print("Done.")


if __name__ == "__main__":
    main()