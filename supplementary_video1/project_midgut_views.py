"""Generate angle-swept Z-projections from an anatomical midgut RGB volume.

Pipeline per view angle:
1. Rotate 3D RGB volume around the X axis.
2. Build depth-dependent fade weights in rotated Z coordinates.
3. Multiply rotated volume by fade weights.
4. Maximum-intensity project along Z.

The default input is expected to be an RGB TIFF volume with shape (Z, Y, X, 3).
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import tifffile
from scipy import ndimage as ndi


DEFAULT_DATADIR = Path(
    "/Users/npmitchell/Dropbox/Soft_Matter_Biophysics/PAPER/gut_chirality/supplementary_videos/"
    "20240531_142410_22C_HandGFPbynGAL4klar_UASmChCAAXH2AviRFP"
)
DEFAULT_INPUT_FILE = "Time_000005_WT_masked_midgutOnly_rgb_anatomical_labeled.tif"
DEFAULT_BATCH_INPUT_FILES = [
    "Time_000005_WT_masked_midgutOnly_rgb_anatomical_labeled.tif",
    "Time_000005_WT_masked_embryo_rgb_anatomical_labeled.tif",
]
DEFAULT_OUTPUT_FILE = "Time_000005_WT_masked_midgutOnly_rgb_anatomical_labeled_xrot_mip_stack.tif"
DEFAULT_ANGLES_FILE = "Time_000005_WT_masked_midgutOnly_rgb_anatomical_labeled_xrot_angles_deg.txt"
DEFAULT_FRAME_DIR = "Time_000005_WT_masked_midgutOnly_rgb_anatomical_labeled_xrot_tiff_frames"
DEFAULT_FRAME_PREFIX = "frame"

DEFAULT_ANGLE_START_DEG = 0.0
DEFAULT_ANGLE_STOP_DEG = 360.0
DEFAULT_ANGLE_STEP_DEG = 2.0

DEFAULT_INTERP_ORDER = 1
DEFAULT_ROTATE_RESHAPE = False

DEFAULT_APPLY_FADE = True
DEFAULT_Z_START_FADE = 220.0
DEFAULT_Z_END_FADE = 280.0
DEFAULT_WRITE_FRAME_TIFFS = True


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Rotate an RGB volume around X across multiple angles, apply optional depth fade "
            "per rotated view, and output Z-MIP image stack."
        )
    )
    parser.add_argument("--datadir", type=Path, default=DEFAULT_DATADIR)
    parser.add_argument(
        "--input-file",
        action="append",
        default=None,
        help=(
            "Input RGB TIFF filename. Provide multiple --input-file flags to process a batch "
            "in one run. If omitted, uses built-in batch defaults (midgut + embryo)."
        ),
    )
    parser.add_argument("--output-file", default=DEFAULT_OUTPUT_FILE)
    parser.add_argument("--angles-file", default=DEFAULT_ANGLES_FILE)
    parser.add_argument(
        "--frame-dir",
        "--png-dir",
        dest="frame_dir",
        default=DEFAULT_FRAME_DIR,
        help="Directory for per-angle TIFF frame exports (legacy alias: --png-dir).",
    )
    parser.add_argument(
        "--frame-prefix",
        "--png-prefix",
        dest="frame_prefix",
        default=DEFAULT_FRAME_PREFIX,
        help="Filename prefix for per-angle TIFF frame exports (legacy alias: --png-prefix).",
    )

    parser.add_argument("--angle-start", type=float, default=DEFAULT_ANGLE_START_DEG)
    parser.add_argument("--angle-stop", type=float, default=DEFAULT_ANGLE_STOP_DEG)
    parser.add_argument("--angle-step", type=float, default=DEFAULT_ANGLE_STEP_DEG)

    parser.add_argument(
        "--interp-order",
        type=int,
        default=DEFAULT_INTERP_ORDER,
        choices=[0, 1, 2, 3, 4, 5],
        help="Interpolation order used during rotate. 0 nearest, 1 linear, 3 cubic.",
    )
    parser.add_argument(
        "--rotate-reshape",
        default=DEFAULT_ROTATE_RESHAPE,
        action=argparse.BooleanOptionalAction,
        help="If true, enlarge output canvas during rotate. Default: False.",
    )

    parser.add_argument(
        "--apply-fade",
        default=DEFAULT_APPLY_FADE,
        action=argparse.BooleanOptionalAction,
        help="Apply depth fade in rotated Z coordinates before MIP. Default: False.",
    )
    parser.add_argument(
        "--z-start-fade",
        type=float,
        default=DEFAULT_Z_START_FADE,
        help="Fade starts at this rotated Z index.",
    )
    parser.add_argument(
        "--z-end-fade",
        type=float,
        default=DEFAULT_Z_END_FADE,
        help="Fade reaches 0 at and beyond this rotated Z index.",
    )
    parser.add_argument(
        "--write-frame-tiffs",
        "--write-png-frames",
        dest="write_frame_tiffs",
        default=DEFAULT_WRITE_FRAME_TIFFS,
        action=argparse.BooleanOptionalAction,
        help=(
            "Write each projected frame as uint16 TIFF while processing "
            "(legacy alias: --write-png-frames). Default: True."
        ),
    )
    return parser.parse_args()


def load_rgb_volume(path: Path) -> np.ndarray:
    """Load TIFF and return RGB volume in (Z, Y, X, 3) order.

    Supports common RGB axis layouts:
    - (Z, Y, X, 3)  channel-last
    - (Z, 3, Y, X)  channel-first per-slice
    - (3, Z, Y, X)  channel-first volume
    """
    arr = tifffile.imread(path)

    if arr.ndim != 4:
        raise ValueError(
            f"Expected 4D RGB volume, got {arr.shape} from {path}"
        )

    if arr.shape[-1] == 3:
        return arr
    if arr.shape[1] == 3:
        return np.transpose(arr, (0, 2, 3, 1))
    if arr.shape[0] == 3:
        return np.transpose(arr, (1, 2, 3, 0))

    raise ValueError(
        "Could not identify RGB channel axis of length 3; "
        f"expected one of (Z,Y,X,3), (Z,3,Y,X), (3,Z,Y,X), got {arr.shape} from {path}"
    )


def build_angles(start: float, stop: float, step: float) -> np.ndarray:
    if step == 0:
        raise ValueError("--angle-step must be non-zero")
    direction = np.sign(stop - start) if stop != start else np.sign(step)
    if direction == 0:
        direction = 1.0
    if np.sign(step) != direction:
        step = abs(step) * direction

    # Include stop within floating tolerance.
    eps = abs(step) * 1e-6 + 1e-12
    return np.arange(start, stop + eps, step, dtype=np.float64)


def _sanitize_for_name(value: float) -> str:
    """Return filesystem-safe compact text for a numeric value."""
    text = f"{value:g}"
    return text.replace("-", "m").replace(".", "p")


def build_fade_tag(apply_fade: bool, z_start_fade: float, z_end_fade: float) -> str:
    if apply_fade:
        return f"fade_s{_sanitize_for_name(z_start_fade)}_e{_sanitize_for_name(z_end_fade)}"
    return "fade_off"


def append_tag_to_filename(path: Path, tag: str) -> Path:
    """Append a tag before file extension."""
    return path.with_name(f"{path.stem}_{tag}{path.suffix}")


def derive_names_from_input(input_file: str) -> tuple[str, str, str]:
    """Derive output filenames from an input filename stem."""
    stem = Path(input_file).stem
    return (
        f"{stem}_xrot_mip_stack.tif",
        f"{stem}_xrot_angles_deg.txt",
        f"{stem}_xrot_tiff_frames",
    )


def depth_fade_weights(nz: int, z_start_fade: float, z_end_fade: float) -> np.ndarray:
    """Piecewise fade in rotated Z coordinates.

    f(z) = 1 for z < z_start
         = 1 - (z - z_start)/(z_end - z_start) for z_start <= z < z_end
         = 0 for z >= z_end

    Edge case: if z_end <= z_start, use hard cutoff at z_start.
    """
    z = np.arange(nz, dtype=np.float32)

    if z_end_fade <= z_start_fade:
        weights = (z < z_start_fade).astype(np.float32)
        return weights

    weights = np.ones_like(z, dtype=np.float32)
    ramp_mask = (z >= z_start_fade) & (z < z_end_fade)
    weights[ramp_mask] = 1.0 - (z[ramp_mask] - z_start_fade) / (z_end_fade - z_start_fade)
    weights[z >= z_end_fade] = 0.0
    return np.clip(weights, 0.0, 1.0)


def rotate_about_x(volume_zyx3: np.ndarray, angle_deg: float, order: int, reshape: bool) -> np.ndarray:
    """Rotate around X by rotating the Z-Y plane (axes 0 and 1) for all channels."""
    # Fast path: avoid interpolation work when requested rotation is effectively identity.
    angle_mod = np.mod(angle_deg, 360.0)
    if np.isclose(angle_mod, 0.0, atol=1e-9) or np.isclose(angle_mod, 360.0, atol=1e-9):
        return volume_zyx3

    # scipy uses CCW-positive in the chosen axes plane.
    rotated = ndi.rotate(
        volume_zyx3,
        angle=angle_deg,
        axes=(1, 0),
        reshape=reshape,
        order=order,
        mode="constant",
        cval=0,
        prefilter=order > 1,
    )

    if np.issubdtype(volume_zyx3.dtype, np.integer):
        info = np.iinfo(volume_zyx3.dtype)
        rotated = np.clip(np.rint(rotated), info.min, info.max).astype(volume_zyx3.dtype)
    return rotated


def write_tiff_frame(path: Path, mip_yx3: np.ndarray) -> None:
    """Write one projected RGB frame as uint16 TIFF, preserving bit depth."""
    frame = np.clip(np.rint(mip_yx3), 0, 65535).astype(np.uint16)
    try:
        tifffile.imwrite(path, frame, photometric="rgb")
    except Exception as exc:
        raise RuntimeError(f"Failed to write TIFF frame: {path}") from exc


def load_tiff_frame(path: Path) -> np.ndarray:
    """Load one projected RGB TIFF frame, preserving stored bit depth."""
    frame = tifffile.imread(path)
    if frame.ndim != 3 or frame.shape[-1] != 3:
        raise ValueError(f"Expected RGB TIFF frame with shape (Y, X, 3), got {frame.shape} from {path}")
    return frame


def project_views(
    volume_zyx3: np.ndarray,
    angles_deg: np.ndarray,
    interp_order: int,
    rotate_reshape: bool,
    apply_fade: bool,
    z_start_fade: float,
    z_end_fade: float,
    write_frame_tiffs: bool,
    frame_dir: Path,
    frame_prefix: str,
) -> np.ndarray:
    views = []
    work = volume_zyx3.astype(np.float32, copy=False)

    if write_frame_tiffs:
        frame_dir.mkdir(parents=True, exist_ok=True)

    for idx, angle in enumerate(angles_deg):
        print(f"[{idx + 1}/{len(angles_deg)}] angle={angle:.3f} deg", flush=True)
        frame_path = frame_dir / f"{frame_prefix}_{idx:04d}_angle_{angle:+08.3f}.tif"

        if write_frame_tiffs and frame_path.exists():
            print(f"... reusing existing TIFF frame: {frame_path.name}", flush=True)
            mip = load_tiff_frame(frame_path)
        else:
            rot = rotate_about_x(work, angle_deg=float(angle), order=interp_order, reshape=rotate_reshape)

            print("... now applying fade" if apply_fade else "skipping fade", flush=True)
            if apply_fade:
                nz = rot.shape[0]
                w = depth_fade_weights(nz, z_start_fade=z_start_fade, z_end_fade=z_end_fade)
                rot = rot * w[:, None, None, None]

            mip = np.max(rot, axis=0)

            if write_frame_tiffs:
                write_tiff_frame(frame_path, mip)

        if mip.shape != volume_zyx3.shape[1:]:
            raise ValueError(
                f"Frame shape mismatch for angle {angle:.3f}: expected {volume_zyx3.shape[1:]}, got {mip.shape}"
            )

        if mip.dtype != np.uint16:
            mip = np.clip(np.rint(mip), 0, 65535).astype(np.uint16)

        views.append(mip)

    out = np.stack(views, axis=0)
    out = np.clip(np.rint(out), 0, 65535).astype(np.uint16)
    return out


def run_one_input(
    args: argparse.Namespace,
    input_file: str,
    output_file: str,
    angles_file: str,
    frame_dir_name: str,
) -> None:
    in_path = args.datadir / input_file
    fade_tag = build_fade_tag(args.apply_fade, args.z_start_fade, args.z_end_fade)
    out_path = append_tag_to_filename(args.datadir / output_file, fade_tag)
    ang_path = args.datadir / angles_file
    frame_base = args.datadir / frame_dir_name
    frame_dir = frame_base.parent / f"{frame_base.name}_{fade_tag}"

    print(f"Input: {in_path}")
    vol = load_rgb_volume(in_path)
    print(f"Loaded shape: {vol.shape}, dtype: {vol.dtype}")

    angles = build_angles(args.angle_start, args.angle_stop, args.angle_step)
    print(f"Angles: {len(angles)} views from {angles[0]:.3f} to {angles[-1]:.3f} deg")
    print(f"Fade enabled: {args.apply_fade}")
    if args.apply_fade:
        print(f"Fade z_start={args.z_start_fade}, z_end={args.z_end_fade}")
    print(f"Write TIFF frames: {args.write_frame_tiffs}")
    if args.write_frame_tiffs:
        print(f"TIFF frame directory: {frame_dir}")
    print(f"Output TIFF: {out_path}")

    stack = project_views(
        vol,
        angles_deg=angles,
        interp_order=args.interp_order,
        rotate_reshape=args.rotate_reshape,
        apply_fade=args.apply_fade,
        z_start_fade=args.z_start_fade,
        z_end_fade=args.z_end_fade,
        write_frame_tiffs=args.write_frame_tiffs,
        frame_dir=frame_dir,
        frame_prefix=args.frame_prefix,
    )

    tifffile.imwrite(out_path, stack, photometric="rgb")
    np.savetxt(ang_path, angles, fmt="%.6f")

    print(f"Wrote MIP stack: {out_path} with shape {stack.shape}")
    print(f"Wrote angle list: {ang_path}")


def main() -> None:
    args = parse_args()

    input_files = args.input_file if args.input_file else list(DEFAULT_BATCH_INPUT_FILES)

    multiple_inputs = len(input_files) > 1
    user_overrode_names = (
        args.output_file != DEFAULT_OUTPUT_FILE
        or args.angles_file != DEFAULT_ANGLES_FILE
        or args.frame_dir != DEFAULT_FRAME_DIR
    )
    if multiple_inputs and user_overrode_names:
        raise ValueError(
            "When passing multiple --input-file values, do not override --output-file/--angles-file/--frame-dir. "
            "Per-input names are auto-derived from each input filename stem."
        )

    for idx, input_file in enumerate(input_files):
        if multiple_inputs:
            output_file, angles_file, frame_dir_name = derive_names_from_input(input_file)
        else:
            output_file, angles_file, frame_dir_name = args.output_file, args.angles_file, args.frame_dir

        if multiple_inputs:
            print("\n" + "=" * 72)
            print(f"Batch item {idx + 1}/{len(input_files)}: {input_file}")
            print("=" * 72)

        run_one_input(
            args,
            input_file=input_file,
            output_file=output_file,
            angles_file=angles_file,
            frame_dir_name=frame_dir_name,
        )


if __name__ == "__main__":
    main()
