"""Rotate masked RGB embryo volumes into an anatomical coordinate system.

The anatomical frame is estimated from the embryo RGB mask:
1. Extract the largest connected component from nonzero RGB voxels.
2. Estimate the embryo long axis from centroids binned along the principal axis.
3. Align that long axis with the longest axis of the data volume.
4. Rotate around that long axis to align the secondary cross-sectional axis with
   the second-longest axis of the data volume.
5. Apply the same affine rotation to both embryo and midgut RGB volumes.

Manual rotation mode reproduces a Fiji macro pipeline step-by-step:

  1. Rotate each XY slice by ``angle1`` degrees         (Fiji: ``run("Rotate...", "angle=<angle1> ...")``)
  2. Reslice from Top : (Z, Y, X, C) → (Y, Z, X, C)   (Fiji: ``run("Reslice [/]...", "start=Top avoid")``)
  3. Rotate each XY slice by ``angle2`` degrees
  4. Reslice from Left: (Z, Y, X, C) → (X, Z, Y, C)   (Fiji: ``run("Reslice [/]...", "start=Left avoid")``)
  5. Rotate each XY slice by ``angle3`` degrees  (folds in the preceding 90° CW step:
     ``angle3 = 90 + <final Fiji angle>`` in CW degrees)
    6. Final reslice from Left (equivalent to Fiji: Image > Stacks > Reslice > Left)

Set ``DEFAULT_FIJI_ANGLES = [angle1, angle2, angle3]`` (degrees) at the top of the file,
or pass ``--fiji-angles A1 A2 A3`` on the command line.

Interpolation order is controlled by ``--interp-order`` (default 0 = nearest-neighbour,
matching Fiji ``interpolation=None``; use 3 for cubic in auto mode).
Canvas enlargement during in-plane rotations can be disabled with ``--no-enlarge``
(default False).

NPMitchell 2026
"""

import argparse
import json
from pathlib import Path
from typing import Optional

import numpy as np
import tifffile
from scipy import ndimage as ndi


DEFAULT_DATADIR = Path(
    "/Users/npmitchell/Dropbox/Soft_Matter_Biophysics/PAPER/gut_chirality/supplementary_videos/"
    "20240531_142410_22C_HandGFPbynGAL4klar_UASmChCAAXH2AviRFP"
)
DEFAULT_TIME_STEM = "Time_000005"
DEFAULT_EMBRYO_FILE = DEFAULT_TIME_STEM + "_WT_masked_embryo_rgb.tif"
DEFAULT_MIDGUT_FILE = DEFAULT_TIME_STEM + "_WT_masked_midgutOnly_rgb.tif"
DEFAULT_EMBRYO_OUT = DEFAULT_TIME_STEM + "_WT_masked_embryo_rgb_anatomical.tif"
DEFAULT_MIDGUT_OUT = DEFAULT_TIME_STEM + "_WT_masked_midgutOnly_rgb_anatomical.tif"
DEFAULT_EMBRYO_CROP_LOG = DEFAULT_TIME_STEM + "_WT_masked_embryo_rgb_anatomical_crop_history.json"
DEFAULT_MIDGUT_CROP_LOG = DEFAULT_TIME_STEM + "_WT_masked_midgutOnly_rgb_anatomical_crop_history.json"
DEFAULT_REPORT_FILE = DEFAULT_TIME_STEM + "_anatomical_rotation_report.txt"
DEFAULT_INTERPOLATION_ORDER = 3 
DEFAULT_ENLARGE = True
DEFAULT_DO_EMBRYO = True
DEFAULT_DO_MIDGUT = True

# Three in-plane rotation angles for the Fiji pipeline steps 1, 3, and 5 (see module docstring).
# angles[2] incorporates the 90° clockwise step: if Fiji used rot90-right then rotate(49),
# the combined angle is 90 + 49 = 139 degrees.
DEFAULT_FIJI_ANGLES = [13.0, -8.0, 140.0]


# ---------------------------------------------------------------------------
# Fiji-equivalent sequential pipeline
# ---------------------------------------------------------------------------

def reslice_top(stack: np.ndarray) -> np.ndarray:
    """Reslice a ZYX[C] stack from the Top.

    Equivalent to Fiji: ``run("Reslice [/]...", "start=Top avoid")``.
    Transposes (Z, Y, X[, C]) → (Y, Z, X[, C]).
    Each output slice is an XZ cross-section of the original volume.
    """
    if stack.ndim == 3:
        return np.transpose(stack, (1, 0, 2))
    elif stack.ndim == 4:
        return np.transpose(stack, (1, 0, 2, 3))
    else:
        raise ValueError(f"reslice_top: unsupported shape {stack.shape}")


def reslice_left(stack: np.ndarray) -> np.ndarray:
    """Reslice a ZYX[C] stack from the Left.

    Equivalent to Fiji: ``run("Reslice [/]...", "start=Left avoid")``.
    Transposes (Z, Y, X[, C]) → (X, Z, Y[, C]).
    Each output slice is a ZY cross-section of the original volume.
    """
    if stack.ndim == 3:
        return np.transpose(stack, (2, 0, 1))
    elif stack.ndim == 4:
        return np.transpose(stack, (2, 0, 1, 3))
    else:
        raise ValueError(f"reslice_left: unsupported shape {stack.shape}")


def rotate_stack_in_plane(
    stack: np.ndarray, angle_deg: float, order: int = 0, enlarge: bool = DEFAULT_ENLARGE
) -> np.ndarray:
    """Rotate each XY slice of a ZYX[C] stack in-plane.

    Equivalent to Fiji: ``run("Rotate...", "angle=<angle_deg> grid=1 ... fill"``.

    Parameters
    ----------
    stack:
        Array with shape (Z, Y, X) or (Z, Y, X, C).
    angle_deg:
        Clockwise rotation angle in degrees (positive = CW), matching Fiji's
        sign convention for ``run("Rotate...", "angle=...")``.  Internally the
        sign is flipped before passing to ``scipy.ndimage.rotate``, which uses
        CCW-positive.
    order:
        Spline interpolation order.  0 = nearest-neighbour (matches Fiji
        ``interpolation=None``).
    enlarge:
        If True, expand the output canvas so no corners are clipped
        (matches Fiji ``enlarge`` option).
    """
    original_dtype = stack.dtype
    rotated = ndi.rotate(
        stack,
        angle=-angle_deg,   # Fiji is CW-positive; scipy ndi.rotate is CCW-positive
        axes=(1, 2),
        reshape=enlarge,
        order=order,
        mode="constant",
        cval=0,
        prefilter=order > 1,
    )
    if np.issubdtype(original_dtype, np.integer):
        info = np.iinfo(original_dtype)
        rotated = np.clip(np.rint(rotated), info.min, info.max).astype(original_dtype)
    return rotated


def rotate_90_right(stack: np.ndarray) -> np.ndarray:
    """Rotate each XY slice 90 degrees clockwise (Fiji: Rotate 90 Degrees Right)."""
    if stack.ndim in (3, 4):
        return np.rot90(stack, k=-1, axes=(1, 2))
    raise ValueError(f"rotate_90_right: unsupported shape {stack.shape}")


def crop_to_nonzero(stack: np.ndarray) -> tuple[np.ndarray, dict]:
    """Crop a ZYX[C] stack to the tight bounding box of nonzero voxels."""
    if stack.ndim not in (3, 4):
        raise ValueError(f"crop_to_nonzero: unsupported shape {stack.shape}")

    spatial_shape = stack.shape[:3]
    if stack.ndim == 4:
        mask = np.any(stack != 0, axis=-1)
    else:
        mask = stack != 0

    if not np.any(mask):
        return stack, {
            "cropped": False,
            "z_start": 0,
            "z_stop": int(spatial_shape[0]),
            "y_start": 0,
            "y_stop": int(spatial_shape[1]),
            "x_start": 0,
            "x_stop": int(spatial_shape[2]),
            "trim_z0": 0,
            "trim_z1": 0,
            "trim_y0": 0,
            "trim_y1": 0,
            "trim_x0": 0,
            "trim_x1": 0,
        }

    coords = np.argwhere(mask)
    mins = coords.min(axis=0)
    maxs = coords.max(axis=0) + 1

    z0, y0, x0 = int(mins[0]), int(mins[1]), int(mins[2])
    z1, y1, x1 = int(maxs[0]), int(maxs[1]), int(maxs[2])

    if stack.ndim == 4:
        cropped = stack[z0:z1, y0:y1, x0:x1, :]
    else:
        cropped = stack[z0:z1, y0:y1, x0:x1]

    return cropped, {
        "cropped": True,
        "z_start": z0,
        "z_stop": z1,
        "y_start": y0,
        "y_stop": y1,
        "x_start": x0,
        "x_stop": x1,
        "trim_z0": z0,
        "trim_z1": int(spatial_shape[0] - z1),
        "trim_y0": y0,
        "trim_y1": int(spatial_shape[1] - y1),
        "trim_x0": x0,
        "trim_x1": int(spatial_shape[2] - x1),
    }


def crop_to_bounds(stack: np.ndarray, crop_info: dict) -> tuple[np.ndarray, dict]:
    """Crop a ZYX[C] stack using explicit bounds from a previous crop log."""
    if stack.ndim not in (3, 4):
        raise ValueError(f"crop_to_bounds: unsupported shape {stack.shape}")

    z0 = int(crop_info["z_start"])
    z1 = int(crop_info["z_stop"])
    y0 = int(crop_info["y_start"])
    y1 = int(crop_info["y_stop"])
    x0 = int(crop_info["x_start"])
    x1 = int(crop_info["x_stop"])

    sz, sy, sx = stack.shape[:3]
    if not (0 <= z0 <= z1 <= sz and 0 <= y0 <= y1 <= sy and 0 <= x0 <= x1 <= sx):
        raise ValueError(
            f"crop_to_bounds: bounds out of range for shape {stack.shape[:3]}: "
            f"z[{z0}:{z1}] y[{y0}:{y1}] x[{x0}:{x1}]"
        )

    if stack.ndim == 4:
        cropped = stack[z0:z1, y0:y1, x0:x1, :]
    else:
        cropped = stack[z0:z1, y0:y1, x0:x1]

    return cropped, {
        "cropped": True,
        "z_start": z0,
        "z_stop": z1,
        "y_start": y0,
        "y_stop": y1,
        "x_start": x0,
        "x_stop": x1,
        "trim_z0": z0,
        "trim_z1": int(sz - z1),
        "trim_y0": y0,
        "trim_y1": int(sy - y1),
        "trim_x0": x0,
        "trim_x1": int(sx - x1),
    }


def apply_fiji_pipeline(
    volume_zyx3: np.ndarray,
    angles: np.ndarray,
    order: int = 0,
    enlarge: bool = True,
    fixed_crop_history: Optional[list[dict]] = None,
) -> tuple[np.ndarray, list[dict]]:
    """Apply the Fiji-equivalent sequential rotation pipeline to a ZYX3 RGB volume.

    Reproduces the following Fiji macro sequence::

        run("Rotate...",      "angle=<angles[0]> grid=1 interpolation=None fill enlarge");
        run("Reslice [/]...", "output=1.000 start=Top avoid");
        run("Rotate...",      "angle=<angles[1]> grid=1 interpolation=None fill enlarge");
        run("Reslice [/]...", "output=1.000 start=Left avoid");
        run("Rotate...",      "angle=<angles[2]> grid=1 interpolation=None fill enlarge");
        run("Reslice [/]...", "output=1.000 start=Left avoid");
        run("Rotate 90 Degrees Right");

    Note: the original Fiji macro had ``run("Rotate 90 Degrees Right")`` followed by
    ``run("Rotate...", "angle=49 ...")``.  In Fiji's CW-positive convention those two
    steps combine to ``90 + 49 = 139`` degrees CW, so the separate rot90 step has been
    folded into ``angles[2]`` (default ``139``).

    Parameters
    ----------
    volume_zyx3:
        Input array with shape (Z, Y, X, 3).
    angles:
        Three in-plane rotation angles ``[angle1, angle2, angle3]`` in degrees.
    order:
        Spline interpolation order (0 = nearest-neighbour, matching Fiji
        ``interpolation=None``).
    enlarge:
        Expand canvas during in-plane rotations (matching Fiji ``enlarge``).
    """
    vol = volume_zyx3
    crop_history: list[dict] = []

    def apply_step(step_index: int, step_label: str, fn) -> None:
        nonlocal vol
        print(f"  [{step_index}/7] {step_label} ...", flush=True)
        shape_before = list(vol.shape)
        vol_after = fn(vol)
        shape_after_op = list(vol_after.shape)

        if fixed_crop_history is None:
            vol_cropped, crop_info = crop_to_nonzero(vol_after)
            crop_source = "nonzero"
        else:
            if len(fixed_crop_history) < step_index:
                raise ValueError(
                    f"fixed_crop_history has {len(fixed_crop_history)} steps, "
                    f"but step {step_index} was requested"
                )
            ref_crop = fixed_crop_history[step_index - 1]
            vol_cropped, crop_info = crop_to_bounds(vol_after, ref_crop)
            crop_source = "embryo_reference"

        shape_after_crop = list(vol_cropped.shape)
        vol = vol_cropped

        crop_history.append(
            {
                "step": step_index,
                "operation": step_label,
                "shape_before": shape_before,
                "shape_after_operation": shape_after_op,
                "shape_after_crop": shape_after_crop,
                "crop_source": crop_source,
                **crop_info,
            }
        )
        print(
            "      crop "
            f"z[{crop_info['z_start']}:{crop_info['z_stop']}] "
            f"y[{crop_info['y_start']}:{crop_info['y_stop']}] "
            f"x[{crop_info['x_start']}:{crop_info['x_stop']}] -> {shape_after_crop[:3]}",
            flush=True,
        )

    apply_step(
        1,
        f"rotate in-plane by {angles[0]:.4g} deg",
        lambda arr: rotate_stack_in_plane(arr, float(angles[0]), order=order, enlarge=enlarge),
    )
    apply_step(2, "reslice from Top", reslice_top)
    apply_step(
        3,
        f"rotate in-plane by {angles[1]:.4g} deg",
        lambda arr: rotate_stack_in_plane(arr, float(angles[1]), order=order, enlarge=enlarge),
    )
    apply_step(4, "reslice from Left", reslice_left)
    apply_step(
        5,
        f"rotate in-plane by {angles[2]:.4g} deg (includes 90-deg-right)",
        lambda arr: rotate_stack_in_plane(arr, float(angles[2]), order=order, enlarge=enlarge),
    )
    apply_step(6, "final reslice from Left", reslice_left)
    apply_step(7, "final Rotate 90 Degrees Right", rotate_90_right)
    return vol, crop_history


def write_crop_history(
    out_path: Path,
    source_file: Path,
    fiji_angles: np.ndarray,
    crop_history: list[dict],
    crop_reference_file: Optional[Path] = None,
) -> None:
    """Write per-step crop metadata so meshes can be transformed consistently."""
    payload = {
        "source_file": str(source_file),
        "fiji_angles_cw_deg": [float(val) for val in fiji_angles],
        "crop_reference_file": None if crop_reference_file is None else str(crop_reference_file),
        "steps": crop_history,
    }
    out_path.write_text(json.dumps(payload, indent=2) + "\n")


def read_crop_history(in_path: Path) -> tuple[list[dict], Optional[np.ndarray]]:
    """Read crop-history JSON and return per-step crop entries and optional angles."""
    payload = json.loads(in_path.read_text())
    steps = payload.get("steps")
    if not isinstance(steps, list) or len(steps) == 0:
        raise ValueError(f"Invalid crop history file (missing non-empty 'steps'): {in_path}")

    angles = payload.get("fiji_angles_cw_deg", None)
    angle_array = None
    if angles is not None:
        angle_array = np.asarray(angles, dtype=np.float64)
    return steps, angle_array


def write_report(report_path: Path, measurement: dict, embryo_file: Path, midgut_file: Path) -> None:
    """Write a text report describing the measured axes and applied rotation."""
    lines = [
        "rotation_mode: Fiji pipeline",
        f"embryo_file: {embryo_file}",
        f"midgut_file: {midgut_file}",
        f"fiji_angles: {np.array2string(measurement['fiji_angles'], precision=6)}",
    ]
    report_path.write_text("\n".join(lines) + "\n")


def load_rgb_volume(path: Path) -> np.ndarray:
    """Load an RGB TIFF volume and validate shape."""
    volume = tifffile.imread(path)
    if volume.ndim != 4 or volume.shape[-1] != 3:
        raise ValueError(f"Expected RGB TIFF with shape (z, y, x, 3), got {volume.shape} from {path}")
    return volume


def build_parser() -> argparse.ArgumentParser:
    """Build CLI parser.

    Supply three in-plane rotation angles via ``--fiji-angles A1 A2 A3``.
    Each is applied at one of the three ``Rotate...`` steps of the Fiji pipeline
    (see module docstring for the full step sequence).
    """
    parser = argparse.ArgumentParser(
        description="Rotate embryo and midgut RGB TIFF volumes into an anatomical frame."
    )
    parser.add_argument("--datadir", type=Path, default=DEFAULT_DATADIR)
    parser.add_argument("--embryo-file", default=DEFAULT_EMBRYO_FILE)
    parser.add_argument("--midgut-file", default=DEFAULT_MIDGUT_FILE)
    parser.add_argument("--embryo-out", default=DEFAULT_EMBRYO_OUT)
    parser.add_argument("--midgut-out", default=DEFAULT_MIDGUT_OUT)
    parser.add_argument("--embryo-crop-log", default=DEFAULT_EMBRYO_CROP_LOG)
    parser.add_argument("--midgut-crop-log", default=DEFAULT_MIDGUT_CROP_LOG)
    parser.add_argument("--report-file", default=DEFAULT_REPORT_FILE)
    parser.add_argument(
        "--fiji-angles",
        type=float,
        nargs=3,
        default=DEFAULT_FIJI_ANGLES,
        metavar=("A1", "A2", "A3"),
        help=(
            "Three in-plane rotation angles [A1 A2 A3] in degrees. "
            "A1: initial XY rotation; A2: after reslice-top; A3: after reslice-left "
            "(folds in the 90-degree-right step). "
            "See module docstring for the full step sequence."
        ),
    )
    parser.add_argument(
        "--mask-threshold",
        type=int,
        default=0,
        help="RGB voxels with any channel above this value are treated as embryo mask.",
    )
    parser.add_argument(
        "--interp-order",
        type=int,
        default=DEFAULT_INTERPOLATION_ORDER,
        choices=[0, 1, 2, 3, 4, 5],
        help=(
            "Spline interpolation order for in-plane rotations. "
            "0=nearest-neighbour, matching Fiji 'interpolation=None'; use 3 for cubic."
        ),
    )
    parser.add_argument(
        "--enlarge",
        default=DEFAULT_ENLARGE,
        action=argparse.BooleanOptionalAction,
        help="Expand canvas during in-plane rotations (Fiji 'enlarge'). Default: False.",
    )
    parser.add_argument(
        "--do-embryo",
        default=DEFAULT_DO_EMBRYO,
        action=argparse.BooleanOptionalAction,
        help="Whether to process/write embryo output. Default: True.",
    )
    parser.add_argument(
        "--do-midgut",
        default=DEFAULT_DO_MIDGUT,
        action=argparse.BooleanOptionalAction,
        help="Whether to process/write midgut output. Default: True.",
    )
    parser.add_argument(
        "--report-only",
        action="store_true",
        help="Measure and report anatomical rotation without writing rotated TIFF volumes.",
    )
    return parser


def main() -> None:
    """Estimate anatomical rotation from embryo RGB and apply it to both RGB volumes."""
    args = build_parser().parse_args()

    embryo_path = args.datadir / args.embryo_file
    midgut_path = args.datadir / args.midgut_file
    embryo_out_path = args.datadir / args.embryo_out
    midgut_out_path = args.datadir / args.midgut_out
    embryo_crop_log_path = args.datadir / args.embryo_crop_log
    midgut_crop_log_path = args.datadir / args.midgut_crop_log
    report_path = args.datadir / args.report_file

    fiji_angles = np.asarray(args.fiji_angles, dtype=np.float64)
    measurement = {"fiji_angles": fiji_angles}

    write_report(report_path, measurement, embryo_path, midgut_path)

    print(f"Rotation source volume: {embryo_path}")
    print(f"fiji_angles: {fiji_angles}")
    print(f"Report written: {report_path}")
    print(f"do_embryo: {args.do_embryo}")
    print(f"do_midgut: {args.do_midgut}")

    if args.report_only:
        return

    if not args.do_embryo and not args.do_midgut:
        print("Nothing to do: both --no-do-embryo and --no-do-midgut are set.")
        return

    embryo_crop_history: Optional[list[dict]] = None

    if args.do_embryo:
        embryo_rgb = load_rgb_volume(embryo_path)
        print(f"Rotating embryo RGB: {embryo_path.name} ...", flush=True)
        embryo_rot, embryo_crop_history = apply_fiji_pipeline(
            embryo_rgb, fiji_angles, order=args.interp_order, enlarge=args.enlarge,
        )
        write_crop_history(
            embryo_crop_log_path,
            embryo_path,
            fiji_angles,
            embryo_crop_history,
            crop_reference_file=None,
        )
        print(f"Wrote embryo crop history: {embryo_crop_log_path}")
        print(f"Writing embryo output ...", flush=True)
        tifffile.imwrite(embryo_out_path, embryo_rot, photometric="rgb")
        print(f"Wrote rotated embryo RGB: {embryo_out_path}")
    elif args.do_midgut:
        print(f"Loading embryo crop history from: {embryo_crop_log_path}")
        embryo_crop_history, log_angles = read_crop_history(embryo_crop_log_path)
        if log_angles is not None and not np.allclose(log_angles, fiji_angles):
            print(
                "Warning: --fiji-angles differs from angles stored in embryo crop log. "
                f"CLI={fiji_angles}, log={log_angles}",
                flush=True,
            )

    if args.do_midgut:
        if embryo_crop_history is None:
            raise ValueError(
                "Midgut processing requires embryo crop history. "
                "Run with --do-embryo first (default), or provide an existing embryo crop log via --embryo-crop-log."
            )

        midgut_rgb = load_rgb_volume(midgut_path)
        print(f"Rotating midgut RGB: {midgut_path.name} ...", flush=True)
        midgut_rot, midgut_crop_history = apply_fiji_pipeline(
            midgut_rgb,
            fiji_angles,
            order=args.interp_order,
            enlarge=args.enlarge,
            fixed_crop_history=embryo_crop_history,
        )
        write_crop_history(
            midgut_crop_log_path,
            midgut_path,
            fiji_angles,
            midgut_crop_history,
            crop_reference_file=embryo_path,
        )
        print(f"Wrote midgut crop history: {midgut_crop_log_path}")
        print(f"Writing midgut output ...", flush=True)
        tifffile.imwrite(midgut_out_path, midgut_rot, photometric="rgb")
        print(f"Wrote rotated midgut RGB: {midgut_out_path}")


if __name__ == "__main__":
    main()