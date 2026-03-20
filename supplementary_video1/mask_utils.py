"""Shared utilities for ilastik-probability-based masking pipelines.

All pure helper functions and the common run_pipeline() entry point used by
mask_embryo.py and mask_midgut.py.

NPMitchell 2026
"""

import itertools
from pathlib import Path
from typing import Optional

import h5py
import numpy as np
import tifffile
from scipy import ndimage as ndi


# ---------------------------------------------------------------------------
# Low-level helpers
# ---------------------------------------------------------------------------

def to_prob_0_1(channel: np.ndarray) -> np.ndarray:
    """Convert a probability channel to float32 in [0, 1]."""
    arr = channel.astype(np.float32, copy=False)
    min_val = float(arr.min())
    max_val = float(arr.max())

    if 0.0 <= min_val and max_val <= 1.0:
        return arr

    if 0.0 <= min_val and max_val <= 65535.0:
        return arr / 65535.0

    raise ValueError(
        f"Probability channel values are outside expected range: min={min_val}, max={max_val}"
    )


def largest_connected_component(binary_mask: np.ndarray) -> np.ndarray:
    """Return a binary mask containing only the single largest connected component."""
    labels, num = ndi.label(binary_mask)
    if num == 0:
        raise ValueError("No connected components found after thresholding. Lower --thres.")

    sizes = np.bincount(labels.ravel())
    sizes[0] = 0
    largest_label = int(np.argmax(sizes))
    return labels == largest_label


def choose_chunks(shape: tuple[int, ...]) -> tuple[int, ...]:
    """Choose moderate HDF5 chunk sizes for 3D/4D volumes."""
    if len(shape) == 4:
        z, y, x, c = shape
        return (max(1, min(z, 16)), max(1, min(y, 128)), max(1, min(x, 128)), c)
    if len(shape) == 3:
        z, y, x = shape
        return (max(1, min(z, 16)), max(1, min(y, 128)), max(1, min(x, 128)))
    return None


def load_data_volume(datadir: Path, data_file: str, data_dataset: str, time_stem: str):
    """Load target data volume from TIFF or H5, with optional auto filename selection."""
    if data_file == "auto":
        tif_path = datadir / f"{time_stem}.tif"
        h5_path = datadir / f"{time_stem}.h5"
        if tif_path.exists():
            data_path = tif_path
        elif h5_path.exists():
            data_path = h5_path
        else:
            raise FileNotFoundError(
                f"Neither TIFF nor H5 found for {time_stem} in {datadir}. "
                f"Looked for {tif_path.name} and {h5_path.name}."
            )
    else:
        data_path = datadir / data_file

    suffix = data_path.suffix.lower()
    if suffix in {".tif", ".tiff"}:
        data = tifffile.imread(data_path)
    elif suffix in {".h5", ".hdf5"}:
        with h5py.File(data_path, "r") as f_data:
            if data_dataset not in f_data:
                raise KeyError(f"Dataset '{data_dataset}' not found in {data_path}")
            data = f_data[data_dataset][()]
    else:
        raise ValueError(f"Unsupported --data-file extension: {data_path.suffix}")

    return data, data_path


def downsample_tiff_to_h5(
    tiff_path: Path,
    h5_path: Path,
    dataset_name: str = "inputData",
    downsample_factor: float = 4.0,
    downsample_axes: str = "xyz",
    compression: Optional[str] = None,
    compression_opts: Optional[int] = None,
) -> Path:
    """Downsample TIFF and write to H5.

    Expected TIFF axis order is ZYX or ZYXC.

    Parameters
    ----------
    tiff_path:
        Input TIFF path.
    h5_path:
        Output H5 path.
    dataset_name:
        Name of output dataset in H5 (default: inputData).
    downsample_factor:
        Downsampling factor applied to selected spatial axes (default: 4.0).
    downsample_axes:
        Spatial axes to downsample. Allowed values are "xyz" (default) and "xy".
    compression:
        HDF5 compression filter. Use None (default) or "none" for no compression.
    compression_opts:
        Compression options for the selected filter (for example gzip level).
    """
    if downsample_factor <= 0:
        raise ValueError("downsample_factor must be > 0")

    axes = downsample_axes.lower()
    if axes not in {"xyz", "xy"}:
        raise ValueError("downsample_axes must be 'xyz' or 'xy'")

    data = tifffile.imread(tiff_path)
    if data.ndim not in (3, 4):
        raise ValueError(f"Expected TIFF with shape ZYX or ZYXC, got {data.shape}")

    zoom_spatial = 1.0 / float(downsample_factor)
    z_zoom = zoom_spatial if axes == "xyz" else 1.0
    if data.ndim == 3:
        zoom = (z_zoom, zoom_spatial, zoom_spatial)
    else:
        zoom = (z_zoom, zoom_spatial, zoom_spatial, 1.0)

    # Use linear interpolation for stable intensity behavior across integer and float inputs.
    ds = ndi.zoom(
        data.astype(np.float32, copy=False),
        zoom=zoom,
        order=1,
        mode="nearest",
        prefilter=False,
    )

    if np.issubdtype(data.dtype, np.integer):
        info = np.iinfo(data.dtype)
        ds = np.clip(np.rint(ds), info.min, info.max).astype(data.dtype)
    else:
        ds = ds.astype(data.dtype, copy=False)

    compression_value: Optional[str]
    if compression is None:
        compression_value = None
    else:
        comp_norm = compression.lower().strip()
        compression_value = None if comp_norm in {"", "none"} else comp_norm

    dset_kwargs = {
        "data": ds,
        "dtype": ds.dtype,
        "chunks": choose_chunks(ds.shape),
    }
    if compression_value is not None:
        dset_kwargs["compression"] = compression_value
        if compression_opts is not None:
            dset_kwargs["compression_opts"] = compression_opts

    with h5py.File(h5_path, "w") as f_out:
        f_out.create_dataset(dataset_name, **dset_kwargs)

    return h5_path


# ---------------------------------------------------------------------------
# Axis reordering
# ---------------------------------------------------------------------------

def reorder_axes_explicit(data: np.ndarray, axis_order: str) -> tuple[np.ndarray, tuple[int, ...]]:
    """Reorder array axes to canonical ZYXC/ZYX using explicit axis-order labels."""
    labels = axis_order.replace(",", "").replace(" ", "").upper()
    if len(labels) != data.ndim:
        raise ValueError(
            f"--axis-order '{axis_order}' has {len(labels)} labels, but data ndim is {data.ndim}."
        )
    if len(set(labels)) != len(labels):
        raise ValueError(f"--axis-order '{axis_order}' has repeated labels.")
    if any(ch not in "ZYXC" for ch in labels):
        raise ValueError(f"--axis-order '{axis_order}' can only use Z,Y,X,C labels.")

    desired = [ax for ax in "ZYXC" if ax in labels]
    perm = [labels.index(ax) for ax in desired]
    return np.transpose(data, axes=perm), tuple(perm)


def best_spatial_permutation(
    spatial_shape: tuple[int, int, int], expected_shape: tuple[int, int, int]
) -> tuple[tuple[int, int, int], tuple[int, int, int], tuple[int, int]]:
    """Find spatial-axis permutation minimizing mismatch to expected ZYX shape."""
    best_perm = (0, 1, 2)
    best_diffs = tuple(int(spatial_shape[i] - expected_shape[i]) for i in best_perm)
    best_score = (max(abs(d) for d in best_diffs), sum(abs(d) for d in best_diffs))

    for perm in itertools.permutations((0, 1, 2)):
        diffs = tuple(int(spatial_shape[perm[i]] - expected_shape[i]) for i in range(3))
        score = (max(abs(d) for d in diffs), sum(abs(d) for d in diffs))
        if score < best_score:
            best_perm = perm
            best_diffs = diffs
            best_score = score

    return best_perm, best_diffs, best_score


def reorder_axes_autodetect(
    data: np.ndarray, expected_spatial_shape: tuple[int, int, int]
) -> tuple[np.ndarray, str, tuple[int, ...]]:
    """Auto-detect channel axis and best spatial permutation to match expected shape."""
    if data.ndim == 3:
        perm, diffs, score = best_spatial_permutation(data.shape, expected_spatial_shape)
        reordered = np.transpose(data, axes=perm)
        return (
            reordered,
            f"autoDetect: 3D spatial perm {perm}, diffs(z,y,x)={diffs}, score={score}",
            tuple(perm),
        )

    if data.ndim == 4:
        likely_channel_axes = [i for i, size in enumerate(data.shape) if size <= 8]
        channel_axes = likely_channel_axes if likely_channel_axes else [3]

        best = None
        for ch_axis in channel_axes:
            moved = np.moveaxis(data, ch_axis, -1)
            move_perm = tuple([i for i in range(data.ndim) if i != ch_axis] + [ch_axis])
            perm, diffs, score = best_spatial_permutation(moved.shape[:3], expected_spatial_shape)
            spatial_perm = (perm[0], perm[1], perm[2], 3)
            reordered = np.transpose(moved, axes=spatial_perm)
            total_perm = tuple(move_perm[i] for i in spatial_perm)
            candidate = (score, ch_axis, perm, diffs, reordered, total_perm)
            if best is None or candidate[0] < best[0]:
                best = candidate

        score, ch_axis, perm, diffs, reordered, total_perm = best
        return (
            reordered,
            "autoDetect: "
            f"channel axis {ch_axis}->3, spatial perm {perm}, "
            f"diffs(z,y,x)={diffs}, score={score}",
            total_perm,
        )

    identity = tuple(range(data.ndim))
    return data, "autoDetect: unchanged (ndim not 3/4)", identity


def reorder_data_axes(
    data: np.ndarray, axis_order: str, expected_spatial_shape: tuple[int, int, int]
) -> tuple[np.ndarray, str, tuple[int, ...]]:
    """Reorder axes either by explicit user order or by auto-detection."""
    if axis_order == "autoDetect":
        return reorder_axes_autodetect(data, expected_spatial_shape)

    reordered, perm = reorder_axes_explicit(data, axis_order)
    return reordered, f"explicit axis order: {axis_order}", perm


# ---------------------------------------------------------------------------
# Permutation / axis utilities
# ---------------------------------------------------------------------------

def invert_permutation(perm: tuple[int, ...]) -> tuple[int, ...]:
    """Return inverse permutation such that arr == arr_perm.transpose(inverse_perm)."""
    return tuple(np.argsort(np.asarray(perm)).tolist())


def infer_native_axes_from_applied_perm(applied_perm: tuple[int, ...]) -> str:
    """Infer native axis labels from permutation that mapped native data to canonical ZYX(C)."""
    canonical_axes = "ZYXC"[: len(applied_perm)]
    inverse_perm = invert_permutation(applied_perm)
    native_axes = [canonical_axes[inverse_perm[i]] for i in range(len(applied_perm))]
    return "".join(native_axes)


def is_imagej_axes_valid(axes: str) -> bool:
    """Return True if axis string is a valid ordered subset of TZCYXS for ImageJ."""
    order = "TZCYXS"
    if any(ax not in order for ax in axes):
        return False
    idx = [order.index(ax) for ax in axes]
    return idx == sorted(idx)


def make_imagej_compatible(
    data: np.ndarray, axes: str
) -> tuple[np.ndarray, str, Optional[tuple[int, ...]]]:
    """Reorder data axes to an ImageJ-valid TZCYXS subset order when needed."""
    if is_imagej_axes_valid(axes):
        return data, axes, None

    target_axes = "".join(ax for ax in "TZCYXS" if ax in axes)
    perm = tuple(axes.index(ax) for ax in target_axes)
    reordered = np.transpose(data, axes=perm)
    return reordered, target_axes, perm


# ---------------------------------------------------------------------------
# Mask alignment and data conversion
# ---------------------------------------------------------------------------

def align_mask_to_target_shape(
    mask: np.ndarray, target_shape: tuple[int, int, int], ssfactor: float
) -> tuple[np.ndarray, tuple[int, int, int]]:
    """Validate mask-target shape agreement and center crop/pad within tolerance."""
    tolerance = int(np.ceil(ssfactor))
    diffs = tuple(int(mask.shape[i] - target_shape[i]) for i in range(3))

    if any(abs(d) > tolerance for d in diffs):
        raise ValueError(
            "Upsampled mask differs too much from target shape. "
            f"mask={mask.shape}, target={target_shape}, diffs(z,y,x)={diffs}, "
            f"allowed per-axis difference <= {tolerance} (from ssfactor={ssfactor})."
        )

    aligned = mask
    for axis in range(3):
        delta = aligned.shape[axis] - target_shape[axis]
        if delta > 0:
            start = delta // 2
            end = start + target_shape[axis]
            slicer = [slice(None)] * 3
            slicer[axis] = slice(start, end)
            aligned = aligned[tuple(slicer)]
        elif delta < 0:
            pad_before = (-delta) // 2
            pad_after = (-delta) - pad_before
            pad_width = [(0, 0), (0, 0), (0, 0)]
            pad_width[axis] = (pad_before, pad_after)
            aligned = np.pad(aligned, pad_width=pad_width, mode="constant", constant_values=0)

    return aligned, diffs


def to_uint16_volume(volume: np.ndarray) -> np.ndarray:
    """Convert a volume to uint16 using round-and-clip in [0, 65535]."""
    return np.clip(np.rint(volume), 0, 65535).astype(np.uint16)


# ---------------------------------------------------------------------------
# RGB colorization
# ---------------------------------------------------------------------------

def colorize_channels_to_rgb_u16(
    masked_data_zyxc_u16: np.ndarray,
    channel_colors: tuple[tuple[float, float, float], ...],
    rgb_scaling: list[float],
    rgb_bg_values: list[float],
) -> np.ndarray:
    """Colorize channels into an RGB ZYX3 uint16 volume using additive mixing.

    Per-channel background is subtracted before brightness scaling. Background
    values are specified in the original uint16 intensity units.
    """
    if masked_data_zyxc_u16.ndim == 3:
        bg_value = rgb_bg_values[0] if rgb_bg_values else 0.0
        scale = rgb_scaling[0] if rgb_scaling else 1.0
        gray = np.clip(masked_data_zyxc_u16.astype(np.float32) / 65535.0 - bg_value / 65535.0, 0.0, 1.0)
        gray = np.clip(gray * scale, 0.0, 1.0)[..., None]
        return np.repeat(gray, 3, axis=-1)

    if masked_data_zyxc_u16.ndim != 4:
        raise ValueError(
            f"Expected ZYXC or ZYX input for RGB colorization, got {masked_data_zyxc_u16.shape}"
        )

    num_channels = masked_data_zyxc_u16.shape[-1]
    colors = np.asarray(channel_colors, dtype=np.float32)
    if colors.ndim != 2 or colors.shape[1] != 3:
        raise ValueError("channel_colors must be a sequence of (R,G,B) triplets")
    if colors.shape[0] < num_channels:
        raise ValueError(
            f"Need at least {num_channels} channel colors, got {colors.shape[0]}"
        )

    data_norm = masked_data_zyxc_u16.astype(np.float32) / 65535.0
    rgb = np.zeros(masked_data_zyxc_u16.shape[:3] + (3,), dtype=np.float32)
    for channel_index in range(num_channels):
        bg_norm = rgb_bg_values[channel_index] / 65535.0 if channel_index < len(rgb_bg_values) else 0.0
        channel_signal = np.clip(data_norm[..., channel_index] - bg_norm, 0.0, 1.0)
        rgb += (
            channel_signal[..., None]
            * rgb_scaling[channel_index]
            * colors[channel_index][None, None, None, :]
        )

    rgb = np.clip(rgb, 0.0, 1.0)
    return np.rint(rgb * 65535.0).astype(np.uint16)


# ---------------------------------------------------------------------------
# Common pipeline
# ---------------------------------------------------------------------------

def run_pipeline(
    args,
    channel_colors: tuple,
    rgb_out_path: Path,
    write_imagej: bool = True,
    write_mask_lowres: bool = True,
    write_rgb: bool = True,
    apply_mode: str = "multiply",
    outside_clip_thres: Optional[float] = None,
    mask_prob_thres: Optional[float] = None,
) -> None:
    """Run the full mask-build → upsample → apply → write pipeline.

    Parameters
    ----------
    args:
        Parsed argparse namespace from either mask_embryo or mask_midgut.
    channel_colors:
        Sequence of (R,G,B) float triplets, one per data channel.
    rgb_out_path:
        Full output path for the RGB composite TIFF.
    write_imagej:
        Save masked data as an ImageJ-compatible TIFF.
    write_mask_lowres:
        Save the low-resolution binary mask.
    write_rgb:
        Colorize channels and save an RGB composite TIFF.
    apply_mode:
        How to apply the mask to data.
        - "multiply": masked_data = data * mask (default)
        - "clip_outside": for mask==0 voxels, clip data to outside_clip_thres
    outside_clip_thres:
        Intensity clipping threshold used only when apply_mode="clip_outside".
    mask_prob_thres:
        Optional probability threshold override for initial binary mask.
        If None, uses args.thres for backward compatibility.
    """
    prob_path = args.datadir / args.prob_file
    out_path = args.datadir / args.out_file

    # --- Load probability volume ---
    with h5py.File(prob_path, "r") as f_prob:
        if args.prob_dataset not in f_prob:
            raise KeyError(f"Dataset '{args.prob_dataset}' not found in {prob_path}")
        prob_data = f_prob[args.prob_dataset][()]

    if prob_data.ndim != 4:
        raise ValueError(f"Expected probability data shape (z,y,x,c), got {prob_data.shape}")
    if not (0 <= args.mask_channel < prob_data.shape[-1]):
        raise ValueError(
            f"mask_channel={args.mask_channel} out of range for {prob_data.shape[-1]} channels"
        )

    prob_ch = to_prob_0_1(prob_data[..., args.mask_channel])

    # --- Build hard binary mask ---
    prob_thres = args.thres if mask_prob_thres is None else mask_prob_thres
    binary = prob_ch >= prob_thres
    largest = largest_connected_component(binary)
    largest = ndi.binary_fill_holes(largest)
    if args.dilate_iters > 0:
        largest = ndi.binary_dilation(largest, iterations=args.dilate_iters)

    mask_lowres = largest.astype(np.float32)

    if args.ssfactor <= 0:
        raise ValueError("--ssfactor must be > 0")

    # --- Upsample mask ---
    if args.ssfactor == 1.0:
        mask = mask_lowres
    else:
        mask = ndi.zoom(
            mask_lowres,
            zoom=(args.ssfactor, args.ssfactor, args.ssfactor),
            order=1,
            mode="nearest",
            prefilter=False,
        ).astype(np.float32)

    # --- Load and reorder data ---
    data, data_path = load_data_volume(
        args.datadir, args.data_file, args.data_dataset, args.time_stem
    )
    data_reordered, axis_order_note, applied_perm = reorder_data_axes(
        data, args.axis_order, mask.shape
    )

    # --- Apply mask ---
    if data_reordered.ndim == 4:
        target_shape = data_reordered.shape[:3]
        mask_aligned, shape_diffs = align_mask_to_target_shape(mask, target_shape, args.ssfactor)
        data_float = data_reordered.astype(np.float32)
        if apply_mode == "multiply":
            masked_data = data_float * mask_aligned[..., None]
        elif apply_mode == "clip_outside":
            if outside_clip_thres is None:
                raise ValueError("outside_clip_thres must be provided when apply_mode='clip_outside'")
            outside = mask_aligned <= 0
            masked_data = data_float.copy()
            masked_data[outside, :] = np.minimum(masked_data[outside, :], outside_clip_thres)
        else:
            raise ValueError(f"Unknown apply_mode: {apply_mode}")
    elif data_reordered.ndim == 3:
        target_shape = data_reordered.shape
        mask_aligned, shape_diffs = align_mask_to_target_shape(mask, target_shape, args.ssfactor)
        data_float = data_reordered.astype(np.float32)
        if apply_mode == "multiply":
            masked_data = data_float * mask_aligned
        elif apply_mode == "clip_outside":
            if outside_clip_thres is None:
                raise ValueError("outside_clip_thres must be provided when apply_mode='clip_outside'")
            outside = mask_aligned <= 0
            masked_data = data_float.copy()
            masked_data[outside] = np.minimum(masked_data[outside], outside_clip_thres)
        else:
            raise ValueError(f"Unknown apply_mode: {apply_mode}")
    else:
        raise ValueError(f"Expected data ndim 3 or 4, got shape {data_reordered.shape}")

    # --- Convert and reorder back to native axes ---
    masked_data_u16 = to_uint16_volume(masked_data)
    inverse_perm = invert_permutation(applied_perm)
    masked_data_u16_native = np.transpose(masked_data_u16, axes=inverse_perm)
    native_axes = infer_native_axes_from_applied_perm(applied_perm)
    imagej_data, imagej_axes, imagej_perm = make_imagej_compatible(masked_data_u16_native, native_axes)

    # --- Write masked data ---
    if write_imagej:
        tifffile.imwrite(
            out_path,
            imagej_data,
            imagej=True,
            metadata={"axes": imagej_axes},
            photometric="minisblack",
        )
    else:
        tifffile.imwrite(out_path, masked_data_u16_native, photometric="minisblack")

    # --- Write low-res mask ---
    lowres_path = out_path.with_name(f"{out_path.stem}_mask_lowres.tif")
    if write_mask_lowres:
        if write_imagej:
            tifffile.imwrite(
                lowres_path,
                mask_lowres,
                imagej=True,
                metadata={"axes": "ZYX"},
                photometric="minisblack",
            )
        else:
            tifffile.imwrite(lowres_path, mask_lowres, photometric="minisblack")

    # --- Write RGB composite ---
    if write_rgb:
        print(f"Colorizing masked data to RGB using channel colors: {channel_colors}")
        rgb_scaling = args.rgb_scaling
        rgb_bg_values = args.rgb_bg_values
        if masked_data_u16.ndim == 4:
            nchan = masked_data_u16.shape[-1]
            ch_sel = list(args.rgb_channels) if args.rgb_channels is not None else list(range(nchan))
            for ci in ch_sel:
                if ci < 0 or ci >= nchan:
                    raise ValueError(f"rgb channel index {ci} out of range 0..{nchan-1}")
            rgb_input = masked_data_u16[..., ch_sel]
            colors_sel = tuple(channel_colors[i] for i in ch_sel)
            scaling_sel = [rgb_scaling[i] if i < len(rgb_scaling) else 1.0 for i in ch_sel]
            bg_sel = [rgb_bg_values[i] if i < len(rgb_bg_values) else 0.0 for i in ch_sel]
        else:
            rgb_input = masked_data_u16
            colors_sel = channel_colors[:1]
            scaling_sel = rgb_scaling[:1] if rgb_scaling else [1.0]
            bg_sel = rgb_bg_values[:1] if rgb_bg_values else [0.0]

        rgb_data_u16 = colorize_channels_to_rgb_u16(rgb_input, colors_sel, scaling_sel, bg_sel)
        tifffile.imwrite(rgb_out_path, rgb_data_u16, photometric="rgb")

    # --- Diagnostics ---
    print(f"Wrote: {out_path}")
    print(f"data_source: {data_path}")
    print(f"data_shape_raw: {data.shape} -> data_shape_reordered: {data_reordered.shape}")
    print(f"applied_perm: {applied_perm} inverse_perm: {inverse_perm}")
    print(f"native_axes: {native_axes}")
    print(f"imagej_axes_used: {imagej_axes}")
    if imagej_perm is not None:
        print(f"imagej_reorder_perm: {imagej_perm}")
    print(f"write_imagej: {write_imagej}")
    print(f"output_shape_native: {masked_data_u16_native.shape}")
    print(f"output_shape_written: {imagej_data.shape}")
    print(f"axis_order: {args.axis_order} ({axis_order_note})")
    print(f"mask_prob_thres: {prob_thres}")
    if args.axis_order == "autoDetect":
        print(f"best_spatial_permutation: {axis_order_note}")
    print(f"shape_diffs(z,y,x): {shape_diffs} (checked against ssfactor={args.ssfactor})")
    print(f"apply_mode: {apply_mode}")
    if apply_mode == "clip_outside":
        print(f"outside_clip_thres: {outside_clip_thres}")
    print(f"mask_lowres: {mask_lowres.shape} dtype={mask_lowres.dtype}")
    print(f"mask_upsampled: {mask_aligned.shape} dtype={mask_aligned.dtype}")
    print(f"masked_data: {masked_data_u16_native.shape} dtype={masked_data_u16_native.dtype}")
    print(f"write_mask_lowres: {write_mask_lowres}")
    if write_mask_lowres:
        print(f"mask_lowres_file: {lowres_path}")
    print(f"write_rgb: {write_rgb}")
    if write_rgb:
        print(f"rgb_file: {rgb_out_path}")
        print("rgb_write_mode: standard TIFF (non-ImageJ; ImageJ does not support uint16 RGB)")
        print(f"channel_colors: {channel_colors}")
        print(f"rgb_scaling: {rgb_scaling}")
        print(f"rgb_bg_values: {rgb_bg_values}")
