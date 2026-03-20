"""Prepare TIFF volumes as H5 for ilastik training.

Loads a TIFF stack, optionally downsamples it, and writes an H5 file with a
single dataset (default name: inputData) that ilastik can open.

Usage:
    python mask_prepare_ilastik.py --tiff-path /path/to/Time_000005.tif
    python mask_prepare_ilastik.py --tiff-path /path/to/Time_000005.tif --downsample-axes xy

NPMitchell 2026
"""

import argparse
from pathlib import Path

from mask_utils import downsample_tiff_to_h5


# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

DEFAULT_DATASET_NAME = "inputData"
DEFAULT_DOWNSAMPLE_FACTOR = 4.0
DEFAULT_DOWNSAMPLE_AXES = "xy"
DEFAULT_COMPRESSION = "none"
DEFAULT_COMPRESSION_OPTS = None
DEFAULT_TIFF_PATH: Path | None = Path("/project/npmitchell/avistrok/2024-06-18_141952_bynGAL4klar_UASmChCAAXUASLamGFP/channel_0_tp_00030_angle_20.ome.tif")
DEFAULT_H5_PATH: Path | None = None


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Convert TIFF to ilastik-ready H5 using shared downsampling utility."
    )
    parser.add_argument(
        "--tiff-path",
        type=Path,
        default=DEFAULT_TIFF_PATH,
        help="Input TIFF path (ZYX or ZYXC).",
    )
    parser.add_argument(
        "--h5-path",
        type=Path,
        default=DEFAULT_H5_PATH,
        help="Output H5 path. Defaults to <tiff-stem>.h5 next to TIFF.",
    )
    parser.add_argument(
        "--dataset-name",
        default=DEFAULT_DATASET_NAME,
        help="Dataset name inside the output H5 (default: inputData).",
    )
    parser.add_argument(
        "--downsample-factor",
        type=float,
        default=DEFAULT_DOWNSAMPLE_FACTOR,
        help="Spatial downsampling factor (default: 4.0).",
    )
    parser.add_argument(
        "--downsample-axes",
        choices=["xyz", "xy"],
        default=DEFAULT_DOWNSAMPLE_AXES,
        help="Axes to downsample: xyz or xy (default: xy).",
    )
    parser.add_argument(
        "--compression",
        default=DEFAULT_COMPRESSION,
        help="HDF5 compression filter (default: none). Use none, gzip, lzf, etc.",
    )
    parser.add_argument(
        "--compression-opts",
        type=int,
        default=DEFAULT_COMPRESSION_OPTS,
        help="Compression options (for example gzip level). Ignored when compression=none.",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()

    if args.tiff_path is None:
        raise ValueError(
            "No TIFF path provided. Set DEFAULT_TIFF_PATH at the top of this script "
            "or pass --tiff-path on the command line."
        )

    tiff_path = args.tiff_path.expanduser().resolve()
    if args.h5_path is None:
        h5_path = tiff_path.with_suffix(".h5")
    else:
        h5_path = args.h5_path.expanduser().resolve()

    out_path = downsample_tiff_to_h5(
        tiff_path=tiff_path,
        h5_path=h5_path,
        dataset_name=args.dataset_name,
        downsample_factor=args.downsample_factor,
        downsample_axes=args.downsample_axes,
        compression=args.compression,
        compression_opts=args.compression_opts,
    )

    print(f"Wrote ilastik training H5: {out_path}")
    print(f"dataset_name: {args.dataset_name}")
    print(f"downsample_factor: {args.downsample_factor}")
    print(f"downsample_axes: {args.downsample_axes}")


if __name__ == "__main__":
    main()
