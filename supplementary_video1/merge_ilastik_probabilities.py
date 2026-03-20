import argparse
from pathlib import Path

import h5py
import numpy as np

"""Script to merge ilastik probability channels into the original inputData H5 files.
Assumes original inputData has 3 channels and ilastik probabilities have 4 channels.
Writes merged data as uint16 with appropriate metadata for downstream use.

NPMitchell 2026
"""

DEFAULT_DATADIR = Path(
    "/Users/npmitchell/Dropbox/Soft_Matter_Biophysics/PAPER/gut_chirality/supplementary_videos/"
    "20240531_142410_22C_HandGFPbynGAL4klar_UASmChCAAXH2AviRFP"
)


def to_uint16_base(data: np.ndarray) -> np.ndarray:
    if data.dtype == np.uint16:
        return data
    if np.issubdtype(data.dtype, np.floating):
        min_val = float(np.min(data))
        max_val = float(np.max(data))
        if 0.0 <= min_val and max_val <= 1.0:
            scaled = np.rint(data * 65535.0)
            return np.clip(scaled, 0, 65535).astype(np.uint16)
        return np.clip(np.rint(data), 0, 65535).astype(np.uint16)
    if np.issubdtype(data.dtype, np.integer):
        return np.clip(data, 0, 65535).astype(np.uint16)
    raise TypeError(f"Unsupported base dtype: {data.dtype}")


def to_uint16_probabilities(data: np.ndarray) -> np.ndarray:
    if data.dtype == np.uint16:
        return data
    if np.issubdtype(data.dtype, np.floating):
        min_val = float(np.min(data))
        max_val = float(np.max(data))
        if 0.0 <= min_val and max_val <= 1.0:
            scaled = np.rint(data * 65535.0)
            return np.clip(scaled, 0, 65535).astype(np.uint16)
        return np.clip(np.rint(data), 0, 65535).astype(np.uint16)
    if np.issubdtype(data.dtype, np.integer):
        return np.clip(data, 0, 65535).astype(np.uint16)
    raise TypeError(f"Unsupported probability dtype: {data.dtype}")


def choose_zyxc_chunks(shape: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    z, y, x, c = shape
    return (
        max(1, min(z, 16)),
        max(1, min(y, 128)),
        max(1, min(x, 128)),
        c,
    )


def merge_channels(datadir: Path, timepoint: str, write_mode: str) -> Path:
    base_path = datadir / f"{timepoint}.h5"
    prob_path = datadir / f"{timepoint}-inputData_Probabilities.h5"
    out_path = base_path if write_mode == "inplace" else datadir / f"{timepoint}_merged.h5"

    if not base_path.exists():
        raise FileNotFoundError(f"Missing base file: {base_path}")
    if not prob_path.exists():
        raise FileNotFoundError(f"Missing probability file: {prob_path}")

    with h5py.File(base_path, "r") as f_base:
        if "inputData" not in f_base:
            raise KeyError(f"Dataset 'inputData' not found in {base_path}")
        base = f_base["inputData"][()]

    with h5py.File(prob_path, "r") as f_prob:
        if "exported_data" not in f_prob:
            raise KeyError(f"Dataset 'exported_data' not found in {prob_path}")
        probs = f_prob["exported_data"][()]

    if base.shape[:3] != probs.shape[:3]:
        raise ValueError(f"Spatial dims differ: {base.shape[:3]} vs {probs.shape[:3]}")

    probs_u16 = to_uint16_probabilities(probs)

    if base.shape[-1] == 3:
        base_u16 = to_uint16_base(base)
        merged = np.concatenate([base_u16, probs_u16], axis=-1)
        source_note = "constructed from original inputData (3 ch) + exported_data (4 ch)"
    elif base.shape[-1] == 7:
        base_u16 = to_uint16_base(base)
        source_note = "copied from already-merged inputData and converted to uint16"
        merged = np.concatenate([base_u16[..., :3], probs_u16], axis=-1)
    else:
        raise ValueError(f"Unexpected inputData channel count: {base.shape[-1]}")

    mode = "r+" if write_mode == "inplace" else "w"
    with h5py.File(out_path, mode) as f_out:
        if write_mode == "inplace" and "inputData" in f_out:
            del f_out["inputData"]
        chunks = choose_zyxc_chunks(merged.shape)
        dset = f_out.create_dataset(
            "inputData",
            data=merged,
            chunks=chunks,
            compression="gzip",
            compression_opts=1,
            shuffle=True,
        )
        dset.attrs["channel_layout"] = "ch0,ch1,ch2 source image; ch3-ch6 ilastik probabilities"
        dset.attrs["source_note"] = source_note
        dset.attrs["DIMENSION_LABELS"] = np.array(["z", "y", "x", "c"], dtype="S1")
        dset.attrs["display_mode"] = "grayscale"

    return out_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Merge Time_XXXXX.h5 inputData channels with ilastik probabilities."
    )
    parser.add_argument(
        "--datadir",
        type=Path,
        default=DEFAULT_DATADIR,
        help=(
            "Full directory path containing the H5 files "
            f"(default: {DEFAULT_DATADIR})."
        ),
    )
    parser.add_argument(
        "--timepoint",
        default="Time_000005",
        help="Timepoint stem (default: Time_000005).",
    )
    parser.add_argument(
        "--write-mode",
        choices=["copy", "inplace"],
        default="copy",
        help="copy writes <timepoint>_merged.h5, inplace overwrites <timepoint>.h5",
    )
    args = parser.parse_args()

    output_path = merge_channels(args.datadir, args.timepoint, args.write_mode)

    with h5py.File(output_path, "r") as f_check:
        ds = f_check["inputData"]
        print(f"Wrote: {output_path}")
        print(f"Shape: {ds.shape} dtype: {ds.dtype}")
        print(f"Chunks: {ds.chunks} compression: {ds.compression}")
        print(f"channel_layout: {ds.attrs.get('channel_layout', 'N/A')}")
        print(f"source_note: {ds.attrs.get('source_note', 'N/A')}")


if __name__ == "__main__":
    main()
