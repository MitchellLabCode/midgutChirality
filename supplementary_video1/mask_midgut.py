"""STEP 2: midgut-only masking pipeline.

Builds a hard binary mask (largest connected component + hole fill) from an
ilastik probability channel (midgut vs background), upsamples it, and applies
it to the original multichannel data.

Usage:
    python mask_midgut.py [options]

NPMitchell 2026
"""

import argparse
from pathlib import Path

import numpy as np

from mask_utils import run_pipeline

# ---------------------------------------------------------------------------
# Defaults — edit these for your dataset
# ---------------------------------------------------------------------------

DEFAULT_DATADIR = Path(
    "/Users/npmitchell/Dropbox/Soft_Matter_Biophysics/PAPER/gut_chirality/supplementary_videos/"
    "20240531_142410_22C_HandGFPbynGAL4klar_UASmChCAAXH2AviRFP"
)

DEFAULT_TIME_STEM = "Time_000005"
DEFAULT_PROB_FILE = DEFAULT_TIME_STEM + "_merged-inputData_Probabilities_midgutOnly.h5"
DEFAULT_DATA_FILE = DEFAULT_TIME_STEM + "_2x.tif"
DEFAULT_OUT_FILE = DEFAULT_TIME_STEM + "_WT_masked_midgutOnly.tif"
DEFAULT_RGB_OUT_FILE = DEFAULT_TIME_STEM + "_WT_masked_midgutOnly_rgb.tif"
DEFAULT_SSFACTOR = 2.5
DEFAULT_AXIS_ORDER = "autoDetect"
DEFAULT_MASK_CHANNEL = 0

# RGB composition — midgut step uses only channels 0 and 2
DEFAULT_RGB_CHANNELS = [0, 2]
DEFAULT_RGB_SCALING = [1.5, 1, 1]
DEFAULT_RGB_BG_VALUES = [280, 1390, 1200]

endoderm = np.array([100, 200, 255]) / 255.0
lvm = np.array([1, 0, 0], dtype=float)
cvm = np.array([255, 200, 100]) / 255.0
DEFAULT_CHANNEL_COLORS = (cvm, lvm, endoderm)

# Write flags
WRITE_IMAGEJ = True
WRITE_MASK_LOWRES = True
WRITE_RGB = True


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="STEP 2 – Build midgut-only mask and apply to data."
    )
    parser.add_argument("--datadir", type=Path, default=DEFAULT_DATADIR)
    parser.add_argument("--prob-file", default=DEFAULT_PROB_FILE)
    parser.add_argument(
        "--data-file",
        default=DEFAULT_DATA_FILE,
        help="Data filename (.tif/.h5). Use 'auto' to detect from --time-stem.",
    )
    parser.add_argument("--out-file", default=DEFAULT_OUT_FILE)
    parser.add_argument("--prob-dataset", default="exported_data")
    parser.add_argument("--data-dataset", default="inputData")
    parser.add_argument("--time-stem", default=DEFAULT_TIME_STEM)
    parser.add_argument("--mask-channel", type=int, default=DEFAULT_MASK_CHANNEL)
    parser.add_argument("--thres", type=float, default=0.3,
                        help="Probability threshold for initial binary mask.")
    parser.add_argument("--dilate-iters", type=int, default=2,
                        help="Binary dilation iterations after hole fill.")
    parser.add_argument("--ssfactor", type=float, default=DEFAULT_SSFACTOR,
                        help="Upsample factor per spatial dimension (z,y,x).")
    parser.add_argument("--axis-order", default=DEFAULT_AXIS_ORDER,
                        help="Axis labels (e.g. ZYXC) or 'autoDetect'.")
    parser.add_argument(
        "--rgb-channels", type=int, nargs="+", default=DEFAULT_RGB_CHANNELS,
        help="0-based channel indices for RGB composite.",
    )
    parser.add_argument(
        "--rgb-scaling", type=float, nargs="+", default=DEFAULT_RGB_SCALING,
        help="Per-channel brightness scale factors for RGB output only.",
    )
    parser.add_argument(
        "--rgb-bg-values", type=float, nargs="+", default=DEFAULT_RGB_BG_VALUES,
        help="Per-channel background values to subtract before RGB scaling, in uint16 intensity units.",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    rgb_out_path = args.datadir / DEFAULT_RGB_OUT_FILE
    run_pipeline(
        args,
        channel_colors=DEFAULT_CHANNEL_COLORS,
        rgb_out_path=rgb_out_path,
        write_imagej=WRITE_IMAGEJ,
        write_mask_lowres=WRITE_MASK_LOWRES,
        write_rgb=WRITE_RGB,
    )


if __name__ == "__main__":
    main()
