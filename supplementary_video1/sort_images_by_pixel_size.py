#!/usr/bin/env python3
"""Inspect and sort images by pixel size.

This script scans an image directory, reads each file's pixel dimensions, and groups
(or optionally copies/moves) files into size-specific buckets.

Expected common sizes in this project:
- 1199x1124
- 1209x1209
"""

from __future__ import annotations

import argparse
import struct
import shutil
from collections import Counter
from pathlib import Path
from typing import Iterable

import tifffile


rootdir = "/Users/npmitchell/Dropbox/Soft_Matter_Biophysics/PAPER/gut_chirality/supplementary_videos/"
DEFAULT_EXTENSIONS = (".png", ".tif", ".tiff")
KNOWN_SIZES = {(1199, 1124), (1209, 1209)}  # (width, height)
DEFAULT_INPUT_DIR = Path(
    rootdir+"/mesh"
)
DEFAULT_OUTPUT_DIR = Path(rootdir+"mesh/sorted_by_size")

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Inspect image dimensions and sort files by pixel size."
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=DEFAULT_INPUT_DIR,
        help="Directory containing the image sequence.",
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Recursively scan subdirectories.",
    )
    parser.add_argument(
        "--mode",
        choices=("report", "copy", "move"),
        default="move",
        help="report: print grouping only; copy/move: place files into grouped folders.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help=(
            "Output root for copy/move mode. Defaults to <input-dir> for move mode "
            "and <input-dir>/sorted_by_size for copy mode."
        ),
    )
    parser.add_argument(
        "--extensions",
        nargs="+",
        default=list(DEFAULT_EXTENSIONS),
        help="File extensions to include (case-insensitive).",
    )
    parser.add_argument(
        "--strict-known-sizes",
        action="store_true",
        help="Exit with code 1 if any image has dimensions outside known sizes.",
    )
    return parser.parse_args()


def iter_image_paths(root: Path, recursive: bool, extensions: Iterable[str]) -> list[Path]:
    exts = {e.lower() if e.startswith(".") else f".{e.lower()}" for e in extensions}
    pattern_iter = root.rglob("*") if recursive else root.glob("*")
    paths = [p for p in pattern_iter if p.is_file() and p.suffix.lower() in exts]
    return sorted(paths)


def shape_to_height_width(shape: tuple[int, ...]) -> tuple[int, int]:
    if len(shape) == 2:
        return int(shape[0]), int(shape[1])
    if len(shape) == 3:
        if shape[-1] in (3, 4):
            return int(shape[0]), int(shape[1])
        if shape[0] in (3, 4):
            return int(shape[1]), int(shape[2])
        return int(shape[-2]), int(shape[-1])
    if len(shape) >= 4:
        return int(shape[-2]), int(shape[-1])
    raise ValueError(f"Unsupported image shape: {shape}")


def read_image_size(path: Path) -> tuple[int, int]:
    """Return image size as (width, height)."""
    suffix = path.suffix.lower()
    if suffix in {".tif", ".tiff"}:
        with tifffile.TiffFile(path) as tif:
            first_page = tif.pages[0]
            height, width = shape_to_height_width(tuple(int(x) for x in first_page.shape))
        return int(width), int(height)

    if suffix == ".png":
        with path.open("rb") as f:
            header = f.read(24)
        if len(header) < 24 or header[:8] != b"\x89PNG\r\n\x1a\n":
            raise ValueError("Invalid PNG header")
        width, height = struct.unpack(">II", header[16:24])
        return int(width), int(height)

    raise ValueError(f"Unsupported extension for size read: {suffix}")


def bucket_name(size_wh: tuple[int, int]) -> str:
    w, h = size_wh
    return f"{w}x{h}"


def ensure_unique_destination(dest: Path) -> Path:
    if not dest.exists():
        return dest

    stem = dest.stem
    suffix = dest.suffix
    parent = dest.parent
    idx = 1
    while True:
        candidate = parent / f"{stem}_{idx}{suffix}"
        if not candidate.exists():
            return candidate
        idx += 1


def main() -> int:
    args = parse_args()

    input_dir = args.input_dir.expanduser().resolve()
    if not input_dir.exists() or not input_dir.is_dir():
        raise FileNotFoundError(f"Input directory does not exist: {input_dir}")

    paths = iter_image_paths(input_dir, args.recursive, args.extensions)
    if not paths:
        print("No matching image files found.")
        return 0

    ext_counts = Counter(p.suffix.lower() for p in paths)

    grouped: dict[tuple[int, int], list[Path]] = {}
    failures: list[tuple[Path, str]] = []

    for path in paths:
        try:
            size = read_image_size(path)
        except Exception as exc:
            failures.append((path, str(exc)))
            continue
        grouped.setdefault(size, []).append(path)

    counts = Counter({bucket_name(k): len(v) for k, v in grouped.items()})

    print(f"Scanned {len(paths)} files from: {input_dir}")
    print("Extension types:")
    for ext, count in sorted(ext_counts.items()):
        print(f"  {ext}: {count}")
    for name, count in sorted(counts.items()):
        print(f"  {name}: {count}")

    if failures:
        print(f"\nFailed to read {len(failures)} files:")
        for path, msg in failures:
            print(f"  {path.name}: {msg}")

    unknown_sizes = [s for s in grouped if s not in KNOWN_SIZES]
    if unknown_sizes:
        print("\nFound sizes outside expected set {1199x1124, 1209x1209}:")
        for s in sorted(unknown_sizes):
            print(f"  {s[0]}x{s[1]}")

    if args.mode == "report":
        print("\nReport mode only: no files were moved or copied.")
        return 1 if (args.strict_known_sizes and unknown_sizes) else 0

    if args.output_dir is not None:
        output_root = args.output_dir.expanduser().resolve()
    else:
        output_root = input_dir if args.mode == "move" else (input_dir / "sorted_by_size")
    output_root.mkdir(parents=True, exist_ok=True)

    op = shutil.copy2 if args.mode == "copy" else shutil.move

    for size, files in grouped.items():
        target_dir = output_root / bucket_name(size)
        target_dir.mkdir(parents=True, exist_ok=True)
        for src in files:
            dest = ensure_unique_destination(target_dir / src.name)
            op(src, dest)

    print(f"\nCompleted {args.mode} into: {output_root}")
    return 1 if (args.strict_known_sizes and unknown_sizes) else 0


if __name__ == "__main__":
    raise SystemExit(main())
