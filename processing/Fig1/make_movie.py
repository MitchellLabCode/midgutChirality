#!/usr/bin/env python3
r"""Assemble PNG frame folders into one MP4, with holds and fades to white.

    python3 make_movie.py OUT.mp4 FOLDER[:HOLD[:RATE]] [...] [--fps 10] [--rate RATE]
                          [--hold 3] [--fade 1] [--size 1920x1080] [--ffmpeg PATH]

--fps is the frame rate of the output file; RATE (per folder, default --rate,
which defaults to --fps) is how many of that folder's PNGs are shown per
second. Each PNG is repeated or skipped as needed, so e.g. --fps 20 with one
folder at RATE 20 and others at 10 plays the first twice as fast.

Each folder's PNGs (sorted by name) play at --fps, are fitted onto a white
canvas of --size (aspect kept, padded white), and end with a freeze of their
last frame for HOLD seconds (per folder, default --hold). Between folders the
picture fades to white over --fade seconds, then fades in from white to the
next folder's first frame over --fade seconds before it plays, so a folder
with a single PNG works as a title or still card.

The colour conversion is done here (BT.709, limited range) and raw yuv420p
is piped to ffmpeg, which only encodes (MPEG-4 Part 2 at -q:v 2), so the
script also works with minimal ffmpeg builds without libx264 or a working
RGB->YUV converter. ffmpeg is taken from --ffmpeg, $FFMPEG, or the PATH.

The Fig. 1 chirality movie (frames from fig1_tantrix_movie.m,
fig1_qmeasure_frame.m and fig1_chiral_movie.m with WT_idx = 4 and
fig_title = 'The resulting shapes display \itperverted\rm chirality'):

    python3 make_movie.py tantrix_qmeasure_chiral3d_s4_s1.mp4 \
        output/tantrix_s4:3:20 output/qmeasure_frame:4:10 output/chiral3d_s4_s1:3:10 \
        --fps 20 --fade 1
"""
import argparse
import glob
import os
import shutil
import subprocess
import sys

import numpy as np
from PIL import Image

FFMPEG_DEFAULT = os.environ.get('FFMPEG') or shutil.which('ffmpeg') or 'ffmpeg'


def fit_on_canvas(img, W, H):
    """Scale img to fit W x H (aspect kept) and centre it on white."""
    w, h = img.size
    s = min(W / w, H / h)
    nw, nh = max(2, round(w * s)), max(2, round(h * s))
    img = img.resize((nw, nh), Image.LANCZOS)
    canvas = Image.new('RGB', (W, H), 'white')
    canvas.paste(img, ((W - nw) // 2, (H - nh) // 2))
    return np.asarray(canvas, dtype=np.float32)


def rgb_to_yuv420_bt709(rgb):
    """float RGB (0..255) -> planar Y, U, V (uint8), BT.709 limited range, 4:2:0."""
    r, g, b = rgb[..., 0] / 255, rgb[..., 1] / 255, rgb[..., 2] / 255
    y = 0.2126 * r + 0.7152 * g + 0.0722 * b
    Y = 16 + 219 * y
    U = 128 + 224 * (b - y) / 1.8556
    V = 128 + 224 * (r - y) / 1.5748
    sub = lambda C: C.reshape(C.shape[0] // 2, 2, C.shape[1] // 2, 2).mean(axis=(1, 3))
    return [np.clip(np.round(p), 0, 255).astype(np.uint8) for p in (Y, sub(U), sub(V))]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('out')
    ap.add_argument('folders', nargs='+')
    ap.add_argument('--fps', type=float, default=10, help='output frame rate')
    ap.add_argument('--rate', type=float, default=None, help='default playback rate of each folder [PNGs/s]')
    ap.add_argument('--hold', type=float, default=3, help='freeze on the last frame of each folder [s]')
    ap.add_argument('--fade', type=float, default=1, help='fade to / from white between folders [s]')
    ap.add_argument('--size', default='1920x1080')
    ap.add_argument('--ffmpeg', default=FFMPEG_DEFAULT)
    a = ap.parse_args()

    W, H = (int(v) for v in a.size.split('x'))
    n_fade = round(a.fade * a.fps)
    rate_default = a.rate if a.rate else a.fps
    segs = []
    for spec in a.folders:
        parts = spec.split(':')
        folder = parts[0]
        hold = float(parts[1]) if len(parts) > 1 and parts[1] else a.hold
        rate = float(parts[2]) if len(parts) > 2 and parts[2] else rate_default
        segs.append((folder, round(hold * a.fps), rate))
    white = np.full((H, W, 3), 255, np.float32)

    cmd = [a.ffmpeg, '-hide_banner', '-loglevel', 'error', '-y',
           '-f', 'rawvideo', '-pix_fmt', 'yuv420p', '-s', f'{W}x{H}', '-r', str(a.fps), '-i', '-',
           '-c:v', 'mpeg4', '-q:v', '2',
           '-colorspace', 'bt709', '-color_primaries', 'bt709', '-color_trc', 'bt709', '-color_range', 'tv',
           a.out]
    enc = subprocess.Popen(cmd, stdin=subprocess.PIPE)

    def write(frame):
        enc.stdin.write(b''.join(p.tobytes() for p in rgb_to_yuv420_bt709(frame)))

    n_out = 0
    for k, (folder, n_hold, rate) in enumerate(segs):
        files = sorted(glob.glob(os.path.join(folder, '*.png')))
        if not files:
            sys.exit(f'no PNGs in {folder}')
        last = len(segs) - 1
        frame = None
        # output frames spanned by this folder's PNGs at the requested rate;
        # PNG i covers output frames [round(i*r), round((i+1)*r)), r = fps/rate
        r = a.fps / rate
        n_play = max(1, round(len(files) * r))
        shown = 0
        for i, f in enumerate(files):
            reps = round((i + 1) * r) - round(i * r)
            if reps <= 0 and i < len(files) - 1:          # faster than the output fps: skip
                continue
            frame = fit_on_canvas(Image.open(f).convert('RGB'), W, H)
            if k > 0 and shown == 0:                      # fade in from white
                for j in range(n_fade):
                    al = (j + 1) / n_fade
                    write(al * frame + (1 - al) * white); n_out += 1
            for _ in range(max(reps, 1 if shown == 0 else reps)):
                write(frame); n_out += 1; shown += 1
        for _ in range(n_hold):                           # hold the last frame
            write(frame); n_out += 1
        if k < last:                                      # fade out to white
            for i in range(n_fade):
                al = (i + 1) / n_fade
                write((1 - al) * frame + al * white); n_out += 1
        print(f'{folder}: ' + (f'{n_fade} fade in + ' if k > 0 else '')
              + f'{len(files)} PNGs at {rate:g}/s ({shown} frames, {shown / a.fps:.1f} s) + {n_hold} hold'
              + (f' + {n_fade} fade out' if k < last else ''))

    enc.stdin.close()
    if enc.wait() != 0:
        sys.exit('ffmpeg failed')
    print(f'{a.out}: {n_out} frames, {n_out / a.fps:.1f} s at {a.fps:g} fps, {W}x{H}')


if __name__ == '__main__':
    main()
