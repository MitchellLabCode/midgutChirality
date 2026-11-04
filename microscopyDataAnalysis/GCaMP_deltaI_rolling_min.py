import numpy as np
import tifffile
import argparse

parser = argparse.ArgumentParser(description="Rolling minimum ΔI for multi-frame TIFF")
parser.add_argument("input", help="Input multi-frame TIFF")
parser.add_argument("-o", "--output", default=None, help="Output TIFF (default: input_deltaI.tif)")
parser.add_argument("-w", "--window-half", type=int, default=5,
                    help="Half-window size (full = 2w+1, default: 5)")
args = parser.parse_args()

base = args.input.replace(".tif", "")
output_di = args.output or f"{base}_deltaI.tif"
window_half = args.window_half

# --- Load ---
stack = tifffile.imread(args.input)
print(f"Raw stack shape: {stack.shape}, dtype {stack.dtype}")
stack = stack.squeeze()
if stack.ndim == 4 and stack.shape[-1] in (3, 4):
    print("RGB(A) input detected, using first channel")
    stack = stack[..., 0]
if stack.ndim != 3:
    raise ValueError(f"Expected a 3D (T, H, W) stack after squeeze, got shape {stack.shape}")
T = stack.shape[0]
window_size = 2 * window_half + 1

if T < window_size:
    print(f"Warning: stack has {T} frames, less than window size {window_size}. "
          f"Using all frames for minimum.")

# --- Rolling minimum and delta ---
stack_f = stack.astype(np.float32)
delta = np.empty_like(stack_f)

for t in range(T):
    start = min(max(t - window_half, 0), max(T - window_size, 0))
    end = start + min(window_size, T)
    rolling_min = np.min(stack_f[start:end], axis=0)

    delta[t] = stack_f[t] - rolling_min

# --- Save ---
tifffile.imwrite(output_di, delta.astype(stack.dtype), imagej=True, metadata={'axes': 'TYX'})
print(f"Saved {T}-frame ΔI stack to {output_di}")
