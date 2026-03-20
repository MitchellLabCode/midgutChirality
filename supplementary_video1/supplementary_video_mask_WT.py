"""Backward-compatibility shim — delegates to mask_embryo.py.

This script is superseded by mask_embryo.py (STEP 1) and mask_midgut.py (STEP 2).
Running it directly is equivalent to running mask_embryo.py with the same arguments.

Script to build a WT embryo mask from an ilastik probability channel, upsample it, and apply it to the original data.
- Loads probability data from an H5 file and extracts a specified channel.
- Thresholds to create a binary mask, keeps the largest connected component, fills interior holes, and optionally dilates it.
- Upsamples the hard binary mask by a specified factor using linear interpolation.
"""

import mask_embryo 
import mask_midgut
from supplementary_video1 import project_midgut_views

if __name__ == "__main__":
    mask_embryo.main()
    mask_midgut.main()
    project_midgut_views.main()
