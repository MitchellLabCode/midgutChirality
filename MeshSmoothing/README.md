Workflow for the smoothing of meshes

# Midgut Surface Mesh Processing

This repository contains Python scripts for post-processing three-dimensional surface meshes of the embryonic *Drosophila* midgut, as used in "Mechanical canalization of chiral gut morphogenesis"

## Overview

The input meshes are generated using `tubi.getMeshes()` from [TubULAR](https://github.com/npmitchell/tubular). These meshes are subsequently processed using PyMeshLab to obtain smooth, closed surfaces with consistent normal vectors for downstream geometric analysis. These scripts go through this process. 

The processing pipeline consists of three sequential steps:

1. **Mesh filtering (`Batch_filtering_and_reconstruction.py`)**  
   We perform multiple rounds of HC Laplacian smoothing to reduce local surface irregularities while preserving the overall geometry. Poisson-disk sampling is then used to obtain a more uniform distribution of surface points (20,000 samples), followed by screened Poisson surface reconstruction (depth = 8, scale = 1.1) to generate a continuous mesh. The reconstructed surface is further smoothed and refined through 10 iterations of Loop subdivision. Finally, remaining holes are closed to obtain a watertight surface.

2. **Genus verification (`Check_genus.py`)**  
   We compute topological measures of the processed meshes and identify any meshes with genus greater than zero. This step is used to verify that the reconstructed surfaces satisfy the genus-zero requirement for subsequent TubULAR analysis. Meshes failing this check are reported for further inspection.

3. **Normal generation (`Generate_normals.py`)**  
   We remove duplicate faces and vertices, then compute surface normals using PyMeshLab to ensure consistent signs. The resulting meshes are saved with vertex normals for subsequent geometric processing.

## Requirements

- Python 3
- [PyMeshLab](https://pymeshlab.readthedocs.io/)

Install PyMeshLab using:

`pip install pymeshlab`

## Usage

1. Generate the initial surface meshes using `tubi.getMeshes()` in TubULAR.
2. Set the input and output directories in `Batch_filtering_and_reconstruction.py` and run the script.
3. Set the directory containing the filtered meshes in `Check_genus.py` and run it to identify any meshes with nonzero genus.
4. 
5. After fix any flagged meshes, set the input and output directories in `Generate_normals.py` and run the script.

The final outputs are processed surface meshes with saved normal vectors, suitable for downstream analysis in TubULAR.

**Note:** The genus-checking script reports meshes with genus greater than zero but does not automatically correct their topology. Closure and genus-zero topology should be done by refinining Ilastik training, and rerunning tubi.getMeshes(), or manually via Meshlab or your mesh-processing software of choice.






