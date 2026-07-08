# HSEDOF-WSI MATLAB Reconstruction

MATLAB implementation for high-speed extended depth-of-field whole-slide
imaging (HSEDOF-WSI).  The repository provides a unified reconstruction entry
point for two methods:

- `RLD`: Richardson-Lucy deconvolution with total variation regularization.
- `WB-ARL`: RLD with a Wiener-Butterworth accelerated mismatched
  back-projector.

The code supports CPU execution, GPU execution, and batched RGB GPU processing.
Small PSF kernels are automatically embedded into the raw image size and
converted to OTFs before frequency-domain convolution.

## Features

- Unified `Main.m` entry point for both RLD and WB-ARL.
- CPU/GPU switch with automatic CPU fallback when no compatible GPU is found.
- Batched RGB GPU mode for synchronized R, G, and B channel computation.
- Precomputed RGB PSF loading and optional optical-model PSF recomputation.
- PSF-to-OTF size conversion for small PSF kernels.
- Shared post-processing pipeline: color transfer, threshold clipping, border
  cropping, image export, and runtime metadata saving.

## Repository Structure

```text
HSEDOF-WSI-main/
|-- Main.m                         # Unified reconstruction script.
|-- Data/
|   |-- RAW_IFS.tif                # Example raw RGB HSEDOF-WSI image.
|   |-- RAW_LBC.tif                # Example raw RGB HSEDOF-WSI image.
|-- Functions/
|   |-- CalculateF.m               # Partially coherent WOTF helper.
|   |-- CalculatePSF.m             # 3D-to-2D projected PSF generation.
|   |-- ColorTransfer.m            # Color transfer / white balance.
|   |-- ThresholdSegmentation.m    # Bright-region clipping and normalization.
|   |-- cropImage.m                # Border cropping.
|   |-- defaultMethodParameters.m  # Default parameters for RLD and WB-ARL.
|   |-- flipPSF.m                  # Shift-aware PSF flipping.
|   |-- initializeGPU.m            # GPU initialization with CPU fallback.
|   |-- loadOrComputePSFs.m        # Load or recompute RGB PSFs.
|   |-- movePSFsToDevice.m         # Optional PSF transfer to GPU.
|   |-- normalizeMethodName.m      # Method-name validation.
|   |-- postprocessRGB.m           # Shared RGB post-processing.
|   |-- printExecutionMode.m       # Runtime mode reporting.
|   |-- psfToOTF2D.m               # PSF-to-OTF size conversion.
|   |-- RLDTV.m                    # Single-channel RLD-TV.
|   |-- RLDTV_RGB.m                # Batched RGB RLD-TV.
|   |-- RLDTV_WBARL.m              # Single-channel WB-ARL.
|   |-- RLDTV_WBARL_RGB.m          # Batched RGB WB-ARL.
|   |-- runRLD.m                   # RLD method wrapper.
|   |-- runWBARL.m                 # WB-ARL method wrapper.
|   |-- splitRGBToDevice.m         # RGB split and optional GPU transfer.
|   `-- WB_back_projector.m        # Wiener-Butterworth back-projector.
|-- PSFR.mat                       # Precomputed red-channel PSF.
|-- PSFG.mat                       # Precomputed green-channel PSF.
|-- PSFB.mat                       # Precomputed blue-channel PSF.
|-- Result/                        # Generated results and run metadata.
|-- README.md
`-- LICENSE
```

## Requirements

- MATLAB R2025a or newer is recommended.
- Image Processing Toolbox.
- Parallel Computing Toolbox and an NVIDIA CUDA GPU are recommended for GPU
  acceleration.

If `requestGPU = true` but MATLAB cannot initialize a compatible GPU, the code
falls back to CPU execution.

## Quick Start

Open MATLAB in the repository root and run:

```matlab
run("Main.m")
```

Configure the reconstruction at the top of `Main.m`:

```matlab
methodName = "WB-ARL";    % "RLD" or "WB-ARL".
requestGPU = true;        % false: force CPU execution.
batchRGBOnGPU = true;     % true: process RGB channels as one GPU batch.
gpuPrecision = 'single';  % 'single' is faster; use 'double' for precision tests.

rawFile = "Raw_LBC.tif";  % "Raw_IFS.tif", "Raw_LBC.tif", or a custom file.
forceRecomputePSF = false;
cropBorder = 50;
thresholdRate = 0.98;
```

## Reconstruction Methods

| Method | Select with | Back-projector | Default iterations |
| --- | --- | --- | --- |
| RLD | `methodName = "RLD"` | Flipped PSF | `100` |
| WB-ARL | `methodName = "WB-ARL"` | Wiener-Butterworth mismatched back-projector | `2` |

`WB-ARL` usually needs fewer iterations, but it has an additional
back-projector generation step.  A single WB-ARL run can therefore take longer
than RLD when comparing only one iteration or when using very small iteration
counts.

## CPU and GPU Modes

The runtime mode is controlled by:

```matlab
requestGPU = true;
batchRGBOnGPU = true;
gpuPrecision = 'single';
```

When `batchRGBOnGPU = true`, the RGB channels are stacked as an
`Ny x Nx x 3` array and processed with page-wise FFT operations.  The channels
are computed in the same GPU operations but are not mixed with each other.

If GPU memory is limited, use:

```matlab
batchRGBOnGPU = false;
```

This keeps GPU acceleration enabled but processes R, G, and B sequentially.
To force CPU execution, set:

```matlab
requestGPU = false;
```

## PSF and OTF Handling

The repository includes precomputed RGB PSFs:

```text
PSFR.mat
PSFG.mat
PSFB.mat
```

They are loaded by default.  To recompute PSFs from the optical model, set:

```matlab
forceRecomputePSF = true;
```

Small PSF kernels are supported.  Before frequency-domain multiplication,
`Functions/psfToOTF2D.m` embeds each PSF into the raw image size and converts
it to an OTF.  This avoids dimension mismatches when the PSF kernel is smaller
than the raw image.

## Default Parameters

Default reconstruction parameters are defined in
`Functions/defaultMethodParameters.m`.

```matlab
% RLD
param.TV = 1e-4;
param.nonnegBeta = 0;
param.niter = 100;

% WB-ARL
param.TV = 1e-5;
param.nonnegBeta = 0;
param.niter = 2;
param.wb.bp_type = 'wiener-butterworth';
param.wb.alpha = 0.1;
param.wb.beta = 0.1;
param.wb.n = 20;
param.wb.resFlag = 1;
param.wb.iRes = [0, 0, 0];
param.wb.verboseFlag = 0;
param.wb.epsValue = 1e-12;
```

For WB-ARL, `param.wb.beta` is the Wiener-Butterworth parameter.  It is
different from `param.nonnegBeta`, which is the lower bound used in the
Richardson-Lucy update.

## Outputs

`Main.m` creates `Result/` if needed and writes:

```text
Result/Result_RLD.tif       # output when methodName = "RLD"
Result/Result_WBARL.tif     # output when methodName = "WB-ARL"
Result/RunInfo_RLD.mat      # runtime metadata for RLD
Result/RunInfo_WBARL.mat    # runtime metadata for WB-ARL
```

The run metadata includes the selected method, elapsed time, GPU mode,
precision, crop indices, and reconstruction parameters.

## Troubleshooting

- If MATLAB reports that no GPU is available, set `requestGPU = false` or
  check the NVIDIA driver and Parallel Computing Toolbox installation.
- If GPU memory is insufficient, keep `requestGPU = true` and set
  `batchRGBOnGPU = false`.
- If a custom PSF is larger than the raw image, recompute the PSF or provide a
  smaller kernel.
- If border artifacts are visible, adjust `cropBorder`.
- If the output is over-clipped or under-clipped, adjust `thresholdRate`.

## Citation

If you use this code in academic work, please cite the corresponding
HSEDOF-WSI paper or project release:

> High-speed extended depth-of-field whole-slide imaging (HSEDOF-WSI):
> bridging resolution, depth, and throughput for advanced digital pathology.

## License

This project is released under the MIT License. 

## Contact

- Runnan Zhang: runnanzhang@njust.edu.com
- Chao Zuo: zuochao@njust.edu.com
