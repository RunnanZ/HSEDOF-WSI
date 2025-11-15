# ⭐ HSEDOF-WSI MATLAB Code Repository
High-speed extended depth-of-field whole-slide imaging (HSEDOF-WSI):
Bridging resolution, depth, and throughput for advanced digital pathology

Version: 1.0

## 📌 Overview

This repository contains the official MATLAB implementation of the algorithms used in:

“High-speed extended depth-of-field whole-slide imaging (HSEDOF-WSI): bridging resolution, depth, and throughput for advanced digital pathology.”

The package reproduces the complete HSEDOF-WSI computational pipeline, including:

- Partially coherent weak object transfer function (WOTF) for annular illumination

- 3D → 2D projected PSF computation

- Richardson–Lucy deconvolution with total variation regularization (RLD–TV)

- Color transfer correction

- Threshold-based artifact suppression

These components enable fast, high-resolution, extended-depth whole-slide imaging reconstruction consistent with the methodology described in the paper.

## 📂 Repository Structure

```text
HSEDOF-WSI/
│
├── Main.m                           # Main workflow script
│
├── Functions/                       # Core algorithm implementations
│   ├── CalculateF.m                 # Partially coherent WOTF calculation
│   ├── CalculatePSF.m               # 3D→2D projected PSF computation
│   ├── RLDTV.m                      # Richardson–Lucy with TV regularization
│   ├── ColorTransfer.m              # Color transfer / white balance
│   └── ThresholdSegmentation.m      # Adaptive artifact suppression
│
├── Data/                            # Example HSEDOF-WSI raw data
│   └── Raw.tif
│
├── Result/                          # Example HSEDOF-WSI deconvolution data
│   └── Result.tif
│
├── PSF_RGB.mat                      # Pre-computed PSF (R/G/B channels)
└── README.md                        # This document


``` 

## 🚀 Quick Start

The code has been tested in **MATLAB 2024b (64-bit)** on **Windows 11 (64-bit)** with the following hardware:

- Intel i9-13900KF CPU  
- NVIDIA GeForce RTX 4060 GPU  
- 64 GB RAM  

Steps:

1. Unpack the package.  
2. Add all subdirectories to your MATLAB path.  
3. Run the scripts whose filenames begin with **“Main”** to process the provided example data.

Additional notes:  
- Raw HSEDOF-WSI example data are stored in the `Data` folder.  
- Pre-computed RGB channel PSFs (`*.mat`) are included in the main directory.  
      Users may also generate their own PSFs using `CalculatePSF.m` or by loading experimentally calibrated PSFs.

For further details, please refer to the individual MATLAB source files.

## 🧩 Function Descriptions

All core algorithms are stored in the functions/ folder.

### 1. CalculateF.m

Computes the partially coherent weak object transfer function (WOTF) under annular illumination.

### 2. CalculatePSF.m

Generates the projected 2D PSF:

Constructs 3D AOTF

Performs inverse FT to obtain APSF

Projects APSF into a 2D PSF

Supports arbitrary ring geometry or calibrated PSFs

### 3. RLDTV.m

Richardson–Lucy deconvolution with total variation regularization:

Noise-suppressing

Edge-preserving

### 4. ColorTransfer.m

Color transfer & white-balance correction based on a reference background image.

### 5. ThresholdSegmentation.m

Performs channel-wise adaptive clipping to suppress oversaturated bright regions and reduce halo-like artifacts.

Provided sample available in /Data/

## 📬 Contact

For questions, please contact:

Runnan Zhang
Email: runnanzhang@njust.edu.com

or

Chao Zuo
Email: zuochao@njust.edu.com
