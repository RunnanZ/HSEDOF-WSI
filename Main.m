% Main  Unified HSEDOF-WSI reconstruction entry point.
%
% Select the reconstruction method and execution device in the configuration
% block below.  The two supported methods are:
%
%   RLD     : Richardson-Lucy deconvolution with a flipped-PSF back-projector.
%   WB-ARL  : RLD with a Wiener-Butterworth unmatched back-projector.

clc
clear
close all

%% User configuration
methodName = "WB-ARL";    % "RLD" or "WB-ARL".
requestGPU = true;        % false: force CPU execution.
batchRGBOnGPU = true;     % true: process RGB channels as one GPU batch.
gpuPrecision = 'single';  % 'single' is faster; use 'double' for precision tests.

rawFile = "Raw_IFS.tif";  % "Raw_IFS.tif" or "Raw_LBC.tif".
forceRecomputePSF = false;
cropBorder = 50;
thresholdRate = 0.98;

%% Paths and execution setup
addpath("Data")
addpath("Functions")

if ~exist("Result", "dir")
    mkdir("Result")
end

methodName = normalizeMethodName(methodName);
useGPU = initializeGPU(requestGPU);

%% Load raw RGB data
Img = im2double(imread(rawFile));
[Ny, Nx, ~] = size(Img);

[ImgR, ImgG, ImgB] = splitRGBToDevice(Img, useGPU, gpuPrecision);

%% Optical and sampling parameters
opt = struct();
opt.Nz = 60;
opt.n_m = 1;
opt.lambdaR = 0.656 / opt.n_m;    % Red wavelength [um].
opt.lambdaG = 0.588 / opt.n_m;    % Green wavelength [um].
opt.lambdaB = 0.486 / opt.n_m;    % Blue wavelength [um].
opt.Mag = 20;                     % Objective magnification.
opt.pixelsize = 6.5 / opt.Mag;    % Effective camera sampling size [um].
opt.z_step = 0.25;                % Axial sampling [um].
opt.NA = 0.8;                     % Objective numerical aperture.
opt.inner_diameter = 0.875;       % Inner diameter of annular illumination.
opt.outer_diameter = 1.0;         % Outer diameter of annular illumination.

%% Load or compute RGB PSFs
[PSFR, PSFG, PSFB] = loadOrComputePSFs( ...
    Nx, Ny, opt, forceRecomputePSF, useGPU, gpuPrecision);
[PSFR, PSFG, PSFB] = movePSFsToDevice(PSFR, PSFG, PSFB, useGPU, gpuPrecision);

%% Method parameters
param = defaultMethodParameters(methodName);

%% Reconstruction
fprintf('Method: %s\n', methodName);
printExecutionMode(useGPU, batchRGBOnGPU);

tStart = tic;
switch char(methodName)
    case 'RLD'
        DeconRGB = runRLD( ...
            ImgR, ImgG, ImgB, PSFR, PSFG, PSFB, param, useGPU, batchRGBOnGPU);
        methodTag = "RLD";
        displayName = "RLD-TV";

    case 'WB-ARL'
        DeconRGB = runWBARL( ...
            ImgR, ImgG, ImgB, PSFR, PSFG, PSFB, param, useGPU, batchRGBOnGPU);
        methodTag = "WBARL";
        displayName = "WB-ARL";

end
elapsedTime = toc(tStart);
fprintf('%s reconstruction finished. Elapsed: %.1f s\n', displayName, elapsedTime);

%% Merge, color-transfer, crop, and save
DeconImg = postprocessRGB(DeconRGB, Img, thresholdRate, cropBorder);
[rawCrop, cropRows, cropCols] = cropImage(Img, cropBorder);

outputFile = "./Result/Result_" + methodTag + ".tif";
imwrite(DeconImg, outputFile);

runInfo = struct();
runInfo.methodName = methodName;
runInfo.elapsedTime = elapsedTime;
runInfo.useGPU = useGPU;
runInfo.batchRGBOnGPU = batchRGBOnGPU;
runInfo.gpuPrecision = gpuPrecision;
runInfo.cropRows = cropRows;
runInfo.cropCols = cropCols;
runInfo.param = param;
save("./Result/RunInfo_" + methodTag + ".mat", "runInfo")

figure
subplot(1, 2, 1), imshow(rawCrop, []), title('Raw HSEDOF data')
subplot(1, 2, 2), imshow(DeconImg, []), title(displayName)
