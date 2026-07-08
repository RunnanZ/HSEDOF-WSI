function [PSFR, PSFG, PSFB] = movePSFsToDevice(PSFR, PSFG, PSFB, useGPU, gpuPrecision)
% movePSFsToDevice  Optionally move RGB PSFs to GPU with selected precision.

if useGPU
    PSFR = gpuArray(cast(PSFR, gpuPrecision));
    PSFG = gpuArray(cast(PSFG, gpuPrecision));
    PSFB = gpuArray(cast(PSFB, gpuPrecision));
end

end
