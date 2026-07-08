function [ImgR, ImgG, ImgB] = splitRGBToDevice(Img, useGPU, gpuPrecision)
% splitRGBToDevice  Split an RGB image and optionally move channels to GPU.

if useGPU
    ImgR = gpuArray(cast(Img(:, :, 1), gpuPrecision));
    ImgG = gpuArray(cast(Img(:, :, 2), gpuPrecision));
    ImgB = gpuArray(cast(Img(:, :, 3), gpuPrecision));
else
    ImgR = Img(:, :, 1);
    ImgG = Img(:, :, 2);
    ImgB = Img(:, :, 3);
end

end
