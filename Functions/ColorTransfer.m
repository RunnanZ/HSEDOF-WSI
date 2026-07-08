% ColorTransfer Perform white balance correction using reference image
%
% This function implements a white balance algorithm that transfers color 
% characteristics from a reference image to a source image in YCbCr color 
% space. The algorithm adjusts the luminance distribution while preserving 
% the chrominance information from the reference.
%
% Inputs:
%   Img - Source image to be white balanced (RGB, MxNx3)
%   Ref - Reference image with desired color characteristics (RGB, MxNx3)
%
% Output:
%   result - White balanced image (RGB, MxNx3)
%
% Algorithm:
%   1. Convert both images from RGB to YCbCr color space
%   2. Extract luminance (Y) and chrominance (Cb, Cr) components
%   3. Calculate mean and standard deviation of luminance channels
%   4. Apply linear transformation to source luminance:
%        Result_Y = (Ref_std/Img_std) * (Img_Y - Img_mean) + Ref_mean
%   5. Combine transformed luminance with reference chrominance
%   6. Convert result back to RGB color space
%
% Method Characteristics:
%   - Operates in YCbCr color space to separate luminance and chrominance
%   - Preserves reference image's color tone (Cb, Cr components)
%   - Adjusts source image's brightness and contrast to match reference
%   - Maintains the structural content of the source image

function result = ColorTransfer(Img, Ref)
%% RGB -> YCbCr
Img_YCbCr = rgb2ycbcr(Img);
Ref_YCbCr = rgb2ycbcr(Ref);
Img_Y = Img_YCbCr(:,:,1);
Ref_Y = Ref_YCbCr(:,:,1);
Ref_Cb = Ref_YCbCr(:,:,2);
Ref_Cr = Ref_YCbCr(:,:,3);

%% Color transmission
Img_Y_mean = mean(Img_Y,'all');
Img_Y_std = std(Img_Y,[],'all');
Ref_Y_mean = mean(Ref_Y,'all'); 
Ref_Y_std = std(Ref_Y,[],'all');

Result_Y = Ref_Y_std / Img_Y_std * (Img_Y - Img_Y_mean) + Ref_Y_mean;
Result_Cb = Ref_Cb;
Result_Cr = Ref_Cr;

%% YCbCr -> RGB
result_YCbCr = zeros(size(Img));
result_YCbCr(:, :, 1) =  Result_Y; result_YCbCr(:, :, 2) =  Result_Cb; result_YCbCr(:, :, 3) =  Result_Cr;
result = ycbcr2rgb(result_YCbCr);

end

