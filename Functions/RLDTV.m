% RLDTV Perform R-L deconvolution with total variation (TV) regularization.
%
% This function perform 2D R-L deconvolution with TV regularization for
% captured raw HSEDOF data
%
% Inputs:
%   raw                 - Captured raw HSEDOF data (blurred)
%   PSF                 - 2D projected Point Spread Function
%   niter               - The number of iteration times
%   beta                - Non-negative constraint value
%   TV                  - TV regularization parameter
%
% Output:
%   ImageEstimate       - 2D deconvolution HSEDOF result



function ImageEstimate = RLDTV(raw, PSF, niter, beta, TV)

raw = max(raw(:))-raw;
ImageEstimate = raw;
psfFlip = flipPSF(PSF);
[Ny Nx] = size(raw);

    for i = 1:niter

    disp(i)
    Hpsf = fftshift(fftn(ifftshift(PSF)));
    HI = (fftshift(fftn(ImageEstimate)));
    Conv = ifftn(ifftshift(Hpsf.*HI));
    DV = raw./Conv;
    HPSFpFlip = fftshift(fftn(ifftshift(psfFlip)));
    HDV = fftshift(fftn(DV));
    DV_Conv = ifftn(ifftshift(HDV.*HPSFpFlip));

    p(:, :, 1)=ImageEstimate(:, [2:Nx, Nx]) - ImageEstimate;
    p(:, :, 2)=ImageEstimate([2:Ny, Ny], :) - ImageEstimate;

    ux = ImageEstimate(:, [2:Nx, Nx]) - ImageEstimate;
    uy = ImageEstimate([2:Ny, Ny], :) - ImageEstimate;

    p = p + 15*cat(3, ux,uy);
    normep = max(1, sqrt(p(:, :, 1).^2 + p(:, :, 2).^2));
    p(:, :, 1) = p(:, :, 1)./normep;
    p(:, :, 2) = p(:, :, 2)./normep;

    div = [p([1:Ny-1], :, 2); zeros(1, Nx)] - [zeros(1, Nx); p([1:Ny-1], :, 2)];
    div = [p(:, [1:Nx-1], 1)  zeros(Ny, 1)] - [zeros(Ny, 1)  p(:, [1:Nx-1], 1)] + div;

    ImageEstimate = DV_Conv.*ImageEstimate./((1-TV.*div));
    ImageEstimate = max(ImageEstimate,beta);     

    end
    
    ImageEstimate = max(ImageEstimate(:))-ImageEstimate;

end

