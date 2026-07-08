function H = psfToOTF2D(psf, outputSize)
% psfToOTF2D  Convert a centered 2D PSF or PSF stack to a shifted OTF.
%
% The reconstruction code stores PSFs with their peak near the spatial
% center.  This helper embeds smaller PSF kernels into the requested image
% size before FFT, so frequency-domain multiplication can be performed
% without requiring the PSF array to have the same size as the raw image.
%
% Inputs:
%   psf        - 2D PSF or H x W x C PSF stack.
%   outputSize - [Ny Nx] size of the raw image.
%
% Output:
%   H          - shifted OTF with size Ny x Nx or Ny x Nx x C.

if numel(outputSize) < 2
    error('psfToOTF2D expects outputSize to contain [Ny Nx].');
end

outNy = outputSize(1);
outNx = outputSize(2);
[psfNy, psfNx, nPages] = size(psf);

if psfNy > outNy || psfNx > outNx
    error('PSF size %s is larger than requested OTF size [%d %d].', ...
          mat2str(size(psf)), outNy, outNx);
end

if psfNy == outNy && psfNx == outNx
    paddedPSF = psf;
else
    paddedPSF = zeros(outNy, outNx, nPages, 'like', psf);
    rowStart = floor((outNy - psfNy) / 2) + 1;
    colStart = floor((outNx - psfNx) / 2) + 1;
    rowIdx = rowStart:(rowStart + psfNy - 1);
    colIdx = colStart:(colStart + psfNx - 1);
    paddedPSF(rowIdx, colIdx, :) = psf;
end

H = fftshift2(fft2(ifftshift2(paddedPSF)));

end

function y = fftshift2(x)
y = fftshift(fftshift(x, 1), 2);
end

function y = ifftshift2(x)
y = ifftshift(ifftshift(x, 1), 2);
end
