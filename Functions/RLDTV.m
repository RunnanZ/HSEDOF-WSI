% RLDTV Perform R-L deconvolution with total variation (TV) regularization.
%
% This function performs 2D R-L deconvolution with TV regularization for
% captured raw HSEDOF data
%
% Inputs:
%   raw                 - Captured raw HSEDOF data (blurred)
%   PSF                 - 2D projected Point Spread Function
%   niter               - The number of iteration times
%   beta                - Non-negative constraint value
%   TV                  - TV regularization parameter
%   showProgress (opt.) - logical, whether to print iteration progress (default: false)
%   label        (opt.) - char/string, label for progress (e.g. 'R','G','B')
%
% Output:
%   ImageEstimate       - 2D deconvolution HSEDOF result

function ImageEstimate = RLDTV(raw, PSF, niter, beta, TV, showProgress, label)

    % -------- handle optional arguments --------
    if nargin < 6 || isempty(showProgress)
        showProgress = false;
    end
    if nargin < 7
        label = '';
    end

    PSF = cast(PSF, 'like', raw);
    if ndims(PSF) == 3 && size(PSF, 3) == 1
        PSF = PSF(:, :, 1);
    end
    if ~ismatrix(raw) || ~ismatrix(PSF)
        error('RLDTV expects 2D raw data and a 2D PSF.');
    end

    raw = max(raw(:)) - raw;
    ImageEstimate = raw;
    psfFlip = flipPSF(PSF);
    [Ny, Nx] = size(raw);

    % FFTs of projectors do not change during iterations.
    Hpsf = psfToOTF2D(PSF, [Ny, Nx]);
    HPSFpFlip = psfToOTF2D(psfFlip, [Ny, Nx]);

    if showProgress
        tStart = tic;
    end

    for i = 1:niter
        % ===== R-L + TV iteration =====
        HI = fftshift(fftn(ImageEstimate));
        Conv = real(ifftn(ifftshift(Hpsf .* HI)));

        DV = raw ./ Conv;
        DV(~isfinite(DV)) = 0;

        HDV = fftshift(fftn(DV));
        DV_Conv = real(ifftn(ifftshift(HDV .* HPSFpFlip)));

        p = zeros(Ny, Nx, 2, 'like', ImageEstimate);
        p(:, :, 1) = ImageEstimate(:, [2:Nx, Nx]) - ImageEstimate;
        p(:, :, 2) = ImageEstimate([2:Ny, Ny], :) - ImageEstimate;

        p = 16 * p;
        normep = max(ones(1, 'like', ImageEstimate), sqrt(p(:, :, 1).^2 + p(:, :, 2).^2));
        p(:, :, 1) = p(:, :, 1) ./ normep;
        p(:, :, 2) = p(:, :, 2) ./ normep;

        zeroRow = zeros(1, Nx, 'like', ImageEstimate);
        zeroCol = zeros(Ny, 1, 'like', ImageEstimate);
        div = [p(1:Ny-1, :, 2); zeroRow] - [zeroRow; p(1:Ny-1, :, 2)];
        div = [p(:, 1:Nx-1, 1), zeroCol] - [zeroCol, p(:, 1:Nx-1, 1)] + div;

        ImageEstimate = DV_Conv .* ImageEstimate ./ (1 - TV .* div);
        ImageEstimate = real(ImageEstimate);
        ImageEstimate(~isfinite(ImageEstimate)) = beta;
        ImageEstimate = max(ImageEstimate, beta);

        % ===== command-line progress (no GUI) =====
        if showProgress
            if isa(ImageEstimate, 'gpuArray')
                wait(gpuDevice);
            end
            elapsed = toc(tStart);
            eta = elapsed / i * (niter - i);
            if isempty(label)
                fprintf('\rRLDTV iter %d/%d | Elapsed %.1fs | ETA %.1fs', ...
                        i, niter, elapsed, eta);
            else
                fprintf('\r[%s] RLDTV iter %d/%d | Elapsed %.1fs | ETA %.1fs', ...
                        label, i, niter, elapsed, eta);
            end
        end
    end

    if showProgress
        fprintf('\n');  
    end

    ImageEstimate = max(ImageEstimate(:)) - ImageEstimate;
    ImageEstimate = real(ImageEstimate);
    ImageEstimate(~isfinite(ImageEstimate)) = 0;

end
