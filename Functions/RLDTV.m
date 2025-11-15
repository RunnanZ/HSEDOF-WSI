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

    raw = max(raw(:)) - raw;
    ImageEstimate = raw;
    psfFlip = flipPSF(PSF);
    [Ny, Nx] = size(raw);

    if showProgress
        tStart = tic;
    end

    for i = 1:niter
        % ===== R-L + TV iteration =====
        Hpsf = fftshift(fftn(ifftshift(PSF)));
        HI = fftshift(fftn(ImageEstimate));
        Conv = ifftn(ifftshift(Hpsf .* HI));

        DV = raw ./ Conv;

        HPSFpFlip = fftshift(fftn(ifftshift(psfFlip)));
        HDV = fftshift(fftn(DV));
        DV_Conv = ifftn(ifftshift(HDV .* HPSFpFlip));

        p(:, :, 1) = ImageEstimate(:, [2:Nx, Nx]) - ImageEstimate;
        p(:, :, 2) = ImageEstimate([2:Ny, Ny], :) - ImageEstimate;

        ux = ImageEstimate(:, [2:Nx, Nx]) - ImageEstimate;
        uy = ImageEstimate([2:Ny, Ny], :) - ImageEstimate;

        p = p + 15 * cat(3, ux, uy);
        normep = max(1, sqrt(p(:, :, 1).^2 + p(:, :, 2).^2));
        p(:, :, 1) = p(:, :, 1) ./ normep;
        p(:, :, 2) = p(:, :, 2) ./ normep;

        div = [p(1:Ny-1, :, 2); zeros(1, Nx)] - [zeros(1, Nx); p(1:Ny-1, :, 2)];
        div = [p(:, 1:Nx-1, 1), zeros(Ny, 1)] - [zeros(Ny, 1), p(:, 1:Nx-1, 1)] + div;

        ImageEstimate = DV_Conv .* ImageEstimate ./ (1 - TV .* div);
        ImageEstimate = max(ImageEstimate, beta);

        % ===== command-line progress (no GUI) =====
        if showProgress
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

end
