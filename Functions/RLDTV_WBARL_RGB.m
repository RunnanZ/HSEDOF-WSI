% RLDTV_WBARL_RGB  Batch RGB RLD-TV with WB mismatched back-projectors.
%
% The WB back-projector is still generated per channel because its cutoff
% depends on each PSF's FWHM.  The iterative deconvolution then runs as
% batched RGB page operations on Ny x Nx x 3 arrays.

function [ImageEstimate, PSF_bp] = RLDTV_WBARL_RGB(rawRGB, PSF_fp, niter, nonnegBeta, TV, wbParam, showProgress, label)

    if nargin < 6 || isempty(wbParam)
        wbParam = struct();
    end
    if nargin < 7 || isempty(showProgress)
        showProgress = false;
    end
    if nargin < 8
        label = 'RGB';
    end

    if ndims(rawRGB) ~= 3 || size(rawRGB, 3) ~= 3
        error('RLDTV_WBARL_RGB expects rawRGB to be an Ny x Nx x 3 RGB array.');
    end
    if ndims(PSF_fp) ~= 3 || size(PSF_fp, 3) ~= 3
        error('RLDTV_WBARL_RGB expects PSF_fp to be an H x W x 3 RGB PSF stack.');
    end

    if ~isa(rawRGB, 'gpuArray')
        rawRGB = double(rawRGB);
    end
    PSF_fp = cast(PSF_fp, 'like', rawRGB);

    psfSum = pageSum(PSF_fp);
    psfSumCheck = psfSum;
    if isa(psfSumCheck, 'gpuArray')
        psfSumCheck = gather(psfSumCheck);
    end
    valid = isfinite(psfSumCheck) & abs(psfSumCheck) > 0;
    if all(valid(:))
        PSF_fp = PSF_fp ./ psfSum;
    else
        warning('HSEDOF:PSFNormalizationSkipped', ...
                'At least one RGB PSF has zero or non-finite sum; PSF normalization was skipped.');
    end

    PSF_bp = getBackProjectorStack(PSF_fp, wbParam, rawRGB);
    epsValue = cast(getOption(wbParam, 'epsValue', 1e-12), 'like', rawRGB);

    rawRGB = pageMax(rawRGB) - rawRGB;
    ImageEstimate = rawRGB;
    [Ny, Nx, Nc] = size(rawRGB);

    Hpsf = psfToOTF2D(PSF_fp, [Ny, Nx]);
    Hbp = psfToOTF2D(PSF_bp, [Ny, Nx]);

    if showProgress
        tStart = tic;
    end

    for iter = 1:niter
        HI = fftshift2(fft2(ImageEstimate));
        Conv = real(ifft2(ifftshift2(Hpsf .* HI)));
        Conv = max(Conv, epsValue);

        DV = rawRGB ./ Conv;
        DV(~isfinite(DV)) = 0;

        HDV = fftshift2(fft2(DV));
        DV_Conv = real(ifft2(ifftshift2(HDV .* Hbp)));

        ux = ImageEstimate(:, [2:Nx, Nx], :) - ImageEstimate;
        uy = ImageEstimate([2:Ny, Ny], :, :) - ImageEstimate;

        px = 16 * ux;
        py = 16 * uy;
        normep = max(ones(1, 'like', ImageEstimate), sqrt(px.^2 + py.^2));
        px = px ./ normep;
        py = py ./ normep;

        zeroRow = zeros(1, Nx, Nc, 'like', ImageEstimate);
        zeroCol = zeros(Ny, 1, Nc, 'like', ImageEstimate);
        div = [py(1:Ny-1, :, :); zeroRow] - [zeroRow; py(1:Ny-1, :, :)];
        div = [px(:, 1:Nx-1, :), zeroCol] - [zeroCol, px(:, 1:Nx-1, :)] + div;

        denom = 1 - TV .* div;
        denom = max(denom, epsValue);

        ImageEstimate = DV_Conv .* ImageEstimate ./ denom;
        ImageEstimate = real(ImageEstimate);
        ImageEstimate(~isfinite(ImageEstimate)) = nonnegBeta;
        ImageEstimate = max(ImageEstimate, nonnegBeta);

        if showProgress
            if isa(ImageEstimate, 'gpuArray')
                wait(gpuDevice);
            end
            elapsed = toc(tStart);
            eta = elapsed / iter * (niter - iter);
            fprintf('\r[%s] RLDTV-WBARL RGB batch iter %d/%d | Elapsed %.1fs | ETA %.1fs', ...
                    label, iter, niter, elapsed, eta);
        end
    end

    if showProgress
        fprintf('\n');
    end

    ImageEstimate = pageMax(ImageEstimate) - ImageEstimate;
    ImageEstimate = real(ImageEstimate);
    ImageEstimate(~isfinite(ImageEstimate)) = 0;
    ImageEstimate = max(ImageEstimate, 0);

end

function PSF_bp = getBackProjectorStack(PSF_fp, wbParam, rawRGB)
    if isfield(wbParam, 'PSF_bp') && ~isempty(wbParam.PSF_bp)
        PSF_bp = cast(wbParam.PSF_bp, 'like', rawRGB);
        if ndims(PSF_bp) ~= 3 || size(PSF_bp, 3) ~= size(PSF_fp, 3)
            error('wbParam.PSF_bp must have the same number of pages as PSF_fp.');
        end
        return
    end

    bp_type = getOption(wbParam, 'bp_type', 'wiener-butterworth');
    wbAlpha = getOption(wbParam, 'alpha', 0.001);
    wbBeta = getOption(wbParam, 'beta', 0.001);
    wbOrder = getOption(wbParam, 'n', 10);
    resFlag = getOption(wbParam, 'resFlag', 1);
    iRes = getOption(wbParam, 'iRes', [0, 0, 0]);
    verboseFlag = getOption(wbParam, 'verboseFlag', 0);

    iRes = iRes(:).';
    if numel(iRes) == 2
        iRes = [iRes, 1];
    elseif numel(iRes) < 2
        iRes = [0, 0, 0];
    elseif numel(iRes) > 3
        iRes = iRes(1:3);
    end

    [Ny, Nx, Nc] = size(PSF_fp);
    PSF_bp = zeros(Ny, Nx, Nc, 'like', rawRGB);
    for c = 1:Nc
        [PSF_bp(:, :, c), ~] = WB_back_projector(PSF_fp(:, :, c), ...
            bp_type, wbAlpha, wbBeta, wbOrder, resFlag, iRes, verboseFlag);
    end
end

function y = fftshift2(x)
    y = fftshift(fftshift(x, 1), 2);
end

function y = ifftshift2(x)
    y = ifftshift(ifftshift(x, 1), 2);
end

function values = pageMax(x)
    [~, ~, pages] = size(x);
    values = reshape(max(reshape(x, [], pages), [], 1), [1, 1, pages]);
end

function values = pageSum(x)
    [~, ~, pages] = size(x);
    values = reshape(sum(reshape(x, [], pages), 1), [1, 1, pages]);
end

function value = getOption(s, fieldName, defaultValue)
    if isstruct(s) && isfield(s, fieldName) && ~isempty(s.(fieldName))
        value = s.(fieldName);
    else
        value = defaultValue;
    end
end
