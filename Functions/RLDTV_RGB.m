% RLDTV_RGB  Batch RGB R-L deconvolution with TV regularization.
%
% rawRGB and PSFRGB are Ny x Nx x 3 arrays.  On a GPU this processes the
% three color pages in the same batched FFT/array operations without mixing
% information across channels.

function ImageEstimate = RLDTV_RGB(rawRGB, PSFRGB, niter, beta, TV, showProgress, label)

    if nargin < 6 || isempty(showProgress)
        showProgress = false;
    end
    if nargin < 7
        label = 'RGB';
    end

    if ndims(rawRGB) ~= 3 || size(rawRGB, 3) ~= 3
        error('RLDTV_RGB expects rawRGB to be an Ny x Nx x 3 RGB array.');
    end
    if ndims(PSFRGB) ~= 3 || size(PSFRGB, 3) ~= 3
        error('RLDTV_RGB expects PSFRGB to be an H x W x 3 RGB PSF stack.');
    end

    PSFRGB = cast(PSFRGB, 'like', rawRGB);

    rawRGB = pageMax(rawRGB) - rawRGB;
    ImageEstimate = rawRGB;
    psfFlip = flipPSFStack(PSFRGB);
    [Ny, Nx, Nc] = size(rawRGB);

    Hpsf = psfToOTF2D(PSFRGB, [Ny, Nx]);
    HpsfFlip = psfToOTF2D(psfFlip, [Ny, Nx]);

    if showProgress
        tStart = tic;
    end

    for iter = 1:niter
        HI = fftshift2(fft2(ImageEstimate));
        Conv = real(ifft2(ifftshift2(Hpsf .* HI)));

        DV = rawRGB ./ Conv;
        DV(~isfinite(DV)) = 0;

        HDV = fftshift2(fft2(DV));
        DV_Conv = real(ifft2(ifftshift2(HDV .* HpsfFlip)));

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

        ImageEstimate = DV_Conv .* ImageEstimate ./ (1 - TV .* div);
        ImageEstimate = real(ImageEstimate);
        ImageEstimate(~isfinite(ImageEstimate)) = beta;
        ImageEstimate = max(ImageEstimate, beta);

        if showProgress
            if isa(ImageEstimate, 'gpuArray')
                wait(gpuDevice);
            end
            elapsed = toc(tStart);
            eta = elapsed / iter * (niter - iter);
            fprintf('\r[%s] RLDTV RGB batch iter %d/%d | Elapsed %.1fs | ETA %.1fs', ...
                    label, iter, niter, elapsed, eta);
        end
    end

    if showProgress
        fprintf('\n');
    end

    ImageEstimate = pageMax(ImageEstimate) - ImageEstimate;
    ImageEstimate = real(ImageEstimate);
    ImageEstimate(~isfinite(ImageEstimate)) = 0;

end

function y = fftshift2(x)
    y = fftshift(fftshift(x, 1), 2);
end

function y = ifftshift2(x)
    y = ifftshift(ifftshift(x, 1), 2);
end

function outPSF = flipPSFStack(inPSF)
    [Ny, Nx, ~] = size(inPSF);
    outPSF = inPSF(end:-1:1, end:-1:1, :);
    dx = abs(mod(Nx, 2) - 1);
    dy = abs(mod(Ny, 2) - 1);
    outPSF = circshift(outPSF, [dy, dx, 0]);
end

function values = pageMax(x)
    [~, ~, pages] = size(x);
    values = reshape(max(reshape(x, [], pages), [], 1), [1, 1, pages]);
end
