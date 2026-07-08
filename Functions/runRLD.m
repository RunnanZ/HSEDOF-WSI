function DeconRGB = runRLD(ImgR, ImgG, ImgB, PSFR, PSFG, PSFB, param, useGPU, batchRGBOnGPU)
% runRLD  Run baseline RLD-TV for RGB data.

if useGPU && batchRGBOnGPU
    ImgRGB = cat(3, ImgR, ImgG, ImgB);
    PSFRGB = cat(3, PSFR, PSFG, PSFB);
    DeconRGB = RLDTV_RGB(ImgRGB, PSFRGB, ...
        param.niter, param.nonnegBeta, param.TV, true, 'RLD');
    wait(gpuDevice);
    DeconRGB = gather(DeconRGB);
    return
end

channels = {'R', 'G', 'B'};
Decon = cell(1, 3);
rawChannels = {ImgR, ImgG, ImgB};
psfChannels = {PSFR, PSFG, PSFB};
for c = 1:3
    fprintf('Processing channel %s with RLD (%d/3)...\n', channels{c}, c);
    Decon{c} = RLDTV(rawChannels{c}, psfChannels{c}, ...
        param.niter, param.nonnegBeta, param.TV, true, channels{c});
end

if useGPU
    wait(gpuDevice);
    DeconRGB = gather(cat(3, Decon{:}));
else
    DeconRGB = cat(3, Decon{:});
end

end
