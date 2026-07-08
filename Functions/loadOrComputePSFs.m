function [PSFR, PSFG, PSFB] = loadOrComputePSFs(Nx, Ny, opt, forceRecomputePSF, useGPU, gpuPrecision)
% loadOrComputePSFs  Load precomputed RGB PSFs or recompute them from optics.

computePSF = forceRecomputePSF || ...
             exist("PSFR.mat", "file") ~= 2 || ...
             exist("PSFG.mat", "file") ~= 2 || ...
             exist("PSFB.mat", "file") ~= 2;

if ~computePSF
    data = load("PSFR.mat", "PSFR");
    PSFR = data.PSFR;
    data = load("PSFG.mat", "PSFG");
    PSFG = data.PSFG;
    data = load("PSFB.mat", "PSFB");
    PSFB = data.PSFB;
    return
end

delta_x = 1 / (opt.pixelsize * Nx);
delta_y = 1 / (opt.pixelsize * Ny);
delta_z = 1 / (opt.z_step * opt.Nz);

fx = (-fix(Nx / 2):fix((Nx - 1) / 2)) * delta_x;
fy = (-fix(Ny / 2):fix((Ny - 1) / 2)) * delta_y;
fz = (-fix(opt.Nz / 2):fix((opt.Nz - 1) / 2)) * delta_z;

if useGPU
    fx = gpuArray(cast(fx, gpuPrecision));
    fy = gpuArray(cast(fy, gpuPrecision));
    fz = gpuArray(cast(fz, gpuPrecision));
end

[rhox, rhoy, eta] = meshgrid(fx, fy, fz);
rhoxy = sqrt(rhox .^ 2 + rhoy .^ 2);
rho = sqrt(rhoxy .^ 2 + eta .^ 2);

PSFR = CalculatePSF(opt.NA, rhoxy, eta, opt.lambdaR, rho, ...
                    opt.outer_diameter, opt.inner_diameter);
PSFG = CalculatePSF(opt.NA, rhoxy, eta, opt.lambdaG, rho, ...
                    opt.outer_diameter, opt.inner_diameter);
PSFB = CalculatePSF(opt.NA, rhoxy, eta, opt.lambdaB, rho, ...
                    opt.outer_diameter, opt.inner_diameter);

PSFR = PSFR / sum(PSFR(:));
PSFG = PSFG / sum(PSFG(:));
PSFB = PSFB / sum(PSFB(:));

if useGPU
    PSFR = gather(PSFR);
    PSFG = gather(PSFG);
    PSFB = gather(PSFB);
end

save("PSFR.mat", "PSFR")
save("PSFG.mat", "PSFG")
save("PSFB.mat", "PSFB")

end
