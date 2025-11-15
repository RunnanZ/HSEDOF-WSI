% CalculatePSF Calculate Point Spread Function for annular illumination
%
% This function computes the 3D Point Spread Function (PSF) for annular
% illumination in bright-field microscopy by calculating the 3D PSF and
% then projecting it onto a 2D plane through summation.
%
% Inputs:
%   NA              - Numerical Aperture of the objective lens
%   rhoxy           - Spatial frequency coordinates in XY plane [1/¦Ìm]
%   eta             - Spatial frequency coordinates in Z direction [1/¦Ìm]
%   lambda          - Wavelength of light [¦Ìm]
%   rho             - 3D spatial frequency coordinates [1/¦Ìm]
%   outer_diameter  - Outer diameter of annular aperture (normalized)
%   inner_diameter  - Inner diameter of annular aperture (normalized)
%
% Output:
%   PSF             - 2D projected Point Spread Function
%
% Algorithm:
%   1. Calculate maximum spatial frequency (rhop) based on NA and lambda
%   2. Compute outer ring AOTF
%   3. Compute inner ring AOTF
%   4. Calculate net annular AOTF as difference: outer - inner
%   5. Perform inverse Fourier transform to get 3D APSF
%   6. Project 3D APSF onto 2D plane by summing along axial dimension


function PSF = CalculatePSF(NA, rhoxy, eta,lambda,rho, outer_diameter, inner_diameter)
    rhop = NA / lambda;
    outer_rhos = outer_diameter * rhop;
    outer_Fp = CalculateF(rhoxy,eta,lambda,rho,outer_rhos,rhop);
    outer_Fm = CalculateF(rhoxy,-eta,lambda,rho,outer_rhos,rhop);
    outer_AOTF = (outer_Fp + outer_Fm) .* lambda ./ (4 * pi);
    
    if inner_diameter == 0
        inner_AOTF = zeros(size(outer_AOTF));
    else
        inner_rhos = inner_diameter * rhop;
        inner_Fp = CalculateF(rhoxy,eta,lambda,rho,inner_rhos,rhop);
        inner_Fm = CalculateF(rhoxy,-eta,lambda,rho,inner_rhos,rhop);
        inner_AOTF = (inner_Fp + inner_Fm) .* lambda ./ (4 * pi);
    end
    
    AOTF = outer_AOTF - inner_AOTF;
    APSF = fftshift(ifftn(ifftshift(AOTF)));
    APSF= max(APSF,0);

    [Ny, Nx, Nz] = size(APSF);
    PSFsum = zeros(Ny,Nx);
    for i = 1:Nz
        PSFsum = PSFsum + APSF(:, :, i);
    end

    PSF = PSFsum(901-101:901+101,901-101:901+101);
    OTF = psf2otf(PSF,[Ny Nx]);
    PSF = otf2psf(OTF,[Ny Nx]);
end