clc
clear 
close all

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: read the image and set parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath("Data\")
addpath("Result\")
addpath("Functions\")

Img = im2double(imread('Raw.tif'));
ImgR = Img(:, :, 1);
ImgG = Img(:, :, 2);
ImgB = Img(:, :, 3);

 
Nx = 1800;
Ny = Nx;
Nz = 60;
n_m = 1;                  % refractive index of objective lens
lambdaR = 0.656;          % typical wavelength of red channel
lambdaG = 0.588;          % typical wavelength of green channel
lambdaB = 0.486;          % typical wavelength of blue channel
lambdaR = lambdaR / n_m; 
lambdaG = lambdaG / n_m; 
lambdaB = lambdaB / n_m; 
Mag = 20;                 % magnification of objective lens
pixelsize = 6.5 / Mag;    % effective spatial sampling size
z_step = 0.25;            % axial sampling
NA = 0.8;                 % numerical aperture of objective lens

inner_diameter = 0.875;   % inner diameter of annular illumination pattern
outer_diameter = 1;       % outer diameter of annular illumination pattern

delta_x = 1 / (pixelsize * Nx);    
delta_y = 1 / (pixelsize * Ny);
delta_z = 1 / (z_step * Nz);

fx = (-fix(Nx / 2): 1: fix((Nx - 1) / 2)) * delta_x;  % x frequency coordinate
fy = (-fix(Ny / 2): 1: fix((Ny - 1) / 2)) * delta_y;  % y frequency coordinate
fz = (-fix(Nz / 2): 1: fix((Nz - 1) / 2)) * delta_z;  % z frequency coordinate

[rhox,rhoy,eta] = meshgrid(fx, fy, fz);
rhoxy = sqrt(rhox .^ 2 + rhoy .^ 2);
rho = sqrt(rhoxy .^ 2 + eta .^ 2);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: Compute PSFs of three RGB channels.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PSFR = CalculatePSF(NA, rhoxy, eta, lambdaR, rho, outer_diameter, inner_diameter);
% PSFG = CalculatePSF(NA, rhoxy, eta, lambdaG, rho, outer_diameter, inner_diameter);
% PSFB = CalculatePSF(NA, rhoxy, eta, lambdaB, rho, outer_diameter, inner_diameter);
% 
% PSFR = PSFR / sum(PSFR(:));
% PSFG = PSFG / sum(PSFG(:));
% PSFB = PSFB / sum(PSFB(:));
% 
% 
% save("PSFR.mat", "PSFR")
% save("PSFG.mat", "PSFG")
% save("PSFB.mat", "PSFB")


% The computed PSFs have already saved for convenient

load("PSFR.mat")
load("PSFG.mat")
load("PSFB.mat")


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: Perform R-L deconvolution with total variation regularization.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TV = 0.0001;        % total variation regularization parameter
beta = 0;      % non-negativity constraint value
niter = 100;   % RLD iteration times

DeconImgR = RLDTV(ImgR, PSFR, niter, beta, TV);
DeconImgG = RLDTV(ImgG, PSFG, niter, beta, TV);
DeconImgB = RLDTV(ImgB, PSFB, niter, beta, TV);



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4: Merge channels and perform color transfer .
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


DeconImg = zeros(Ny, Nx, 3);
DeconImg(:, :, 1) = DeconImgR/max(DeconImgR(:));
DeconImg(:, :, 2) = DeconImgG/max(DeconImgG(:));
DeconImg(:, :, 3) = DeconImgB/max(DeconImgB(:));
DeconImg = DeconImg/max(DeconImg(:));

DeconImg = ColorTransfer(DeconImg, Img);



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 5: Process and save images using threshold segmentation to handle edge artifacts.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Rate = 0.94;
DeconImg = ThresholdSegmentation(DeconImg, Rate);
DeconImg = DeconImg(101: 1700, 101: 1700, :);
imwrite(DeconImg, './Result/Result.tif');

figure,
subplot(1,2,1),imshow(Img(101: 1700, 101: 1700, :),[]);title('Raw HSEDOF data')
subplot(1,2,2),imshow(DeconImg,[]);title('Deconvolution HSEDOF data')


