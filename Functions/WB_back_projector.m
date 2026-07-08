function [PSF_bp, OTF_bp] = WB_back_projector(PSF_fp, bp_type, alpha, beta, n, resFlag, iRes, verboseFlag)
%This code is based on the work by Min Guo (April 12, 2019).
% Generate Backprojector: PSF and OTF
% % % Output
% PSF_bp: Back projector in spatial domain
% OTF_bp: Back projector in Fourier domain
% % % Input
% PSF_fp: Forward projector
% bp_type: 'wiener', 'wiener-butterworth'
% alpha: 0.0001 ~ 0.001
%      1: use OTF value of PSF_bp at resolution limint; 
%      else: alpha = input alpha;
% beta: 0.001 ~ 0.01
%      1: use OTF value of PSF_bp at resolution limint;
%      else: beta = input beta;
% n: 4 ~ 15
%      order of the Butterworth filter
% resFlag: 
%      0: use PSF_fp FWHM/root(2) as resolution limit (for iSIM);
%      1: use PSF_fp FWHM as resoltuion limit; 
%      2: use input values (iRes) as resoltuion limit;
% iRes: 1 x 3 array
%      input resolution limit in 3 dimensions in terms of pixels;
% verboseFlag:
%      0: hide log and intermediate results 
%      1: show log and intermediate results 

if(nargin == 1)
    bp_type = 'wiener-butterworth';
    alpha = 0.001;
    beta = 1;
    n =10;
    resFlag = 1;
    iRes = [0,0,0];
    verboseFlag = 0;
end
if(nargin == 2)
    alpha = 0.001;
    beta = 1;
    n =10;
    resFlag = 1;
    iRes = [0,0,0];
    verboseFlag = 0;
end
if(nargin ==3)
    beta = 1;
    n =10;
    resFlag = 1;
    iRes = [0,0,0];
    verboseFlag = 0;
end
if(nargin ==4)
    n =10;
    resFlag = 1;
    iRes = [0,0,0];
    verboseFlag = 0;
end
if(nargin ==5)
    resFlag = 1;
    iRes = [0,0,0];
    verboseFlag = 0;
end
if(nargin ==7)
    verboseFlag = 0;
end
% input PSF size and center
[Sx, Sy, Sz] = size(PSF_fp);
PSF_fp = single(PSF_fp);
PSF_fwhm = PSF_fp;
if isa(PSF_fwhm, 'gpuArray')
    PSF_fwhm = gather(PSF_fwhm);
end
Scx = (Sx+1)/2; 
Scy = (Sy+1)/2;
Scz = (Sz+1)/2;
Sox = round((Sx+1)/2); 
Soy = round((Sy+1)/2);
Soz = round((Sz+1)/2);
if(verboseFlag)
    disp(['Back projector type:' bp_type]);
end

%%%% Calculate PSF and OTF size
[FWHMx, FWHMy, FWHMz] = fwhm_PSF(PSF_fwhm);
if(verboseFlag)
    disp(['Forward projector FWHMs:' num2str(FWHMx) ' x ' num2str(FWHMy) ' x ' num2str(FWHMz)]);
end


% normalize flipped PSF: traditional back projector
flippedPSF = flipPSF(PSF_fp);
OTF_flip = fftn(ifftshift(flippedPSF));
OTF_abs = fftshift(abs(OTF_flip));
OTFmax = max(OTF_abs(:)); % find maximum value and position
M =  OTFmax(1);
OTF_abs_norm = OTF_abs/M;

% set resolution cutoff
switch(resFlag)
    case 0 % Set resolution as 1/root(2) of PSF_fp FWHM: iSIM case
        resx = FWHMx/2^0.5;resy = FWHMy/2^0.5;resz = FWHMz/2^0.5; 
    case 1 % Set resolution as PSF_fp FWHM
        resx = FWHMx;resy = FWHMy;resz = FWHMz; 
    case 2 % Set resolution based input values
        resx = iRes(1);resy = iRes(2);resz = iRes(3);
    otherwise
        error('Processing terminated, please set resFlag as 0, 1, or 2')  
end
% pixel size in Fourier domain
px = 1/Sx; py = 1/Sy; pz = 1/Sz; 
% frequency cutoff in terms of pixels
tx = 1/resx/px; ty = 1/resy/py; tz = 1/resz/pz;
if(verboseFlag)
    disp(['Resolution cutoff in spatial domain:' num2str(resx) ' x ' num2str(resy) ' x ' num2str(resz)]);
    disp(['Resolution cutoff in Fourier domain:' num2str(tx) ' x ' num2str(ty) ' x ' num2str(tz)]);
end

%%% Check cutoff gains of traditional back projector
tplane = squeeze(max(OTF_abs_norm,[],3));
tline = max(tplane,[],2);
to1 = max(round(Scx -tx), 1); to2 = min(round(Scx+tx), Sx); 
beta_fpx = (tline(to1) + tline(to2))/2; % OTF frequency intensity at cutoff:x

tplane = squeeze(max(OTF_abs_norm,[],3));
tline = max(tplane,[],1);
to1 = max(round(Scy -ty), 1); to2 = min(round(Scy+ty), Sy); 
beta_fpy = (tline(to1) + tline(to2))/2; % OTF frequency intensity at cutoff:y

tplane = squeeze(max(OTF_abs_norm,[],1));
tline = max(tplane,[],1);
to1 = max(round(Scz -tz), 1); to2 = min(round(Scz+tz), Sz); 
beta_fpz = (tline(to1) + tline(to2))/2; % OTF frequency intensity at cutoff:z

beta_fp = (beta_fpx + beta_fpy + beta_fpz)/3;
if(verboseFlag)
    disp(['Cutoff gain of forward projector:' num2str(gatherIfGPU(beta_fpx)) ' x ' num2str(gatherIfGPU(beta_fpy))...
    ' x ' num2str(gatherIfGPU(beta_fpz)) ', Average = ' num2str(gatherIfGPU(beta_fp))]);
end      

% % % parameter for wiener filter;
if(alpha==1)
    alpha = beta_fp;
    if(verboseFlag)
        disp(['Wiener parameter adjusted as traditional BP cutoff gain: alpha = ' num2str(gatherIfGPU(alpha))]);
    end
else
    if(verboseFlag)
        disp(['Wiener parameter set as input: alpha = ' num2str(gatherIfGPU(alpha))]);
    end
end
if(beta==1)
    beta = beta_fp;
    if(verboseFlag)
        disp(['Cutoff gain adjusted as traditional BP cutoff gain: beta = ' num2str(gatherIfGPU(beta))]);
    end
else
    if(verboseFlag)
        disp(['Cutoff gain set as input: beta = ' num2str(gatherIfGPU(beta))]);
    end
end
% % order of Butterworth filter
% pn = 2*n;
if(verboseFlag)
    disp(['Butterworth order (slope parameter) set as: n = ' num2str(n)]);
end

switch bp_type       
    case 'wiener'
        OTF_flip_norm = OTF_flip/M; % Normalized OTF_flip
        OTF_bp = OTF_flip_norm ./(abs(OTF_flip_norm).^2+alpha); % Wiener filter
        PSF_bp = fftshift(real(ifftn(OTF_bp)));
        
    case 'wiener-butterworth'
        % *** OTF_wiener-butterworth = Winer .* 1/sqrt(1+ee*(kx/kcx)^pn)
        % *** beta = beta_wienerx * 1/sqrt(1+ee) --> ee = beta_wienerx/beta^2 - 1;
        % % create Wiener filter
        OTF_flip_norm = OTF_flip/M;
        OTF_Wiener = OTF_flip_norm ./(abs(OTF_flip_norm).^2+alpha);
        % cutoff gain for winer filter
        OTF_Wiener_abs = fftshift(abs(OTF_Wiener));
        tplane = abs(squeeze(OTF_Wiener_abs(:,:,Soz))); % central slice
        tline = max(tplane,[],2);
        to1 = max(round(Scx -tx), 1); to2 = min(round(Scx+tx), Sx);
        beta_wienerx = (tline(to1) + tline(to2))/2; % OTF frequency intensity at cutoff:x
        if(verboseFlag)
            disp(['Wiener cutoff gain: beta_wienerx = ' num2str(gatherIfGPU(beta_wienerx))]);
        end
        % % create Butteworth Filter
        kcx = tx; % width of Butterworth Filter
        kcy = ty; % width of Butterworth Filter
        kcz = tz; % width of Butterworth Filter
        ee = beta_wienerx/beta^2 - 1;
        kx = reshape(cast(((1:Sx) - Scx) ./ kcx, 'like', PSF_fp), [Sx, 1, 1]);
        ky = reshape(cast(((1:Sy) - Scy) ./ kcy, 'like', PSF_fp), [1, Sy, 1]);
        kz = reshape(cast(((1:Sz) - Scz) ./ kcz, 'like', PSF_fp), [1, 1, Sz]);
        mask = kx.^2 + ky.^2 + kz.^2;
        mask = 1 ./ sqrt(1 + ee .* mask.^n); % w^n = (kx/kcx)^pn
        mask = ifftshift(mask); % Butterworth Filter
        % % % % % create Wiener-Butteworth Filter
        OTF_bp = mask.*OTF_Wiener;% final OTF_bp cutfoff gain: beta
        PSF_bp = fftshift(real(ifftn(OTF_bp)));
    otherwise
        error('bp_type does not match any back-projector type')
end
if(verboseFlag)
    line1 = gatherIfGPU(squeeze(flippedPSF(:,Soy,Soz)));
    line2 = gatherIfGPU(squeeze(PSF_bp(:,Soy,Soz)));
    figure,subplot(1,2,1);
    plot(1:Sx,line1/max(line1(:)),1:Sx,line2/max(line2(:)),'LineWidth',2);
    xlabel('Pixel Position');
    ylabel('Normalized Value');
    legend('Trad. bp',[ bp_type ' bp']);
    title('Back Projector Profiles(Spacial Domain)');
 
    line1 = gatherIfGPU(squeeze(abs(OTF_flip(1:Sox,1,1))));
    line2 = gatherIfGPU(squeeze(abs(OTF_bp(1:Sox,1,1))));
    line3 = line1.*line2;

    Xline1=line1/max(line1(1));%    tra
    Xline2=line2/max(line2(1));%    bp
    Xline3=line3/max(line3(1));%    multi

    XLOGline1=log(Xline1+1);
    XLOGline1=XLOGline1/max(XLOGline1(1));
    XLOGline2=log(Xline2+1);
    XLOGline2=XLOGline2/max(XLOGline2(1));
    XLOGline3=log(Xline3+1);
    XLOGline3=XLOGline3/max(XLOGline3(1));

    Xmatrix = [Xline1, Xline2, Xline3];
    XLOGmatrix = [XLOGline1, XLOGline2, XLOGline3];

    subplot(1,2,2);
    plot(0:Sox-1,line1/max(line1(1)),0:Sox-1,line2/max(line2(1)), ...
        0:Sox-1,line3/max(line3(1)),'LineWidth',2);
    hold on, plot(ones(1,13)*tx,0:0.1:1.2,'--r.');
    xlabel('Pixel Position');
    ylabel('Normalized Value');
    legend('Trad. bp',[ bp_type ' bp'], 'multi. bp','resolution limit');
    title('Back Projectors Profiles(Frequency Domain)');
    disp('Back projector generated!!!');
end
    
 

% % % Functions
function [FWHMx,FWHMy,FWHMz] = fwhm_PSF(PSF, pixelSize, cFlag, fitFlag)
% Feed back the full width at half maximun of the input PSF
% fwhm.m and mygaussfit.m are needed
% cFlag
%       0: use maximum's position as PSF center position
%       1: use matrix's center position as PSF center position
% fitFlag
%       0: no fitting before calculate FWHM
%       1: spine fitting before calculate FWHM
%       2: gaussian fitting before calculate FWHM
% 
if(nargin == 1)
    pixelSize = 1;
    cFlag = 0;
    fitFlag = 0;
end

if(nargin == 2)
    cFlag = 0;
    fitFlag = 0;
end

if(nargin == 3)
    fitFlag = 0;
end

% PSF = PSF - mean(PSF(:));
[Sx,Sy,Sz] = size(PSF);
if((Sx ==1)||(Sy==1)) % 1D input
    x = 1:max(Sx,Sy);
    x = x';
    y = PSF(:);
    FWHMx = fwhm(x, y);
    FWHMy = 0;
    FWHMz = 0;
    else if(Sz == 1) % 2D input
        if(cFlag)  
            indx = floor((Sx+1)/2);
            indy = floor((Sy+1)/2);
        else
            [~, ind] = max(PSF(:)); % find maximum value and position 
            [indx,indy] = ind2sub([Sx,Sy],ind(1));
        end
     
        x = 1:Sx;
        x = x';
        y = PSF(:,indy);
        y = y(:);
        if(fitFlag==1)
            xq = 1:0.1:Sx;
            yq = interp1(x, y, xq, 'spline');
            FWHMx = fwhm(xq, yq);
        elseif(fitFlag==2)
            [sig,~,~] = mygaussfit(x,y);
            FWHMx = sig*2.3548;
        else
            FWHMx = fwhm(x, y);
        end
       
        
        x = 1:Sy;
        x = x';
        y = PSF(indx,:);
        y = y(:);
        if(fitFlag==1)
            xq = 1:0.1:Sx;
            yq = interp1(x, y, xq, 'spline');
            FWHMy = fwhm(xq, yq);
        elseif(fitFlag==2)
            [sig,~,~] = mygaussfit(x,y);
            FWHMy = sig*2.3548;
        else
            FWHMy = fwhm(x, y);
        end
     
        FWHMz = 0;
     else % 3D input
         if(cFlag)  
            indx = floor((Sx+1)/2);
            indy = floor((Sy+1)/2);
            indz = floor((Sz+1)/2);
        else
            [~, ind] = max(PSF(:)); % find maximum value and position 
            [indx,indy,indz] = ind2sub([Sx,Sy,Sz],ind(1));
        end
        
        
        x = 1:Sx;
        x = x';
        y = PSF(:,indy,indz);
        y = y(:);
        if(fitFlag==1)
            xq = 1:0.1:Sx;
            yq = interp1(x, y, xq, 'spline');
            FWHMx = fwhm(xq, yq);
        elseif(fitFlag==2)
            [sig,~,~] = mygaussfit(x,y);
            FWHMx = sig*2.3548;
        else
            FWHMx = fwhm(x, y);
        end
        x = 1:Sy;
        x = x';
        y = PSF(indx,:,indz);
        y = y(:);
        if(fitFlag==1)
            xq = 1:0.1:Sy;
            yq = interp1(x, y, xq, 'spline');
            FWHMy = fwhm(xq, yq);
        elseif(fitFlag==2)
            [sig,~,~] = mygaussfit(x,y);
            FWHMy = sig*2.3548;
        else
            FWHMy = fwhm(x, y);
        end
        
        x = 1:Sz;
        x = x';
        y = PSF(indx,indy,:);
        y = y(:);
        if(fitFlag==1)
            xq = 1:0.1:Sz;
            yq = interp1(x, y, xq, 'spline');
            FWHMz = fwhm(xq, yq);
        elseif(fitFlag==2)
            [sig,~,~] = mygaussfit(x,y);
            FWHMz = sig*2.3548;
        else
            FWHMz = fwhm(x, y);
        end
%         FWHMz = fwhm(x, y);
    end
end

FWHMx = FWHMx*pixelSize;
FWHMy = FWHMy*pixelSize;
FWHMz = FWHMz*pixelSize;



function width = fwhm(x,y)
% Full-Width at Half-Maximum (FWHM) of the waveform y(x)
% and its polarity.
% The FWHM result in 'width' will be in units of 'x'
%
% Rev 1.2, April 2006 (Patrick Egan)

y = y / max(y);
N = length(y);
lev50 = 0.5;
if y(1) < lev50                  % find index of center (max or min) of pulse
    [garbage,centerindex]=max(y);
    Pol = +1;
%     disp('Pulse Polarity = Positive')
else
    [garbage,centerindex]=min(y);
    Pol = -1;
%     disp('Pulse Polarity = Negative')
end
i = 2;
while sign(y(i)-lev50) == sign(y(i-1)-lev50)
    i = i+1;
end                                   %first crossing is between v(i-1) & v(i)
interp = (lev50-y(i-1)) / (y(i)-y(i-1));
tlead = x(i-1) + interp*(x(i)-x(i-1));
i = centerindex+1;                    %start search for next crossing at center
while ((sign(y(i)-lev50) == sign(y(i-1)-lev50)) & (i <= N-1))
    i = i+1;
end
if i ~= N
    Ptype = 1;  
%     disp('Pulse is Impulse or Rectangular with 2 edges')
    interp = (lev50-y(i-1)) / (y(i)-y(i-1));
    ttrail = x(i-1) + interp*(x(i)-x(i-1));
    width = ttrail - tlead;
else
    Ptype = 2; 
%     disp('Step-Like Pulse, no second edge')
    ttrail = NaN;
    width = NaN;
end


function outPSF = flipPSF(inPSF)
% function outPSF = flipPSF(inPSF)
% outPSF(i,j,k) = inPSF(m-i+1,n-j+1,l-k+1);
outPSF = inPSF(end:-1:1, end:-1:1, end:-1:1);


function value = gatherIfGPU(value)
if isa(value, 'gpuArray')
    value = gather(value);
end

