function param = defaultMethodParameters(methodName)
% defaultMethodParameters  Return default reconstruction parameters.

methodName = normalizeMethodName(methodName);
param = struct();

switch methodName
    case "RLD"
        param.TV = 1e-4;
        param.nonnegBeta = 0;
        param.niter = 100;

    case "WB-ARL"
        param.TV = 1e-5;
        param.nonnegBeta = 0;
        param.niter = 2;
        param.wb = struct();
        param.wb.bp_type = 'wiener-butterworth';
        param.wb.alpha = 0.1;
        param.wb.beta = 0.1;
        param.wb.n = 20;
        param.wb.resFlag = 1;      % 1: use FWHM of the forward PSF.
        param.wb.iRes = [0, 0, 0]; % Used only when resFlag = 2.
        param.wb.verboseFlag = 0;
        param.wb.epsValue = 1e-12;
end

end
