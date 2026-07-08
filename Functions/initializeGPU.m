function useGPU = initializeGPU(requestGPU)
% initializeGPU  Try to initialize a CUDA GPU for heavy array computation.
%
% Set requestGPU to false to force the original CPU execution path.

if nargin < 1 || isempty(requestGPU)
    requestGPU = true;
end

useGPU = false;
if ~requestGPU
    fprintf('GPU disabled; using CPU arrays.\n');
    return
end

try
    dev = gpuDevice();
    useGPU = true;
    fprintf('Using GPU: %s (available memory %.2f GB).\n', ...
            dev.Name, dev.AvailableMemory / 1024^3);
catch ME
    warning('HSEDOF:GPUUnavailable', ...
            'GPU unavailable; falling back to CPU. %s', ME.message);
end

end
