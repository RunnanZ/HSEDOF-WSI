function printExecutionMode(useGPU, batchRGBOnGPU)
% printExecutionMode  Print the selected CPU/GPU execution mode.

if useGPU
    fprintf('Device: GPU');
    if batchRGBOnGPU
        fprintf(' with batched RGB processing\n');
    else
        fprintf(' with sequential RGB processing\n');
    end
else
    fprintf('Device: CPU\n');
end

end
