function result = postprocessRGB(DeconRGB, RefRGB, thresholdRate, cropBorder)
% postprocessRGB  Normalize, color-transfer, threshold, and crop RGB output.

result = zeros(size(RefRGB));
for c = 1:3
    channel = DeconRGB(:, :, c);
    result(:, :, c) = channel / max(max(channel(:)), eps);
end

result = result / max(max(result(:)), eps);
result = ColorTransfer(result, RefRGB);
result = ThresholdSegmentation(result, thresholdRate);
result = cropImage(result, cropBorder);

end
