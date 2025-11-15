% ThresholdSegmentation  Apply channel-wise adaptive clipping to suppress 
%                        over-saturated regions in RGB images.
%
% This function performs a simple yet effective dynamic-threshold 
% segmentation/clipping operation on an RGB image. For each channel (R, G, B),
% the algorithm determines a cutoff threshold as:
%
%           threshold = Rate × max(channel intensity)
%
% Pixels exceeding this threshold are clipped to the threshold value.
% The entire image is then normalized to [0, 1] to maintain consistency.
%
% Inputs:
%   Img       - Input RGB image, can be uint8 / uint16 / double, size H×W×3.
%   Rate      - Ratio (0–1) that defines how aggressively to clip bright pixels.
%
% Output:
%   result    - Output RGB image after channel-wise threshold clipping and
%             global normalization to [0, 1].
%
% Algorithm:
%   1. Normalize input image globally to [0, 1].
%   2. Split image into R/G/B channels.
%   3. For each channel:
%            Compute threshold = Rate × channel_max
%            Replace values greater than threshold with the threshold.
%   4. Recombine the clipped channels.
%   5. Normalize the final result to [0, 1].
%
% This approach is useful for suppressing extremely bright pixels while
% preserving overall color information, commonly used in illumination
% correction, color normalization, and preprocessing before segmentation.




function result = ThresholdSegmentation(Img, Rate)

    % Global normalization
    M = max(Img, [], 'all');
    m = min(Img, [], 'all');
    Img = (Img - m) / (M - m);

    % Split RGB channels
    R = Img(:, :, 1);
    G = Img(:, :, 2);
    B = Img(:, :, 3);

    % Compute adaptive thresholds
    ThresholdR = max(R, [], 'all') * Rate;
    ThresholdG = max(G, [], 'all') * Rate;
    ThresholdB = max(B, [], 'all') * Rate;

    % Apply clipping
    R(R > ThresholdR) = ThresholdR;
    G(G > ThresholdG) = ThresholdG;
    B(B > ThresholdB) = ThresholdB;

    % Reconstruct clipped image
    result = cat(3, R, G, B);

    % Normalize final output
    result = result / max(result(:));
end
