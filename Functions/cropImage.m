function [cropped, rows, cols] = cropImage(img, cropBorder)
% cropImage  Crop a fixed border when the image is large enough.

[Ny, Nx, ~] = size(img);
if cropBorder > 0 && Ny > 2 * cropBorder && Nx > 2 * cropBorder
    rows = (cropBorder + 1):(Ny - cropBorder);
    cols = (cropBorder + 1):(Nx - cropBorder);
else
    rows = 1:Ny;
    cols = 1:Nx;
end

cropped = img(rows, cols, :);

end
