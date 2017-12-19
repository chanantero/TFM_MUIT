function im = changePixelColor(im, colorRGB, ind)
% Input arguments:
% - im. RGB image. (numRows x numCols x 3) array
% - colorRGB. 3 element vector, RGB triplet.
% - ind. (numRows x numCols) logical matrix. Each element corresponds to a
% pixel. If ind(i,j) == true, the (i,j)-th pixel of the image im will be
% set to the color in the colorRGB variable. ind can also be an vector of integers, being each elent the index of one
% of the pixels that must be set to colorRGB.

[numRows, numCols, ~] = size(im);

if islogical(ind)
    ind = find(ind);
end

N = numel(ind);
[rows, cols] = ind2sub([numRows, numCols], N);

indR = sub2ind([numRows, numCols, 3], rows, cols, ones(N, 1));
indG = sub2ind([numRows, numCols, 3], rows, cols, 2*ones(N, 1));
indB = sub2ind([numRows, numCols, 3], rows, cols, 3*ones(N, 1));

im(indR) = colorRGB(1);
im(indG) = colorRGB(2);
im(indB) = colorRGB(3);
    
end