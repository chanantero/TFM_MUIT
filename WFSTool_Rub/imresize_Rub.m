function resIm = imresize_Rub(im, newSize)
% Versión arcaica de la versión imresize hecha para ordenadores que no
% tienen la Image Processing Toolbox

width = size(im, 2);
height = size(im, 1);
numChann = size(im, 3); % Number of color channels

newHeight = round(newSize(1));
newWidth = round(newSize(2));

sampWidthVector = round(linspace(1, width, newWidth));
sampHeightVector = round(linspace(1, height, newHeight));


resIm = cast(zeros(newHeight, newWidth, numChann), class(im));


for c = 1:numChann
    resIm(:, :, c) = im(sampHeightVector, sampWidthVector, c);
end

end