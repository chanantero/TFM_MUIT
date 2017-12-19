function [ extendedColormap ] = extendColormap( inputColormap, numOutputColors )

numColInput = size(inputColormap, 1);

colorMap_R = interp1(1:numColInput, inputColormap(:,1), linspace(1, numColInput, numOutputColors)');
colorMap_G = interp1(1:numColInput, inputColormap(:,2), linspace(1, numColInput, numOutputColors)');
colorMap_B = interp1(1:numColInput, inputColormap(:,3), linspace(1, numColInput, numOutputColors)');

extendedColormap = [colorMap_R, colorMap_G, colorMap_B];


end

