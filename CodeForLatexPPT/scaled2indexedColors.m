function indices = scaled2indexedColors(colormapLength, CLim, C)
cmin = CLim(1);
cmax = CLim(2);

indices = fix((C-cmin)/(cmax-cmin)*colormapLength)+1;
%Clamp values outside the range [1 m]
indices(indices<1) = 1;
indices(indices>colormapLength) = colormapLength;
end