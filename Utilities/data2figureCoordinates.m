function figCoord = data2figureCoordinates( ax, dataPoints )
% Input arguments:
% - ax. Axes object.
% - dataPoints. (numPoints x 2) or (numPoints x 3) matrix. For 2D axes,
% (numPoints x 2) is enough. Each row is the x and y coordinates. For 3D
% view, (numPoints x 3) matrix is necessary. If (numPoints x 2) matrix is
% provided, the function assumes the third column (z coordinate) to be 0.
% This 3D functionality might not be implementated yet in this version of
% the function.
% The function works with normalized units

[numPoints, numCol] = size(dataPoints);

assert(numCol == 2 || numCol == 3, 'axs2figureCoordinates:wrongNumberOfColumns', 'The number of columns should be 2 or three');

if numCol == 2
    dataPoints = [dataPoints, zeros(numPoints, 1)];
end

dataXlength = abs(ax.XLim(2) - ax.XLim(1));
dataYlength = abs(ax.YLim(2) - ax.YLim(1));
axPos = ax.Position;

x = dataPoints(:, 1);
y = dataPoints(:, 2);
figCoordX = axPos(1) + axPos(3)*(x - ax.XLim(1))/dataXlength;
figCoordY = axPos(2) + axPos(4)*(y - ax.YLim(1))/dataYlength;

figCoord = [figCoordX, figCoordY];

end

