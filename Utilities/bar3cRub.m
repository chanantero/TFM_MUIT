function sur = bar3cRub(data, xEdges, yEdges, varargin)

p = inputParser;
addOptional(p, 'axes', []);
parse(p, varargin{:})
predefAx = ~ismember('axes', p.UsingDefaults);

[M, N] = size(data);

xEdges = xEdges(:);
xEdges = stretchArray(xEdges, 2);

yEdges = yEdges(:);
yEdges = stretchArray(yEdges, 2);

dataStretch = stretchArray(data, [2 2]); % Own function. Look for it in the Matlab folder of Master 2º curso
Z = zeros(2*M+2, 2*N+2);
Z(2:end-1, 2:end-1) = dataStretch;
CData = Z;

% Modify CData so the lateral side of the bars have the right color
increaseX = find(diff([zeros(1, N); data; zeros(1, N)], 1, 1) > 0);
wallsXDir_x = 1:2:size(CData, 1)-1;
wallsXDir_y = 2:2:size(CData, 2);
[X, Y] = ndgrid(wallsXDir_x, wallsXDir_y);
indWallXmat = sub2ind(size(CData), X, Y);
CData(indWallXmat(increaseX)) = data(increaseX);

increaseY = find(diff([zeros(M, 1), data, zeros(M, 1)], 1, 2) > 0);
wallsYDir_y = 1:2:size(CData, 2)-1;
wallsYDir_x = 2:2:size(CData, 1);
[X, Y] = ndgrid(wallsYDir_x, wallsYDir_y);
indWallYmat = sub2ind(size(CData), X, Y);
CData(indWallYmat(increaseY)) = data(increaseY);

if predefAx
    ax = p.Results.axes;
else
    ax = axes(figure);
end
sur = surf(ax, xEdges, yEdges, Z', CData');
sur.CData(:);

end