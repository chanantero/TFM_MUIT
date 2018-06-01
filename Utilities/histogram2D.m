function ax = histogram2D( A, xDim, xValues, varargin )

p = inputParser;

addOptional(p, 'xEdges', [])
addOptional(p, 'yEdges', [])

parse(p, varargin{:});

xEdges = p.Results.xEdges;
yEdges = p.Results.yEdges;

if isempty(yEdges);
    yEdges = linspace(min(A(:)), max(A(:)), 40);
end

if isempty(xEdges);
    minX = min(xValues(:));
    maxX = max(xValues(:));
    xEdges = linspace(minX, maxX, 40);
end

if ndims(xValues) == ndims(A) && all(size(xValues) == size(A))
    xValuesMat = xValues;
else    
    xValuesMat = pointWiseExtend(permute(xValues(:), [2:xDim, 1]), A);
end

N = histcounts2(xValuesMat, A, xEdges, yEdges);

Nnorm = N./repmat(sum(N, 2), [1, size(N, 2)]);

ax = axes(figure);
C = zeros(length(xEdges), length(yEdges));
C(1:end-1, 1:end-1) = Nnorm;
pcolor(ax, xEdges, yEdges, C');

end

