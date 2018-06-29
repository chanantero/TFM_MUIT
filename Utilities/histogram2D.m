function ax = histogram2D( A, xDim, xValues, varargin )

p = inputParser;

addOptional(p, 'xEdges', [])
addOptional(p, 'yEdges', [])
addParameter(p, 'visual3D', false)

parse(p, varargin{:});

xEdges = p.Results.xEdges;
yEdges = p.Results.yEdges;
visual3D = p.Results.visual3D;

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
    xValuesMat = pointWiseExtend(permute(xValues(:), [2:xDim, 1, xDim+1]), A); % xDim+1 is put there because if xDim == 1, then permute throws an error
end

N = histcounts2(xValuesMat, A, xEdges, yEdges);

Nnorm = N./repmat(sum(N, 2), [1, size(N, 2)]);

ax = axes(figure);
if ~visual3D
    C = zeros(length(xEdges), length(yEdges));
    C(1:end-1, 1:end-1) = Nnorm;
    pcolor(ax, xEdges, yEdges, C');
else
    bar3cRub(Nnorm, xEdges, yEdges, ax);
end

end


