function y = shiftArray(x, shiftPos, varargin)
% It shifts array 'x', 'shiftPos' positions, along dimension dim

p = inputParser;
addOptional(p, 'dim', [])
addOptional(p, 'newLength', [])
parse(p, varargin{:})

dim = p.Results.dim;
newLength = p.Results.newLength;

if isempty(dim)
    dim = find(size(x) > 1, 1, 'first');
end

N = size(x, dim);

if isempty(newLength)
    newLength = N;
    ySize = size(x);
else
    ySize = size(x);
    ySize(dim) = newLength;
end

numDims = ndims(x);
order = [dim, 1:dim-1, dim+1:numDims];
newSize = ySize(order);
yPerm = zeros(newSize);
xPerm = permute(x, order);

firstYind = max(1, 1 + shiftPos);
lastYind = min(newLength, N + shiftPos);

firstXind = max(1, 1 - shiftPos);
lastXind = min(N, newLength - shiftPos);

yPerm(firstYind:lastYind, :) = xPerm(firstXind:lastXind, :);

y = ipermute(yPerm, order);

end

