function [ result ] = extendMultiply( varargin )
% Sometimes it is tedious to multiply arrays when we need to repeat its
% value along dimensions. For example, if we want to multiply a column
% vector x by the values [1 2 3], we would need to apply repmat to the
% operation:
% repmat(x, 1, 3).*repmat([1 2 3], numel(x), 1)
% This function does that for you

numArrays = nargin;
arrays = varargin;

sizes = cell(numArrays, 1);
numDims = zeros(numArrays, 1);
for a = 1:numArrays
    sizes{a} = size(arrays{a});
    numDims(a) = ndims(arrays{a});
end

% Extend sizes to the maximum number of dimensions and put them together in
% a matrix
finalNumDim = max(numDims); 
sizesMat = ones(numArrays, finalNumDim);
for a = 1:numArrays
    sizesMat(a, 1:numDims(a)) = sizes{a};
end

% Operation can only be done if in each dimension there are not more than
% two different sizes: 1 and another one. In any other case, the operation
% cannot be done.
for dim = 1:finalNumDim
    if numel(unique(sizesMat(:, dim))) > 2
        error('Operation cannot be done')
    end
end

% Now that we know the operation can be done, we extend the arrays.
% First, we find the size of the result.
finalSize = max(sizesMat, [], 1);

% Then, we extend each array
extArrays = cell(numArrays, 1);
for a = 1:numArrays
    repVec = finalSize;
    repVec(sizesMat(a, :) == finalSize) = 1;
    extArrays{a} = repmat(arrays{a}, repVec);
end

% Concatenate arrays along a superior dimension
catExtArrays = cat(finalNumDim + 1, extArrays{:});

% Multiply
result = prod(catExtArrays, finalNumDim + 1);

end

