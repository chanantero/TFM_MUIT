function [ result ] = extendMultiply( varargin )
% Sometimes it is tedious to multiply arrays when we need to repeat its
% value along dimensions. For example, if we want to multiply a column
% vector x by the values [1 2 3], we would need to apply repmat to the
% operation:
% repmat(x, 1, 3).*repmat([1 2 3], numel(x), 1)
% This function does that for you

numArrays = nargin;
extArrays = cell(numArrays, 1);
[extArrays{:}] = pointWiseExtend(varargin);

% Concatenate arrays along a superior dimension
catExtArrays = cat(finalNumDim + 1, extArrays{:});

% Multiply
result = prod(catExtArrays, finalNumDim + 1);

end

