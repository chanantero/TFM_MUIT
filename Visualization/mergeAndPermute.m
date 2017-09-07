function N = mergeAndPermute( M, modVec, varargin )
% An array A is modified by changing dimensions and reducing certain
% dimensions by merging them in ndgrid format in specified order
% One mathematical way of expressing it is by a cell vector, with as many
% dimensions as the resulting array B. Each cell contains a vector of
% positive integers (the values are unique). That vector contains the
% dimensions of the original array A that have been merged in ndgrid format
% in the dimension order specified by that vector

p = inputParser;

addOptional(p, 'inverse', false)
addOptional(p, 'sizeA', size(M))

parse(p, varargin{:});
inverse = p.Results.inverse;
sizeA = p.Results.sizeA;

if inverse
    B = M;
    
    ind = mergeAndPermute_ind( sizeA, modVec );
    A = zeros(sizeA);
    A(ind) = B;
    
    N = A;
else
    A = M;
    
    sizeA = size(A);
    ind = mergeAndPermute_ind( sizeA, modVec );
    B = A(ind);
    
    N = B;
end


end


function ind = mergeAndPermute_ind( sizeA, modVec )
% An array A is modified by changing dimensions and reducing certain
% dimensions by merging them in ndgrid format in specified order
% One mathematical way of expressing it is by a cell vector, with as many
% dimensions as the resulting array B. Each cell contains a vector of
% positive integers (the values are unique). That vector contains the
% dimensions of the original array A that have been merged in ndgrid format
% in the dimension order specified by that vector.
numDimA = numel(sizeA);

numDimB = numel(modVec);
sizeB = zeros(1, numDimB);
for d = 1:numDimB
    sizeB(d) = prod(sizeA(modVec{d}));
end

indArray = cell(numDimA, 1);
for d = 1:numDimB
    % Create a matrix that has the size of the matrix B. The variation will
    % be along the dimension of B to which it corresponds
    
    % mergedDim contains the merged dimensions of A
    mergedDim = modVec{d};
    numMergedDim = numel(mergedDim);

    % inds is a cell vector with as many elements as merged dimensions.
    % Each cell contains a double vector with the indexes from one to the
    % size of the dimension in A
    inds = cell(numMergedDim, 1);
    for k = 1:numMergedDim
        inds{k} = 1:sizeA(mergedDim(k));
    end
    
    % M is a cell vector with the ndgrid matrices of inds
    M = cell(numMergedDim, 1);
    [M{:}] = ndgrid(inds{:});
    
    ind = cell(numDimB, 1);
    for k = 1:numDimB
        ind{k} = 1:sizeB(k);
    end
    
    for k = 1:numMergedDim
        % ind is a cell vector with the default indices of B. We are only
        % interested in the vector of the d-th dimension. It is substituted
        % for the vectorized version of M{k}
        ind{d} = M{k}(:);
        aux = cell(numDimB, 1);
        [aux{:}] = ndgrid(ind{:});
        indArray{mergedDim(k)} = aux{d};
    end
    
end

ind = sub2ind(sizeA, indArray{:});
end

