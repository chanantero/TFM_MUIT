function reducArray = reduceDimensionalityGrid(dataArray, points, dimensToReduce, oper)
% Points is in grid format
numDim = numel(points);
dimensToReduce = dimensToReduce(:);

permOrder = [dimensToReduce; find(~ismember((1:numDim)', dimensToReduce))];
permDataArray = permute(dataArray, permOrder);

% Sizes
sizeArray = cellfun(@(p) numel(p), points);
sizeReducedArray = sizeArray; sizeReducedArray(dimensToReduce) = 1;
permSizeRed = sizeReducedArray(permOrder);

% Subindices
subs = cell(1, numDim);
[subs{:}] = ind2sub(permSizeRed, (1:prod(permSizeRed))');
subsCellArray = num2cell(cell2mat(subs));

numRedDims = numel(dimensToReduce);
fullIndRedDims = cell(numRedDims, 1);
for k = 1:numRedDims
    fullIndRedDims{k} = 1:sizeArray(dimensToReduce(k));
end
redPointMats = cell(numRedDims, 1);
[redPointMats{:}] = ndgrid(points{dimensToReduce});

reducArray = cell(permSizeRed);
for k = 1:prod(permSizeRed)
    currInd = subsCellArray(k, :);
    currInd(1:numRedDims) = fullIndRedDims;
    subArray = permDataArray(currInd{:});
    
    % Make wathever you want to do with the subArray
    reducArray{k} = oper(subArray, redPointMats);
end

% Permute to turn it back the the original dimension order
[~, permInv] = sort(permOrder);
reducArray = permute(reducArray, permInv);
end