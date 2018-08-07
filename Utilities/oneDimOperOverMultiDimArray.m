function [ result ] = oneDimOperOverMultiDimArray( operation, array, dim )
% - operation. Function Handle. The output will be converted to a column
% vector
memoryFriendly = false;

numDims = ndims(array);
sizes = size(array);

order = [dim, 1:dim-1, dim+1:numDims];
try
    auxArray = permute(array, order);
catch ME
    if strcmp(ME.identifier, 'MATLAB:nomem')
        memoryFriendly = true;
    else
        rethrow(ME);
    end
end

if memoryFriendly
    sizesCell = sizes;
    sizesCell(dim) = 1;
    
    resultCell = cell(sizesCell);
    inds = cell(numDims, 1);
    for k = 1:numel(resultCell)
        [inds{:}] = ind2sub(sizesCell, k);
        inds{dim} = 1:sizes(dim);
        res = operation(squeeze(array(inds{:})));
        resultCell{k} = ipermute(res(:), order); 
    end
    result = cell2mat(resultCell);
else
    sizesCell = [1, sizes(order(2:end))];

    resultCell = cell(sizesCell);
    for k = 1:numel(resultCell)
        res = operation(auxArray(:, k));
        resultCell{k} = res(:);
    end
    result = cell2mat(resultCell);
    result = ipermute(result, order);
end

end

