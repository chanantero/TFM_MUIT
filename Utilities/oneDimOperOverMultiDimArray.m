function [ result ] = oneDimOperOverMultiDimArray( operation, array, dim )
% - operation. Function Handle. The output will be converted to a column
% vector

numDims = ndims(array);
sizes = size(array);

order = [dim, 1:dim-1, dim+1:numDims];
auxArray = permute(array, order);

sizesCell = [1, sizes(order(2:end))];

resultCell = cell(sizesCell);
for k = 1:numel(resultCell)
    res = operation(auxArray(:, k));
    resultCell{k} = res(:);
end
result = cell2mat(resultCell);
result = ipermute(result, order);

end

