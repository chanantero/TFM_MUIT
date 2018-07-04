function [ result ] = MdimOperOverNdimArray( operation, array, dims )
% - operation. Function Handle.
% To be completed

numDims = ndims(array);
sizes = size(array);

order = [dim, 1:dim-1, dim+1:numDims];
auxArray = permute(array, order);

sizesCell = [1, sizes(order(2:end))];

resultCell = cell(sizesCell);
for k = 1:numel(resultCell)
    resultCell{k} = operation(auxArray(:, k));
end
result = cell2mat(resultCell);
result = ipermute(result, order);

end
