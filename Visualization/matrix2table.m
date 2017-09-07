function T = matrix2table( parameters, dataMatrix, dimensionOrder, names )

% Check inputs: how many parameters and data matrices are there?
numParameters = numel(parameters);
numDataMatrices = numel(dataMatrix);
assert(numParameters + numDataMatrices == numel(names), 'matrix2table:wrongInput', 'Size of input names is not correct')

% Reorder dimensions
parameters = parameters(dimensionOrder);
if ~isscalar(dimensionOrder)
    dataMatrix = cellfun(@(dataMat) permute(dataMat, dimensionOrder), dataMatrix, 'UniformOutput', false);
end

% Generate parameter columns
mats = cell(1, numParameters);
[mats{:}] = ndgrid(parameters{:});
paramColumns = cellfun(@(mat) mat(:), mats, 'UniformOutput', false);

% Generate data columns
dataColumns = cellfun(@(mat) mat(:), dataMatrix, 'UniformOutput', false);

% Identify column names and reorder them and the parameter columns
parametersNames = names(1:numParameters);
parametersNames = parametersNames(flip(dimensionOrder));
dataNames = names(numParameters+1:end);
names = [parametersNames(:); dataNames(:)];
paramColumns = flip(paramColumns);
   
T = table(paramColumns{:}, dataColumns{:}, 'VariableNames', names);

end

