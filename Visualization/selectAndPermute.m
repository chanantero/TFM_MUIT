function [ selInd, dimOrder, parameterMats_perm, newParameters ] = selectAndPermute( parameters, particularizedValues, completeDimensions, dimOrder )

numDims = numel(parameters);

% Which parameters are selected as independent? selDims
notCompDims = find(~ismember(1:numDims, completeDimensions));

% Generate indices that, applied to the array, select the dimensions correspondent to the selected
% paramteters, particularizing the rest of parameters to the value that
% best approximates the chosen values.
numDims = numel(parameters);
sizes = cellfun(@(param) numel(param), parameters);
selInd = cell(1, numDims);
selInd(completeDimensions) = arrayfun(@(selDim) 1:sizes(selDim) , completeDimensions, 'UniformOutput', false);
selInd(notCompDims) = cellfun(@(param, selValue) findIndBestApprox(param, selValue), parameters(notCompDims), particularizedValues(notCompDims), 'UniformOutput', false);

% Generate dimension order for the permutation so the selected ones are the first ones, and the rest
% are the latest
if(nargin < 4)
    dimOrder = [completeDimensions, notCompDims];
end
% We have what we wanted: selInd and dimOrder

% New parameters
newParameters = arrayfun(@(k) parameters{k}(selInd{k}), 1:numel(parameters), 'UniformOutput', false);
newParameters = newParameters(dimOrder);

% Calculate the parameter arrays
parameterMats_perm = cell(1, numDims);
[parameterMats_perm{:}] = ndgrid(parameters{dimOrder});
parameterMats_perm = cellfun(@(paramMat) paramMat(selInd{dimOrder}), parameterMats_perm, 'UniformOutput', false);

end

