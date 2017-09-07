function [ selIndices ] = findIndBestApprox( values, desiredValues )
% Find indices that approximate better the desiredValues by values
values = values(:);
desiredValues = desiredValues(:);

[valMat, desirValMat] = ndgrid(desiredValues, values);
[~, selIndices] = min(abs(valMat - desirValMat), [], 2);

end

